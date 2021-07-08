library(dplyr)
library(DESeq2)
sys = import('sys')
idmap = import('process/idmap')
util = import('../rev_effect/util')

batch_anchors = function(has_batch=TRUE) {
    if (!has_batch)
        return(tibble(sample_id="invalid", anchor=NA))
    batch_anchors = c(
        "SU_BT549_08_S15"=1, "SU_BT549_22_S29"=1, # batch 1, stat1+dmso
        "SU_CH210215_17_S35"=1, # batch 4, stat1+dmso
        "SU_CH210215_15_S33"=2, # batch 4, cgas+il6
        "SU_cGAS_KO_IL6_1_S21"=2, "SU_cGAS_KO_IL6_2_S22"=2, # batch 3, cgas+il6
        "SU_CH210215_13_S31"=3, # batch 4, wt+il6Ab
        "SU_WT_anti_IL6_48h_S9"=3 # batch 3, wt+il16Ab
    )
    tibble(sample_id = names(batch_anchors),
           anchor = factor(batch_anchors))
}

deseq_one = function(rec, name=NULL, contrast_only=FALSE) {
    message(name)
    ba = batch_anchors(grepl("batch", rec$design))
    keep = as.data.frame(rec$samples) %>% mutate(in_comparison=1) %>%
        right_join(as.data.frame(colData(eset))) %>%
        filter(!is.na(in_comparison) | sample_id %in% ba$sample_id) %>%
        left_join(ba) %>%
        mutate(anchor = factor(ifelse(is.na(anchor), "0", anchor)))
    eset2 = eset[,keep$sample_id]
    colData(eset2) = DataFrame(keep)

    for (v in c("batch", "genotype", "treatment")) {
        colData(eset2)[[v]] = droplevels(factor(colData(eset2)[[v]]))
        if (v %in% names(rec$samples))
            colData(eset2)[[v]] = relevel(colData(eset2)[[v]], rev(rec$samples[[v]])[1])
    }
    for (v in c("genotype", "treatment")) {
        colData(eset2)[[v]][is.na(eset2$in_comparison)] = levels(colData(eset2)[[v]])[1]
        colData(eset2)[[v]] = droplevels(factor(colData(eset2)[[v]]))
    }
    design(eset2) = as.formula(sub("batch", "batch + anchor", rec$design, fixed=TRUE))

    mm = model.matrix(design(eset2), colData(eset2))
    cmp = cbind(keep[c("batch", "genotype", "treatment")], mm) %>%
        as.data.frame() %>% tibble::rownames_to_column("sample_id") %>% as_tibble()
    if (contrast_only)
        return(cmp)

    res = DESeq2::DESeq(eset2) %>%
        DESeq2::results(name=rec$extract) %>%
        as.data.frame() %>% tibble::rownames_to_column("ensembl_gene_id") %>%
        as_tibble() %>%
        mutate(gene_name = idmap$gene(ensembl_gene_id, to="external_gene_name")) %>%
        select(ensembl_gene_id, gene_name, everything()) %>%
        arrange(padj, pvalue)

    list(compare=cmp, genes=res)
}

args = sys$cmd$parse(
    opt('e', 'eset', 'rds', '../data/rnaseq.rds'),
    opt('o', 'outfile', 'rds', 'deseq.rds'),
    arg('setfiles', 'rds', arity='*', '../data/genesets/human/MSigDB_Hallmark_2020.rds')
)

cfg = yaml::read_yaml("comps.yaml")$comparisons

sets = sapply(args$setfiles, readRDS, simplify=FALSE)
names(sets) = basename(tools::file_path_sans_ext(names(sets)))

eset = readRDS(args$eset)
eset = eset[, eset$time == "48"]
eset$treatment = sub("none", "dmso", eset$treatment)
eset$treatment = relevel(factor(eset$treatment), "dmso")

res = tibble(name=names(cfg), rec=unname(cfg)) %>%
    head(2) %>% #DEBUG
    rowwise() %>%
        mutate(res = list(deseq_one(rec, name)),
               compare = list(res$compare),
               genes = list(res$genes)) %>%
    ungroup() %>%
    select(-res, -rec)



for (ns in names(sets))
    res[[ns]] = lapply(res$genes, util$test_gsets, sets=sets[[ns]])

saveRDS(res, file=args$outfile)
