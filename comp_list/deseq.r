library(dplyr)
library(DESeq2)
sys = import('sys')
idmap = import('process/idmap')
util = import('../rev_effect/util')

deseq_one = function(rec, name=NULL) {
    message(name)
    eset2 = eset
    for (i in seq_along(rec$samples)) {
        key = names(rec$samples)[i]
        value = rec$samples[[i]]
        eset2 = eset2[,colData(eset2)[[key]] %in% value]
        colData(eset2)[[key]] = factor(colData(eset2)[[key]], levels=rev(value))
    }
    design(eset2) = as.formula(rec$design)

    mm = model.matrix(design(eset2), colData(eset2))
    cmp = cbind(colData(eset2)[c("batch", "genotype", "treatment")], mm) %>%
        as.data.frame() %>% tibble::rownames_to_column("sample_id") %>% as_tibble()

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
eset$cond = sub("\\+?(rev|dmso)", "",
                         paste(eset$genotype, eset$treatment, sep="+"))
eset$cond = relevel(factor(eset$cond), "wt")
eset$rev = ifelse(grepl("rev", eset$treatment), 1, 0)


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
