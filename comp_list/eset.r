library(dplyr)
library(DESeq2)
sys = import('sys')

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

make_eset = function(rec) {
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
        cur = colData(eset2)[[v]]
        cmp_terms = unique(cur[!is.na(eset2$in_comparison)])
        colData(eset2)[[v]][!cur %in% cmp_terms] = levels(cur)[1]
        colData(eset2)[[v]] = droplevels(factor(cur))
    }
    design(eset2) = as.formula(sub("batch", "batch + anchor", rec$design, fixed=TRUE))

    mm = model.matrix(design(eset2), colData(eset2))
    cmp = cbind(keep[c("batch", "genotype", "treatment")], mm) %>%
        as.data.frame() %>% tibble::rownames_to_column("sample_id") %>% as_tibble()

    attr(eset2, "extract") = rec$extract
    attr(eset2, "cmp") = cmp
    eset2
}

sys$run({
    args = sys$cmd$parse(
        opt('c', 'config', 'yaml', 'comps.yaml'),
        opt('k', 'key', 'key in yaml', 'BT549 WT reversine vs BT549 cGAS KO reversine'),
        opt('i', 'infile', 'rds', '../data/rnaseq.rds'),
        opt('o', 'outfile', 'rds', 'deseq.rds')
    )

    cfg = yaml::read_yaml("comps.yaml")
    rec = cfg$comparisons[[args$key]]

    eset = readRDS(args$infile)
    eset = eset[, eset$time == "48" | eset$treatment == "ifng"]
    eset$treatment = sub("none", "dmso", eset$treatment)
    eset$treatment = relevel(factor(eset$treatment), "dmso")

    eset2 = make_eset(rec)
    saveRDS(eset2, file=args$outfile)
})
