library(dplyr)
library(depmap)
library("ExperimentHub")
sys = import('sys')
gset = import('genesets')
util = import('./cgasdep')

sys$run({
    args = sys$cmd$parse(
        opt('t', 'tissue', 'pan|BRCA', 'BRCA'),
        opt('d', 'dep', 'rnai|ko', 'ko'),
        opt('s', 'set', 'gene set identifier', 'CIN'),
        opt('o', 'outfile', 'rds', 'cgasdep/MSigDB_Hallmark_2020/BRCA-ko.rds')
    )

    meta = depmap::depmap_metadata()
    dmat = switch(args$dep, rnai=depmap::depmap_rnai(), ko=depmap::depmap_crispr()) %>%
        util$long2wide("dependency", meta)
    gex_mat = depmap::depmap_TPM() %>% util$long2wide("rna_expression", meta)

    if (args$tissue == "BRCA") {
        brca_lines = meta %>% filter(primary_disease == "Breast Cancer") %>% pull(cell_line)
        gex_mat = gex_mat[,colnames(gex_mat) %in% brca_lines]
    }

    smat = switch(args$set,
        genes = dmat,
        CIN = util$cin(meta, dmat),
        GSVA::gsva(dmat, gset$get_human(args$set))
    )

    dep_genes = c("CGAS", "MB21D1", "IL6", "IL6R", "IL6ST", "RELB", "RELA")

    res = util$do_test(meta, smat, gex_mat, dgenes=dep_genes)
    saveRDS(res, file=args$outfile)
})
