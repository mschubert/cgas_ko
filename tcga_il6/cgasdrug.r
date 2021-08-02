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
        opt('o', 'outfile', 'rds', 'cgasdep/MSigDB_Hallmark_2020/BRCA-drug.rds')
    )

    meta = depmap::depmap_metadata()
    gex_mat = depmap::depmap_TPM() %>% util$long2wide("rna_expression", meta)

    dresp = readr::read_csv("secondary-screen-dose-response-curve-parameters.csv")
    dmat = util$long2wide(dplyr::rename(dresp, gene_name=name, cell_line=ccle_name), "auc")

    if (args$tissue == "BRCA") {
        brca_lines = meta %>% filter(primary_disease == "Breast Cancer") %>% pull(cell_line)
        gex_mat = gex_mat[,colnames(gex_mat) %in% brca_lines]
    }

    smat = switch(args$set,
        genes = gex_mat[c("CGAS", "MB21D1", "IL6", "IL6R", "IL6ST", "RELB", "RELA"),],
        CIN = util$cin(meta, gex_mat),
        GSVA::gsva(gex_mat, gset$get_human(args$set))
    )

    res = util$do_test(meta, dresp, smat)
    saveRDS(res, file=args$outfile)
})
