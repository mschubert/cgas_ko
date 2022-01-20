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
    prot_mat = depmap::depmap_proteomic() %>% util$long2wide("protein_expression", meta)
    rownames(prot_mat) = paste0("prot.", rownames(prot_mat))
    gpmat = narray::stack(gex_mat, prot_mat, along=1)

#    dresp = readr::read_csv("secondary-screen-dose-response-curve-parameters.csv")
#    dmat = util$long2wide(dplyr::rename(dresp, gene_name=name, cell_line=ccle_name), "auc")

    dresp = readxl::read_xlsx("41586_2020_3114_MOESM8_ESM.xlsx")
    dmat = t(as.matrix(setNames(dresp$`Reversine AUC`, dresp$CCLE_ID)))
    rownames(dmat) = "Reversine"

    if (args$tissue == "BRCA") {
        brca_lines = meta %>% filter(primary_disease == "Breast Cancer") %>% pull(cell_line)
        gex_mat = gex_mat[,colnames(gex_mat) %in% brca_lines]
        gpmat = gpmat[,colnames(gpmat) %in% brca_lines]
    }

    smat = switch(args$set,
        genes = gpmat[c("CGAS", "IL6", "IL6R", "IL6ST", "RELB", "RELA",
                        "prot.MB21D1", "prot.IL6ST", "prot.RELA", "prot.RELB"),],
        CIN = util$cin(meta, gex_mat),
        GSVA::gsva(gex_mat, gset$get_human(args$set))
    )

    res = util$do_test(meta, dmat, smat) %>%
        dplyr::rename(mcol=dcol, dcol=mcol)
    saveRDS(res, file=args$outfile)
})
