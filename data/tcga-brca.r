library(dplyr)
sys = import('sys')
tcga = import('data/tcga')
gset = import('genesets')
idmap = import('process/idmap')

load_expr = function(x) {
    ensg = c("ENSG00000160712", "ENSG00000136244", "ENSG00000134352", "ENSG00000164430")
    re = tcga$rna_seq(x, trans="vst")[ensg,] %>%
        tcga$filter(cancer=TRUE, primary=TRUE) %>% t()
    tibble(Sample=rownames(re), IL6R=re[,1], IL6=re[,2], IL6ST=re[,3], CGAS=re[,4])
}

brca_gsva = function() {
    gex = tcga$rna_seq("BRCA", trans="vst")
    rownames(gex) = idmap$gene(rownames(gex), to="hgnc_symbol")
    sets = gset$get_human(c("MSigDB_Hallmark_2020", "CIN")) %>% unname() %>% do.call(c, .)
    t(GSVA::gsva(gex, sets))
}

load_thorsson = function(cohorts="BRCA") {
    immune_df = tcga$immune() %>%
        filter(cohort %in% cohorts) %>%
        mutate(Sample = paste0(barcode, "-01A")) %>%
#        select(-cohort, -OS, -`OS Time`, -PFI, -`PFI Time`) %>%
        as.data.frame()
    immune_df
}

brca_meta = function() {
    meta = tcga$rna_seq("BRCA", annot=TRUE) %>%
        SummarizedExperiment::colData() %>%
        as.data.frame() %>% as_tibble() %>%
        filter(shortLetterCode == "TP") %>% # can compare vs NT, TM here
        transmute(Sample = sample,
                  ER = case_when(
                      subtype_ER.Status == "Positive" ~ 1,
                      subtype_ER.Status == "Negative" ~ 0,
                      TRUE ~ NA_real_
                  ),
                  PR = case_when(
                      subtype_PR.Status == "Positive" ~ 1,
                      subtype_PR.Status == "Negative" ~ 0,
                      TRUE ~ NA_real_
                  ),
                  HER2 = case_when(
                      subtype_HER2.Final.Status == "Positive" ~ 1,
                      subtype_HER2.Final.Status == "Negative" ~ 0,
                      TRUE ~ NA_real_
                  ),
                  ER_PR = pmax(ER, PR),
                  TNBC = pmin(ER, PR, HER2),
                  tumor_stage, age_at_diagnosis
        )
}

load_cohort = function(cohort) {
    scores = brca_gsva() %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Sample") %>%
        as_tibble() %>%
        mutate(Sample = gsub(".", "-", Sample, fixed=TRUE))

    brca = tcga$aneuploidy(cohort) %>%
        inner_join(load_expr(cohort)) %>%
        inner_join(tcga$purity()) %>%
        left_join(brca_meta()) %>%
        left_join(load_thorsson() %>% select(Sample, `Immune Subtype`, `OS Time`, `OS`)) %>%
        left_join(scores) %>%
        mutate(aneuploidy = aneup_log2seg / estimate) # cancer aneup - stroma
}

sys$run({
    args = sys$cmd$parse(
        opt('o', 'outfile', 'rds', 'tcga-brca.rds')
    )

    brca = load_cohort("BRCA")
    saveRDS(brca, file=args$outfile)
})
