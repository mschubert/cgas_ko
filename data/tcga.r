library(dplyr)
sys = import('sys')
tcga = import('data/tcga')
idmap = import('process/idmap')

gex_vst = function(cohort, genes) {
    re = tcga$rna_seq(cohort, trans="vst")[genes,]
    rownames(re) = names(genes)
    re
}

brca_meta = function() {
    query = TCGAbiolinks::GDCquery(project="TCGA-BRCA", data.category="Clinical",
                                   data.type="Clinical Supplement", data.format="BCR Biotab")
    TCGAbiolinks::GDCdownload(query)
    clin = TCGAbiolinks::GDCprepare(query)

    clin$clinical_patient_brca %>%
        transmute(patient = bcr_patient_barcode,
                  `ER/PR` = case_when(
                      er_status_by_ihc == "Positive" | pr_status_by_ihc == "Positive" ~ "Positive",
                      er_status_by_ihc == "Unknown" | pr_status_by_ihc == "Unknown" ~ "Unknown",
                      TRUE ~ "Negative"),
                  HER2 = case_when(
                      her2_status_by_ihc %in% c("Negative", "Positive") ~ her2_status_by_ihc,
                      TRUE ~ "Unknown")) %>%
        filter(grepl("^TCGA-", patient))
}

sys$run({
    args = sys$cmd$parse(
        opt('o', 'outfile', 'rds', 'tcga.rds')
    )

#    incl = c("BRCA", "LUAD", "LUSC", "OV", "COAD", "SKCM")
    incl = c("BLCA", "BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "SKCM", "READ", "UCEC")
    genes = c(IL6R="ENSG00000160712", IL6="ENSG00000136244", IL6ST="ENSG00000134352", CGAS="ENSG00000164430")
    hms = c("E2F Targets", "Interferon Gamma Response", "IL-6/JAK/STAT3 Signaling")
    CIN = c("CIN70_Carter2006")

    pur = tcga$purity_aran2015() %>% select(Sample, cohort, purity=estimate) %>%
        filter(cohort %in% incl, substr(Sample, 14, 16) == "01A")

    meta = tcga$clinical() %>%
        transmute(patient = submitter_id,
                  sex = gender,
                  stage = case_when(
                      grepl("iv", tumor_stage) ~ "iv",
                      grepl("iii", tumor_stage) ~ "iii",
                      grepl("i", tumor_stage) ~ "i/ii",
                      TRUE ~ NA_character_),
                  age_days = age_at_diagnosis,
                  vital_status = c(alive=0L, dead=1L)[vital_status],
                  os_years = pmax(0,
                                  days_to_death,
                                  days_to_last_known_disease_status,
                                  days_to_last_follow_up,
                                  na.rm=TRUE) / 365,
                  vital_status = ifelse(os_years > 5, 0L, vital_status),
                  os_years = pmin(os_years, 5L))

    scores = c(lapply(incl, gex_vst, genes=genes),
               lapply(incl, function(i) tcga$gsva(i, "MSigDB_Hallmark_2020")[hms,]),
               lapply(incl, function(i) tcga$gsva(i, "CIN")[CIN,,drop=FALSE])) %>%
        narray::stack(along=2) %>% t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Sample") %>%
        as_tibble()

    aneup = tcga$aneuploidy() %>%
        select(Sample, aneuploidy=aneup_log2seg)

    dset = pur %>%
        mutate(patient = substr(Sample, 1, 12)) %>%
        inner_join(meta) %>%
        left_join(aneup) %>%
        left_join(scores) %>%
        left_join(brca_meta()) %>%
        select(-patient)

    saveRDS(dset, file=args$outfile)
})
