library(dplyr)
sys = import('sys')
tcga = import('data/tcga')
idmap = import('process/idmap')

gex_tmm = function(cohort, genes) {
    re = tcga$rna_seq(cohort, trans="tmm")[genes,]
    rownames(re) = names(genes)
    re
}

sys$run({
    args = sys$cmd$parse(
        opt('o', 'outfile', 'rds', 'tcga-pan.rds')
    )

    incl = c("BRCA", "LUAD", "LUSC", "OV", "COAD", "SKCM")
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
                  vital_status = factor(vital_status, levels=c("alive", "dead")),
                  os_days = pmax(0,
                                 days_to_death,
                                 days_to_last_known_disease_status,
                                 days_to_last_follow_up,
                                 na.rm=TRUE))

    scores = c(lapply(incl, gex_tmm, genes=genes),
               lapply(incl, function(i) tcga$gsva(i, "MSigDB_Hallmark_2020")[hms,]),
               lapply(incl, function(i) tcga$gsva(i, "CIN")[CIN,,drop=FALSE])) %>%
        narray::stack(along=2) %>% t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Sample") %>%
        as_tibble()

    aneup = tcga$aneuploidy() %>%
       tcga$filter(along="Sample")

    dset = pur %>%
        mutate(patient = substr(Sample, 1, 12)) %>%
        inner_join(meta) %>%
        left_join(aneup) %>%
        left_join(scores) %>%
        select(-patient)

    saveRDS(dset, file=args$outfile)
})
