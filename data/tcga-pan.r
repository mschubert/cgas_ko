library(dplyr)
sys = import('sys')
tcga = import('data/tcga')
idmap = import('process/idmap')

sys$run({
    args = sys$cmd$parse(
        opt('o', 'outfile', 'rds', 'tcga-pan.rds')
    )

    pur = tcga$purity_estimate() %>% select(Sample, cohort, purity) # aran2015 different(!)
    incl = sort(unique(pur$cohort))

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
                  os_days = pmax(days_to_death,
                                 days_to_last_known_disease_status,
                                 days_to_last_follow_up,
                                 na.rm=TRUE))

    scores = c(lapply(incl, tcga$gsva, setname="MSigDB_Hallmark_2020"),
               lapply(incl, tcga$gsva, setname="CIN")) %>%
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
        left_join(scores %>% select(Sample,
                                    `CIN70_Carter2006`,
                                    `Interferon Gamma Response`,
                                    `IL-6/JAK/STAT3 Signaling`)) %>%
        select(-patient)

    saveRDS(dset, file=args$outfile)
})
