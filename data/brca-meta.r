library(dplyr)
library(TCGAbiolinks)

query = GDCquery(project = "TCGA-BRCA",
                 data.category = "Clinical",
                 data.type = "Clinical Supplement",
                 data.format = "BCR Biotab")
GDCdownload(query)
clin = GDCprepare(query)

res = clin$clinical_patient_brca %>%
    transmute(patient = bcr_patient_barcode,
              ER = case_when(
                  er_status_by_ihc %in% c("Negative", "Positive") ~ er_status_by_ihc,
                  TRUE ~ "Unknown"),
              PR = case_when(
                  er_status_by_ihc %in% c("Negative", "Positive") ~ pr_status_by_ihc,
                  TRUE ~ "Unknown"),
              `ER/PR` = case_when(
                  ER == "Positive" | PR == "Positive" ~ "Positive",
                  ER == "Unknown" | PR == "Unknown" ~ "Unknown",
                  TRUE ~ "Negative"),
              HER2 = case_when(
                  her2_status_by_ihc %in% c("Negative", "Positive") ~ her2_status_by_ihc,
                  TRUE ~ "Unknown")) %>%
    filter(grepl("^TCGA-", patient))

saveRDS(res, file="brca-meta.rds")
