library(dplyr)
library(ggplot2)
tcga = import('data/tcga')

load_expr = function(x) {
    ensg = c("ENSG00000160712", "ENSG00000136244", "ENSG00000134352", "ENSG00000164430")
    re = tcga$rna_seq(x, trans="vst")[ensg,] %>%
        tcga$filter(cancer=TRUE, primary=TRUE) %>% t()
    tibble(Sample=rownames(re), IL6R=re[,1], IL6=re[,2], IL6ST=re[,3], CGAS=re[,4])
}

#todo: er/pr/her2 status? (@brca only script)

load_thorsson = function() {
    immune_df = tcga$immune() %>%
        filter(cohort %in% cohorts) %>%
#        select(-cohort, -OS, -`OS Time`, -PFI, -`PFI Time`) %>%
        as.data.frame()
#    immune = immune_df[-1]
#    rownames(immune) = paste0(immune_df$barcode, "-01A")
#    colnames(immune) = make.names(colnames(immune))
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
                  TNBC = pmin(ER, PR, HER2)
        )
}

cohorts = c("BRCA", "LUAD", "COAD", "OV")
expr = sapply(cohorts, load_expr, simplify=FALSE) %>% bind_rows(.id="cohort")
pur = tcga$purity()
#pur = tcga$purity_estimate()
aneup = lapply(cohorts, tcga$aneuploidy) %>%
    bind_rows() %>%
    inner_join(expr) %>%
    inner_join(pur) %>%
    left_join(brca_meta()) %>%
    mutate(aneuploidy = aneup_log2seg / estimate) # cancer aneup - stroma
immune = load_thorsson()
#TODO: add different measures of aneup, cor with IL6/R/ST,CGAS

brca = aneup %>% filter(cohort == "BRCA")
lm(IL6 ~ estimate + aneuploidy, data=brca) %>% broom::tidy()
lm(IL6 ~ estimate * aneuploidy, data=brca) %>% broom::tidy()
lm(IL6R ~ estimate + aneuploidy, data=brca) %>% broom::tidy()
lm(CGAS ~ aneuploidy, data=brca) %>% broom::tidy()
lm(CGAS ~ estimate + aneuploidy, data=brca) %>% broom::tidy()

tnbc = aneup %>% filter(TNBC == 1)
lm(IL6 ~ estimate + aneuploidy, data=tnbc) %>% broom::tidy()
lm(IL6R ~ estimate + aneuploidy, data=tnbc) %>% broom::tidy()
lm(CGAS ~ aneuploidy, data=tnbc) %>% broom::tidy()
lm(CGAS ~ estimate + aneuploidy, data=tnbc) %>% broom::tidy()

aneup %>%
    filter(cohort == "BRCA", estimate >= 0.8,
           aneuploidy < quantile(aneuploidy, 0.25, na.rm=TRUE)) %>%
    lm(IL6 ~ estimate, data=.) %>% broom::tidy()
aneup %>%
    filter(cohort == "BRCA", estimate >= 0.8,
           aneuploidy > quantile(aneuploidy, 0.75, na.rm=TRUE)) %>%
    lm(IL6 ~ estimate, data=.) %>% broom::tidy()


pdf("il6r-aneup.pdf")
ggplot(aneup, aes(x=aneuploidy, y=IL6)) +
    geom_point() +
    geom_smooth(method="lm") +
    facet_wrap(~  cohort)
ggplot(aneup, aes(x=aneuploidy, y=IL6R)) +
    geom_point() +
    geom_smooth(method="lm") +
    facet_wrap(~  cohort)
ggplot(aneup, aes(x=aneuploidy, y=IL6ST)) +
    geom_point() +
    geom_smooth(method="lm") +
    facet_wrap(~  cohort)
ggplot(aneup, aes(x=aneuploidy, y=CGAS)) +
    geom_point() +
    geom_smooth(method="lm") +
    facet_wrap(~  cohort)
dev.off()
