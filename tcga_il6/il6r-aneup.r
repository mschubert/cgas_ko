library(dplyr)
library(ggplot2)
tcga = import('data/tcga')
gset = import('genesets')
idmap = import('process/idmap')

load_expr = function(x) {
    ensg = c("ENSG00000160712", "ENSG00000136244", "ENSG00000134352", "ENSG00000164430")
    re = tcga$rna_seq(x, trans="vst")[ensg,] %>%
        tcga$filter(cancer=TRUE, primary=TRUE) %>% t()
    tibble(Sample=rownames(re), IL6R=re[,1], IL6=re[,2], IL6ST=re[,3], CGAS=re[,4])
}

#todo: er/pr/her2 status? (@brca only script)

brca_gsva = function() {
    gex = tcga$rna_seq("BRCA", trans="vst")
    rownames(gex) = idmap$gene(rownames(gex), to="hgnc_symbol")
    sets = gset$get_human(c("MSigDB_Hallmark_2020", "CIN")) %>% unname() %>% do.call(c, .)
    t(GSVA::gsva(gex, sets))
}
scores = brca_gsva() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Sample") %>%
    as_tibble() %>%
    mutate(Sample = gsub(".", "-", Sample, fixed=TRUE))

load_thorsson = function() {
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
                  TNBC = pmin(ER, PR, HER2)
        )
}

cohorts = c("BRCA", "LUAD", "COAD", "OV")
expr = sapply(cohorts, load_expr, simplify=FALSE) %>% bind_rows(.id="cohort")
pur = tcga$purity()
immune = load_thorsson()
aneup = lapply(cohorts, tcga$aneuploidy) %>%
    bind_rows() %>%
    inner_join(expr) %>%
    inner_join(pur) %>%
    left_join(brca_meta()) %>%
    left_join(immune %>% select(Sample, `Immune Subtype`)) %>%
    left_join(scores) %>%
    mutate(aneuploidy = aneup_log2seg / estimate) # cancer aneup - stroma

brca = aneup %>% filter(cohort == "BRCA") %>%
    mutate(cgas_quart = cut(CGAS, breaks=quantile(CGAS, c(0,0.25,0.75,1)), labels=c("low", NA, "high")))
coad = aneup %>% filter(cohort == "COAD")
lm(IL6 ~ estimate + aneuploidy, data=brca) %>% broom::tidy()
lm(IL6 ~ estimate * aneuploidy, data=brca) %>% broom::tidy()
lm(IL6 ~ cgas_quart, data=brca) %>% broom::tidy()
lm(IL6 ~ CGAS, data=brca) %>% broom::tidy()
lm(IL6 ~ estimate + CGAS, data=brca) %>% broom::tidy()
lm(IL6R ~ estimate + aneuploidy, data=brca) %>% broom::tidy()
lm(CGAS ~ aneuploidy, data=brca) %>% broom::tidy()
lm(CGAS ~ estimate + aneuploidy, data=brca) %>% broom::tidy()
lm(CGAS ~ estimate + CIN70_Carter2006, data=brca) %>% broom::tidy()
lm(CGAS ~ estimate + Buccitelli_up, data=brca) %>% broom::tidy()
lm(`IL-6/JAK/STAT3 Signaling` ~ estimate + IL6R, data=brca) %>% broom::tidy() # sign: cgas, il6, il6r
lm(`Interferon Gamma Response` ~ estimate + CGAS, data=brca) %>% broom::tidy() # sign: only cgas
lm(IL6R ~ estimate + `Interferon Gamma Response`, data=brca) %>% broom::tidy()
# ifn pos + cgas top/bottom quartile?
ifn = brca %>% filter(`Interferon Gamma Response` > 0) %>%
    mutate(cgas_quart = cut(CGAS, breaks=quantile(CGAS, c(0,0.25,0.75,1)), labels=c("low", NA, "high")))
lm(IL6 ~ estimate + `Interferon Gamma Response`, data=ifn) %>% broom::tidy()
lm(IL6R ~ cgas_quart, data=ifn) %>% broom::tidy()
no_ifn = brca %>% filter(`Interferon Gamma Response` < 0) %>%
    mutate(cgas_quart = cut(CGAS, breaks=quantile(CGAS, c(0,0.25,0.75,1)), labels=c("low", NA, "high")))
lm(IL6 ~ estimate + `Interferon Gamma Response`, data=no_ifn) %>% broom::tidy()
lm(IL6R ~ cgas_quart, data=no_ifn) %>% broom::tidy()


# we have:
# * cgas up with aneup, CIN70, Buccitelli_up, IFNg resp, IL6 sig
# * IL6/R up with cgas, IL6 sig but not with IFNg resp
# -> does [aneup/CIN] surv [with inflamm env] depend on cgas/IL6 sig?

ggplot(brca, aes(x=CGAS, y=`Interferon Gamma Response`, color=factor(TNBC))) +
    geom_point() +
    geom_smooth(method="lm", se=FALSE)

tnbc = brca %>% filter(TNBC == 1)
lm(IL6 ~ estimate + aneuploidy, data=tnbc) %>% broom::tidy()
lm(IL6 ~ estimate * aneuploidy, data=tnbc) %>% broom::tidy()
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
