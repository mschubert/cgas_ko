library(dplyr)
library(ggplot2)
tcga = import('data/tcga')

load_expr = function(x) {
    ensg = c("ENSG00000160712", "ENSG00000136244", "ENSG00000134352", "ENSG00000164430")
    re = tcga$rna_seq(x, trans="vst")[ensg,] %>%
        tcga$filter(cancer=TRUE, primary=TRUE) %>% t()
    tibble(Sample=rownames(re), IL6R=re[,1], IL6=re[,2], IL6ST=re[,3], CGAS=re[,4])
}

cohorts = c("BRCA", "LUAD", "COAD", "OV")
expr = sapply(cohorts, load_expr, simplify=FALSE) %>% bind_rows(.id="cohort")
meta = readxl::read_xlsx(".....")
pur = tcga$purity()
aneup = lapply(cohorts, tcga$aneuploidy) %>%
    bind_rows() %>%
    inner_join(expr) %>%
    left_join(pur)
#TODO: add different measures of aneup, cor with IL6/R/ST,CGAS

lm(IL6 ~ estimate * aneuploidy, data=aneup %>% filter(cohort == "BRCA")) %>% broom::tidy()
lm(IL6R ~ estimate * aneuploidy, data=aneup %>% filter(cohort == "BRCA")) %>% broom::tidy()
lm(IL6ST ~ estimate * aneuploidy, data=aneup %>% filter(cohort == "BRCA")) %>% broom::tidy()
lm(CGAS ~ estimate + aneuploidy, data=aneup %>% filter(cohort == "BRCA")) %>% broom::tidy()

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
