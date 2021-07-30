library(dplyr)
library(ggplot2)
tcga = import('data/tcga')

load_expr = function(x) {
    re = t(tcga$rna_seq(x, trans="vst")[c("ENSG00000160712", "ENSG00000136244"),])
    tibble(Sample=rownames(re), IL6R=re[,1], IL6=re[,2])
}

cohorts = c("BRCA", "LUAD", "COAD", "OV")
expr = sapply(cohorts, load_expr, simplify=FALSE) %>% bind_rows(.id="cohort")
aneup = lapply(cohorts, tcga$aneuploidy) %>%
    bind_rows() %>%
    inner_join(expr)

ggplot(aneup, aes(x=aneuploidy, y=IL6)) +
    geom_point() +
    geom_smooth(method="lm") +
    facet_wrap(~  cohort)
