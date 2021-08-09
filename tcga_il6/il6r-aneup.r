library(dplyr)
library(ggplot2)
library(patchwork)
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
                  TNBC = pmin(ER, PR, HER2),
                  tumor_stage, age_at_diagnosis
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
    left_join(immune %>% select(Sample, `Immune Subtype`, `OS Time`, `OS`)) %>%
    left_join(scores) %>%
    mutate(aneuploidy = aneup_log2seg / estimate) # cancer aneup - stroma
brca = aneup %>% filter(cohort == "BRCA") %>%
    mutate(cgas_quart = cut(CGAS, breaks=quantile(CGAS, c(0,0.25,0.75,1)), labels=c("low", NA, "high")))
coad = aneup %>% filter(cohort == "COAD")

scatter_with_correction = function(df, x, ys, cor) {
    df_one_y = function(y) {
        m = lm(as.formula(paste(y, "~", cor)), data=df)
        broom::tidy(m)
        df$cor = df[[y]] - predict(m, newdata=df[cor])

        df %>%
            select(Sample, !!! rlang::syms(c(x, sample=y, cor, cancer="cor"))) %>%
            tidyr::gather("type", "expression", -estimate, -rlang::sym(x), -Sample) %>%
            mutate(type = factor(type, levels=c("sample", "cancer")))
    }

    comb = sapply(ys, df_one_y, simplify=FALSE) %>% dplyr::bind_rows(.id="gene") %>%
        left_join(brca %>% select(Sample, HER2, ER_PR) %>% distinct())
    comb$gene = factor(comb$gene, levels=ys)

    ggplot(comb, aes(x=!! rlang::sym(x), y=!! rlang::sym("expression"))) +
        geom_point(aes(shape=factor(HER2), fill=factor(ER_PR)), alpha=0.1) +
        geom_smooth(method="lm", se=FALSE) +
        ggpmisc::stat_fit_glance(method="lm", geom="text_npc", size=3.5, label.x=0.1,
                                 aes(label = paste("p", signif(..p.value.., digits = 2)))) +
        scale_shape_manual(values=c("0"=21, "1"=24), na.value=22) +
        scale_fill_manual(values=c("0"="red", "1"="blue"), na.value="black") +
        guides(shape = guide_legend("HER2 amp", override.aes = list(size=2, alpha=0.5)),
               fill = guide_legend("ER/PR status", override.aes = list(size=2, alpha=0.5, shape=21))) +
        theme_classic() +
        theme(strip.text = element_text(size=10)) +
        facet_grid(gene ~ type, scales="free", switch="y")
}

p11 = scatter_with_correction(brca, "CGAS", c("aneuploidy", "CIN70_Carter2006"), "estimate") +
    ggtitle("CIN/Aneuploidy with CGAS")
p12 = scatter_with_correction(brca, "CGAS", c("IL6", "IL6R"), "estimate") +
    ggtitle("IL6/IL6R with CGAS")
p21 = scatter_with_correction(brca, "Interferon Gamma Response", c("CGAS", "IL6", "IL6R"), "estimate") +
    ggtitle("Interferon Reponse")
p22 = scatter_with_correction(brca, "IL-6/JAK/STAT3 Signaling", c("CGAS", "IL6", "IL6R"), "estimate") +
    ggtitle("IL6-STAT3 axis")
pdf("tcga-cgas.pdf", 10, 11)
(p11 | p12) / (p21 | p22) + plot_layout(heights=c(2,3), guides="collect")
dev.off()

# todo: add residuals for IFN, IL6 that are not explained by purity
m1 = lm(`Interferon Gamma Response` ~ estimate, data=brca)
brca$ifn_cor = brca$`Interferon Gamma Response` - predict(m1, newdata=data.frame(brca["estimate"]))
m2 = lm(`IL-6/JAK/STAT3 Signaling` ~ estimate, data=brca)
brca$il6_cor = brca$`IL-6/JAK/STAT3 Signaling` - predict(m2, newdata=data.frame(brca["estimate"]))
brca$il6_bias = brca$il6_cor - brca$ifn_cor

library(ggpmisc)
psc = ggplot(brca, aes(x=`Interferon Gamma Response`, y=il6_cor)) +
    geom_vline(xintercept=0.25, color="grey", linetype="dashed", size=2) +
    geom_segment(y=0, yend=0, x=-5, xend=0.25, color="grey", linetype="dashed", size=2) +
    xlim(c(min(brca$`Interferon Gamma Response`, na.rm=TRUE), NA)) +
    annotate("text_npc", npcx=0.1, npcy=0.9, label="IL6-driven", color="grey", size=6) +
    annotate("text_npc", npcx=0.3, npcy=0.05, label="Immune cold", color="grey", size=6) +
    annotate("text_npc", npcx=0.95, npcy=0.05, label="Ifn-driven", color="grey", size=6) +
    geom_point(aes(shape=factor(HER2), fill=factor(ER_PR), size=CIN70_Carter2006), alpha=0.3) +
    scale_shape_manual(values=c("0"=21, "1"=24), na.value=22) +
    scale_fill_manual(values=c("0"="red", "1"="blue"), na.value="black") +
    scale_size_continuous(range=c(0.3,4)) +
    guides(shape = guide_legend("HER2 amp", override.aes = list(size=2, alpha=0.5)),
           fill = guide_legend("ER/PR status", override.aes = list(size=2, alpha=0.5, shape=21))) +
    theme_classic() +
    labs(x = "Interferon Gamma Response (sample)",
         y = "IL-6/JAK/STAT3 Signaling (cancer)")
pdx = ggplot(brca, aes(x=`Interferon Gamma Response`, fill=factor(ER_PR))) +
    geom_density(color="#ffffff00", alpha=0.15) +
    scale_fill_manual(values=c("0"="red", "1"="blue")) + theme_void() + guides(fill=FALSE) +
    coord_cartesian(expand=c(0,0))
pdy = ggplot(brca, aes(x=il6_cor, fill=factor(ER_PR))) +
    geom_density(color="#ffffff00", alpha=0.15) +
    scale_fill_manual(values=c("0"="red", "1"="blue")) + theme_void() + guides(fill=FALSE) +
    coord_flip(expand=c(0,0))
asm1 = pdx + plot_spacer() + psc + pdy + plot_layout(ncol=2, widths=c(10,1), heights=c(1,10), guides="collect")

library(survminer)
library(survival)
x = brca %>%
#    filter(CIN70_Carter2006 > 0) %>%
#    filter(aneuploidy > quantile(aneuploidy, 0.2, na.rm=T)) %>%
    mutate(OS_years = `OS Time` / 365,
           inflamm = `IL-6/JAK/STAT3 Signaling` + `Interferon Gamma Response`,
           il6_bias = `IL-6/JAK/STAT3 Signaling` - `Interferon Gamma Response`,
           cin = case_when(
                CIN70_Carter2006 > 0 ~ "high",
                CIN70_Carter2006 < 0 ~ "low",
                TRUE ~ NA_character_
            ),
           cgas = case_when(
                CGAS > quantile(CGAS, 0.75) ~ "high",
                CGAS < quantile(CGAS, 0.25) ~ "low",
                TRUE ~ NA_character_
            ),
           il6 = case_when(
#             il6_cor > quantile(il6_cor, 0.75, na.rm=T) ~ "il6",
#             il6_cor < quantile(il6_cor, 0.25, na.rm=T) ~ "ifn",
            `IL-6/JAK/STAT3 Signaling` > 0 & `Interferon Gamma Response` < 0.0 ~ "il6",
#                (`IL-6/JAK/STAT3 Signaling` - `Interferon Gamma Response`) > 0.2 ~ "il6",
            `IL-6/JAK/STAT3 Signaling` < 0 & `Interferon Gamma Response` > 0.0 ~ "ifn",
#                (`IL-6/JAK/STAT3 Signaling` - `Interferon Gamma Response`) < -0.2 ~ "ifn",
            TRUE ~ NA_character_
            ),
            qq = case_when(
                il6_cor > 0 & `Interferon Gamma Response`> 0.25 ~ "both",
                il6_cor > 0 & `Interferon Gamma Response` < 0.25 ~ "il6",
                il6_cor < 0 & `Interferon Gamma Response` > 0.25 ~ "both",
                il6_cor < 0 & `Interferon Gamma Response` < 0.25 ~ "neither",
                TRUE ~ NA_character_
            )) # last 2 -0.25 gives stronger sep, but interpretable?
m1 = coxph(Surv(OS_years, OS) ~ age_at_diagnosis + estimate + qq, data=x %>% filter(CIN70_Carter2006 > 0)) %>% broom::tidy(); m1
m2 = coxph(Surv(OS_years, OS) ~ age_at_diagnosis + estimate + qq, data=x %>% filter(CIN70_Carter2006 < 0)) %>% broom::tidy()

pal = c("blue", "#ad07e3", "#ababab")
lab = c("Ifn-driven", "IL6-driven", "Immune cold")
fit = survfit(Surv(OS_years, OS) ~ qq, data=x %>% filter(CIN70_Carter2006 > 0)); surv_pvalue(fit)
x %>% filter(CIN70_Carter2006 > 0) %>% pull(qq) %>% table()
ps1 = ggsurvplot(fit, data=x, pval=TRUE, xlim=c(0,5), break.time.by=2.5, palette=pal, legend.labs=lab)$plot +
    ylim(c(0.7,1)) + xlab("Overall survival (years)") + ggtitle("CIN") +
    annotate("text_npc", npcx=0.1, npcy=0.1,
             label=sprintf("p %.2g", m1 %>% filter(term=="qqil6") %>% pull(p.value)))
fit = survfit(Surv(OS_years, OS) ~ qq, data=x %>% filter(CIN70_Carter2006 < 0)); surv_pvalue(fit)
x %>% filter(CIN70_Carter2006 < 0) %>% pull(qq) %>% table()
ps2 = ggsurvplot(fit, data=x, pval=TRUE, xlim=c(0,5), break.time.by=2.5, palette=pal, legend.labs=lab)$plot +
    ylim(c(0.7,1)) + xlab("Overall survival (years)") + ggtitle("Non-CIN") +
    annotate("text_npc", npcx=0.1, npcy=0.1,
             label=sprintf("p %.2g", m2 %>% filter(term=="qqil6") %>% pull(p.value)))
asm2 = ps1 + ps2 + plot_layout(guides="collect") & theme(legend.direction = "vertical")
# p=0.1 for CIN70>0, il6

pdf("tcga-cgas2.pdf", 8, 10)
asm1 / asm2 + plot_layout(heights=c(2,0.85))
dev.off()

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
