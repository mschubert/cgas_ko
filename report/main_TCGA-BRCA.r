library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpmisc)
library(survminer)
library(survival)

brca = readRDS("../data/tcga-brca.rds")

# todo: add residuals for IFN, IL6 that are not explained by purity
m1 = lm(`Interferon Gamma Response` ~ estimate, data=brca)
brca$ifn_cor = brca$`Interferon Gamma Response` - predict(m1, newdata=data.frame(brca["estimate"]))
m2 = lm(`IL-6/JAK/STAT3 Signaling` ~ estimate, data=brca)
brca$il6_cor = brca$`IL-6/JAK/STAT3 Signaling` - predict(m2, newdata=data.frame(brca["estimate"]))
brca$il6_bias = brca$il6_cor - brca$ifn_cor

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

x = brca %>%
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
            `IL-6/JAK/STAT3 Signaling` > 0 & `Interferon Gamma Response` < 0.0 ~ "il6",
            `IL-6/JAK/STAT3 Signaling` < 0 & `Interferon Gamma Response` > 0.0 ~ "ifn",
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

pdf("main_TCGA-BRCA.pdf", 8, 10)
asm1 / asm2 + plot_layout(heights=c(2,0.85))
dev.off()
