library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpmisc)
library(survminer)
library(survival)

brca = readRDS("../data/tcga-brca.rds") %>%
    mutate(HER2 = case_when(
        HER2 == 0 ~ "normal",
        HER2 == 1 ~ "amplified",
        TRUE ~ "unknown"
    )) %>%
    mutate(ER_PR = case_when(
        ER_PR == 0 ~ "negative",
        ER_PR == 1 ~ "positive",
        TRUE ~ "unknown"
    ))

tt = theme(
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    legend.title = element_text(size=14),
    legend.text = element_text(size=12),
    plot.title = element_text(size=14)
)

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
    annotate("text_npc", npcx=0.1, npcy=0.9, label="IL6-driven", color="black", size=6) +
    annotate("text_npc", npcx=0.3, npcy=0.05, label="Immune cold", color="black", size=6) +
    annotate("text_npc", npcx=0.95, npcy=0.05, label="Ifn-driven", color="black", size=6) +
    geom_point(aes(shape=factor(HER2), fill=factor(ER_PR), size=CIN70_Carter2006), alpha=0.3) +
    scale_shape_manual(values=c(normal=21, amplified=24, unknown=22)) +
    scale_fill_manual(values=c(negative="red", positive="blue", unknown="#ffffff00")) +
    scale_size_continuous(range=c(0.3,4)) +
    guides(shape = guide_legend("HER2 status", override.aes = list(size=3, alpha=0.5)),
           fill = guide_legend("ER/PR status", override.aes = list(size=3, alpha=0.5, shape=21))) +
    theme_classic() + tt +
    labs(size = "CIN70",
         x = "Interferon Gamma Response (sample)",
         y = "IL-6/JAK/STAT3 Signaling (cancer)")
pdx = ggplot(brca %>% filter(ER_PR != "unknown"), aes(x=`Interferon Gamma Response`, fill=factor(ER_PR))) +
    geom_density(color="#ffffff00", alpha=0.15) +
    scale_fill_manual(values=c(negative="red", positive="blue")) +
    theme_void() + guides(fill="none") +
    coord_cartesian(expand=c(0,0))
pdy = ggplot(brca %>% filter(ER_PR != "unknown"), aes(x=il6_cor, fill=factor(ER_PR))) +
    geom_density(color="#ffffff00", alpha=0.15) +
    scale_fill_manual(values=c(negative="red", positive="blue")) +
    theme_void() + guides(fill="none") +
    coord_flip(expand=c(0,0))
asm1 = pdx + plot_spacer() + psc + pdy + plot_layout(ncol=2, widths=c(10,1), heights=c(1,10), guides="collect")

x = brca %>%
    mutate(OS_years = `OS Time` / 365,
           qq = case_when(
               `Interferon Gamma Response` > 0.25 ~ "both",
               il6_cor > 0 & `Interferon Gamma Response` < 0.25 ~ "il6",
               il6_cor < 0 & `Interferon Gamma Response` < 0.25 ~ "neither",
               TRUE ~ NA_character_
           )) %>%
    mutate(OS = ifelse(OS_years>5, 0, OS),
           OS_years = pmin(OS_years, 5))

m1 = coxph(Surv(OS_years, OS) ~ age_at_diagnosis + estimate + qq,
           data=x %>% filter(CIN70_Carter2006 > 0)) %>% broom::tidy(); m1
m2 = coxph(Surv(OS_years, OS) ~ age_at_diagnosis + estimate + qq,
           data=x %>% filter(CIN70_Carter2006 < 0)) %>% broom::tidy(); m2
m1p = coxph(Surv(OS_years, OS) ~ age_at_diagnosis + estimate + `E2F Targets` + qq,
            data=x %>% filter(CIN70_Carter2006 > 0)) %>% broom::tidy(); m1p
m2p = coxph(Surv(OS_years, OS) ~ age_at_diagnosis + estimate + `E2F Targets` + qq,
            data=x %>% filter(CIN70_Carter2006 < 0)) %>% broom::tidy(); m2p
p_il6 = . %>% filter(term=="qqil6") %>% pull(p.value)

pal = c("blue", "#ad07e3", "#ababab")
lab = c("Ifn-driven", "IL6-driven", "Immune cold")
fit = survfit(Surv(OS_years, OS) ~ qq, data=x %>% filter(CIN70_Carter2006 > 0)); surv_pvalue(fit)
x %>% filter(CIN70_Carter2006 > 0) %>% pull(qq) %>% table()
ps1 = ggsurvplot(fit, data=x, pval=TRUE, xlim=c(0,5), break.time.by=2.5, palette=pal, legend.labs=lab)$plot +
    ylim(c(0.7,1)) + xlab("Overall survival (years)") + ggtitle("CIN70 high") +
    annotate("text_npc", npcx=0.05, npcy=0.07,
             label=sprintf("p=%.2g\np=%.2g (Proliferation corr.)", p_il6(m1), p_il6(m1p)))
fit = survfit(Surv(OS_years, OS) ~ qq, data=x %>% filter(CIN70_Carter2006 < 0)); surv_pvalue(fit)
x %>% filter(CIN70_Carter2006 < 0) %>% pull(qq) %>% table()
ps2 = ggsurvplot(fit, data=x, pval=TRUE, xlim=c(0,5), break.time.by=2.5, palette=pal, legend.labs=lab)$plot +
    ylim(c(0.7,1)) + xlab("Overall survival (years)") + ggtitle("CIN70 low") +
    annotate("text_npc", npcx=0.05, npcy=0.07,
             label=sprintf("p=%.2g\np=%.2g (Proliferation corr.)", p_il6(m2), p_il6(m2p)))
asm2 = (ps1 + ps2 + plot_layout(guides="collect")) & theme(legend.direction = "vertical") & tt

pdf("main_TCGA-BRCA.pdf", 8.35, 10)
asm1 / asm2 + plot_layout(heights=c(2,0.85))
dev.off()
