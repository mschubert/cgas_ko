library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpmisc)
library(survminer)
library(survival)
sys = import('sys')

tt = theme(
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    legend.title = element_text(size=14),
    legend.text = element_text(size=12),
    plot.title = element_text(size=14)
)

plot_partition = function(brca) {
    psc = ggplot(brca, aes(x=`Interferon Gamma Response`, y=il6_cor)) +
        geom_vline(xintercept=0.25, color="grey", linetype="dashed", size=2) +
        geom_segment(y=0, yend=0, x=-5, xend=0.25, color="grey", linetype="dashed", size=2) +
        xlim(c(min(brca$`Interferon Gamma Response`, na.rm=TRUE), NA)) +
        annotate("text_npc", npcx=0.1, npcy=0.9, label="IL6-driven", color="black", size=6) +
        annotate("text_npc", npcx=0.3, npcy=0.05, label="Immune cold", color="black", size=6) +
        annotate("text_npc", npcx=0.95, npcy=0.05, label="Ifn-driven", color="black", size=6) +
        geom_point(aes(shape=factor(HER2), fill=factor(`ER/PR`), size=aneuploidy), alpha=0.3) +
        scale_shape_manual(values=c(Negative=21, Positive=24, Unknown=22)) +
        scale_fill_manual(values=c(Negative="red", Positive="blue", Unknown="#ffffff00")) +
        scale_size_continuous(range=c(0.2,5)) +
        guides(shape = guide_legend("HER2 status", override.aes = list(size=3, alpha=0.5)),
               fill = guide_legend("ER/PR status", override.aes = list(size=3, alpha=0.5, shape=21))) +
        theme_classic() + tt +
        labs(size = "Aneuploidy",
             x = "Interferon Gamma Response (sample)",
             y = "IL-6/JAK/STAT3 Signaling (cancer)")

    pdx = ggplot(brca %>% filter(`ER/PR` != "Unknown"), aes(x=`Interferon Gamma Response`, fill=factor(`ER/PR`))) +
        geom_density(color="#ffffff00", alpha=0.15) +
        scale_fill_manual(values=c(Negative="red", Positive="blue")) +
        theme_void() + guides(fill="none") +
        coord_cartesian(expand=c(0,0))
    pdy = ggplot(brca %>% filter(`ER/PR` != "Unknown"), aes(x=il6_cor, fill=factor(`ER/PR`))) +
        geom_density(color="#ffffff00", alpha=0.15) +
        scale_fill_manual(values=c(Negative="red", Positive="blue")) +
        theme_void() + guides(fill="none") +
        coord_flip(expand=c(0,0))

    pdx + plot_spacer() + psc + pdy + plot_layout(ncol=2, widths=c(10,1), heights=c(1,10), guides="collect")
}

plot_surv = function(x) {
    surv_one = function(dx) {
        pal = c("blue", "#ad07e3", "#ababab")
        lab = c("Ifn-driven", "IL6-driven", "Immune cold")

        m1 = coxph(Surv(OS_years, OS) ~ age_days + estimate + qq, data=dx) %>% broom::tidy(); m1
        m1p = coxph(Surv(OS_years, OS) ~ age_days + estimate + `E2F Targets` + qq, data=dx) %>% broom::tidy(); m1p
        p_il6 = . %>% filter(term=="qqil6") %>% pull(p.value)

        fit = survfit(Surv(OS_years, OS) ~ qq, data=dx); surv_pvalue(fit)
        dx %>% pull(qq) %>% table()

        ggsurvplot(fit, data=dx, pval=TRUE, xlim=c(0,5), break.time.by=2.5, palette=pal, legend.labs=lab)$plot +
            ylim(c(0.7,1)) + xlab("Overall survival (years)") +
            annotate("text_npc", npcx=0.05, npcy=0.07,
                     label=sprintf("p=%.2g\np=%.2g (Proliferation corr.)", p_il6(m1), p_il6(m1p)))
    }

    survs = list(
        surv_one(x %>% filter(aneuploidy > 0.2)) + ggtitle("Aneuploidy high"),
        surv_one(x %>% filter(aneuploidy < 0.2)) + ggtitle("Aneuploidy low"),
        surv_one(x %>% filter(CIN70_Carter2006 > 0)) + ggtitle("CIN70 high"),
        surv_one(x %>% filter(CIN70_Carter2006 < 0)) + ggtitle("CIN70 low")
    )
    (wrap_plots(survs) + plot_layout(guides="collect")) & theme(legend.direction = "vertical") & tt
}

sys$run({
    brca = readRDS("../data/tcga-pan.rds") %>%
        filter(cohort == "BRCA",
               substr(Sample, 14,16) == "01A") %>%
        mutate(patient = substr(Sample, 1, 12),
               OS = as.integer(vital_status) - 1,
               OS_years = os_days / 365,
               OS = ifelse(OS_years>5, 0, OS),
               OS_years = pmin(OS_years, 5)) %>%
        left_join(readRDS("../data/brca-meta.rds")) %>%
        select(Sample, estimate=purity, aneuploidy=aneup_log2seg, `ER/PR`, HER2,
               OS, OS_years, age_days, CIN70_Carter2006, `E2F Targets`,
               `Interferon Gamma Response`, `IL-6/JAK/STAT3 Signaling`)

    m1 = lm(`Interferon Gamma Response` ~ estimate, data=brca)
    brca$ifn_cor = brca$`Interferon Gamma Response` - predict(m1, newdata=data.frame(brca["estimate"]))
    m2 = lm(`IL-6/JAK/STAT3 Signaling` ~ estimate, data=brca)
    brca$il6_cor = brca$`IL-6/JAK/STAT3 Signaling` - predict(m2, newdata=data.frame(brca["estimate"]))

    brca = brca %>%
        mutate(qq = case_when(
                   `Interferon Gamma Response` > 0.25 ~ "both",
                   il6_cor > 0 & `Interferon Gamma Response` < 0.25 ~ "il6",
                   il6_cor < 0 & `Interferon Gamma Response` < 0.25 ~ "neither",
                   TRUE ~ NA_character_))

    dx = brca # wtf, ggsurvplot
    asm1 = plot_partition(brca)
    asm2 = plot_surv(brca)

    pdf("main_TCGA-BRCA.pdf", 18, 8)
    print((asm1 | wrap_elements(asm2)) + plot_layout(widths=c(1,1.2)))
    dev.off()
})
