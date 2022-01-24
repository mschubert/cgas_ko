library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpmisc)
library(survminer)
library(survival)
plt = import('plot')

tt = theme(
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    legend.title = element_text(size=14),
    legend.text = element_text(size=12),
    plot.title = element_text(size=14)
)

pan = readRDS("../data/tcga.rds")

cors = pan %>%
    group_by(cohort) %>%
        summarize(mod_ifn = list(lm(`Interferon Gamma Response` ~ purity)),
                  mod_il6 = list(lm(`IL-6/JAK/STAT3 Signaling` ~ purity))) %>%
    ungroup()

dset_surv = left_join(pan, cors) %>%
    rowwise() %>%
        mutate(ifn_cor = `Interferon Gamma Response` - predict(mod_ifn, newdata=data.frame(purity=purity)),
               il6_cor = `IL-6/JAK/STAT3 Signaling` - predict(mod_il6, newdata=data.frame(purity=purity))) %>%
    ungroup() %>%
    select(-mod_ifn, -mod_il6) %>%
    mutate(il6_bias = il6_cor - `Interferon Gamma Response`)

dset_expr = dset_surv %>%
    select(Sample, cohort, aneuploidy, CGAS, IL6, IL6R) %>%
    filter(!is.na(aneuploidy), cohort != "BRCA") %>%
    mutate(aneuploidy = cut(aneuploidy, c(-Inf,0.1,Inf), labels=c("low", "high"))) %>%
    tidyr::gather(field, value, -aneuploidy,  -Sample, -cohort) %>%
    group_by(cohort) %>%
        mutate(n_cohort = n()) %>%
    group_by(cohort, aneuploidy) %>%
        filter(n() > 0.05*n_cohort) %>%
    group_by(field) %>%
        mutate(value = scale(value)) %>%
    ungroup()

plot_cohort = function(.cohort) {
    if (.cohort %in% c("OV", "SKCM", "COAD"))
        dset_surv$il6_cor = dset_surv$il6_bias

    dc = dset_surv %>% filter(cohort == .cohort) %>%
        mutate(il6_cor2 = cut(il6_cor, labels=c("low", "middle", "high"),
               c(-Inf, quantile(il6_cor, 0.25, na.rm=T), quantile(il6_cor, 0.75, na.rm=T), Inf)))

    top = dset_expr %>%
        filter(cohort == .cohort) %>%
        ggplot(aes(x=field, y=value, color=aneuploidy, fill=field)) +
            geom_hline(yintercept=0, linetype="dotted", size=2, alpha=0.5) +
            geom_boxplot(outlier.shape=NA, alpha=0.9, position=position_dodge2(preserve="single")) +
            scale_x_discrete(drop=FALSE) +
            theme_minimal() +
            coord_cartesian(ylim=c(-2.5,2.5))

    p_il6 = . %>% filter(term=="il6_cor") %>% pull(p.value)

    dca = dc %>% filter(aneuploidy > 0.1)
    m1 = coxph(Surv(os_years, vital_status) ~ age_days + purity + il6_cor, data=dca) %>% broom::tidy()
    fit = survfit(Surv(os_years, vital_status) ~ il6_cor2, data=dca)
    surv1 = ggsurvplot(fit, data=dca, xlim=c(0,5), break.time.by=2.5) +
        ylab("OS aneuploid")
    p1 = annotate("text_npc", npcx=0.1, npcy=0.1, label=sprintf("p=%.2g", p_il6(m1)))

    dce = dc %>% filter(aneuploidy < 0.1)
    m2 = coxph(Surv(os_years, vital_status) ~ age_days + purity + il6_cor, data=dce) %>% broom::tidy()
    fit = survfit(Surv(os_years, vital_status) ~ il6_cor2, data=dce)
    surv2 = ggsurvplot(fit, data=dce, xlim=c(0,5), break.time.by=2.5) +
        ylab("OS euploid")
    p2 = annotate("text_npc", npcx=0.1, npcy=0.1, label=sprintf("p=%.2g", p_il6(m2)))

    (plt$text(.cohort, size=5) / top / (surv1$plot + p1) / (surv2$plot + p2)) +
        plot_layout(heights=c(0.08,1,1,1))
}

plots = lapply(c("LUAD", "LUSC", "OV", "SKCM", "COAD"), plot_cohort)

pdf("main_TCGA-pan.pdf", 16, 9)
wrap_plots(plots, nrow=1) + plot_layout(guides="collect") & theme(legend.direction = "vertical")
dev.off()
