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
    filter(!is.na(aneuploidy)) %>%
    tidyr::gather(field, value, -aneuploidy,  -Sample, -cohort) %>%
    group_by(field) %>%
        mutate(value = scale(value)) %>%
    ungroup()

plot_cohort = function(.cohort, aneup_thresh=0.2) {
    if (identical(.cohort, c("OV")))
        dset_surv$il6_cor = dset_surv$il6_bias

    dc = dset_surv %>% filter(cohort %in% .cohort) %>%
        mutate(il6_cor2 = cut(il6_cor, labels=c("low", "high"),
               c(-Inf, quantile(il6_cor, 0.25, na.rm=T), Inf)))

    pe = dset_expr %>%
        filter(cohort %in% .cohort) %>%
        ggplot(aes(x=field, y=value, fill=field)) +
            geom_hline(yintercept=0, linetype="dotted", size=2, alpha=0.5) +
            geom_boxplot(outlier.shape=NA, alpha=0.9, position=position_dodge2(preserve="single")) +
            scale_x_discrete(drop=FALSE) +
            theme_minimal() + theme(axis.text.x=element_blank()) +
            coord_cartesian(ylim=c(-2.5,2.5)) +
            labs(fill="Gene")

    p_il6 = function(mod) {
        re = broom::tidy(mod) %>% filter(term=="il6_cor") %>% pull(p.value)
        p = ifelse(re<0.05, sprintf("%.2g", re), "n.s.")
        sprintf("n=%i\np=%s", mod$n, p)
    }

    dca = dc %>% filter(aneuploidy > aneup_thresh)
    m1 = coxph(Surv(os_years, vital_status) ~ age_days + purity + il6_cor, data=dca)
    fit = survfit(Surv(os_years, vital_status) ~ il6_cor2, data=dca)
    surv1 = ggsurvplot(fit, data=dca, xlim=c(0,5), break.time.by=5, censor.size=3,
                       legend.title="IL6 (cancer)", legend.labs=c("low","high"))
    p1 = list(scale_color_manual(values=c(low="grey", high="#ad07e3"), drop=FALSE),
              scale_y_continuous(breaks=c(0,1)),
              annotate("text_npc", npcx=0.1, npcy=0.1, label=p_il6(m1)))

    dce = dc %>% filter(aneuploidy < aneup_thresh)
    m2 = coxph(Surv(os_years, vital_status) ~ age_days + purity + il6_cor, data=dce)
    fit = survfit(Surv(os_years, vital_status) ~ il6_cor2, data=dce)
    surv2 = ggsurvplot(fit, data=dce, xlim=c(0,5), break.time.by=5, censor.size=3,
                       legend.title="IL6 (cancer)", legend.labs=c("low","high"))
    p2 = list(scale_color_manual(values=c(low="grey", high="#ad07e3"), drop=FALSE),
              scale_y_continuous(breaks=c(0,1)),
              annotate("text_npc", npcx=0.1, npcy=0.1, label=p_il6(m2)))

    list(plt$text(paste(.cohort, collapse="/"), size=5, angle=90), pe, surv1$plot + p1, surv2$plot + p2)
}

cohorts = list(c("LUAD", "LUSC"), "OV", "SKCM")
plots = lapply(cohorts, plot_cohort)
wp = c(list(plot_spacer(), plt$text("Expression\nz-score") + coord_cartesian(clip="off"),
            plt$text("Aneuploidy high"), plt$text("Aneuploidy low")),
       do.call(c, plots))
asm = wrap_plots(wp, ncol=4, widths=c(1,3,5,5), heights=c(1.5,rep(5,length(cohorts)))) + plot_layout(guides="collect") &
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.direction = "vertical",
          legend.title = element_text(size=12),
          legend.text = element_text(size=10))

pdf("main_TCGA-pan.pdf", 6.5, 6.5)
print(asm)
dev.off()
