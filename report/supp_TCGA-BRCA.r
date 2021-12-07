library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')

scatter_with_correction = function(df, x, ys, cor, cor_value=1) {
    df_one_y = function(y) {
        m = lm(as.formula(paste(y, "~", cor)), data=df)
        broom::tidy(m)
        pred_cor = predict(m, newdata=df[cor])
        pred_intcp = predict(m, newdata=tibble(!! rlang::sym(cor) := cor_value))
        df$cor = df[[y]] - (pred_cor - pred_intcp)

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

sys$run({
    brca = readRDS("../data/tcga-brca.rds")

    p11 = scatter_with_correction(brca, "CGAS", c("aneuploidy", "CIN70_Carter2006"), "estimate") +
        ggtitle("CIN/Aneuploidy with CGAS")
    p12 = scatter_with_correction(brca, "CGAS", c("IL6", "IL6R"), "estimate") +
        ggtitle("IL6/IL6R with CGAS")
    p21 = scatter_with_correction(brca, "Interferon Gamma Response", c("CGAS", "IL6", "IL6R"), "estimate") +
        ggtitle("Interferon Reponse")
    p22 = scatter_with_correction(brca, "IL-6/JAK/STAT3 Signaling", c("CGAS", "IL6", "IL6R"), "estimate") +
        ggtitle("IL6-STAT3 axis")

    asm = ((p11 | p12) / (p21 | p22)) +
        plot_annotation(tag_levels="a") +
        plot_layout(heights=c(2,3), guides="collect") &
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf("supp_TCGA-BRCA.pdf", 10, 11)
    print(asm)
    dev.off()
})
