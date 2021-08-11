library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
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

load_thorsson = function(cohorts="BRCA") {
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

load_brca = function() {
    scores = brca_gsva() %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Sample") %>%
        as_tibble() %>%
        mutate(Sample = gsub(".", "-", Sample, fixed=TRUE))

    brca = tcga$aneuploidy("BRCA") %>%
        inner_join(load_expr("BRCA")) %>%
        inner_join(tcga$purity()) %>%
        left_join(brca_meta()) %>%
        left_join(load_thorsson() %>% select(Sample, `Immune Subtype`, `OS Time`, `OS`)) %>%
        left_join(scores) %>%
        mutate(aneuploidy = aneup_log2seg / estimate) # cancer aneup - stroma
}

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
    brca = load_brca()

    p11 = scatter_with_correction(brca, "CGAS", c("aneuploidy", "CIN70_Carter2006"), "estimate") +
        ggtitle("CIN/Aneuploidy with CGAS")
    p12 = scatter_with_correction(brca, "CGAS", c("IL6", "IL6R"), "estimate") +
        ggtitle("IL6/IL6R with CGAS")
    p21 = scatter_with_correction(brca, "Interferon Gamma Response", c("CGAS", "IL6", "IL6R"), "estimate") +
        ggtitle("Interferon Reponse")
    p22 = scatter_with_correction(brca, "IL-6/JAK/STAT3 Signaling", c("CGAS", "IL6", "IL6R"), "estimate") +
        ggtitle("IL6-STAT3 axis")

    pdf("supp_TCGA-BRCA.pdf", 10, 11)
    print((p11 | p12) / (p21 | p22) + plot_layout(heights=c(2,3), guides="collect"))
    dev.off()

    saveRDS(brca, file="supp_TCGA-BRCA.rds")
})
