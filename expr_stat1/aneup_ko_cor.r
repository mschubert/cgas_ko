library(dplyr)
library(ggplot2)
io = import('io')
sys = import('sys')
idmap = import('process/idmap')
gset = import('data/genesets')

colors = setNames(c("#a50f15", "#006d2c", "#045a8d", "#cccccc", "#8a8a8a"),
                  c("compensated", "hyper-dereg", "inconsistent", "no change", "only 1 dset"))

merge_one = function(subs, wald=1.5) {
    both = list(ko_diff=stat1[[subs]], aneup=aneup) %>%
        bind_rows(.id="dset") %>%
        transmute(dset = dset,
                  gene_name = gene_name,
                  padj = padj,
                  stat = log2FoldChange / lfcSE) %>%
        tidyr::pivot_wider(names_from="dset", values_from=c("stat", "padj")) %>%
        mutate(type = case_when(
                is.na(stat_ko_diff) | is.na(stat_aneup) ~ "only 1 dset",
                abs(stat_ko_diff) < wald | abs(stat_aneup) < wald ~ "no change",
                stat_ko_diff < 0 & stat_aneup < 0 ~ "compensated",
                stat_ko_diff > 0 & stat_aneup > 0 ~ "hyper-dereg",
                TRUE ~ "inconsistent"
            ),
            type = factor(type, levels=names(colors))) %>%
        na.omit()
}

merge_cgas = function(subs, wald=1.5) {
    both = list(ko_diff=stat1[[subs]], aneup=stat1$rev24_cgas_over_wt) %>%
        bind_rows(.id="dset") %>%
        transmute(dset = dset,
                  gene_name = gene_name,
                  padj = padj,
                  stat = log2FoldChange / lfcSE) %>%
        tidyr::pivot_wider(names_from="dset", values_from=c("stat", "padj")) %>%
        mutate(type = case_when(
                is.na(stat_ko_diff) | is.na(stat_aneup) ~ "only 1 dset",
                abs(stat_ko_diff) < wald | abs(stat_aneup) < wald ~ "no change",
                stat_ko_diff < 0 & stat_aneup < 0 ~ "compensated",
                stat_ko_diff > 0 & stat_aneup > 0 ~ "hyper-dereg",
                TRUE ~ "inconsistent"
            ),
            type = factor(type, levels=names(colors))) %>%
        na.omit()
}

plot_one = function(merged) {
    merged = merged %>%
        mutate(label = ifelse(padj_ko_diff < 0.1 & padj_aneup < 0.1, gene_name, NA)) %>%
        filter(padj_ko_diff < 0.1 | padj_aneup < 0.1)
    ggplot(merged, aes(x=stat_ko_diff, y=stat_aneup)) +
        geom_point(color="grey") +
        scale_color_manual(values=colors) +
        ggrepel::geom_label_repel(aes(label=label, color=type), size=3,
                                  na.rm=TRUE, segment.alpha=0.3, fill="#ffffffc0",
                                  label.padding=0.1, max.iter=1e4, min.segment.length=0) +
        geom_hline(yintercept=0, linetype="dashed", size=2, alpha=0.3) +
        geom_vline(xintercept=0, linetype="dashed", size=2, alpha=0.3)
}

args = sys$cmd$parse(
    opt('d', 'diff_expr', 'rds', 'diff_expr.rds'),
    opt('a', 'diff_aneup', 'RData', '../expr_diff/de_Mad2PB.RData'),
    opt('p', 'plotfile', 'pdf', 'aneup_ko_cor.pdf'))

stat1 = readRDS(args$diff_expr)
aneup = io$load(args$diff_aneup)$aneuploidy
aneup$gene_name = toupper(aneup$gene_name)
#aneup$gene_name = idmap$orthologue(aneup$gene_name, dset="mmusculus_gene_ensembl",
#                                   from="mgi_symbol", to="hgnc_symbol")

pdf(args$plotfile, 16, 14)
plot_one(merge_one("rev24_cgas_over_wt")) + ggtitle("rev24_cgas_over_wt")
plot_one(merge_one("rev24_stat1_over_wt")) + ggtitle("rev24_stat1_over_wt")
plot_one(merge_one("rev48_stat1_over_wt")) + ggtitle("rev48_stat1_over_wt")
plot_one(merge_cgas("rev24_stat1_over_wt")) + ggtitle("rev24 cgas vs stat1 KO")
plot_one(merge_cgas("rev48_stat1_over_wt")) + ggtitle("rev24 cgas vs 48 stat1 KO")
dev.off()
