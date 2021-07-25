library(dplyr)
library(ggplot2)
library(DESeq2)
sys = import('sys')
plt = import('plot')
gset = import('genesets')

args = sys$cmd$parse(
    opt('i', 'infile', 'rds', '../comp_list/summary.rds'),
    opt('p', 'plotfile', 'pdf', '1-rev_cgas_wt.pdf')
)

lookup = c(
    "BT549 cGAS KO reversine vs BT549 WT reversine" = "diff",
    "BT549 cGAS KO Reversine vs DMSO" = "y",
    "BT549 WT Reversine vs DMSO" = "x"
)

res = readRDS(args$infile) %>%
    filter(comparison %in% names(lookup)) %>%
    mutate(comparison = lookup[comparison])

## purity adjustment (better to do it with genes?)
#res$genes[[3]] = res$genes[[3]] %>%
#    left_join(res$genes[[2]] %>% select(ensembl_gene_id, wt_lfc=log2FoldChange)) %>%
#    mutate(log2FoldChange = log2FoldChange - (1-0.83) * wt_lfc)

#de_genes = list(cgas_dep = res$cgas_dep$ensembl_gene_id[res$cgas_dep$padj < 0.1],
#                cgas_indep = res$cgas_indep$ensembl_gene_id[res$cgas_indep$padj < 0.1])
#plt$venn(de_genes)

cmp_sets = res %>%
    select(cond=comparison, MSigDB_Hallmark_2020, DoRothEA, CIN) %>%
    tidyr::gather("collection", "stats", -cond) %>%
    tidyr::unnest() %>%
#    filter(collection != "DoRothEA" | size >= 50) %>%
    filter(size >= 50) %>%
    mutate(label = sub(" (a)", "", label, fixed=TRUE))
cmp_sets = cmp_sets %>%
    left_join(cmp_sets %>% filter(cond == "x") %>% select(label, wt_lfc=log2FoldChange)) %>%
    mutate(log2FoldChange = ifelse(cond == "y", log2FoldChange - (1-0.83) * wt_lfc, log2FoldChange)) %>%
    select(-wt_lfc)

show_lab = c(
    "Interferon Gamma Response", "Interferon Alpha Response", "TNF-alpha Signaling via NF-kB",
    "E2F Targets", "G2-M Checkpoint", "Myc Targets V1",
    "STAT1", "RELA", "E2F4", "NFKB1", "STAT3", "JUN", "FOS", "CREB1"
)
sres2 = cmp_sets %>%
    filter(cond != "diff") %>%
    select(cond, collection, label, size, log2FoldChange) %>%
    tidyr::pivot_wider(names_from=c(cond), values_from=c(log2FoldChange)) %>%
    left_join(cmp_sets %>% filter(cond == "diff") %>% select(label, p.value, adj.p)) %>%
    mutate(label = ifelse(adj.p < 1e-10 | label %in% show_lab, label, NA),
           min_p = cut(adj.p, breaks=c(0, 1e-10, 0.05, Inf), labels=c("<1e-10", "<0.01", "n.s.")))

arr = c("Interferon Gamma Response", "Interferon Alpha Response", "STAT1",
        "TNF-alpha Signaling via NF-kB", "Myc Targets V1")
arts = c("cGas-independent", "cGas-dependent")
arrws = sres2 %>%
    filter(label %in% arr) %>%
    transmute(y, from1=0, to1=y, from2=y, to2=x) %>%
    tidyr::pivot_longer(!y, names_to = c(".value", "set"), names_pattern = "([^0-9]+)(.)") %>%
    mutate(set=arts[as.integer(set)])

cairo_pdf(args$plotfile, 9, 7)
ggplot(sres2, aes(x=x, y=y)) +
    geom_hline(yintercept=0, color="darkgrey", linetype="dotted") +
    geom_vline(xintercept=0, color="darkgrey", linetype="dotted") +
    geom_abline(slope=1, size=2, color="grey", linetype="dashed") +
    geom_point(aes(size=size, fill=collection, alpha=min_p), shape=21) +
    scale_alpha_manual(values=c("<1e-10"=0.9, "<0.01"=0.4, "n.s."=0.1)) +
    scale_size_area() +
    geom_segment(data=arrws, aes(x=from, xend=to, y=y, yend=y, color=set),
                 arrow = arrow(length = unit(0.01, "npc"), type="closed")) +
    scale_color_brewer(palette="Dark2") +
    ggrepel::geom_label_repel(aes(label=label), size=3, max.iter=1e5, label.size=NA,
        min.segment.length=0, max.overlaps=Inf, segment.alpha=0.3, fill="#ffffffc0",
        label.padding=unit(0.2, "lines")) +
#    geom_text(data=data.frame(x=0.5, y=0.5, txt="more cGas-independent ⇖    ⇘ more cGas-dependent  "),
#              aes(x=x, y=y, label=txt), inherit.aes=FALSE, size=3.5, fontface="bold", color="#757575") +
    theme_classic() +
    labs(x = "log2 FC wt: rev vs. DMSO",
         y = "log2 FC cGas KO: rev vs. DMSO",
         size = "Set genes",
         color = "FC origin",
         fill = "Set type",
         alpha = "FDR Δ")
dev.off()
