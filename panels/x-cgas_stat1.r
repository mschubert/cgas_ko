library(dplyr)
library(ggplot2)
library(DESeq2)
sys = import('sys')
plt = import('plot')
gset = import('genesets')
idmap = import('process/idmap')

args = sys$cmd$parse(
    opt('i', 'infile', 'rds', '../data/rnaseq.rds'),
    opt('p', 'plotfile', 'pdf', 'x-rev_cgas_stat1.pdf')
)

eset = readRDS(args$infile)
eset = eset[, eset$genotype %in% c("cgas", "stat1") & eset$batch != "4" &
              eset$treatment %in% c("dmso", "rev") & eset$time == "48"]

eset$genotype = droplevels(eset$genotype)
eset$treatment = droplevels(eset$treatment)
eset$cgas_ko = (eset$genotype == "stat1" & eset$treatment == "rev") + 0
eset$stat1_ko = (eset$genotype == "cgas" & eset$treatment == "rev") + 0
design(eset) = ~ genotype + cgas_ko + stat1_ko #genotype:treatment
model.matrix(design(eset), colData(eset))

mod = DESeq2::DESeq(eset)
resultsNames(mod)
res = list(cgas_ko = DESeq2::results(mod, name="cgas_ko"),
           stat1_ko = DESeq2::results(mod, name="stat1_ko")) %>%
    lapply(function(x) as.data.frame(x) %>% tibble::rownames_to_column("ensembl_gene_id") %>%
           mutate(gene_name = idmap$gene(ensembl_gene_id, to="hgnc_symbol")) %>%
           as_tibble() %>% select(ensembl_gene_id, gene_name, everything()) %>%
           arrange(padj, pvalue) %>% filter(!is.na(padj)))

de_genes = list(cgas_ko = res$cgas_ko$gene_name[res$cgas_ko$padj < 0.1],
                stat1_ko = res$stat1_ko$gene_name[res$stat1_ko$padj < 0.1])
#plt$venn(de_genes)

snames = c("MSigDB_Hallmark_2020", "DoRothEA")
sets = gset$get_human(snames, conf="A") %>%
    lapply(gset$filter, min=50)

#volc w/ wt rev eff?

sres = tidyr::crossing(tibble(cond=c("cgas_ko", "stat1_ko")), tibble(set=snames)) %>%
    rowwise() %>%
    mutate(res = list(gset$test_lm(res[[cond]], sets[[set]], stat="log2FoldChange"))) %>%
    tidyr::unnest() %>%
    mutate(label = sub(" (a)", "", label, fixed=TRUE))

show_lab = c(
    "Interferon Gamma Response", "Interferon Alpha Response", "TNF-alpha Signaling via NF-kB",
    "E2F Targets", "G2-M Checkpoint", "Apoptosis", "Myc Targets V1", "UV Response Up",
    "Inflammatory Response", "KRAS Signaling Up", "IL-6/JAK/STAT3 Signaling",
    "Epithelial Mesenchymal Transition", "Allograft Rejection", "Complement", "Estrogen Response Late",
    "STAT1", "RELA", "E2F4", "NFKB1", "STAT3", "JUN", "FOS", "CREB1", "SPI1", "SP3"
)
sres2 = sres %>%
    select(cond, set, label, size, estimate, adj.p) %>%
    tidyr::pivot_wider(names_from=c(cond), values_from=c(estimate, adj.p)) %>%
    mutate(label = ifelse(label %in% show_lab, label, NA),
           min_p = apply(cbind(adj.p_cgas_ko, adj.p_stat1_ko), 1, function(x) min(x, na.rm=TRUE)),
           min_p = cut(min_p, breaks=c(0, 1e-10, 0.01, Inf), labels=c("<1e-10", "<0.01", "n.s.")))

cairo_pdf(args$plotfile, 9, 7)
ggplot(sres2, aes(x=estimate_cgas_ko, y=estimate_stat1_ko)) +
    geom_hline(yintercept=0, color="darkgrey", linetype="dotted") +
    geom_vline(xintercept=0, color="darkgrey", linetype="dotted") +
    geom_abline(slope=1, size=2, color="grey", linetype="dashed", inherit.aes=F) +
    geom_point(aes(size=size, color=set, alpha=min_p)) +
    scale_alpha_manual(values=c("<1e-10"=0.9, "<0.01"=0.4, "n.s."=0.1)) +
    scale_size_area() +
    ggrepel::geom_label_repel(aes(label=label), size=3, max.iter=1e5, label.size=NA,
        min.segment.length=0, max.overlaps=Inf, segment.alpha=0.3, fill="#ffffffc0",
        label.padding=unit(0.2, "lines")) +
#    geom_text(data=data.frame(x=0.5, y=0.5, txt="more cGas-independent ⇖    ⇘ more cGas-dependent  "),
#              aes(x=x, y=y, label=txt), inherit.aes=FALSE, size=3.5, fontface="bold", color="#757575") +
    theme_classic() +
    labs(#x = "log2 FC wt: rev vs. DMSO",
         #y = "log2 FC cGas KO: rev vs. DMSO",
         size = "Set genes",
         color = "Set type",
         alpha = "FDR")
dev.off()
