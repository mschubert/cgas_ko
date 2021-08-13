library(dplyr)
library(ggplot2)
library(DESeq2)
library(pheatmap)
sys = import('sys')
gset = import('genesets')
idmap = import('process/idmap')

args = sys$cmd$parse(
    opt('e', 'eset', 'rds', '../data/rnaseq.rds'),
    opt('p', 'plotfile', 'pdf', 'cytokine.pdf')
)

eset = readRDS(args$eset)
#eset = eset[,eset$genotype %in% c("wt", "cgas", "stat1", "stat3", "rela", "relb") &
eset = eset[,eset$genotype %in% c("wt", "cgas", "stat1", "stat3", "rela", "relb", "cgas+stat1") &
            eset$treatment %in% c("dmso", "rev") & eset$time == "48"]

wt_rev = readRDS("../comp_list/summary.rds") %>%
    filter(comparison == "BT549 WT Reversine vs DMSO") %>%
    pull(genes) %>% `[[`(1) %>%
    filter(padj < 0.1)

#go = gset$get_human("GO_Biological_Process_2021")
dor = gset$get_human("DoRothEA")[c("STAT1 (a)", "STAT3 (a)", "RELA (a)", "NFKB1 (a)")]
hm = gset$get_human("MSigDB_Hallmark_2020")[c("Interferon Alpha Response",
        "Interferon Gamma Response", "TNF-alpha Signaling via NF-kB")]

# missing?: CASP8, CDKN2A, FASL
hlg = c("IL6", "TRAF2", "TNFRSF10A", "FOS", "JUN", "IRF1", "OAS3", "FAS", "NFKB1", "ISG15", "ISG20", "STAT3",
        "TP53", "VEGFA", "VEGFB", "VEGFC", "KIF2A", "CXCL10", "IRF7", "MMP2", "SLC25A27", "CASP4", "CASP8",
        "CCL2", "CCL20", "TNFSF10", "CREM", "BCL6", "IRF2", "HDAC1", "MYC", "CCND1", "BCL2L1", "IDO1", "ISG20",
        "OAS1", "CXCL3", "CASP8", "CDKN1A", "TLR3", "CD274", "CXCL2", "EGFR", "FGF2", "HIF1A", "CDK6", "CSF1",
        "IL6ST", "CCND3", "HDAC7", "NFKBIA", "NFKBIA", "IFNAR1", "TRIM25", "TWIST1")


design(eset) = ~1
vs = assay(DESeq2::varianceStabilizingTransformation(eset))
rownames(vs) = idmap$gene(rownames(vs), to="external_gene_name")
colnames(vs) = paste(eset$batch, eset$genotype, eset$treatment, eset$replicate)

genes = intersect(unique(unlist(dor), unlist(hm)), rownames(vs))
#vs = vs[intersect(rownames(vs), genes),]
vs = vs[intersect(genes, wt_rev$label),]
vs_scaled = vs %>% narray::map(along=2, scale) #%>% narray::map(along=1, scale)
#rownames(vs_scaled)[! rownames(vs_scaled) %in% hlg] = ""

pdf(args$plotfile)
on.exit(dev.off)

pheatmap(t(vs), main="48 h")
pheatmap(t(vs_scaled), main="z-score genes WT DMSO>REV 48h FDR<0.05")



um2 = uwot::umap(vs, n_components=2)
um3 = uwot::umap(vs, n_components=5)
#knn_lou = igraph::cluster_louvain(scran::buildSNNGraph(t(um2)))$membership
knn_lou = igraph::cluster_louvain(scran::buildSNNGraph(t(vs)))$membership
df = tibble(
    gene = rownames(vs),
    um1 = um2[,1],
    um2 = um2[,2],
    clust = factor(knn_lou)
)
ggplot(df, aes(x=um1, y=um2, color=clust)) +
    geom_point() +
    ggrepel::geom_text_repel(label = ifelse(df$gene %in% hlg, df$gene, NA), color="black") +
    scale_color_brewer(palette="Paired")



koset = split(df$gene, df$clust)
scores = GSVA::gsva(vs_scaled, koset)
pheatmap(t(scores))
