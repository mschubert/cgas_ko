library(dplyr)
library(ggplot2)
library(DESeq2)
sys = import('sys')
idmap = import('process/idmap')

plot_pca = function(eset, vst) {
    pcadata = DESeq2::plotPCA(vst, intgroup=c("batch", "genotype", "treatment","time"), returnData=TRUE)
    pcavar = round(100 * attr(pcadata, "percentVar"))

    ggplot(pcadata, aes(PC1, PC2)) +
        geom_point(aes(shape=factor(batch), fill=genotype, color=treatment), size=3) +
#        scale_shape_manual(values=setNames(c(21, 24, 25), c("parental", "gains", "losses"))) +
        xlab(paste0("PC1: ", pcavar[1], "% variance")) +
        ylab(paste0("PC2: ", pcavar[2], "% variance")) +
        coord_fixed() +
        ggrepel::geom_text_repel(aes(label=paste(genotype, treatment, sep=":")), size=2) +
        theme_classic()
}

#plot_dist = function(eset, vst) {
#    meta = as.data.frame(SummarizedExperiment::colData(eset)) %>%
#        mutate()
#    sample_dist = as.matrix(dist(t(SummarizedExperiment::assay(vst))))
#
#}

plot_kos = function(eset) {
    genes = c("CGAS", "STAT1", "STAT3", "IL6", "IFNA", "RELB", "CCL5", "CXCL10", "MYC", "BIRC5")
    nr = counts(eset, normalized=TRUE)
    rownames(nr) = idmap$gene(rownames(nr), to="hgnc_symbol")
    names(dimnames(nr)) = c("gene", "sample_id")
    dset = reshape2::melt(nr, value.name="reads") %>%
        as_tibble() %>%
        filter(gene %in% genes) %>%
        inner_join(colData(eset) %>% as.data.frame()) %>%
        mutate(cond = paste(genotype, sub("\\+$", "", sub("\\+?dmso|none|rev", "", treatment))),
               rev = ifelse(grepl("rev", treatment), "rev", ifelse(grepl("dmso", treatment), "dmso", "none")),
               rev = factor(rev, levels=c("rev", "dmso", "none")),
               treatment = gsub("\\+| ", "", sub("dmso|rev|none", "", treatment)),
               time = ifelse(time == 48, "48", NA))
    ggplot(dset, aes(x=treatment, y=reads)) +
        geom_point(aes(color=rev, shape=factor(replicate)), size=2,
                   position=position_jitter(w=0.15,h=0), alpha=0.5) +
        geom_text(aes(label=time), size=1.5, position=position_jitter(w=0.15,h=0)) +
        facet_grid(gene ~ genotype, scales="free", space="free_x") +
        scale_y_log10() +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        theme_light()
}

sys$run({
    args = sys$cmd$parse(
        opt('s', 'samples', 'tsv', 'rnaseq.tsv'),
        opt('o', 'outfile', 'rds', 'rnaseq.rds'),
        opt('p', 'plotfile', 'pdf', 'rnaseq.pdf')
    )

    sept2019 = readr::read_tsv("rnaseq_new/count_matrix_known_barcodes_STL_and_USS_genes_2019.txt.gz")
    oct2020 = readr::read_tsv("rnaseq_new/count_matrix_known_barcodes_STL_and_USS_genes_2020.txt.gz")

    batch1 = data.matrix(sept2019[-1])
    rownames(batch1) = sub("\\.[0-9]+$", "", sept2019$gene_id)
    batch2 = data.matrix(oct2020[-1])
    rownames(batch2) = sub("\\.[0-9]+$", "", oct2020$gene_id)

    # 18k genes unique to batch1, 9k unique to batch2, 21k common
    reads = na.omit(narray::stack(batch1, batch2, along=2))
    reads = reads[rowSums(reads) > ncol(reads),]
    samples = readr::read_tsv("rnaseq.tsv") %>%
        mutate(batch = factor(batch),
               genotype = relevel(factor(genotype), "wt"),
               treatment = relevel(factor(treatment), "dmso"))
    narray::intersect(reads, samples$sample_id, along=2)
    eset = DESeq2::DESeqDataSetFromMatrix(reads, samples, ~1) %>%
        DESeq2::estimateSizeFactors()

    vst = DESeq2::varianceStabilizingTransformation(eset)
    pdf(args$plotfile, 10, 8)
    print(plot_pca(eset, vst))
#    print(plot_dist(eset, vst))
    print(plot_kos(eset))
    dev.off()

    saveRDS(eset, file=args$outfile)
})
