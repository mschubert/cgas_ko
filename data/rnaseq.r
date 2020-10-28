library(dplyr)
library(ggplot2)
sys = import('sys')
idmap = import('process/idmap')

plot_pca = function(eset) {
    vst = DESeq2::varianceStabilizingTransformation(eset)
    pcadata = DESeq2::plotPCA(vst, intgroup=c("batch", "genotype", "treatment"), returnData=TRUE)
    pcavar = round(100 * attr(pcadata, "percentVar"))

    ggplot(pcadata, aes(PC1, PC2)) +
        geom_point(aes(shape=factor(batch), fill=genotype, color=treatment), size=3) +
#        scale_shape_manual(values=setNames(c(21, 24, 25), c("parental", "gains", "losses"))) +
        xlab(paste0("PC1: ", pcavar[1], "% variance")) +
        ylab(paste0("PC2: ", pcavar[2], "% variance")) +
        coord_fixed() +
        ggrepel::geom_text_repel(aes(label=paste(genotype, treatment, sep=":"))) +
        theme_classic()
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
    samples = readr::read_tsv("rnaseq.tsv")
    eset = DESeq2::DESeqDataSetFromMatrix(reads, samples, ~1)

    p = plot_pca(eset)
    pdf(args$plotfile, 10, 8)
    print(p)
    dev.off()

    saveRDS(eset, file=args$outfile)
})
