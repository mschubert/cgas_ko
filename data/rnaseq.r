library(dplyr)
library(ggplot2)
library(DESeq2)
sys = import('sys')
idmap = import('process/idmap')

plot_pca = function(eset, vst) {
    pcadata = DESeq2::plotPCA(vst,
        intgroup=c("batch", "genotype", "treatment", "drug", "rev", "time", "replicate"), returnData=TRUE)
    pcavar = round(100 * attr(pcadata, "percentVar"))

    revs = as.data.frame(pcadata) %>%
        filter(time != "0") %>%
        group_by(genotype, drug, time, replicate) %>%
            filter(n() == 2 & all(c("dmso", "rev") %in% rev)) %>%
            summarize(revPC1 = PC1[rev == "rev"],
                      revPC2 = PC2[rev == "rev"],
                      PC1 = PC1[rev == "dmso"],
                      PC2 = PC2[rev == "dmso"])

    norev = as.data.frame(pcadata) %>%
        filter(rev == "dmso")

    shapes = c(dmso=21, ifng=24, ifna=25, il6=22, ask=23, il6Ab=23, jnk=1)
    colors = c(wt="#33a02c", stat1="#fb9a99", cgas="#1f78b4", stat3="#e31a1c",
               `cgas+stat1`="#cab2d6", `cgas+stat3`="#6a3d9a", relb="#b15928")

    ggplot(norev, aes(PC1, PC2)) +
        geom_point(aes(shape=drug, fill=genotype), size=3, alpha=0.8) +
        geom_segment(data=revs, aes(xend=revPC1, yend=revPC2, color=genotype),
                     arrow=arrow(length=unit(1,"mm"), type="closed"), alpha=0.8) +
        geom_text(aes(label=time), size=1.5) +
        xlab(paste0("PC1: ", pcavar[1], "% variance")) +
        ylab(paste0("PC2: ", pcavar[2], "% variance")) +
        coord_fixed() +
        scale_shape_manual(values=shapes) +
        scale_fill_manual(values=colors) +
        scale_color_manual(values=colors) +
        guides(fill = guide_legend(override.aes=list(shape=21))) +
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
        mutate(time = ifelse(time == "24", NA, time))
    ggplot(dset, aes(x=drug, y=reads)) +
        geom_point(aes(color=rev, shape=factor(replicate)), size=2,
                   position=position_jitter(w=0.15,h=0), alpha=0.5) +
        geom_text(aes(label=time), size=1.5, position=position_jitter(w=0.15,h=0)) +
        facet_grid(gene ~ genotype, scales="free", space="free_x") +
        scale_y_log10() +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        theme_light()
}

#' Read a count table into a matrix
load_table = function(fname) {
    reads = readr::read_tsv(fname)
    rmat = data.matrix(reads[-1])
    rownames(rmat) = sub("\\.[0-9]+$", "", reads$gene_id)
    rmat
}

sys$run({
    args = sys$cmd$parse(
        opt('c', 'config', 'yaml', '../config.yaml'),
        opt('s', 'samples', 'tsv', 'rnaseq.tsv'),
        opt('o', 'outfile', 'rds', 'rnaseq.rds'),
        opt('p', 'plotfile', 'pdf', 'rnaseq.pdf'),
        arg('infiles', 'count tables', arity='*', list.files("count_tables",
            "STL_and_USS_genes.*\\.txt\\.gz", recursive=TRUE, full.names=TRUE))
    )

    cfg = yaml::read_yaml(args$config)
    rtabs = lapply(args$infiles, load_table)

    # 18k genes unique to batch1, 9k unique to batch2, 21k common
    reads = na.omit(narray::stack(rtabs, along=2))
    reads = reads[rowSums(reads) > ncol(reads),]
    samples = readr::read_tsv("rnaseq.tsv") %>%
        mutate(batch = factor(batch),
               genotype = relevel(factor(genotype), "wt"),
               treatment = relevel(factor(treatment), "dmso"))
    samples$rev = relevel(factor(ifelse(grepl("rev", samples$treatment), "rev", "dmso")), "dmso")
    samples$drug = sub("\\+?rev\\+?", "", samples$treatment)
    samples$drug[samples$drug == ""] = "dmso"
    samples$drug = relevel(factor(samples$drug), "dmso")
    narray::intersect(reads, samples$sample_id, along=2)
    eset = DESeq2::DESeqDataSetFromMatrix(reads, samples, ~1) %>%
        DESeq2::estimateSizeFactors()

    vst = DESeq2::varianceStabilizingTransformation(eset)
    pdf(args$plotfile, 10, 8)
    print(plot_pca(eset, vst) + ggtitle("PCA all samples"))
    for (batch in levels(colData(eset)$batch)) {
        b = colData(eset)$batch == batch
        print(plot_pca(eset[,b], vst[,b]) + ggtitle(paste("PCA batch", batch)))
    }

    eset = eset[,! colnames(eset) %in% names(cfg$rnaseq_exclude)]
    vst = vst[,! colnames(vst) %in% names(cfg$rnaseq_exclude)]
    print(plot_pca(eset, vst) + ggtitle("PCA all samples (no outliers)"))
    for (batch in levels(colData(eset)$batch)) {
        b = colData(eset)$batch == batch
        print(plot_pca(eset[,b], vst[,b]) + ggtitle(paste("PCA batch", batch, "(no outliers)")))
    }

    all48 = colData(eset)$time == "48"
    print(plot_pca(eset[,all48], vst[,all48]) + ggtitle("PCA 48h samples (no outliers)"))
    rev48 = colData(eset)$time == "48" & colData(eset)$drug != "ifna"
    print(plot_pca(eset[,rev48], vst[,rev48]) + ggtitle("PCA rev 48h samples (no outliers)"))

#    print(plot_dist(eset, vst))
    print(plot_kos(eset))
    dev.off()

    saveRDS(eset, file=args$outfile)
})
