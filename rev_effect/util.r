library(dplyr)
library(ggplot2)
library(patchwork)
plt = import('plot')
idmap = import('process/idmap')
gset = import('data/genesets')

plot_pcs = function(idx, pca, x, y, hl=c()) {
    imp = summary(pca)$importance[2,c(x,y)] * 100
    pcs = names(imp)
    df = cbind(idx, pca$x) %>% mutate(ins=sample %in% hl)
    ggplot(df, aes_string(x=pcs[1], y=pcs[2])) +
        geom_point(aes(size=aneuploidy, shape=tissue, fill=type, color=ins), stroke=1) +
        scale_shape_manual(values=c(22:30)) +
        scale_color_manual(values=c("#ffffff00", "black")) +
        ggrepel::geom_text_repel(aes(label=sample), color="black") +
        labs(x = sprintf("%s (%.1f%%)", pcs[1], imp[1]),
             y = sprintf("%s (%.1f%%)", pcs[2], imp[2]),
             title = "PCA plot") +
        guides(fill = guide_legend(override.aes=list(shape=21)))
}

extract_coef = function(res, coef, type="apeglm") {
    DESeq2::lfcShrink(res, coef=coef, type=type) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("label") %>%
        tbl_df() %>%
        arrange(pvalue)
}

#' Extract a data.frame with DE genes from a DESeq2 object
#'
#' @param res  DESeq2 object
#' @param rn   Name of the coefficient to extract
#' @return     A tibble with DE stats incl. 'ensembl_gene_id', 'label'
extract_result = function(res, rn) {
    rdf = DESeq2::results(res, name=rn) %>%
        as.data.frame() %>%
        as_data_frame(rownames="ensembl_gene_id") %>%
        arrange(padj) %>%
        mutate(label = idmap$gene(ensembl_gene_id, to="hgnc_symbol"),
               label = ifelse(is.na(label), ensembl_gene_id, label)) %>%
        select(ensembl_gene_id, label, everything())
}

#' Volcano plots for genes and gene sets
#'
#' @param genes  Association data.frame for genes
#' @param ...    Association data.frame for any gene sets (list of tibbles)
#' @param title  Title for the combined plot
#' @param hl     Character vector of genes or gene sets to highlight
#' @return       A ggplot2 object with multiple volcano plots
plot_volcanos = function(genes, ..., title="untitled", hl=c()) {
    gene = genes %>%
        mutate(size = sqrt(baseMean) / 5,
               circle = label %in% hl) %>%
        plt$p_effect("padj", "log2FoldChange") %>%
        plt$volcano(repel=TRUE) + ggtitle("genes")
    set = mapply(plot_gset, res=list(...), title=names(list(...)), SIMPLIFY=FALSE)
    plt$text(title) /
        patchwork::wrap_plots(c(list(genes=gene), set), nrow=1) +
        patchwork::plot_layout(heights=c(1,10))
}

#' Test genes and gene sets for a given DESeq2 data set
#'
#' @param eset    DESeq2 object
#' @param design  Design formula for differential expression
#' @param sets    A list of gene set collections (lists of character vectors)
deseq_genes_and_sets = function(eset, design, sets) {
    DESeq2::design(eset) = design
    for (v in all.vars(design))
        if (is.factor(colData(eset)[[v]]))
            colData(eset)[[v]] = droplevels(colData(eset)[[v]])

    res = DESeq2::DESeq(eset)
    res = tibble(clone=DESeq2::resultsNames(res)[-1]) %>%
        mutate(genes = lapply(clone, extract_result, res=res),
               clone = sub("clone_(.*)_vs_p53kd", "\\1", clone))
    for (ns in names(sets))
        res[[ns]] = lapply(res$genes, test_gsets, sets=sets[[ns]])
    res
}

plot_volcano = function(res, highlight=NULL) {
    if ("baseMean" %in% colnames(res))
        df = res %>%
            mutate(circle = label %in% highlight,
                   size = log10(baseMean + 1)) %>%
            plt$p_effect("padj", "log2FoldChange", thresh=0.1)
    else
        df = res %>%
            mutate(circle = label %in% highlight,
                   size = mean_expr) %>%
            plt$p_effect("adj.p", "estimate", thresh=0.1)

    plt$volcano(df, p=0.1, base.size=5, label_top=30, repel=TRUE)
}

do_wald = function(eset, fml, ex=NULL) {
    design(eset) = fml
    res = DESeq2::estimateDispersions(eset) %>%
        DESeq2::nbinomWaldTest(maxit=1000)
    if (length(ex) == 0)
        ex = setdiff(DESeq2::resultsNames(res), "Intercept")
    else
        ex = grep(ex, DESeq2::resultsNames(res), value=TRUE)
    res = sapply(ex, extract_coef, res=res, simplify=FALSE)
    if (length(res) == 1)
        res = res[[1]]
    res
}

do_lrt = function(eset, fml, red) {
    design(eset) = fml
    DESeq2::estimateDispersions(eset) %>%
        DESeq2::nbinomLRT(reduced=red, maxit=1000) %>%
        DESeq2::results() %>%
        as.data.frame() %>%
        tibble::rownames_to_column("label") %>%
        arrange(padj, pvalue)
}

test_gset = function(res, set) {
    if ("log2FoldChange" %in% colnames(res))
        cur = res %>% mutate(stat = log2FoldChange / lfcSE)
    else
        cur = res %>% mutate(stat = statistic)
    fdata = mutate(cur, in_set = label %in% set)
    mod = try(lm(stat ~ in_set, data=fdata))
    if (class(mod) == "try-error")
        return()
    broom::tidy(mod) %>%
        filter(term == "in_setTRUE") %>%
        select(-term) %>%
        mutate(size = sum(fdata$in_set, na.rm=TRUE))
}

test_gsets = function(res, sets) {
    sets = gset$filter(sets, valid=na.omit(res$label), min=4)
    result = lapply(sets, test_gset, res=res) %>%
        setNames(names(sets)) %>%
        dplyr::bind_rows(.id="label") %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)
}

plot_gset = function(res, title=NULL, highlight=NULL, fdr=0.1, base.size=0.1,
                     label_top=30, repel=TRUE) {
    p = res %>%
        plt$p_effect("adj.p", thresh=fdr) %>%
        plt$volcano(p=fdr, base.size=base.size, label_top=label_top,
                    repel=repel, text.size=2)
    if (!is.null(title))
        p = p + ggtitle(title)
    plt$try(p)
}
