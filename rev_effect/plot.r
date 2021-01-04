library(dplyr)
library(DESeq2)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
util = import('./util')

#' Compare separation statistics between two samples
#'
#' @param res         A data.frame with columns specified in `ref`, `cmp` and `col`
#' @param ref         Condition to be plotted on the x axis (character)
#' @param cmp         Condition to be plotted on the y axis (character)
#' @param col         Column in data.frame `res` to use (default: 'genes')
#' @param thresh      Do not plot points below this value to prevent crowding (default: 2)
#' @param label_base  How many points to label, scaled by label length (default: 200)
#' @param hl          Character vector of filled+border circles (default: none)
#' @param hollow      Character vector of circle-only points (default: none)
#' @return            A ggplot2 object with a scatter plot of x/y separation
compare = function(res, ref, cmp, col="genes", thresh=2, label_base=200, hl=c(), hollow=c()) {
    sym = list(x=rlang::sym(ref), y=rlang::sym(cmp), col=rlang::sym(col))
    refdf = res %>% filter(cond == ref) %>% pull(!! sym$col) %>% `[[`(1)
    cmpdf = res %>% filter(cond == cmp) %>% pull(!! sym$col) %>% `[[`(1)
    if (! "stat" %in% colnames(refdf)) {
        refdf = refdf %>% mutate(stat = statistic)
        cmpdf = cmpdf %>% mutate(stat = statistic)
    }

    avg_chr_label = mean(nchar(refdf$label), na.rm=TRUE)
    topf = function(dist, frac)
        rank(-dist, ties.method="random") < frac*label_base/sqrt(avg_chr_label)
    colors = setNames(c("#a50f1567", "#006d2c67", "#045a8d67", "#78787822", "#4d004b67"),
                      c("down", "up", "opposite", "no change", "condition specific"))

    rdf = inner_join(refdf %>% transmute(label, !! sym$x := stat),
                     cmpdf %>% transmute(label, !! sym$y := stat)) %>%
        filter((!! sym$x)^2 + (!! sym$y)^2 > thresh^2) %>%
        mutate(quad = factor(case_when(
                  abs(!! sym$x) > thresh*2 & abs(!! sym$y) < thresh ~ "condition specific",
                  abs(!! sym$x) < thresh & abs(!! sym$y) > thresh*2 ~ "condition specific",
                  abs(!! sym$x) < thresh | abs(!! sym$y) < thresh ~ "no change",
                  sign(!! sym$x) > 0 & sign(!! sym$y) > 0 ~ "up",
                  sign(!! sym$x) < 0 & sign(!! sym$y) < 0 ~ "down",
                  TRUE ~ "opposite"
               ), levels=names(colors)),
               dist = sqrt(scale(!! sym$x, center=FALSE)^2 +
                           scale(!! sym$y, center=FALSE)^2),
               showlab = ifelse(topf(dist,0.33), TRUE, FALSE)) %>% # third of labels for top
        group_by(quad) %>%
            mutate(showlab = ifelse(topf(dist,0.66*0.25), TRUE, showlab)) %>% # others per quad
        ungroup() %>%
        mutate(shape = ifelse(label %in% hollow, 1, 16),
               label = ifelse(showlab | label %in% hl, label, NA))

    ggplot(rdf, aes_string(x=ref, y=cmp)) +
        geom_vline(xintercept=0, size=1, color="grey", linetype="dashed") +
        geom_hline(yintercept=0, size=1, color="grey", linetype="dashed") +
        geom_point(aes(color=quad, shape=shape), size=2) +
        scale_shape_identity() +
        ggrepel::geom_label_repel(aes(label=label, color=quad), size=2,
            label.size=NA, alpha=1, segment.alpha=0.3, min.segment.length=0,
            fill="#ffffffc0", label.padding=0.1, max.iter=1e5) +
        theme_classic() +
        scale_color_manual(name="change", values=colors, drop=FALSE)
}

#' Volcano plot of differential expression for a sample
#'
#' @param res    A data.frame with columns specified in `ref`, `cmp` and `col`
#' @param rname  Condition to be plotted (character)
#' @param col    Column in data.frame `res` to use (default: 'genes')
#' @param hl     Character vector of filled+border circles (default: none)
#' @return       A ggplot2 object of volcano plot
volcano = function(res, rname, col, hl=c()) {
    rdf = res %>%
        filter(cond == rname) %>%
        pull(!! rlang::sym(col)) %>%
        `[[`(1) %>%
        mutate(circle = label %in% hl)

    if ("baseMean" %in% colnames(rdf))
        rdf %>%
            mutate(size = log10(baseMean + 1)) %>%
            plt$p_effect("padj", "log2FoldChange", thresh=0.1) %>%
            plt$volcano(p=0.1, repel=TRUE) + ggtitle(rname)
    else
        rdf %>%
            plt$p_effect("adj.p", "estimate", thresh=0.1) %>%
            plt$volcano(p=0.1, repel=TRUE, base.size=0.1, text.size=2.5) + ggtitle(rname)
}

#' Plot grid of volcano (diagonal) or condition compare (off-diagonal)
#'
#' @param res    A data.frame with columns specified in `ref`, `cmp` and `col`
#' @param col    Column in data.frame `res` to use (default: 'genes')
#' @return       A patchwork grid of volcano plots or correlation plots
plot_matrix = function(res, col="genes") {
    volcano_or_compare = function(ref, cmp, col="genes", thresh=2, label_base=200, hl=c(), hollow=c()) {
        if (ref == cmp)
            volcano(res, ref, col, hl)
        else
            compare(res, ref, cmp, col, thresh, label_base, hl, hollow)
    }

    cmp = tidyr::crossing(tibble(cmp = res$cond),
                          tibble(ref = res$cond)) %>%
        rowwise() %>%
            mutate(plot = list(volcano_or_compare(ref, cmp, col))) %>%
        ungroup()

    plt$text(col, size=12) /
        wrap_plots(cmp$plot, ncol=nrow(res), guides="collect") +
        plot_layout(heights=c(1,100))
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'rds', 'deseq.rds'),
        opt('p', 'plotfile', 'pdf', 'rev_effect.pdf')
    )

    res = readRDS(args$infile)
    res$cond = make.names(res$cond)
    plots = sapply(colnames(res)[-1], plot_matrix, res=res, simplify=FALSE)

    pdf(args$plotfile, 70, 70)
    for (p in plots)
        print(p)
    dev.off()
})
