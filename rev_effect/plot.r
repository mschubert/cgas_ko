library(dplyr)
library(DESeq2)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
util = import('./util')

volcano_or_compare = function(ref, cmp) {
    if (ref == cmp)
        return(volcano(ref))

    sr = rlang::sym(ref)
    sc = rlang::sym(cmp)
    rdf = inner_join(res[[ref]] %>% transmute(label, !! sr := stat),
                     res[[cmp]] %>% transmute(label, !! sc := stat)) %>%
        filter((!! sr)^2 + (!! sc)^2 > 2.5^2)

    ggplot(rdf, aes_string(x=cmp, y=ref)) +
        geom_point() +
        theme_classic() #todo: colors for up, down, inconsistent, select labels, etc
}

volcano = function(rname, hl=c()) {
    res[[rname]] %>%
        mutate(size = log10(baseMean + 1),
               circle = label %in% hl) %>%
        plt$p_effect("padj", "log2FoldChange") %>%
        plt$volcano(repel=TRUE) + ggtitle(rname)
}

args = sys$cmd$parse(
    opt('e', 'eset', 'rds', '../data/rnaseq.rds'),
    opt('o', 'outfile', 'rds', 'deseq.rds'),
    opt('p', 'plotfile', 'pdf', 'deseq.pdf'),
    arg('setfiles', 'rds', arity='*', '../data/genesets/human/CH.HALLMARK.rds')
)

sets = sapply(args$setfiles, readRDS)
names(sets) = basename(tools::file_path_sans_ext(names(sets)))

eset = readRDS(args$eset)
colData(eset)$cond = ifelse(grepl("il6", colData(eset)$treatment), "wt+il6",
                            as.character(colData(eset)$genotype))
colData(eset)$cond = relevel(factor(colData(eset)$cond), "wt")
colData(eset)$rev = ifelse(grepl("rev", colData(eset)$treatment), 1, 0)

design(eset) = ~ cond:rev
res = DESeq(eset)
res = resultsNames(res) %>%
    setdiff("Intercept") %>%
    sapply(util$extract_result, res=res, simplify=FALSE)
res = tibble(cond = sub("\\.rev", "", sub("^cond", "", names(res))),
             genes = unname(res))


cmp = tidyr::crossing(tibble(ref = res$cond),
                      tibble(cmp = res$cond)) %>%
    rowwise() %>%
        mutate(plot = list(volcano_or_compare(ref, cmp))) %>%
    ungroup()

plots = wrap_plots(cmp$plot, ncol=length(res))

pdf("rev_effect.pdf", 50, 50)
print(plots)
dev.off()
