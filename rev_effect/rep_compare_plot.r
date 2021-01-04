library(dplyr)
library(DESeq2)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('./plot')

do_compare = function(res24, res48, .cond) {
    subs = . %>%
        filter(cond == .cond) %>%
        mutate(cond = paste(make.names(cond), rep, sep="."))
    subs24 = subs(res24)
    subs48 = subs(res48)
    conds = subs24$cond

    plots = list(
        plt$compare(subs24, conds[1], conds[2], "genes") + ggtitle("genes 24h"),
        plt$compare(subs24, conds[1], conds[2], "MSigDB_Hallmark_2020") + ggtitle("Hallmarks 24h"),
        plt$compare(subs48, conds[1], conds[2], "genes") + ggtitle("genes 48h"),
        plt$compare(subs48, conds[1], conds[2], "MSigDB_Hallmark_2020") + ggtitle("Hallmarks 48h")
    )
    patchwork::wrap_plots(plots)
}

sys$run({
    args = sys$cmd$parse(
        opt('s', 'short', 'rds', 'rep_compare-24h.rds'),
        opt('l', 'long', 'rds', 'rep_compare-48h.rds'),
        opt('p', 'plotfile', 'pdf', 'rep_compare_plot.pdf')
    )

    res24 = readRDS(args$short)
    res48 = readRDS(args$long)

    conds = intersect(res24$cond, res48$cond)
    pdf(args$plotfile, 12, 12)
    for (cond in conds) {
        print(do_compare(res24, res48, cond))
    }
    dev.off()
})
