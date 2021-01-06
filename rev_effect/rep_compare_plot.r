library(dplyr)
library(DESeq2)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
plt2 = import('./plot')

do_compare = function(res24, res48, .cond) {
    subs = . %>%
        filter(cond == .cond) %>%
        mutate(cond = paste(make.names(cond), rep, sep="."))
    subs24 = subs(res24)
    subs48 = subs(res48)
    conds = subs24$cond

    do_plot = function(df, col, time)
        plt2$compare(df, ref=conds[1], cmp=conds[2], col=col, thresh=1.5) +
            labs(title=col, subtitle=time)

    plots = list(
        do_plot(subs24, "genes", "24"),
        do_plot(subs24, "MSigDB_Hallmark_2020", "24"),
        do_plot(subs24, "GO_Biological_Process_2020", "24"),
        do_plot(subs48, "genes", "48"),
        do_plot(subs48, "MSigDB_Hallmark_2020", "48"),
        do_plot(subs48, "GO_Biological_Process_2020", "48")
    )
    plt$text(.cond, size=10) / wrap_plots(plots, guides="collect", nrow=2) +
        plot_layout(ncol=1, heights=c(1,20))
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
    pdf(args$plotfile, 18, 12)
    for (cond in conds) {
        print(do_compare(res24, res48, cond))
    }
    dev.off()
})
