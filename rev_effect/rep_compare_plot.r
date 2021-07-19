library(dplyr)
library(DESeq2)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
plt2 = import('./plot')

do_compare = function(res24, res48, .rep, .cond) {
    message(.cond)
    subs = . %>%
        filter(cond == .cond,
               rep %in% c("1", .rep)) %>%
        mutate(cond = paste(make.names(cond), rep, sep="."))
    subs24 = subs(res24)
    subs48 = subs(res48)
    conds = unique(c(subs24$cond, subs48$cond))

    do_plot = function(df, col, time) {
        if (nrow(df) == 2)
            plt2$compare(df, ref=conds[1], cmp=conds[2], col=col, thresh=1.5) +
                labs(title=col, subtitle=time)
        else
            plt$text(sprintf("%s : %i samples", col, nrow(df)))
    }

    plots = list(
        do_plot(subs24, "genes", "24"),
        do_plot(subs24, "MSigDB_Hallmark_2020", "24"),
        do_plot(subs24, "GO_Biological_Process_2021", "24"),
        do_plot(subs48, "genes", "48"),
        do_plot(subs48, "MSigDB_Hallmark_2020", "48"),
        do_plot(subs48, "GO_Biological_Process_2021", "48")
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

    idx = bind_rows(res24[c("cond", "rep")], res48[c("cond", "rep")]) %>%
        distinct() %>%
        filter(rep != "1") %>%
        arrange(cond, rep) %>%
        rowwise() %>%
        mutate(plot = list(do_compare(res24, res48, rep, cond)))

    pdf(args$plotfile, 18, 12)
    for (i in seq_len(nrow(idx)))
        print(idx$plot[[i]])
    dev.off()
})
