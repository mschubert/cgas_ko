library(dplyr)
library(DESeq2)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')

args = sys$cmd$parse(
    opt('s', 'short', 'rds', 'rep_compare-24h.rds'),
    opt('l', 'long', 'rds', 'rep_compare-48h.rds'),
    opt('p', 'plotfile', 'pdf', 'inflamm_resp_overview.pdf')
)

cmp = c("Interferon Gamma Response",
        "Inflammatory Response",
        "TNF-alpha Signaling via NF-kB")

res = list("24h" = readRDS(args$short), "48h" = readRDS(args$long)) %>%
    bind_rows(.id="time") %>%
    select(cond, time, rep, MSigDB_Hallmark_2020) %>%
    tidyr::unnest("MSigDB_Hallmark_2020") %>%
    filter(label %in% cmp)

pdf(args$plotfile, 8, 12)
ggplot(res, aes(y=label, fill=rep, x=statistic)) +
    geom_col(position = position_dodge2(preserve = "single")) +
    geom_vline(xintercept=0) +
    facet_grid(cond ~ time)
dev.off()
