library(dplyr)
library(ggplot2)
library(patchwork)
plt = import('plot')

wt = "../comp_list/cmp/BT549 WT Reversine vs DMSO.xlsx"
cgas = "../comp_list/cmp/BT549 cGAS KO Reversine vs DMSO.xlsx"
sheets = c("genes", "MSigDB_Hallmark_2020", "DoRothEA")

res = list(
    wt = sapply(sheets, function(s) readxl::read_excel(wt, sheet=s)),
    cgas = sapply(sheets, function(s) readxl::read_excel(cgas, sheet=s))
)

wt = list(
    plt$volcano(res$wt[[1]], p=0.25, label_top=40) + labs(title="BT549 wild-type", subtitle="Genes", y="FDR"),
    plt$volcano(res$wt[[2]], p=0.05, label_top=20) + labs(subtitle="MSigDB Hallmarks", y="FDR"),
    plt$volcano(res$wt[[3]], p=0.05, label_top=20) + labs(subtitle="DoRothEA TF regulons", y="FDR")
)
cgas = list(
    plt$volcano(res$cgas[[1]], p=0.25, label_top=40) + labs(title="BT549 cGAS KO", subtitle="Genes", y="FDR"),
    plt$volcano(res$cgas[[2]], p=0.05, label_top=20) + labs(subtitle="MSigDB Hallmarks", y="FDR"),
    plt$volcano(res$cgas[[3]], p=0.05, label_top=20) + labs(subtitle="DoRothEA TF regulons", y="FDR")
)
volcs = wrap_plots(wt) / wrap_plots(cgas)

pdf("supp_rnaseq.pdf", 14, 14)
print(volcs)
dev.off()
