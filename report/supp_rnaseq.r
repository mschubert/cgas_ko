library(dplyr)
library(ggplot2)
library(patchwork)
plt = import('plot')

wt = "../comp_list/cmp/BT549 WT Reversine vs DMSO.xlsx"
cgas = "../comp_list/cmp/BT549 cGAS KO Reversine vs DMSO.xlsx"
cmp = "../comp_list/cmp/BT549 cGAS KO reversine vs BT549 WT reversine.xlsx"
sheets = c("genes", "MSigDB_Hallmark_2020", "DoRothEA")

res = list(
    wt = sapply(sheets, function(s) readxl::read_excel(wt, sheet=s)),
    cgas = sapply(sheets, function(s) readxl::read_excel(cgas, sheet=s)),
    cmp = sapply(sheets, function(s) readxl::read_excel(cmp, sheet=s))
)

wt = list(
    plt$volcano(res$wt[[1]], p=0.25, label_top=40) +
        labs(title="BT549 wild-type: Reversine vs. DMSO", subtitle="Genes", y="Adjusted p-value (FDR)"),
    plt$volcano(res$wt[[2]], p=0.05, label_top=20) + labs(subtitle="MSigDB Hallmarks", y="Adjusted p-value (FDR)"),
    plt$volcano(res$wt[[3]], p=0.05, label_top=20) + labs(subtitle="DoRothEA TF regulons", y="Adjusted p-value (FDR)")
)
cgas = list(
    plt$volcano(res$cgas[[1]], p=0.25, label_top=40) +
        labs(title="BT549 cGAS KO: Reversine vs. DMSO", subtitle="Genes", y="Adjusted p-value (FDR)"),
    plt$volcano(res$cgas[[2]], p=0.05, label_top=20) + labs(subtitle="MSigDB Hallmarks", y="Adjusted p-value (FDR)"),
    plt$volcano(res$cgas[[3]], p=0.05, label_top=20) + labs(subtitle="DoRothEA TF regulons", y="Adjusted p-value (FDR)")
)
cmp = list(
    plt$volcano(res$cmp[[1]], p=0.25, label_top=40) +
        labs(title="BT549 Reversine: cGAS KO vs. WT", subtitle="Genes", y="Adjusted p-value (FDR)"),
    plt$volcano(res$cmp[[2]], p=0.05, label_top=20) + labs(subtitle="MSigDB Hallmarks", y="Adjusted p-value (FDR)"),
    plt$volcano(res$cmp[[3]], p=0.05, label_top=20) + labs(subtitle="DoRothEA TF regulons", y="Adjusted p-value (FDR)")
)
volcs = wrap_plots(wt) / wrap_plots(cgas) / wrap_plots(cmp)

pdf("supp_rnaseq.pdf", 17, 21)
print(volcs)
dev.off()
