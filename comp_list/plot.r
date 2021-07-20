library(dplyr)
library(DESeq2)
sys = import('sys')
plt = import('plot')

args = sys$cmd$parse(
    opt('i', 'infile', 'xlsx', 'cmp/BT549 WT reversine vs BT549 cGAS KO reversine.xlsx'),
    opt('p', 'plotfile', 'pdf', 'cmp/BT549 WT reversine vs BT549 cGAS KO reversine.pdf')
)

sheets = readxl::excel_sheets(args$infile)
res = sapply(sheets, readxl::read_xlsx, path=args$infile, simplify=FALSE)

pdf(12, 10, file=args$plotfile)
for (i in seq_len(length(res)-1))
    print(plt$volcano(res[[i]], label_top=30, repel=TRUE) + ggtitle(names(res)[i]))

gridExtra::grid.arrange(top="Fitting coefficients", gridExtra::tableGrob(res$cmp))
dev.off()
