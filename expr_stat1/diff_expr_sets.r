library(dplyr)
io = import('io')
sys = import('sys')
gset = import('data/genesets')
util = import('../expr_diff/util')

args = sys$cmd$parse(
    opt('d', 'diff_expr', 'rds', 'diff_expr.rds'),
    opt('s', 'setfile', 'RData', '../data/genesets/human/KEA_2015.RData'),
    opt('p', 'plotfile', 'pdf', 'sets/KEA_2015.pdf'))

res = readRDS(args$diff_expr)
sets = io$load(args$setfile) %>%
    gset$filter(min=5, valid=res[[1]]$gene_name)

pdf(args$plotfile)
for (rname in names(res)) {
    message(rname)
    print(util$plot_gset(res[[rname]], sets) + ggtitle(rname))
}
dev.off()
