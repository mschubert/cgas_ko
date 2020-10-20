library(dplyr)
io = import('io')
sys = import('sys')
gset = import('data/genesets')
util = import('../expr_diff/util')

args = sys$cmd$parse(
    opt('d', 'diff_expr', 'rds', 'diff_expr.rds'),
    opt('s', 'setfile', 'RData', '../data/genesets/human/CH.HALLMARK.RData'),
    opt('p', 'plotfile', 'pdf', 'stat1ko_diff_inflamm_resp.pdf'))

res = readRDS(args$diff_expr)
sets = io$load(args$setfile)
sets = sets[grepl("INTERFERON|NFKB", names(sets))]

assocs = expand.grid(rname = names(res), sname=names(sets), stringsAsFactors=FALSE) %>%
    mutate(assocs = purrr::map2(rname, sname, function(r, s) util$test_gset(res[[r]], sets[[s]])),
           rname = factor(rname, levels=names(res)),
           sname = factor(sname, levels=names(sets))) %>%
    tidyr::unnest(cols="assocs")

p = ggplot(assocs, aes(x=sname, y=statistic)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~ rname, ncol=1)

pdf(args$plotfile, 8, 12)
print(p)
dev.off()
