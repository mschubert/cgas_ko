library(dplyr)
library(ggplot2)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'RData', '../expr_diff/eset_MILE.RData'),
    opt('p', 'plotfile', 'pdf', 'mile_stat-ets-erg.pdf'))

eset = io$load(args$infile)
dset = cbind(eset$meta, t(eset$expr[c("ETS1","ERG","STAT1"),])) %>%
    select(annot, lineage, type, ETS1, ERG, STAT1) %>%
    tidyr::gather("gene", "expr", -(annot:type))

pdf(args$plotfile, 14, 10)
ggplot(dset %>% filter(grepl("(mature B-ALL)|(ALL.*MLL)", annot)), aes(x=annot, y=expr)) +
#    ggbeeswarm::geom_beeswarm() +
    geom_boxplot() +
    facet_grid(. ~ gene) +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

ggplot(dset %>% filter(type == "B_like"), aes(x=annot, y=expr)) +
#    ggbeeswarm::geom_beeswarm() +
    geom_boxplot() +
    facet_grid(. ~ gene) +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

ggplot(dset, aes(x=annot, y=expr)) +
#    ggbeeswarm::geom_beeswarm() +
    geom_boxplot() +
    facet_grid(. ~ gene) +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
dev.off()
