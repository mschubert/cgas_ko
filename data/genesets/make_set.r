library(dplyr)
sys = import('sys')
gset = import('genesets')

args = sys$cmd$parse(
    opt('g', 'geneset', 'Identifier of the gene set', 'KEA_2015'),
    opt('h', 'human', 'save to rds', 'human/KEA_2015.rds'),
    opt('m', 'mouse', 'save to rds', 'mouse/KEA_2015.rds')
)

sets = gset$get_human(args$geneset)
mouse = gset$hu2mouse(sets)

saveRDS(sets, file=args$human)
saveRDS(mouse, file=args$mouse)
