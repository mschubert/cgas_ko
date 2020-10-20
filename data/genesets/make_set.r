sys = import('sys')
idmap = import('process/idmap')
enr = import('tools/enrichr')
msdb = import('tools/msigdb')

args = sys$cmd$parse(
    opt('g', 'geneset', 'Identifier of the gene set', 'KEA_2015'),
    opt('h', 'human', 'save to RData', 'human/KEA_2015.RData'),
    opt('m', 'mouse', 'save to RData', 'mouse/KEA_2015.RData'))

if (args$geneset %in% enr$dbs()$name) {
    sets = enr$genes(args$geneset)
} else if (args$geneset %in% msdb$dbs()) {
    sets = msdb$genes(args$geneset)
} else
    stop("invalid gene set: ", args$geneset)

mouse = stack(sets)
mouse$values = unname(idmap$orthologue(mouse$values, from="external_gene_name", to="mgi_symbol"))
mouse = unstack(na.omit(mouse))

save(sets, file=args$human)
save(mouse, file=args$mouse)
