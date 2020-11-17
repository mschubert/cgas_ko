sys = import('sys')
idmap = import('process/idmap')
gset = import('data/genesets')
enr = import('tools/enrichr')
msdb = import('tools/msigdb')

args = sys$cmd$parse(
    opt('g', 'geneset', 'Identifier of the gene set', 'KEA_2015'),
    opt('h', 'human', 'save to RData', 'human/KEA_2015.rds'),
    opt('m', 'mouse', 'save to RData', 'mouse/KEA_2015.rds'))

if (args$geneset == "GO_Biological_Process_2020") {
    sets = gset$go(as_list=TRUE)
} else if (args$geneset == "DoRothEA") {
    sets = dorothea::dorothea_hs %>%
        filter(mor == 1) %>%
        group_by(tf) %>%
            filter(case_when(
                sum(confidence %in% c("A")) >= 20 ~ confidence %in% c("A"),
                sum(confidence %in% c("A", "B")) >= 20 ~ confidence %in% c("A", "B"),
                sum(confidence %in% c("A", "B", "C")) >= 20 ~ confidence %in% c("A", "B", "C"),
                sum(confidence %in% c("A", "B", "C", "D")) >= 20 ~ confidence %in% c("A", "B", "C", "D"),
                TRUE ~ FALSE
            )) %>%
            mutate(tf = sprintf("%s (%s)", tf, tolower(tail(sort(confidence), 1)))) %>%
        ungroup() %>%
        select(target, tf) %>%
        unstack()
} else if (args$geneset %in% enr$dbs()$name) {
    sets = enr$genes(args$geneset)
} else if (args$geneset %in% msdb$dbs()) {
    sets = msdb$genes(args$geneset)
} else
    stop("invalid gene set: ", args$geneset)

mouse = stack(sets)
mouse$values = unname(idmap$orthologue(mouse$values, from="external_gene_name", to="mgi_symbol"))
mouse = unstack(na.omit(mouse))

saveRDS(sets, file=args$human)
saveRDS(mouse, file=args$mouse)
