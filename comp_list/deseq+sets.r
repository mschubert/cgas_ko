library(dplyr)
library(DESeq2)
sys = import('sys')
idmap = import('process/idmap')
gset = import('genesets')
util = import('../rev_effect/util')

deseq_and_sets = function(eset, sets) {
    genes = DESeq2::DESeq(eset) %>%
        DESeq2::results(name=attr(eset, "extract")) %>%
        as.data.frame() %>% tibble::rownames_to_column("ensembl_gene_id") %>%
        as_tibble() %>%
        mutate(label = idmap$gene(ensembl_gene_id, to="external_gene_name")) %>%
        select(ensembl_gene_id, label, everything()) %>%
        arrange(padj, pvalue)

    sets = lapply(sets, util$test_gsets, res=genes)
    c(list(genes=genes), sets)
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'rds', 'eset/BT549 WT reversine vs BT549 cGAS KO reversine.rds'),
        opt('o', 'outfile', 'xlsx', 'cmp/BT549 WT reversine vs BT549 cGAS KO reversine.xlsx')
    )

    eset = readRDS(args$infile)
    sets = gset$get_human(c(
        "MSigDB_Hallmark_2020",
        "GO_Biological_Process_2021",
        "KEGG_2021_Human"
    ))

    res = deseq_and_sets(eset, sets)
    cmp = attr(eset, "cmp")
    res = c(res, list(cmp=cmp))

    writexl::write_xlsx(res, path=args$outfile)
})
