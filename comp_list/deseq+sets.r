library(dplyr)
library(DESeq2)
sys = import('sys')
idmap = import('process/idmap')
gset = import('genesets')

deseq_and_sets = function(eset, sets) {
    genes = DESeq2::DESeq(eset) %>%
        DESeq2::results(name=attr(eset, "extract")) %>%
        as.data.frame() %>% tibble::rownames_to_column("ensembl_gene_id") %>%
        as_tibble() %>%
        mutate(label = idmap$gene(ensembl_gene_id, to="external_gene_name")) %>%
        select(ensembl_gene_id, label, everything()) %>%
        arrange(padj, pvalue)

    sres = lapply(sets, function(s) gset$test_lm(genes, s, add_means="log2FoldChange"))
    c(list(genes=genes), sres)
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'rds', 'eset/BT549 cGAS KO Reversine vs DMSO.rds'),
        opt('o', 'outfile', 'xlsx', 'cmp/BT549 cGAS KO Reversine vs DMSO.xlsx')
    )

    eset = readRDS(args$infile)
    sets = gset$get_human(c(
        "MSigDB_Hallmark_2020",
        "GO_Biological_Process_2021",
#        "KEGG_2021_Human",
        "DoRothEA",
        "CIN"
    ))

    res = deseq_and_sets(eset, sets)
    cmp = attr(eset, "cmp")
    res = c(res, list(cmp=cmp))

    writexl::write_xlsx(res, path=args$outfile)
})
