library(dplyr)
library(DESeq2)
sys = import('sys')
idmap = import('process/idmap')
util = import('../rev_effect/util')

deseq_and_sets = function(eset, sets) {
    mm = model.matrix(design(eset), colData(eset))
    cmp = cbind(keep[c("batch", "genotype", "treatment")], mm) %>%
        as.data.frame() %>% tibble::rownames_to_column("sample_id") %>% as_tibble()

    genes = DESeq2::DESeq(eset) %>%
        DESeq2::results(name=rec$extract) %>%
        as.data.frame() %>% tibble::rownames_to_column("ensembl_gene_id") %>%
        as_tibble() %>%
        mutate(gene_name = idmap$gene(ensembl_gene_id, to="external_gene_name")) %>%
        select(ensembl_gene_id, gene_name, everything()) %>%
        arrange(padj, pvalue)

    sets = lapply(sets, util$test_gsets, genes=genes)
    c(list(genes=genes), sets, design=cmp)
}

sys$run({
    args = sys$cmd$parse(
        opt('e', 'eset', 'rds', '../data/rnaseq.rds'),
        opt('o', 'outfile', 'xlsx', 'deseq.rds'),
        arg('setfiles', 'rds', arity='*', '../data/genesets/human/MSigDB_Hallmark_2020.rds')
    )

    eset = readRDS(args$infile)

    sets = sapply(args$setfiles, readRDS, simplify=FALSE)
    names(sets) = basename(tools::file_path_sans_ext(names(sets)))

    res = deseq_and_sets(eset, sets)

    writexl::write.xlsx(res, path=args$outfile)
})
