library(dplyr)
library(ggplot2)
library(depmap)
library("ExperimentHub")
gset = import('genesets')
seq = import('seq')

long2wide = function(meta, df, field) {
    wide = df %>%
        group_by(cell_line, gene_name) %>%
            summarize(value = mean(!! rlang::sym(field))) %>%
        ungroup() %>%
        tidyr::spread("gene_name", value)
    mat = t(data.matrix(wide[-1]))
    colnames(mat) = wide$cell_line
    mat
}

test_one = function(mat, mcol, dmat, dcol) {
    dset = meta %>%
        inner_join(tibble(cell_line=colnames(mat), gene=mat[mcol,])) %>%
        inner_join(tibble(cell_line=colnames(dmat), dep=dmat[dcol,]))

    lm(dep ~ gene, data=dset) %>% broom::tidy() %>%
        filter(term != "(Intercept)") %>%
        arrange(p.value)
}

sys$run({
    eh = ExperimentHub()
    query(eh, "depmap")

    rnai = depmap::depmap_rnai()
    ko = depmap::depmap_crispr()
    gex = depmap::depmap_TPM()
    aneup = readr::read_csv("aneuploidy_scores.csv") %>%
        dplyr::rename(depmap_id = DepMap_ID)

    aneup2 = readr::read_csv("CCLE_segment_cn.csv") %>%
        mutate(width = End - Start, ploidy = Segment_Mean*2) %>%
        seq$aneuploidy(seqnames="Chromosome", sample="DepMap_ID") %>%
        select(depmap_id=DepMap_ID, aneup=aneuploidy)

    hms = gset$get_human("MSigDB_Hallmark_2020")
    gex_mat = long2wide(gex, "rna_expression")
    scores = GSVA::gsva(gex_mat, hms) #, parallel.sz=30)
    sdf = as_tibble(t(scores)) %>%
        mutate(cell_line = colnames(scores))

    rnaimat = long2wide(rnai, "dependency")
    komat = long2wide(ko, "dependency")

    # RPPA: STAT3_pY705 STAT3_Caution NF-kB-p65_pS536_Caution

    meta = depmap::depmap_metadata() %>%
        left_join(aneup) %>%
        left_join(aneup2)

    res = tidyr::crossing(tibble(set=rownames(scores)), tibble(dep=c("CGAS", "IL6R", "IL6ST"))) %>%
        rowwise() %>%
            mutate(res = list(test_one(scores, set, komat, dep))) %>%
        ungroup() %>%
        tidyr::unnest() %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)

    ggplot(df, aes(x=rna_expression, y=dependency)) +
        geom_point(aes(size=`Aneuploidy score`, color=lineage_sub_subtype)) +
        ggrepel::geom_text_repel(aes(label=cell_line), size=3)
})
