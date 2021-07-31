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

test_one = function(mat, mcol, dmat, dcol, quartile=FALSE) {
    dset = meta %>%
        inner_join(tibble(cell_line=colnames(mat), gene=mat[mcol,])) %>%
        inner_join(tibble(cell_line=colnames(dmat), dep=dmat[dcol,]))

    if (quartile) {
        val = rep(NA, nrow(dset))
        val[dset$gene < quantile(dset$gene, 0.25)] = 0
        val[dset$gene > quantile(dset$gene, 0.75)] = 1
        dset$gene = val
    }

    lm(dep ~ gene, data=dset) %>% broom::tidy() %>%
        filter(term != "(Intercept)") %>%
        arrange(p.value)
}

cin = function(meta, gex_mat) {
    seg_aneup = readr::read_csv("CCLE_segment_cn.csv") %>%
        mutate(width = End - Start, ploidy = Segment_Mean*2) %>%
        seq$aneuploidy(seqnames="Chromosome", sample="DepMap_ID") %>%
        select(depmap_id=DepMap_ID, seg_aneup=aneuploidy)
    aneup = readr::read_csv("aneuploidy_scores.csv") %>%
        dplyr::rename(depmap_id = DepMap_ID) %>%
        full_join(seg_aneup, by="depmap_id") %>%
        left_join(meta %>% select(depmap_id, cell_line)) %>%
        filter(!is.na(cell_line))

    sets = gset$get_human("CIN")
    scores = GSVA::gsva(gex_mat, sets)

    narray::intersect(aneup$cell_line, scores, along=2)
    amat = aneup %>% select(-cell_line, -depmap_id) %>% data.matrix() %>% t()
    rbind(amat, scores)
}

sys$run({
    eh = ExperimentHub()
    query(eh, "depmap")

    rnai = depmap::depmap_rnai()
    ko = depmap::depmap_crispr()
    gex = depmap::depmap_TPM()

    hms = gset$get_human("CIN")
    gex_mat = long2wide(gex, "rna_expression")
    scores = GSVA::gsva(gex_mat, hms) #, parallel.sz=30)

    rnaimat = long2wide(rnai, "dependency")
    komat = long2wide(ko, "dependency")

    # RPPA: STAT3_pY705 STAT3_Caution NF-kB-p65_pS536_Caution

    meta = depmap::depmap_metadata() %>%
        left_join(aneup) %>%
        left_join(aneup2)

    cin_mat = cin(meta, gex_mat)
    hm_mat = GSVA::gsva(gex_mat, gset$get_human("MSigDB_Hallmark_2021"))

    # tryouts
#    scores = gex_mat[c("CGAS", "IL6", "IL6R", "IL6ST", "MAD2L1", "BUB1", "NFKB1", "RELA", "RELB"),]
#    meta = meta %>% filter(primary_disease == "Breast Cancer")

    res = tidyr::crossing(tibble(set=rownames(scores)), tibble(dep=c("CGAS", "IL6", "IL6R", "IL6ST", "RELB", "RELA"))) %>%
        rowwise() %>%
            mutate(res = list(test_one(scores, set, komat, dep, quartile=F))) %>%
        ungroup() %>%
        tidyr::unnest() %>%
    filter(term == "gene") %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)

    # usable findings:
    #  NC NFkB reg pos more dependent on IL6ST pan-can 0.01 (0.1 pan-cov; reg 0.09 BRCA; reg pos 0.254 BRCA & 0.07 quartile)
    #    HET70 BRCA dep on RELB 0.03 (quartile 0.04); pancan on RELA (0.007) CGAS (0.02)
    #  BRCA IL6R expr dep on CGAS 0.052 (IL6ST 0.3) quartile: 0.07/0.07

    ggplot(df, aes(x=rna_expression, y=dependency)) +
        geom_point(aes(size=`Aneuploidy score`, color=lineage_sub_subtype)) +
        ggrepel::geom_text_repel(aes(label=cell_line), size=3)
})
