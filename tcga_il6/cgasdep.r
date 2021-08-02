library(dplyr)
library(ggplot2)
library(depmap)
library("ExperimentHub")
sys = import('sys')
gset = import('genesets')
seq = import('seq')

long2wide = function(df, field, meta) {
    wide = df %>%
        group_by(cell_line, gene_name) %>%
            summarize(value = mean(!! rlang::sym(field))) %>%
        ungroup() %>%
        tidyr::spread("gene_name", value)
    mat = t(data.matrix(wide[-1]))
    colnames(mat) = wide$cell_line
    mat
}

do_test = function(meta, mat, dmat, mgenes=rownames(mat), dgenes=rownames(dmat), quartile=FALSE) {
    test_one = function(mcol, dcol) {
        dset = meta %>%
            inner_join(tibble(cell_line=colnames(mat), x=mat[mcol,])) %>%
            inner_join(tibble(cell_line=colnames(dmat), y=dmat[dcol,]))

        if (quartile) {
            val = rep(NA, nrow(dset))
            val[dset$x < quantile(dset$x, 0.25)] = 0
            val[dset$x > quantile(dset$x, 0.75)] = 1
            dset$x = val
        }

        lm(y ~ x, data=dset) %>% broom::tidy() %>%
            filter(term == "x") %>%
            select(-term)
    }

    mcols = intersect(mgenes, rownames(mat))
    dcols = intersect(dgenes, rownames(dmat))
    tidyr::crossing(tibble(mcol=mcols), tibble(dcol=dcols)) %>%
        rowwise() %>%
            mutate(res = list(test_one(mcol, dcol))) %>%
        ungroup() %>%
        tidyr::unnest() %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)
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
    args = sys$cmd$parse(
        opt('t', 'tissue', 'pan|BRCA', 'BRCA'),
        opt('d', 'dep', 'rnai|ko', 'ko'),
        opt('s', 'set', 'gene set identifier', 'CIN'),
        opt('o', 'outfile', 'rds', 'cgasdep/MSigDB_Hallmark_2020/BRCA-ko.rds')
    )

    meta = depmap::depmap_metadata()
    dmat = switch(args$dep, rnai=depmap::depmap_rnai(), ko=depmap::depmap_crispr()) %>%
        long2wide("dependency", meta)
    gex_mat = depmap::depmap_TPM() %>% long2wide("rna_expression", meta)
    # RPPA: STAT3_pY705 STAT3_Caution NF-kB-p65_pS536_Caution

    if (args$tissue == "BRCA") {
        brca_lines = meta %>% filter(primary_disease == "Breast Cancer") %>% pull(cell_line)
        gex_mat = gex_mat[,colnames(gex_mat) %in% brca_lines]
    }

    smat = switch(args$set,
        genes = gex_mat,
        CIN = cin(meta, gex_mat),
        GSVA::gsva(gex_mat, gset$get_human(args$set))
    )

    dep_genes = c("CGAS", "MB21D1", "IL6", "IL6R", "IL6ST", "RELB", "RELA")

    res = do_test(meta, smat, dmat, dgenes=dep_genes)
    saveRDS(res, file=args$outfile)

    # usable findings:
    #  NC NFkB reg pos more dependent on IL6ST pan-can 0.01 (0.1 pan-cov; reg 0.09 BRCA; reg pos 0.254 BRCA & 0.07 quartile)
    #    HET70 BRCA dep on RELB 0.03 (quartile 0.04); pancan on RELA (0.007) CGAS (0.02)
    #  BRCA IL6R expr dep on CGAS 0.052 (IL6ST 0.3) quartile: 0.07/0.07
    #  pancan Ploidy<RELB 0.03, AneupS<RELB 0.04
})
