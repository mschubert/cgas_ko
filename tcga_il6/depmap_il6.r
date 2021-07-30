library(dplyr)
library(ggplot2)
library(depmap)
library("ExperimentHub")
gset = import('genesets')

eh = ExperimentHub()
query(eh, "depmap")

meta = depmap::depmap_metadata()
rnai = depmap::depmap_rnai()
ko = depmap::depmap_crispr()
gex = depmap::depmap_TPM()
aneup = readr::read_csv("aneuploidy_scores.csv") %>%
    dplyr::rename(depmap_id = DepMap_ID)

hms = gset$get_human("MSigDB_Hallmark_2020")
gex_wide = gex %>%
    group_by(cell_line, gene_name) %>%
        summarize(rna_expression = mean(rna_expression)) %>%
    ungroup() %>%
    tidyr::spread("gene_name", "rna_expression")
gex_mat = data.matrix(gex_wide[-1])
rownames(gex_mat) = gex_wide$cell_line
scores = GSVA::gsva(t(gex_mat), hms) #, parallel.sz=30)
sdf = as_tibble(t(scores)) %>%
    mutate(cell_line = colnames(scores))

# RPPA: STAT3_pY705 STAT3_Caution NF-kB-p65_pS536_Caution

#df = inner_join(rnai %>% filter(gene_name == "MB21D1") %>% select(cell_line, dependency),
df = inner_join(ko %>% filter(gene_name == "IL6R") %>% select(cell_line, dependency),
                gex %>% filter(gene_name == "BUB1") %>% select(cell_line, rna_expression)) %>%
    inner_join(meta) %>% #filter(primary_disease == "Breast Cancer")) %>%
    inner_join(aneup) %>%
    left_join(sdf)

lm(dependency ~ rna_expression, data=df) %>% broom::tidy() %>%
    arrange(p.value)

ggplot(df, aes(x=rna_expression, y=dependency)) +
    geom_point(aes(size=`Aneuploidy score`))
