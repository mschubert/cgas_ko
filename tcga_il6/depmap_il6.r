library(dplyr)
library(ggplot2)
library(depmap)
library("ExperimentHub")

eh = ExperimentHub()
query(eh, "depmap")

meta = depmap::depmap_metadata()
rnai = depmap::depmap_rnai()
ko = depmap::depmap_crispr()
gex = depmap::depmap_TPM()

df = inner_join(rnai %>% filter(gene_name == "MB21D1") %>% select(cell_line, dependency),
                gex %>% filter(gene_name == "CGAS") %>% select(cell_line, rna_expression)) %>%
    inner_join(meta)

lm(dependency ~ primary_disease + rna_expression, data=df) %>% broom::tidy() %>%
    arrange(p.value)
