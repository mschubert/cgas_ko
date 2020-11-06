library(dplyr)
library(DESeq2)
sys = import('sys')
util = import('./util')

args = sys$cmd$parse(
    opt('e', 'eset', 'rds', '../data/rnaseq.rds'),
    opt('o', 'outfile', 'rds', 'deseq.rds'),
    arg('setfiles', 'rds', arity='*', '../data/genesets/human/CH.HALLMARK.rds')
)

sets = sapply(args$setfiles, readRDS)
names(sets) = basename(tools::file_path_sans_ext(names(sets)))

eset = readRDS(args$eset)
colData(eset)$cond = ifelse(grepl("il6", colData(eset)$treatment), "wt+il6",
                            as.character(colData(eset)$genotype))
colData(eset)$cond = relevel(factor(colData(eset)$cond), "wt")
colData(eset)$rev = ifelse(grepl("rev", colData(eset)$treatment), 1, 0)

design(eset) = ~ cond:rev
res = DESeq(eset)
res = resultsNames(res) %>%
    setdiff("Intercept") %>%
    sapply(util$extract_result, res=res, simplify=FALSE)
res = tibble(cond = sub("\\.rev", "", sub("^cond", "", names(res))),
             genes = unname(res))

for (ns in names(sets))
    res[[ns]] = lapply(res$genes, util$test_gsets, sets=sets[[ns]])

saveRDS(res, file=args$outfile)
