library(dplyr)
library(DESeq2)
sys = import('sys')
util = import('./util')

test_rev = function(eset, cond) {
    eset2 = eset[, colData(eset)$cond == cond]
    design(eset2) = ~ replicate + rev
    res = DESeq(eset2) %>% util$extract_result("rev")
}

args = sys$cmd$parse(
    opt('e', 'eset', 'rds', '../data/rnaseq.rds'),
    opt('o', 'outfile', 'rds', 'deseq.rds'),
    arg('setfiles', 'rds', arity='*', '../data/genesets/human/CH.HALLMARK.rds')
)

sets = sapply(args$setfiles, readRDS)
names(sets) = basename(tools::file_path_sans_ext(names(sets)))

eset = readRDS(args$eset)
eset = eset[,! colData(eset)$treatment %in% c("none", "ifna", "ifng") & colData(eset)$time == 24]
colData(eset)$treatment = sub("none", "dmso", colData(eset)$treatment)
colData(eset)$cond = sub("\\+?(rev|dmso)", "",
                         paste(colData(eset)$genotype, colData(eset)$treatment, sep="+"))
colData(eset)$cond = relevel(factor(colData(eset)$cond), "wt")
colData(eset)$rev = ifelse(grepl("rev", colData(eset)$treatment), 1, 0)

cond = levels(colData(eset)$cond)
res = sapply(cond, test_rev, eset=eset, simplify=FALSE)
res = tibble(cond = cond, genes = unname(res))

for (ns in names(sets))
    res[[ns]] = lapply(res$genes, util$test_gsets, sets=sets[[ns]])

saveRDS(res, file=args$outfile)
