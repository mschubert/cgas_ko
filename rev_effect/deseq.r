library(dplyr)
library(DESeq2)
sys = import('sys')
util = import('./util')

test_rev = function(eset, cond) {
    eset2 = eset[, colData(eset)$cond == cond]
    print(colData(eset2))
    design(eset2) = ~  rev
#    design(eset2) = ~ replicate + rev
    res = DESeq(eset2) %>% util$extract_result("rev")
}

args = sys$cmd$parse(
    opt('e', 'eset', 'rds', '../data/rnaseq.rds'),
    opt('t', 'time', '24|48', '24'),
    opt('o', 'outfile', 'rds', 'deseq.rds'),
    arg('setfiles', 'rds', arity='*', '../data/genesets/human/MSigDB_Hallmark_2020.rds')
)

sets = sapply(args$setfiles, readRDS, simplify=FALSE)
names(sets) = basename(tools::file_path_sans_ext(names(sets)))

eset = readRDS(args$eset)
eset = eset[, colData(eset)$time == args$time & colData(eset)$treatment != "ifna"]
colData(eset)$treatment = sub("none", "dmso", colData(eset)$treatment)
colData(eset)$cond = sub("\\+?(rev|dmso)", "",
                         paste(colData(eset)$genotype, colData(eset)$treatment, sep="+"))
colData(eset)$cond = relevel(factor(colData(eset)$cond), "wt")
colData(eset)$rev = ifelse(grepl("rev", colData(eset)$treatment), 1, 0)

cond = levels(eset$cond)
cond = cond[table(eset$cond) == 4] # cgas[24h], il6Ab[48h] only one rep
res = sapply(cond, test_rev, eset=eset, simplify=FALSE)
res = tibble(cond = cond, genes = unname(res))

for (ns in names(sets))
    res[[ns]] = lapply(res$genes, util$test_gsets, sets=sets[[ns]])

saveRDS(res, file=args$outfile)
