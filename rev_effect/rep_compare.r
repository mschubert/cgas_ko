library(dplyr)
library(DESeq2)
sys = import('sys')
util = import('./util')

test_rev = function(eset, cond) {
    eset2 = eset[, colData(eset)$cond == cond]
    print(colData(eset2))
    design(eset2) = ~  replicate:rev
    res = DESeq(eset2)
    re = list(rep1 = util$extract_result(res, "replicate1.rev"),
         rep2 = util$extract_result(res, "replicate2.rev"))
}

args = sys$cmd$parse(
    opt('e', 'eset', 'rds', '../data/rnaseq.rds'),
    opt('t', 'time', 'hours', '48'),
    opt('o', 'outfile', 'rds', 'rep_compare.rds'),
    arg('setfiles', 'rds', arity='*', '../data/genesets/human/MSigDB_Hallmark_2020.rds')
)

sets = sapply(args$setfiles, readRDS, simplify=FALSE)
names(sets) = basename(tools::file_path_sans_ext(names(sets)))

eset = readRDS(args$eset)
eset = eset[, colData(eset)$time == args$time & colData(eset)$treatment != "ifna" &
              ! grepl("il6Ab", colData(eset)$treatment)] # only one rep
colData(eset)$treatment = sub("none", "dmso", colData(eset)$treatment)
colData(eset)$cond = sub("\\+?(rev|dmso)", "",
                         paste(colData(eset)$genotype, colData(eset)$treatment, sep="+"))
colData(eset)$cond = relevel(factor(colData(eset)$cond), "wt")
colData(eset)$rev = ifelse(grepl("rev", colData(eset)$treatment), 1, 0)
colData(eset)$replicate = factor(colData(eset)$replicate)

cond = levels(colData(eset)$cond)
res = sapply(cond, test_rev, eset=eset, simplify=FALSE) %>%
    lapply(. %>% bind_rows(.id="rep")) %>%
    bind_rows(.id="cond") %>%
    group_by(cond, rep) %>%
    tidyr::nest(.key="genes")

for (ns in names(sets))
    res[[ns]] = lapply(res$genes, util$test_gsets, sets=sets[[ns]])

saveRDS(res, file=args$outfile)