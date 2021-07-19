library(dplyr)
library(DESeq2)
sys = import('sys')
util = import('./util')

test_rev = function(eset, cond) {
    message("* condition: ", cond)
    eset2 = eset[, eset$cond == cond]
    eset2$replicate = factor(eset2$replicate) # be sure to drop levels
    design(eset2) = ~ replicate:rev
    res = DESeq(eset2)
    rns = grep("replicate.\\.rev", resultsNames(res), value=TRUE)
    lapply(rns, util$extract_result, res=res) %>%
        setNames(sub("replicate([0-9]+)\\.rev", "\\1", rns))
}

args = sys$cmd$parse(
    opt('e', 'eset', 'rds', '../data/rnaseq.rds'),
    opt('t', 'time', 'hours', '48'),
    opt('o', 'outfile', 'rds', 'rep_compare-48h.rds'),
    arg('setfiles', 'rds', arity='*', '../data/genesets/human/MSigDB_Hallmark_2020.rds')
)

sets = sapply(args$setfiles, readRDS, simplify=FALSE)
names(sets) = basename(tools::file_path_sans_ext(names(sets)))

eset = readRDS(args$eset)
eset = eset[, eset$time == args$time & eset$treatment != "ifna"]
eset$treatment = sub("none", "dmso", eset$treatment)
eset$cond = sub("\\+?(rev|dmso)", "", paste(eset$genotype, eset$treatment, sep="+"))
eset$rev = ifelse(grepl("rev", eset$treatment), 1, 0)

conds = as.data.frame(colData(eset)) %>%
    group_by(cond) %>%
    filter(n_distinct(replicate) > 1) %>%
    pull(cond) %>% unique()

res = sapply(conds, test_rev, eset=eset, simplify=FALSE) %>%
    lapply(. %>% bind_rows(.id="rep")) %>%
    bind_rows(.id="cond") %>%
    group_by(cond, rep) %>%
    tidyr::nest(.key="genes")

for (ns in names(sets))
    res[[ns]] = lapply(res$genes, util$test_gsets, sets=sets[[ns]])

saveRDS(res, file=args$outfile)
