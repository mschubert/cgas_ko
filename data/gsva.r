library(dplyr)
sys = import('sys')
idmap = import('process/idmap')

args = sys$cmd$parse(
    opt('e', 'eset', 'rds', 'rnaseq.rds'),
    opt('s', 'setfile', 'rds', 'genesets/human/MSigDB_Hallmark_2020.rds'),
    opt('o', 'outfile', 'rds', 'gsva/MSigDB_Hallmark_2020.rds')
)

eset = readRDS(args$eset)
sets = readRDS(args$setfile)

vst = SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(eset))
rownames(vst) = idmap$gene(rownames(vst), to="hgnc_symbol")

scores = GSVA::gsva(expr=vst, gset.idx.list=sets, parallel.sz=0, min.sz=3)

saveRDS(t(scores), file=args$outfile)
