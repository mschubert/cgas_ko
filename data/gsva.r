library(dplyr)
sys = import('sys')

args = sys$cmd$parse(
    opt('e', 'eset', 'rds', 'rnaseq.rds'),
    opt('s', 'setfile', 'rds', 'genesets/human/MSigDB_Hallmark_2020.rds'),
    opt('o', 'outfile', 'rds', 'gsva/MSigDB_Hallmark_2020.rds')
)
