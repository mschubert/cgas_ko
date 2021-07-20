library(dplyr)
library(ggplot2)
sys = import('sys')

args = sys$cmd$parse(
    opt('c', 'config', 'yaml', 'comps.yaml'),
    opt('p', 'plotfile', 'pdf', 'summary.pdf'),
    arg('infiles', '', arity='*', list.files("cmp", "\\.xlsx$", full.names=TRUE))
)
