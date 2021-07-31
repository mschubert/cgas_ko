library(ggplot2)
library(dplyr)
sys = import('sys')
plt = import('plot')

args = sys$cmd$parse(
    opt('i', 'infile', 'rds', 'cgasdep/CIN/BRCA-ko.rds'),
    opt('p', 'plotfile', 'pdf', 'cgasdep/CIN/BRCA-ko.pdf')
)

res = readRDS(args$infile) %>%
    split(.$dep)

plots = lapply(res, plt$volcano, y="p.value")

pdf(args$plotfile, 10, 8)
for (i in seq_along(plots))
    print(plots[[i]] + ggtitle(names(plots)[i]))
dev.off()
