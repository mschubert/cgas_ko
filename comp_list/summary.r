library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')

read_xlsx = function(fpath) {
    sheets = setdiff(readxl::excel_sheets(fpath), "cmp")
    re = sapply(sheets, readxl::read_xlsx, path=fpath, simplify=FALSE)
    do.call(tibble, lapply(re, list))
}

plot_page = function(comparison, ...) {
    args = list(...)
    plots = lapply(args, plt$volcano, label_top=30)
    plt$text(comparison, size=7) / wrap_plots(plots, nrow=1) + plot_layout(heights=c(1,20))
}

args = sys$cmd$parse(
    opt('c', 'config', 'yaml', 'comps.yaml'),
    opt('o', 'outfile', 'rds', 'summary.rds'),
    opt('p', 'plotfile', 'pdf', 'summary.pdf'),
    arg('infiles', '', arity='*', list.files("cmp", "\\.xlsx$", full.names=TRUE))
)

res = lapply(args$infiles, read_xlsx) %>%
    setNames(tools::file_path_sans_ext(basename(args$infiles))) %>%
    bind_rows(.id="comparison") %>%
    mutate(plot = purrr::pmap(., plot_page))

pdf(args$plotfile, 24, 9)
for (i in seq_len(nrow(res)))
    plot(res$plot[[i]])
dev.off()

saveRDS(res %>% select(-plot), file=args$outfile)
