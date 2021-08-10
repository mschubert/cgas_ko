library(dplyr)
library(DESeq2)
.b = import('base')
sys = import('sys')
plt = import('plot')
gset = import('genesets')

#todo: move this to a "scales" util
.reverselog_trans = function(base=exp(1)) {
    scales::trans_new(paste0("log-", format(base)),
                      function(x) -log(x, base),
                      function(x) base^-x,
                      scales::log_breaks(base = base),
                      domain = c(1e-100, Inf))
}
.scientific_10 = function(x) {
    fmt = ifelse(x < 0.01, scales::scientific_format()(x), x)
    parse(text=gsub("1e", "10^", fmt))
}
.simplify = function(df, y, p=0.25) {
    set.seed(123456)
    idx = which(df[[y]] >= .b$minN(df[[y]][df[[y]] > p], 100))
    prob = 1 - df[[y]][idx]+.Machine$double.eps*2
#    prob[df$circle[idx]] = 1
    keep = sample(idx, size=500, replace=FALSE, prob=prob)
    df[[y]][setdiff(idx, keep)] = NA
    df$label[idx] = NA
    df
}

args = sys$cmd$parse(
#    opt('e', 'eset', 'rds', '../comp_list/eset/BT549 WT Reversine vs DMSO.rds'),
    opt('i', 'infile', 'xlsx', '../comp_list/cmp/BT549 WT Reversine vs DMSO.xlsx'),
    opt('p', 'plotfile', 'pdf', 'main_rev-IL6.pdf')
)

is_cyt = gset$get_human("GO_Molecular_Function_2021")[["cytokine activity (GO:0005125)"]]

cgas_indep = readxl::read_xlsx("../comp_list/cmp/BT549 cGAS KO Reversine vs DMSO.xlsx") %>%
    filter(padj < 0.25) %>% pull(label)

res = readxl::read_xlsx(args$infile) %>%
    mutate(
        cyt = label %in% is_cyt,
        use_lab = ifelse(cyt & padj < 0.25, label, NA),
        cgasdep = ! label %in% cgas_indep,
        class = case_when(
            padj > 0.25 ~ "n.s.",
            log2FoldChange > 0 ~ "up",
            log2FoldChange < 0  ~"down"
    )) %>% .simplify("padj")

plt$volcano(res, p=0.2, label_top=10, pos_label_bias=0.5)
breaks_with_thresh = function(...) c(0.25, scales::log_breaks(base=10)(res$padj, 4))

p = ggplot(res, aes(x=log2FoldChange, y=padj)) +
    geom_hline(yintercept=0.25, color="grey", linetype="dashed") +
    geom_vline(xintercept=0, color="#858585") +
    geom_point(aes(size=sqrt(baseMean), color=class, fill=class, alpha=cyt, shape=cgasdep)) +
    ggrepel::geom_text_repel(aes(label=use_lab)) +
    scale_y_continuous(trans = .reverselog_trans(base=10),
                       labels = .scientific_10,
                       breaks = breaks_with_thresh) +
    scale_color_manual(values=c("n.s."="#c8c8c8", "up"="#00971e", "down"="#e10000")) +
    scale_fill_manual(values=c("n.s."="#c8c8c8", "up"="#00971e", "down"="#e10000")) +
    scale_alpha_manual(values=c("TRUE"=0.9, "FALSE"=0.1)) +
    scale_size_continuous(breaks=c(0, 50, 200)) +
    scale_shape_manual(values=c("TRUE"=21, "FALSE"=1)) +
    theme_classic() +
    guides(color=FALSE, fill=FALSE, shape=guide_legend(override.aes=list(fill="black"))) +
    labs(y = "adjusted p-value (FDR)",
         size = "Expression level",
         shape = "cGAS-dependent",
         alpha = "Cytokine")

pdf(args$plotfile, 5.8, 6)
print(p)
dev.off()
