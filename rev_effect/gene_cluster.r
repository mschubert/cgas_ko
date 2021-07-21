library(dplyr)
library(DESeq2)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')

args = sys$cmd$parse(
#    opt('s', 'short', 'rds', 'rep_compare-24h.rds'),
    opt('l', 'long', 'rds', 'rep_compare-48h.rds'),
    opt('p', 'plotfile', 'pdf', 'inflamm_resp_overview.pdf')
)

dset = readRDS("deseq-48h.rds")

mat = setNames(dset$genes, dset$cond) %>%
    bind_rows(.id="cond") %>%
    mutate(signed_p = sign(stat) * pmax(0, 1 - 1e3*pvalue)) %>%
    narray::construct(signed_p ~ ensembl_gene_id + cond)

mat[is.na(mat)] = 0
mat = mat[!rowSums(mat == 0) == ncol(mat),]
dim(mat)

um1 = uwot::umap(mat, n_components=1)
um2 = uwot::umap(mat, n_components=2)
#um2 = uwot::umap(mat[,c("wt", "cgas", "stat1", "stat3")], n_components=2)
um1t = uwot::umap(t(mat), n_components=1, n_neighbors=5)

df = tibble(ensembl_gene_id = rownames(mat),
            umap1 = um2[,1],
            umap2 = um2[,2],
            cl = factor(igraph::cluster_louvain(scran::buildSNNGraph(t(mat), k=10))$membership))

ggplot(df, aes(x=umap1, y=umap2)) +
    geom_point(aes(color=cl))

df2 = reshape2::melt(mat) %>%
    as_tibble() %>%
    inner_join(df) %>%
    mutate(ensembl_gene_id = factor(ensembl_gene_id, levels=rownames(mat)[rank(um1, ties="random")]),
           cond = factor(cond, levels=colnames(mat)[rank(um1t, ties="random")]))
ggplot(df2, aes(x=cond, y=ensembl_gene_id)) +
   geom_tile(aes(fill=value)) +
   scale_fill_distiller(palette="RdBu", direction=1)
ggplot(df2, aes(x=cond, y=value, group=ensembl_gene_id)) +
    geom_line() +
    facet_wrap(~cl)
