library(dplyr)
library(ggplot2)
theme_set(cowplot::theme_cowplot())
io = import('io')
st = import('stats')
plt = import('plot')
enr = import('tools/enrichr')
idmap = import('process/idmap')

# response signatures: Stat1 ChIP (enrichr), Ifn response
expr = io$load("../expr_diff/eset_Mad2PB.RData")
#genes = c("Stat1", "Pias1", "Mb21d1", "Ifng", "Ifngr1", "Trp53", "Wrap53")
cis_genes = c("Stat1", "Pias1", "Ifng", "Trp53", "Nfkb1", "Nfkbib", "Nfkbiz",
              "Nfkbil1", "Nfkbia", "Tbk1", "Mb21d1", "Notch1", "Xrcc6", "Pten",
              "Tlr11", "Tlr4", "Tlr13", "Rel", "Reln")
expr_genes = c(cis_genes, "Ets1", "Erg", "Stat3")
ee = io$load("../expr_diff/eset_Mad2PB.RData")$vs
stat1 = ee["Stat1",]
rest = ee[-which(rownames(ee) == "Stat1"),]
cc = data.frame(
    gene = rownames(rest),
    cor = narray::map(rest, along=2, function(x) cor(x, stat1))
) %>% arrange(-cor)

encc = io$load("../data/genesets/mouse/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.RData")
chea = io$load("../data/genesets/mouse/ChEA_2016.RData")
go = io$load("../data/genesets/mouse/GO_Biological_Process_2018.RData")
hm = io$load("../data/genesets/mouse/CH.HALLMARK.RData")
mile = io$load("../data/genesets/mouse/MILE_regulons.RData")
sets = c(go["regulation of complement activation (GO:0030449)"],
         chea["STAT1_20625510_ChIP-Seq_HELA_Human"],
         chea["STAT1_17558387_ChIP-Seq_HELA_Human"],
         encc["STAT3_CHEA"],
         encc["STAT3_ENCODE"],
         hm["HALLMARK_INTERFERON_GAMMA_RESPONSE"],
         STAT1_complement = list(intersect(
            go[["regulation of complement activation (GO:0030449)"]],
            chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]])),
         STAT1_mhc = list(intersect(
            go[["antigen receptor-mediated signaling pathway (GO:0050851)"]],
            chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]])),
         STAT1_apop = list(intersect(
            go[["apoptotic process (GO:0006915)"]],
            chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]])),
         STAT1_cor = list(intersect(
            head(cc$gene, 1000),
            chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]]))
)
sets$HALLMARK_INTERFERON_GAMMA_RESPONSE = setdiff(sets$HALLMARK_INTERFERON_GAMMA_RESPONSE,
                                                  chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]])
# + MHC, IFN(g), apop(?) -> do GO enrichment, then choose (or: clust?)
gsva = GSVA::gsva(expr$vs, sets)
scores = rbind(expr$vs[expr_genes,], gsva)

# expression changes with insertions (in different subtypes)
cis = io$load("../cis_analysis/poisson.RData")$samples %>%
    filter(external_gene_name %in% cis_genes, sample %in% colnames(expr$vs)) %>%
    narray::construct(n_ins ~ sample + external_gene_name, data=., fill=0) %>%
    apply(1, function(x) paste(names(x)[x != 0], collapse="+"))
annot = as.data.frame(SummarizedExperiment::colData(expr$eset))
annot$ins = "none"
annot$ins[match(names(cis), annot$sample)] = cis
annot$ins = relevel(factor(annot$ins), "none")

both = cbind(annot, t(scores)) %>%
    filter(type != "unknown") %>%
    mutate(Stat1_act = STAT1_cor > 0)

p = ggplot(both, aes(x=Stat1, y=STAT1_cor)) +
    geom_point(aes(size=aneuploidy, color=type, shape=ins)) +
    geom_text_repel(aes(label=sample)) +
    geom_hline(yintercept=0, linetype="dashed")

x = t.test(both %>% filter(type=="Other", Stat1_act) %>% pull(aneuploidy),
           both %>% filter(type=="Other", !Stat1_act) %>% pull(aneuploidy)) %>% broom::tidy()

b2 = both %>%
    mutate(status = case_when(
        type == "T-cell" & !Stat1_act ~ "-",
        type == "Myeloid" & Stat1_act ~ "+",
        type == "Other" & Stat1_act ~ "+",
        type == "Other" & !Stat1_act ~ "-"
    ), status = factor(status, levels=c("+", "-"))) %>% filter(!is.na(status))
p2 = ggplot(b2, aes(x=status, y=aneuploidy, color=type)) +
    geom_boxplot() +
    ggbeeswarm::geom_quasirandom(aes(shape=ins), size=5) +
    facet_wrap(~ type, scales="free_x", nrow=1)

# add cor plots for stat1 act <> stat1 subsets (complement, mhc, apop)
# + <> aneup
# pt size: stat1 expr?
cors = both %>% 
    tidyr::gather("subs", "expr", STAT1_complement:STAT1_apop) %>%
    dplyr::rename(IFNg_response = HALLMARK_INTERFERON_GAMMA_RESPONSE)
#    tidyr::gather("subs", "expr", aneuploidy, STAT1_complement:STAT1_apop)

stats = split(cors, cors$subs) %>%
    lapply(function(x) broom::tidy(lm(expr ~ STAT1_cor, data=x))) %>%
    dplyr::bind_rows(.id="with") %>%
    filter(term == "STAT1_cor") %>%
    select(-term)
stats2 = cors %>%
    filter(type == "Other") %>%
    split(.$subs) %>%
    lapply(function(x) broom::tidy(lm(expr ~ STAT1_cor, data=x))) %>%
    dplyr::bind_rows(.id="with") %>%
    filter(term == "STAT1_cor") %>%
    select(-term)
stats3 = broom::tidy(lm(aneuploidy ~ STAT1_cor, data=cors)) %>%
    filter(term == "STAT1_cor") %>%
    select(-term)
stats4 = broom::tidy(lm(aneuploidy ~ STAT1_cor, data=filter(cors, type=="Other"))) %>%
    filter(term == "STAT1_cor") %>%
    select(-term)

both$ins[both$ins == "none"] = NA
cors$High_Erg = cors$Erg > 6 & cors$Ets1 < 10 & cors$type == "Other"
p1 = ggplot(cors, aes(x=STAT1_cor, y=aneuploidy)) +
    geom_vline(xintercept=0, linetype="dashed") +
    geom_boxplot(aes(group=Stat1_act), outlier.shape=NA, color="grey", varwidth=TRUE) +
    geom_boxplot(aes(group=Stat1_act, color=type), outlier.shape=NA, color="grey", varwidth=TRUE) +
    geom_point(aes(fill=type, size=IFNg_response, shape=High_Erg)) +
    ggrepel::geom_text_repel(data=both, aes(label=ins), segment.alpha=0.7, min.segment.length=0) +
    guides(color = guide_legend(title="Cancer type")) +
    labs(x = "Stat1 activity (inferred)",
         y = "Aneuploidy") +
    scale_shape_manual(values=c(21,22))
cors2 = cors %>%
#    filter(subs != "STAT1_mhc") %>%
    mutate(subs = factor(subs))
#levels(cors2$subs) = c("Apoptotic process\n(GO:0006915)",
#                       "Regulation of complement\nactivation (GO:0030449)")
p2 = ggplot(cors2, aes(x=STAT1_cor, y=expr, color=type)) +
    geom_point(aes(size=aneuploidy), alpha=0.8) +
    geom_smooth(method='lm', color="black", se=FALSE) +
    geom_smooth(method='lm', aes(color=type), se=FALSE) +
    coord_fixed() +
    facet_wrap(~ subs) +
    guides(color = guide_legend(title="Cancer type"),
           size = guide_legend(title="Aneuploidy")) +
    labs(x = "Stat1 activity (inferred)",
         y = "GO signature score")

library(gridExtra)
pdf("stat1_tf_activity.pdf", 10, 8)
print(p1)
print(p2)
grid.arrange(grid::textGrob("Aneuploidy Stat1 low vs high, all types"), tableGrob(stats3))
grid.arrange(grid::textGrob("Aneuploidy Stat1 low vs high, B-like only"), tableGrob(stats4))
grid.arrange(grid::textGrob("Stat1+GO line fit, all types"), tableGrob(stats))
grid.arrange(grid::textGrob("Stat1+GO line fit, B-like only"), tableGrob(stats2))
dev.off()
