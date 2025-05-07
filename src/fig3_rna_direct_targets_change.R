# Figure 3 - Check expression changes at 'direct' vs 'indirect' targets
# last updated 2025-01-03

# Direct = gene overlaps a Men1/MLL1/H3K4me3 region 
# Indirect = !direct. 
library(tidyverse)
library(tidyHeatmap)
library(plyranges)
library(ComplexHeatmap)
library(plotgardener)

dmso_col <- 'grey40'
sndx_col <- '#d95f02'

# read in peak data
menin <- read_tsv('data/derived_data/cnr/menin_diff_menpeaks_res.tsv') |> as_granges()
mll   <- read_tsv('data/derived_data/cnr/mll_diff_menpeaks_res.tsv') |> as_granges()
k4me3 <- read_tsv('data/derived_data/cnr/k4me3_diff_menpeaks_res.tsv') |> as_granges()
# read in expresion data from drug
rna <- read_tsv('data/derived_data/rna/d4_sndx_all_genes.tsv') |> as_granges()
seqlevelsStyle(rna) <- 'UCSC'
# get set of peaks that are down in menin/mll/h3k4me3
all_down <- menin |> filter(padj < 0.05) |> 
  join_overlap_inner(mll |> filter(padj < 0.05), suffix = c('.menin', '.mll')) |> 
  join_overlap_inner(k4me3 |> filter(padj < 0.05) )

direct_target <- rna |> 
  filter_by_overlaps(all_down) |> 
  mutate(group = 'Direct')

indirect_target <- rna |> 
  filter_by_non_overlaps(all_down) |> 
  mutate(group = 'Indirect')


all_targets <- bind_ranges(direct_target, indirect_target)
targ_plt <- all_targets |> 
  as_tibble() |> 
  filter(padj < 0.05) |> 
  mutate(group = factor(group, levels = c('Indirect', 'Direct'))) |> 
  ggplot(aes(x = group, y = log2FoldChange, fill = group)) + 
  geom_boxplot() + 
  theme_minimal() + 
  geom_hline(aes(yintercept = 0), linetype = 'dashed', color = 'grey30')+
  xlab('') +
  ylab('Log2 Fold Change SNDX-5613/DMSO') + 
  scale_fill_manual(values = c(dmso_col, sndx_col)) + 
  theme(legend.position = 'none',
        axis.text = element_text(size = 8, family = 'sans'), 
        axis.title = element_text(size = 9, family = 'sans') )

tab <- all_targets |> 
  as_tibble() |> 
  filter(padj < 0.05) |> 
  count(group)

pval <- wilcox.test(log2FoldChange ~ group, data = all_targets)$p.value

pdf('plots/cnr/fig3_rna_fc_at_menintargets.pdf', width = 8.5, height = 11)
pageCreate(height = 11, width = 8.5, showGuide = F)
plotGG(targ_plt, x = 0.5, y = 0.5, width = 2, height=3)
plotText(print(as.data.frame(tab)), x = 2.5, y= 4)
plotText(paste0('P-value: ', print(pval)), x = 4, y = 0.5) 
dev.off()

write_tsv(as.data.frame(all_targets), file = 'data/derived_data/rna/rna_targets_at_men1_targets.tsv')

