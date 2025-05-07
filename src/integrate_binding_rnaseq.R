# combined heatmaps for menin/mll/k4me3/rnaseq of nearest genes
# figure 3

library(tidyverse)
library(tidyHeatmap)
library(plyranges)
library(ComplexHeatmap)

menin <- read_tsv('data/derived_data/cnr/menin_diff_menpeaks_res.tsv') |> as_granges()
mll   <- read_tsv('data/derived_data/cnr/mll_diff_menpeaks_res.tsv') |> as_granges()
k4me3 <- read_tsv('data/derived_data/cnr/k4me3_diff_menpeaks_res.tsv') |> as_granges()

rna <- read_tsv('data/derived_data/rna/d4_sndx_all_genes.tsv') |> as_granges()
seqlevelsStyle(rna) <- 'UCSC'

cd <- menin |> filter(padj < 0.05) |> 
  join_overlap_inner(mll , suffix = c('.menin', '.mll')) |> 
  join_overlap_inner(k4me3) |> 
  join_nearest(rna, suffix = c('.k4me3', '.rna')) 

all_down <- menin |> filter(padj < 0.05) |> 
  join_overlap_inner(mll , suffix = c('.menin', '.mll')) |> 
  join_overlap_inner(k4me3)

pal <- circlize::colorRamp2(c(-3,0,3), c('blue', 'black', 'yellow'))


  
notrna <- cd |> 
  as_tibble() |> 
  dplyr::select(log2FoldChange.menin, log2FoldChange.mll, log2FoldChange.k4me3, log2FoldChange.rna, gene_name) |> 
  pivot_longer(names_to = 'sample', values_to = 'score', -gene_name) |> 
  filter(sample != 'log2FoldChange.rna') 
  group_by(gene_name, sample) |> 
  mutate(score = mean(score, na.rm = T)) 
  
rna_d <- cd |> 
  as_tibble() |> 
  dplyr::select(log2FoldChange.menin, log2FoldChange.mll, log2FoldChange.k4me3, log2FoldChange.rna, gene_name) |> 
  pivot_longer(names_to = 'sample', values_to = 'score', -gene_name) |> 
  filter(sample == 'log2FoldChange.rna') |> 
  group_by(gene_name, sample) |> 
  mutate(score = mean(score, na.rm = T)) 
  


mini <- cd |> filter(padj.rna < 0.05) |> filter(padj.menin < 0.01) |> filter(padj.mll < 0.05) |> filter(padj.k4me3 < 0.1)
mini <- mini |> filter(log2FoldChange.rna < -0.6)
fs <- mini |>  as_tibble() |> 
  dplyr::select(log2FoldChange.menin, log2FoldChange.mll, log2FoldChange.k4me3, log2FoldChange.rna, gene_name) |> 
  pivot_longer(names_to = 'sample', values_to = 'score', -gene_name) |> 
  group_by(gene_name, sample) |> 
  mutate(score = mean(score, na.rm = T))  |> unique() |> 
  pivot_wider(names_from = sample, values_from = score)

mat_rn <- fs$gene_name
mat <- as.matrix(fs[,2:ncol(fs)])
rownames(mat) <- mat_rn
Heatmap(mat, show_row_names = F, col = pal)

fs_best <- fs$gene_name

cd |> filter(gene_name %in% fs_best) |> 
  as_tibble() |> 
  ggplot(aes(x = log2FoldChange.rna, y = -log10(padj.rna))) + geom_point()

direct_target <- rna |> 
  filter_by_overlaps(all_down) |> 
  mutate(group = 'direct')
  
indirect_target <- rna |> 
  filter_by_non_overlaps(all_down) |> 
  mutate(group = 'indirect')


all_targets <- bind_ranges(direct_target, indirect_target)
all_targets |> 
  as_tibble() |> 
  filter(padj < 0.05) |> 
  ggplot(aes(x = group, y = log2FoldChange)) + 
  geom_boxplot() 


all_targets |> 
  as_tibble() |> 
  filter(padj < 0.05) |> 
  count(group)
