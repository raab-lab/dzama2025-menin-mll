# plot menin/mll fc vs h3k4me3 fc color by expression (dmso v sndx) or (wt v mut)
# author: Jesse Raab
# Date: 2025-01-21 # cleaning up final versions

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(csaw)
library(ChIPpeakAnno)
library(plyranges)
library(ggplot2)
library(tidyverse)
library(plotgardener)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
chrs <- paste0("chr", c(1:22, "X", "Y"))
fdr_cut <- 0.01
lfc_cut <- 0


load('data/processed_data/cnr/menin_peaks_diff_all.Rda')
rna <- read_tsv('data/processed_data/rna/hlf_sndx_v_dmso.tsv') |> as_granges()
seqlevelsStyle(rna) <- 'UCSC'

table(k4me3_res$padj< 0.05, k4me3_res$log2FoldChange < 0)
table(men_res$padj < 0.05, men_res$log2FoldChange < 0)
table(mll_res$padj < 0.05, mll_res$log2FoldChange < 0)



# Plot differential occupancy MA plots
k4me3_res |> 
  as_tibble() |> 
  ggplot(aes( x= log2(baseMean), y = log2FoldChange, color = padj < 0.05) ) + geom_point()

men_res |> 
  as_tibble() |> 
  ggplot(aes( x= log2(baseMean), y = log2FoldChange, color = padj < 0.05) ) + geom_point()


mll_res |> 
  as_tibble() |> 
  ggplot(aes( x= log2(baseMean), y = log2FoldChange, color = padj < 0.05) ) + geom_point()


k4me3_plot <- k4me3_res |> 
  as_tibble() |> 
  ggplot(aes( x= log2FoldChange, y = -log10(padj), color = padj < 0.05) ) + geom_point() + 
  theme_bw() + 
  scale_color_manual(values = c( 'grey20', 'red2')) + 
  theme( panel.grid.minor =element_blank()) + 
  geom_vline(aes(xintercept = 0), linetype = 'dashed', color= 'grey10') + 
  ggtitle('H3K4me3 SNDX/DMSO')

mll_plot<- mll_res |> 
  as_tibble() |> 
  ggplot(aes( x= log2FoldChange, y = -log10(padj), color = padj < 0.05) ) + geom_point() + 
  theme_bw() + 
  scale_color_manual(values = c( 'grey20', 'red2')) + 
  theme(panel.grid.minor =element_blank()) + 
  geom_vline(aes(xintercept = 0), linetype = 'dashed', color= 'grey10') +
  ggtitle('MLL1 SNDX/DMSO')

men_plot <- men_res |> 
  as_tibble() |> 
  ggplot(aes( x= log2FoldChange, y = -log10(padj), color = padj < 0.05) ) + geom_point() + 
  theme_bw() + 
  scale_color_manual(values = c( 'grey20', 'red2')) + 
  theme( panel.grid.minor =element_blank()) + 
  geom_vline(aes(xintercept = 0), linetype = 'dashed', color= 'grey10') + 
  ggtitle('Menin SNDX/DMSO')



pdf('plots/figs/cnr/differential_peak_volcano_2F.pdf', width = 8.5, height = 11)
pageCreate(height = 11, width = 8.5, showGuide = F)
legend <- cowplot::get_legend(k4me3_plot)
plotGG(men_plot + theme(legend.position = 'none'),
       x = 0.5, y =0.5, width = 2.5, height = 2.5)
plotGG(mll_plot + theme(legend.position = 'none')+ ylab(''), 
       x = 3.1, y = 0.5, width = 2.5, height = 2.5)
plotGG(k4me3_plot + theme(legend.position = 'none') +ylab(''), 
       x = 0.5, y = 4, width = 2.5, height = 2.5)

plotGG(legend, x = 0.5, y = 8, width = 1, height = 1, just = c('top', 'left')) 
dev.off()




men_res_anno <- men_res |> join_nearest(rna, suffix = c('.men', '.rna') )


men_res_anno |> join_overlap_inner(k4me3_res) |> 
  
  as_tibble() |> 
  filter(abs(log2FoldChange.rna) > 0.6) |> 
  ggplot(aes(x = log2FoldChange.men, y = log2FoldChange, color = log2FoldChange.rna)) + 
  geom_point() + 
  scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint =0)
