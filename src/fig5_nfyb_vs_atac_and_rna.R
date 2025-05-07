# FIgure 5 plots comparing nfyb atac and rna
# This plots Figure 5 and 5S - scatter/boxplots of atac/rna vs nfyb and distributions
# Author: Jesse Raab
# Date: 2025-02-24

library(tidyverse)
library(DESeq2)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene) 
library(plyranges)
library(plotgardener)
library(grid)
library(gridExtra)
library(ChIPseeker)

################################################################################
fix_category <- function(tib) {
  tib |> 
    mutate(new_cat = case_when(grepl ('Promoter', annotation) ~ 'Promoter', 
                               grepl ('Intron', annotation) ~ 'Intron', 
                               grepl ('Exon', annotation) ~ 'Exon', 
                               .default = as.character(annotation) ) )
  
  
  
}
################################################################################
outdir <- 'data/derived_data/integrated/'
dir.create(outdir, recursive = T)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
load('data/derived_data/cnr/nfyb_analysis/nfyb_analysis_obj.Rda')
nfyb_dar <- read_tsv('data/derived_data/cnr/nfyb_analysis/nfyb_diff_at_atac.tsv')
nfyb_dor <- read_tsv('data/derived_data/cnr/nfyb_analysis/nfyb_diff_at_nfyb_peaks.tsv')
ex <- read_tsv('data/derived_data/rna/d4_sndx_all_genes.tsv') |> as_granges()
seqlevelsStyle(ex) <- 'UCSC'

# all I really want from expression is padj, fold change, and basemean, and gene
ex <- ex |> as_tibble() |> dplyr::select(log2FoldChange, padj, baseMean, gene_name)

#load accessibility regions
load('data/derived_data/atac/atac_de.Rda')
# re-shrink the atac
resultsNames(des)
dar <- lfcShrink(des, type = 'apeglm', coef = 2, format = 'GRanges')

# Combine NFYB, ATAC, RNA datasets
nfyb_dar_anno_all <- nfyb_dor |> 
  as_granges() |> 
  annotatePeak(TxDb = txdb, tssRegion = c(-2000, 200), annoDb='org.Hs.eg.db'  ) 
nfyb_dar_anno_all <- nfyb_dar_anno_all |> as_granges() |> fix_category()
nfyb_dar_anno_all <- nfyb_dar_anno_all |> mutate(group = 
                                                   case_when( log2FoldChange < 0 & padj < 0.1 ~ 'Down', 
                                                              log2FoldChange > 0 & padj < 0.1 ~ 'Up', 
                                                              padj > 0.1 ~ 'Stable') )


# Add ATAC data to above - #2025-02-26 - changed to dor instead of dar regions - so join_overlap_left
nfyb_dar_anno_all <- nfyb_dar_anno_all |> 
  filter(!is.na(padj)) |> 
  join_nearest(dar, suffix = c('.nfyb', '.atac')) 

# Add expression data to NFYB 
nfyb_dar_anno_all <- nfyb_dar_anno_all |> 
  as_tibble() |> 
  left_join(ex |> as_tibble(), by = c('SYMBOL'='gene_name'), relationship = 'many-to-many') 


# minimize the data so we only have columns we car about
nfyb_dar_anno_all <- nfyb_dar_anno_all |> 
  dplyr::select('seqnames', 'start', 'end', 
                'baseMean.nfyb', 'log2FoldChange.nfyb', 'padj.nfyb', 
                'baseMean.atac', 'log2FoldChange.atac', 'padj.atac', 
                'baseMean', 'log2FoldChange', 'padj',
                'new_cat', 'group', 'SYMBOL', 'distanceToTSS')

# Rename a couple columns for ease

colnames(nfyb_dar_anno_all)[10:12] <- c('baseMean.rna', 'log2FoldChange.rna', 'padj.rna')

plt_rna_by_nfyb <- nfyb_dar_anno_all |> 
  ggplot(aes(x = group, y = log2FoldChange.rna)) + 
  geom_boxplot(width = 0.3) + 
  theme_bw() + 
  theme(axis.text = element_text(size = 6, family = 'sans', color = 'black'), 
        axis.title = element_text(size = 7, family = 'sans', color = 'black'),
        )  +
  xlab('') + ylab('Log2 Fold Change RNA') + 
  geom_hline(aes(yintercept = 0), linetype = 'dashed')

plt_atac_by_nfyb <- nfyb_dar_anno_all |> 
  ggplot(aes(x = group, y = log2FoldChange.atac)) + 
  geom_boxplot(width = 0.3) + 
  theme_bw() + 
  theme(axis.text = element_text(size = 6, family = 'sans'), 
        axis.title = element_text(size = 7, family = 'sans'),
  )  +
  xlab('') + ylab('Log2 Fold Change RNA') + 
  geom_hline(aes(yintercept = 0), linetype = 'dashed')



# scatter plot version
plt_rna_by_nfyb_scatter <- nfyb_dar_anno_all |> 
  ggplot(aes(x = log2FoldChange.nfyb, y = log2FoldChange.rna)) + 
  geom_point() +  
  geom_smooth(method = 'lm' ) + 
  theme_bw() + 
  theme(axis.text = element_text(size = 6, family = 'sans'), 
        axis.title = element_text(size = 7, family = 'sans'),
  )  + 
  xlab('Log2 Fold Change NFYB') + ylab('Log2 Fold Change RNA')

cor_rna_nfyb <- cor.test(nfyb_dar_anno_all$log2FoldChange.nfyb, nfyb_dar_anno_all$log2FoldChange.rna)
cor_rna_nfyb$p.value
cor_rna_nfyb$estimate

plt_atac_by_nfyb_scatter <- nfyb_dar_anno_all |> 
  ggplot(aes(x = log2FoldChange.nfyb, y = log2FoldChange.atac)) + 
  geom_point() +  
  geom_smooth(method = 'lm' ) + 
  theme_bw() + 
  theme(axis.text = element_text(size = 6, family = 'sans', color = 'black'), 
        axis.title = element_text(size = 7, family = 'sans', color = 'black'),
  )  + 
  xlab('Log2 Fold Change NFYB') + ylab('Log2 Fold Change ATAC')
  
cor_atac_nfyb <- cor.test(nfyb_dar_anno_all$log2FoldChange.nfyb, nfyb_dar_anno_all$log2FoldChange.atac)
cor_atac_nfyb$p.value
cor_atac_nfyb$estimate

# Distribution Plot
nfyb_distribution_plot <- nfyb_dar_anno_all |> 
  as_tibble() |>  
  group_by(group, new_cat) |> 
  dplyr::summarise(cat_total = dplyr::n() ) |> 
  mutate(group_total = sum(cat_total), 
         frac = cat_total/group_total) |> 
  mutate(new_cat = factor(new_cat, levels = rev(unique(new_cat) ) )) |>  
  filter(!is.na(group)) |> 
  ggplot(aes(  x= group, y = frac, fill = new_cat)) +
  geom_col(color = 'grey30', position = 'stack') + 
  scale_fill_brewer(palette =   'Paired') + 
  theme_minimal() + 
  theme(axis.text = element_text(size = 6, family = 'sans', color = 'black'), 
        axis.title = element_text(size = 7, family = 'sans', color = 'black') ) + 
  xlab('')+ ylab('Fraction of Peaks') +
  geom_hline(yintercept = -0.01, linetype = 'solid', color = 'grey10')
tab_n_dist <-nfyb_dar_anno_all |> dplyr::count(group) |> tableGrob()


plt_volcano <- nfyb_dar_anno_all |> 
  ggplot(aes(x = log2FoldChange.nfyb, y = -log10(padj.nfyb), color = padj.nfyb < 0.1)) + 
  geom_point() + 
  theme_bw() + 
  ylab('Log10 padj NFYB') + xlab('Log2 FC NFYB DMSO/SNDX')+
  theme(axis.text = element_text(size = 6, family = 'sans', color = 'black'), 
             axis.title = element_text(size = 7, family = 'sans', color = 'black') )  +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', color = 'grey30')+
  scale_color_manual(values = c('grey50','red3'))

pdf('plots/cnr/fig5_nfyb_overivew.pdf', width= 8.5, height = 11) 
pageCreate(showGuides=F)
plotGG(nfyb_distribution_plot, height = 2.5, width = 4, x = 0.5, y = 0.5)
plotGG(plt_atac_by_nfyb_scatter + ggtitle('ATAC'), height = 2, width = 2, x = 0.5, y = 3.0)
plotGG(plt_rna_by_nfyb_scatter + ggtitle('RNA'), height = 2, width = 2, x = 3, y = 3.0)
plotText(label = paste0("r=", round(cor_atac_nfyb$estimate, 2)), x = 1.0, y = 5.15, just = 'left', fontsize = 6) 
plotText(label = paste0("r=", round(cor_rna_nfyb$estimate, 2)), x = 3, y = 5.15, just = 'left', fontsize = 6) 
plotText(label = paste0("p=", sprintf("%.2e", cor_atac_nfyb$p.value)), x = 1, y = 5.25, just = 'left', fontsize = 6)
plotText(label = paste0("p=", sprintf("%.2e", cor_rna_nfyb$p.value)), x = 3, y = 5.25, just = 'left', fontsize = 6) 

plotGG(plt_atac_by_nfyb + ggtitle('ATAC'), height = 2, width = 2, x = 0.5, y = 6)
plotGG(plt_rna_by_nfyb + ggtitle('RNA'), height = 2, width = 2, x= 3, y = 6)
plotGG(plt_volcano+ theme(legend.position = 'none'), height  = 2, width = 2, x = 0.5, y = 8.5)
plotText(label = paste0('Down:', length(nfyb_dar_anno_all |> filter(group == 'Down') |> pull(group))), x = 0.8, y = 9.1, fontsize = 6)
plotText(label = paste0('Up:', length(nfyb_dar_anno_all |> filter(group == 'Up') |> pull(group))), x = 2, y = 9.1, fontsize = 6, just = 'right')

plotGG(tab_n_dist, x = 6, y = 6, height = 1, width = 1)

dev.off()# Plot these in plotgardener



write_bed(nfyb_dar_anno_all |> filter(group == 'Up'), file = 'data/processed_data/integrated/nfyb_dor_anno_up.bed')
write_bed(nfyb_dar_anno_all |> filter(group == 'Down'), file = 'data/processed_data/integrated/nfyb_dor_anno_down.bed')
write_tsv(nfyb_dar_anno_all, file = 'data/processed_data/integrated/nfyb_dor_anno_all.tsv')

