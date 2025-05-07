# basic RNAseq plots - figure 3

library(DESeq2)
library(tidyverse)
library(plotgardener)
library(UpSetR)

hlf_sndx <- read_tsv('data/processed_data/rna/hlf_sndx_v_dmso.tsv')
hlf_ko   <- read_tsv('data/processed_data/rna/hlf_menko_v_wt.tsv')

plc_sndx <- read_tsv('data/processed_data/rna/plc_sndx_v_dmso.tsv')
plc_ko   <- read_tsv('data/processed_data/rna/plc_ko_v_wt.tsv')

plc_ko |> ggplot(aes(x = log2FoldChange, y = -log10(padj))) + geom_point()

# just so I can check counts  
load('data/processed_data/rna/des_obj_jr.Rda')

plotCounts(des_hlf_men, gene = 'ENSG00000025039', intgroup = 'treatment')
colData(des_hlf_men)

# plot scatter of shrunken fold change
hlf_comb <- left_join(hlf_sndx, hlf_ko, by = c('rowname', 'gene_name'), suffix = c('.drug', '.ko') )
hlf_scatter <- hlf_comb |> 
  ggplot(aes(x = log2FoldChange.drug, y = log2FoldChange.ko))+
  geom_point() +
  geom_smooth(method = 'lm') + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 6, family = 'sans'),
        axis.title = element_text(size = 8, family = 'sans')) + 
  xlab('Log2 Fold Change SNDX/DMSO') + 
  ylab("Log2 Fold Change NTG/MEN1-KO")
hlf_val <- cor.test(hlf_comb$log2FoldChange.drug, hlf_comb$log2FoldChange.ko)$p.value
hlf_cor <- cor.test(hlf_comb$log2FoldChange.drug, hlf_comb$log2FoldChange.ko)$estimate

plc_comb <- left_join(plc_sndx, plc_ko, by = c('rowname', 'gene_name'), suffix = c('.drug', '.ko'))
plc_scatter <- plc_comb |> 
  ggplot(aes(x = log2FoldChange.drug, y = log2FoldChange.ko))+
  geom_point() +
  geom_smooth(method = 'lm') + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 6, family = 'sans'),
        axis.title = element_text(size = 8, family = 'sans')) + 
  xlab('Log2 Fold Change SNDX/DMSO') + 
  ylab("Log2 Fold Change NTG/MEN1-KO")

plc_val <- cor.test(plc_comb$log2FoldChange.drug, plc_comb$log2FoldChange.ko)$p.value
plc_cor <- cor.test(plc_comb$log2FoldChange.drug, plc_comb$log2FoldChange.ko)$estimate

plc_val
plc_cor

# Print out
pdf('plots/rna/final_scatter_fig3.pdf', height = 11, width = 8.5)
pageCreate(showGuide = F)
plotGG(hlf_scatter, height = 2, width = 2, y = 0.5, x = 0.5) 
plotGG(plc_scatter, height = 2, width = 2, y = 0.5, x = 3) 
plotText(label = paste0('r = ', round(hlf_cor, 2)), fontsize = 7, x = 1.15, y = 0.8)
plotText(label = paste0('r = ', round(plc_cor, 2)), fontsize = 7, x = 3.55, y = 0.8)
dev.off()
 
# upset plots for up and down in each condition
hlf_drug_up <- hlf_sndx |> filter(padj < 0.05 & log2FoldChange > 1) 
hlf_drug_down <- hlf_sndx |> filter(padj < 0.05 & log2FoldChange < 1)

hlf_ko_up <- hlf_ko |> filter(padj < 0.05 & log2FoldChange > 1) 
hlf_ko_down <- hlf_ko |> filter(padj < 0.05 & log2FoldChange < 1 ) 

plc_drug_up <- plc_sndx |> filter(padj < 0.05 & log2FoldChange > 1) 
plc_drug_down <- plc_sndx |> filter(padj < 0.05  & log2FoldChange < 1 ) 

plc_ko_up <- plc_ko |> filter(padj < 0.05 & log2FoldChange > 1) 
plc_ko_down <- plc_ko |> filter(padj < 0.05 & log2FoldChange < 1) 


# up 
up_list <- list(hlf_drug_up$rowname, hlf_ko_up$rowname, plc_drug_up$rowname, plc_ko_up$rowname) 
names(up_list) <- c('HLF SNDX/DMSO','HLF KO/WT', 'PLC SNDX/DMSO' , 'PLC KO/WT')
down_list <- list(hlf_drug_down$rowname, hlf_ko_down$rowname, plc_drug_down$rowname, plc_ko_down$rowname)
names(down_list) <- c('HLF SNDX/DMSO','HLF KO/WT', 'PLC SNDX/DMSO' , 'PLC KO/WT')

up_plot <- upset(fromList(up_list))
down_plot <- upset(fromList(down_list))

# Print out
pdf('plots/rna/final_scatter_fig3.pdf', height = 11, width = 8.5)
pageCreate(showGuide = F)
plotGG(hlf_scatter, height = 2, width = 2, y = 0.5, x = 0.5) 
plotGG(plc_scatter, height = 2, width = 2, y = 0.5, x = 3) 
plotText(label = paste0('r = ', round(hlf_cor, 2)), fontsize = 7, x = 1.15, y = 0.8)
plotText(label = paste0('r = ', round(plc_cor, 2)), fontsize = 7, x = 3.55, y = 0.8)
plotGG(up_plot, height = 2.5, width = 4.5, x = 0.5, y = 3)
plotGG(down_plot, height = 2.5, width = 4.5, x = 0.5, y = 6)
dev.off()


