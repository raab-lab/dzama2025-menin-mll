# Comparison plots for each cell line and for drug vs ko etc

suppressPackageStartupMessages({
  library(tidyverse)
  library(plotgardener)
  
})
paper_theme <- function() {
  theme(plot.title = element_text(family ="sans", color = 'black', 
                                  face = "bold", size = 10),
        axis.title = element_text(family = "sans", color = 'black', 
                                    size = 8),
        axis.text = element_text(size = 6, color = 'black', family = 'sans'), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = 'dashed', color = 'grey90'),
        panel.border = element_rect(color = 'grey10', fill = NA),
        panel.background = element_blank() )
}


hlf_men_ko_v_wt <- read_tsv('data/processed_data/rna/hlf_menko_v_wt.tsv')
hlf_drug_v_veh   <- read_tsv('data/processed_data/rna/hlf_sndx_v_dmso.tsv')

plc_men_ko_v_wt <- read_tsv('data/processed_data/rna/plc_ko_v_wt.tsv')
plc_drug_v_veh  <- read_tsv('data/processed_data/rna/plc_sndx_v_dmso.tsv')


hlf_men_ko_v_wt


# Make data frames for comparison by ggplot
hlf_gen_v_drug <- left_join(hlf_men_ko_v_wt, hlf_drug_v_veh, by = c('rowname', 'seqnames', 'start', 'end', 'strand', 'gene_name'), 
                            suffix = c('.gene', '.drug') )

hlf_gen_v_drug_plot <- hlf_gen_v_drug |> 
  filter(baseMean.gene > 20) |> 
  ggplot(aes(y = log2FoldChange.gene, x = log2FoldChange.drug)) + 
  geom_point() + 
  geom_smooth(method = 'lm') + 
  paper_theme()+
  ggtitle('HLF MEN1 KO vs SNDX')
hlf_gen_v_drug_cor <-cor.test(hlf_gen_v_drug$log2FoldChange.gene, hlf_gen_v_drug$log2FoldChange.drug)


# plc men ko vs drug
plc_ko_v_drug <- left_join(plc_men_ko_v_wt, plc_drug_v_veh, by = c('rowname', 'seqnames', 'start', 'end', 'strand', 'gene_name'),
                           suffix = c('.ko', '.drug'))
plc_plot <- plc_ko_v_drug |> 
  filter(baseMean.ko > 20) |> 
  ggplot(aes(y = log2FoldChange.ko, x = log2FoldChange.drug)) +
  geom_point() + 
  geom_smooth(method = 'lm') + 
  paper_theme()+ 
  ggtitle('PLCPRF5 MEN1 KO vs SNDX')
plc_cor <- cor.test(plc_ko_v_drug$log2FoldChange.ko, plc_ko_v_drug$log2FoldChange.drug)

##### Plot these as a supplemental page
pdf('plots/rna/scatter_comparisons.pdf', height = 11, width = 8.5)
pageCreate(width= 8.5, height = 11, default.units = 'inches', showGuides = F)
plotGG(hlf_gen_v_drug_plot, x = 0.5, y = 0.5, width = 3, height = 3)
plotText(label = paste0('r: ', round(hlf_gen_v_drug_cor$estimate, 2)), x = 3.25, y = 2.75, just = c('right', 'bottom'))
plotText(label = paste0('pval: ', round(hlf_gen_v_drug_cor$p.value)), x = 3.25, y = 3, just = c('right', 'bottom'))

plotGG(plc_plot, x = 0.5, y = 7.1, width = 3, height= 3)
plotText(label = paste0('r: ', round(plc_cor$estimate, 2)), x = 3.25, y = 9.3, just = c('right', 'bottom'))
plotText(label = paste0('pval: ', round(plc_cor$p.value, 2)), x = 3.25, y = 9.5, just = c('right', 'bottom'))


dev.off()
