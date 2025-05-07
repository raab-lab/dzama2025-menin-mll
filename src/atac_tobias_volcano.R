# library plot output of bind detect for figures
library(tidyverse)
library(plotgardener)

bd <- read_tsv('data/external_data/output_hocomoco_motifs/bindetect_results.txt')
bd <- bd |> 
  mutate(short_name = str_extract(motif_id, pattern = "^[^.]*") )


label_df <- bd |> filter(grepl(short_name, pattern = 'TEAD|NFY|SP5|IRF')) 


  
plt <- bd |> 
  ggplot(aes(x = SNDX_footprints_DMSO_footprints_change, 
             y = -log10(SNDX_footprints_DMSO_footprints_pvalue))) + 
  geom_point() + 
  ggrepel::geom_label_repel(aes(label = short_name), data = label_df)


plt <- plt  + theme_bw() +  
  ylab('-log10 Pvalue SNDX/DMSO') + 
  xlab( 'TOBIAS Footprint Change SNDX/DMSO')

pdf('plots/atac/fig4_tobias_volanco.pdf', width = 8.5, height = 11)
pageCreate(height = 11, width = 8.5, showGuide = F)
plotGG(plt, x = 0.5, y =0.5, width =4 , height = 4)

dev.off() 
