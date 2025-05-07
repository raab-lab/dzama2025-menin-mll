# Plot homoer motifs
# homer run separately
suppressPackageStartupMessages({
  library(motifStack)
  library(universalmotif)
  library(plotgardener)
  library(monaLisa)
  library(ggseqlogo)
  })
source('src/paper_themes.R')

up_motifs <-read_tsv(file = 'data/derived_data/homer/sndx_up_vs_unchanged/knownResults.txt', 
                     col_types = c('c', 'c', 'd', 'd', 'd', 'd', 'd', 'd', 'd')) |> 
  janitor::clean_names() |> 
  mutate(percent_of_target_sequences_with_motif = as.numeric(str_replace(percent_of_target_sequences_with_motif, "%", ""))) |> 
  mutate(percent_of_background_sequences_with_motif = as.numeric(str_replace(percent_of_background_sequences_with_motif, "%", "") ) ) |> 
  mutate(group = 'Gained')

down_motifs <-read_tsv(file = 'data/derived_data/homer/sndx_down_vs_unchanged/knownResults.txt', 
                     col_types = c('c', 'c', 'd', 'd', 'd', 'd', 'd', 'd', 'd')) |> 
  janitor::clean_names() |> 
  mutate(percent_of_target_sequences_with_motif = as.numeric(str_replace(percent_of_target_sequences_with_motif, "%", ""))) |> 
  mutate(percent_of_background_sequences_with_motif = as.numeric(str_replace(percent_of_background_sequences_with_motif, "%", "") ) ) |> 
  mutate(group = 'Lost')

motif_all <- bind_rows(up_motifs, down_motifs)
motif_all<- motif_all|> 
  separate(motif_name, into = c('motif_name', 'part2'), sep = "\\/") |> 
  mutate(motif_name = str_replace(motif_name, "\\([^()]*\\)", '') ) |> 
  arrange(log_p_value) 

motif_plot <- motif_all |> 
  ggplot(aes (x =log2(percent_of_target_sequences_with_motif/percent_of_background_sequences_with_motif), y = -log_p_value)) + 
  geom_point(size = 1) +  
  ggrepel::geom_text_repel(aes(label = motif_name),data = motif_all |> filter(-log10(p_value) > 5 ),
                           max.overlaps = 30, 
                           size =2.5)   + 
  facet_wrap(~group, scales = 'free') + 
  theme_bw()  +
  ylab('-p-value')  + 
  xlab('Target/Background') +
  paper_theme() + 
  theme(strip.text = element_text(size = 10, color = 'grey10'),
        strip.background = element_rect(fill = 'grey90') )
motif_plot
# Plot logos from these
read_homer_folder <- function(path){
  # returns list of universal motif objects
  files <- list.files(path, pattern = '*.motif', full.names = T)
  motifs <- unlist(lapply(files, read_homer))
  return(motifs)
} 



known_down <- read_homer_folder('data/derived_data/homer/sndx_down_vs_unchanged/knownResults')
known_up <- read_homer_folder('data/derived_data/homer/sndx_up_vs_unchanged//knownResults') 

down_1 <- ggseqlogo(known_down[[1]]@motif) + ggtitle(known_down[[1]]@name)
down_2 <- ggseqlogo(known_down[[2]]@motif) + ggtitle(known_down[[2]]@name)
up_1 <- ggseqlogo(known_up[[1]]@motif) + ggtitle(known_up[[1]]@name)
up_2 <-ggseqlogo(known_up[[2]]@motif) + ggtitle(known_up[[2]]@name) 

pdf('plots/atac/atac_homer_motif_plots.pdf', height = 11, width = 8.5) 
pageCreate(showGuides =F)
plotGG(motif_plot, x = 0.5, y = 0.5, width = 4, height = 2.5)
plotGG(down_1, x =0.5, y = 4.5, width = 3, height = 2.5, just = c('left', 'top') )
plotGG(down_2, x = 0.5, y = 6.75, width = 3, height = 2.5, just = c('left', 'top')) 
plotGG(up_1, x = 4.5, y =4.5, width = 3, height = 2.5, just = c('left', 'top'))
plotGG(up_2, x = 4.5, y = 6.75, width = 3, height = 2.5, just = c('left', 'top') )
plotText(label = 'Lost', fontsize = 18, rot = 0, x = 0.5, y = 4.0, height = 0.25, width = 2, just = c('center', 'center') ) 
plotText(label = 'Gained', fontsize = 18, rot = 0, x = 4.5, y = 4.0 , height = 0.25, width = 2, just = c('center', 'center') ) 
dev.off()

# Write out the processed files
write_tsv(motif_all, file = 'data/derived_data/atac/homer_motifs_processed.tsv')
