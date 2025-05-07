# Plot homer data from NFYB motifs
# compare overlap atac with those that do not

# Author: Jesse Raab
# Date: 2025-01-27
library(tidyverse)
library(plotgardener)

overlap_motifs <-read_tsv(file = 'data/derived_data/homer/nfyb_atac_overlap/knownResults.txt', 
                     col_types = c('c', 'c', 'd', 'd', 'd', 'd', 'd', 'd', 'd')) |> 
  janitor::clean_names() |> 
  mutate(percent_of_target_sequences_with_motif = as.numeric(str_replace(percent_of_target_sequences_with_motif, "%", ""))) |> 
  mutate(percent_of_background_sequences_with_motif = as.numeric(str_replace(percent_of_background_sequences_with_motif, "%", "") ) ) |> 
  mutate(group = 'Overlap') |> 
  mutate(enrichment = percent_of_target_sequences_with_motif/percent_of_background_sequences_with_motif) |> 
  dplyr::select(motif_name, log_p_value, enrichment)

alone_motifs <-read_tsv(file = 'data/derived_data/homer/nfyb_atac_nooverlap/knownResults.txt', 
                          col_types = c('c', 'c', 'd', 'd', 'd', 'd', 'd', 'd', 'd')) |> 
  janitor::clean_names() |> 
  mutate(percent_of_target_sequences_with_motif = as.numeric(str_replace(percent_of_target_sequences_with_motif, "%", ""))) |> 
  mutate(percent_of_background_sequences_with_motif = as.numeric(str_replace(percent_of_background_sequences_with_motif, "%", "") ) ) |> 
  mutate(group = 'Alone') |> 
  mutate(enrichment = percent_of_target_sequences_with_motif/percent_of_background_sequences_with_motif) |> 
  dplyr::select(motif_name, log_p_value, enrichment)

comb_motifs <- full_join(overlap_motifs, alone_motifs, by = 'motif_name', suffix = c('.overlap', '.alone'))
nfyb_motif_plot <- comb_motifs |> 
  ggplot(aes(x = enrichment.overlap, y = enrichment.alone)) + 
  geom_point() +
  ggrepel::geom_label_repel(aes(x = enrichment.overlap, y = enrichment.alone), 
                            label = 'NFY', data = comb_motifs |> filter(motif_name == 'NFY(CCAAT)/Promoter/Homer'),
                            size = 2.4)+
  theme_bw() + 
  labs(x ="Motif Enrichment at NFYB\n with ATAC sites", 
       y = "Motif Enrichment at NFYB\n without ATAC peak") + 
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size =8))

pdf('plots/cnr/fig5_nfyb_motif_enrichment.pdf', height = 11, width = 8.5)
pageCreate(height = 11, width = 8.5, showGuides = F)
plotGG(nfyb_motif_plot, width = 2.5, height = 2.5, x = 1, y = 1)
dev.off()

# Create a new column to identify points below the line - this ended up bein gunintersting -only ZNF528 and not much known 
# Assuming the rightmost line has slope=1 and intercept=-1 (adjust these values as needed)
# code from claude 3.5
# is_below_line <- function(x, y, slope, intercept) {
#   y < (slope * x + intercept)
# }
# 
# is_above_line <- function(x, y, slope, intercept) {
#   y > (slope * x + intercept)
# }
# comb_motifs$below_line <- is_below_line(
#   x = comb_motifs$enrichment.overlap,
#   y = comb_motifs$enrichment.alone,
#   slope = 1,
#   intercept = -1
# )
# comb_motifs$above_line <- is_above_line(
#   x = comb_motifs$enrichment.overlap,
#   y = comb_motifs$enrichment.alone,
#   slope = 1,
#   intercept = 1
# )
# 
# comb_motifs$label_points <-max(comb_motifs$above_line, comb_motifs$below_line) 
# comb_motifs |> 
#   ggplot(aes(x = enrichment.overlap, y = enrichment.alone)) + 
#   geom_point() 