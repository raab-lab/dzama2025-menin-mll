# CRISPR plotting
# PLots for Figure 1
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(ComplexHeatmap)
  library(depmap)
  library(circlize)
})
output_dir <- 'plots/crispr/'
if(!dir.exists(output_dir)) {dir.create(output_dir)}
# read in the data from previous processing
all_combined <- read_tsv('data/derived_data/crispr/tables/all_scores_combined.rra.tsv')
all_wide     <- read_tsv('data/derived_data/crispr/tables/all_scores_combined_wide.rra.tsv')

# Bring some depmap data for filtering
depmap_scores <- depmap::depmap_crispr()
 
depmap_genes <- depmap::gene_summary_22Q1()
common_essential <- depmap_genes |> 
  filter(dataset == 'Chronos_Combined') |> 
  filter(common_essential == TRUE)
common_essential

strongly_selective <- depmap_genes |> 
  filter(dataset == 'Chronos_Combined') |> 
  filter(strongly_selective == TRUE)

# Need a filter for my significant genes
all_combined_sig <- all_combined |> 
  filter(FDR < 0.01) |> 
  filter(!id %in% common_essential$gene_name) # lots of genes so worth limiting

epi_selective <- strongly_selective |> filter(gene_name %in% all_combined_sig$id)
epi_selective

all_wide_sig <-  all_combined|>  pivot_wider(id_cols = id, names_from = comp, values_from = c(Score, FDR) ) 
all_wide_sig <- all_wide_sig |> filter(id %in% all_combined_sig$id)
all_wide_sig_ids <- all_wide_sig$id
score_mat <- as.matrix(all_wide_sig[,2:8])
rownames(score_mat) <- all_wide_sig_ids
Heatmap(score_mat)

hist(all_combined$FDR)

# Define concordance score to find genes most essential
all_ranked <- all_combined |> 
  group_by(id) |> 
  summarise(mean_score  = mean(Score), 
            total_score = sum(Score), 
            max_sore = max(Score), 
            min_score = min(Score)) |> 
  filter(!id %in% common_essential$gene_name) |> 
  arrange(total_score)

all_ranked
top_n <- all_ranked |> head(30)

# Plot top 20 candidates
# if you want to order them by hand
top_heatmap <- all_combined |> filter(id %in% top_n$id) |> 
  separate(comp, into  = c('cell_line', 'treatment', 'vs', 't0', 'expt'), sep = '_', remove = F) |> 
  filter(!cell_line == 'Huh7') |> 
  mutate(label = paste(cell_line, treatment, sep = '_') ) |> 
  ggplot(aes(x = label, y = id, fill = Score)) + 
  geom_tile(color = 'grey50') + 
  scale_fill_gradientn(colors = c('blue', 'white'))  + 
  theme_minimal() + coord_fixed() + 
  theme(axis.text.x = element_text(angle = 90)) + 
  xlab('') + ylab('')
  

ggsave(filename = file.path(output_dir, 'top_30_tileplot.pdf'), height = 6, width = 2.5)

# If you want a clustered version
col_fun = colorRamp2(c(-5, 0, 5), c("blue", "black", "yellow"))
top_20_mat <- score_mat[rownames(score_mat) %in% top_n$id, ]
pdf(file.path(output_dir, 'top_30_heatmap.pdf'), height = 9, width = 4.5)
Heatmap(top_20_mat, col = col_fun)
dev.off()
#HLF

# MENIN COMPLEX 
men <- c('MEN1', 'ASH2L', 'WDR5', 'RBBP5', 'DPY30', 'KMT2A', 'PSIP1', 'MLL1')
#Plot HLF and PLCPRF5 2D vs 3D 
hlf_scatter <- all_wide |> 
  ggplot(aes(x = Score_HLF_2D_vs_T0_818, y = Score_HLF_3D_vs_T0_818)) + 
  geom_point() + 
  ggrepel::geom_label_repel(aes(label = id), data = all_wide |> filter(id %in% men),   min.segment.length = 2) + 
  geom_point(color = 'red2', data = all_wide |> filter(id %in% men), size = 2) + 
  geom_hline(yintercept = c(-1,1), linetype = 'dashed', color = 'darkred') + 
  geom_vline(xintercept = c(-1,1), linetype = 'dashed', color = 'darkred') + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16) ) + 
  xlab('2D RRA Score') + ylab('3D RRA Score')


#PLCPRF5
plcprf5_scatter <- all_wide |> 
  ggplot(aes(x = Score_PLCPRF5_2D_vs_T0_844, y = Score_PLCPRF5_3D_vs_T0_844)) + 
  geom_point() + 
  ggrepel::geom_label_repel(aes(label = id), data = all_wide |> filter(id %in% men), min.segment.length = 2) + 
  geom_point(color = 'red2', data = all_wide |> filter(id %in% men), size = 2) + 
  geom_hline(yintercept = c(-1,1), linetype = 'dashed', color = 'darkred') + 
  geom_vline(xintercept = c(-1,1), linetype = 'dashed', color = 'darkred') + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16) ) + 
  xlab('2D RRA Score') + ylab('3D RRA Score')


#HepG2
hepg2_scatter <- all_wide |> 
  ggplot(aes(x = Score_HepG2_2D_vs_T0_844, y = Score_HepG2_3D_vs_T0_844)) + 
  geom_point() + 
  ggrepel::geom_label_repel(aes(label = id), data = all_wide |> filter(id %in% men), min.segment.length = 2) + 
  geom_point(color = 'red2', data = all_wide |> filter(id %in% men), size = 2) + 
  geom_hline(yintercept = c(-1,1), linetype = 'dashed', color = 'darkred') + 
  geom_vline(xintercept = c(-1,1), linetype = 'dashed', color = 'darkred') + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14) ) + 
  xlab('2D RRA Score') + ylab('3D RRA Score')

#Huh7
huh7_scatter <- all_wide |> 
  ggplot(aes(x = Score_Huh7_2D_vs_T0_885, y = Score_Huh7_3D_vs_T0_885)) + 
  geom_point() + 
  ggrepel::geom_label_repel(aes(label = id), data = all_wide |> filter(id %in% men), min.segment.length = 2) + 
  geom_point(color = 'red2', data = all_wide |> filter(id %in% men), size = 2) + 
  geom_hline(yintercept = c(-1,1), linetype = 'dashed', color = 'darkred') + 
  geom_vline(xintercept = c(-1,1), linetype = 'dashed', color = 'darkred') + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14) ) + 
  xlab('2D RRA Score') + ylab('3D RRA Score')

# save plots
ggsave(filename = file.path(output_dir, 'hlf_2d_vs_3d_scatter.pdf'), hlf_scatter, height = 2, width = 2, scale = 2)
ggsave(filename = file.path(output_dir, 'plcprf5_2d_vs_3d_scatter.pdf'), plcprf5_scatter, height = 2, width = 2, scale = 2)
ggsave(filename = file.path(output_dir, 'hepg2_vs_3d_scatter.pdf'), hepg2_scatter, height = 2, width = 2, scale = 2)
ggsave(filename = file.path(output_dir, 'huh7_2d_vs_3d_scatter.pdf'), huh7_scatter, height = 2, width = 2, scale = 2)
