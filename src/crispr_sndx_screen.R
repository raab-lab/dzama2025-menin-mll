suppressPackageStartupMessages({
  library(tidyverse)
  library(MAGeCKFlute)
  library(plotgardener)
  library(grid)
})

sndx_v_dmso <- ReadRRA('data/source_data/crispr/results/HLF_SNDX_vs_DMSO_885/HLF_SNDX_vs_DMSO_885.gene_summary.txt')
sndx_v_t0   <- ReadRRA('data/source_data/crispr/results/HLF_SNDX_vs_T0_885/HLF_SNDX_vs_T0_885.gene_summary.txt')
dmso_v_t0   <- ReadRRA('data/source_data/crispr/results/HLF_DMSO_vs_T0_885/HLF_DMSO_vs_T0_885.gene_summary.txt')

sndx_v_dmso |> filter(FDR < 0.1) |> arrange(FDR)
sndx_v_t0 |> filter(FDR < 0.05) |> arrange(FDR)

ggplot(sndx_v_dmso, aes(x = Score, y = -log10(FDR)) ) + 
  geom_point()


t0_comb <- left_join(dmso_v_t0, sndx_v_t0, by = 'id', suffix = c('.dmso', '.sndx'))

sndx_down <- sndx_v_dmso |> 
  filter(FDR < 0.1) |> 
  pull(id)

t0_comb |> 
  ggplot(aes(x = Score.dmso, y = Score.sndx)) + geom_point() + 
  geom_abline(slope = 1, intercept = c(-1,0,1) ) +
  geom_point(data = t0_comb |> filter(id %in% sndx_down), color = 'red2', size = 2)

# code from claude 3.5
is_below_line <- function(x, y, slope, intercept) {
  y < (slope * x + intercept) 
}

is_above_line <- function(x, y, slope, intercept) {
  y > (slope * x + intercept) 
}

# Create a new column to identify points below the line
# Assuming the rightmost line has slope=1 and intercept=-1 (adjust these values as needed)
t0_comb$below_line <- is_below_line(
  x = t0_comb$Score.dmso,
  y = t0_comb$Score.sndx,
  slope = 1,
  intercept = -1
)

t0_comb$above_line <- is_above_line(
  x = t0_comb$Score.dmso,
  y = t0_comb$Score.sndx,
  slope = 1,
  intercept = 1
)

keepers <- t0_comb |> filter(below_line == T) |> pull(id)
keepers <- c(keepers, 'CBX4')

# Create the plot
dvs_plot <- ggplot(t0_comb, aes(x = Score.dmso, y = Score.sndx)) +
  geom_point(aes(color = below_line)) +
  geom_abline(slope = 1, intercept = 1, linetype = 'dashed', color = 'red2') +
  geom_abline(slope = 1, intercept = -1, linetype = 'dashed', color = 'red2') +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', color = 'grey80') + 
  geom_vline(aes(xintercept = 0), linetype = 'dashed', color = 'grey80') + 
  scale_color_manual(values = c("black", "red3")) +
  theme_minimal() +
  theme(legend.position = "none")+ 
  ggrepel::geom_label_repel(
    aes(label = id),
    size = 5,
    box.padding = 0.5,
    point.padding = 0.2,
    segment.color = "grey50",
    show.legend = FALSE, 
    max.overlaps = 100, 
    data = t0_comb |> filter(id %in% keepers) 
  ) +
  xlab('Score DMSO') + ylab('Score SNDX')

pdf('plots/crispr/fig5_scatter.pdf', width = 8.5, height = 11)
pageCreate(height = 11, width = 8.5, showGuide = F)
plotGG(dvs_plot, x = 0.5, y = 0.5, width = 5, height=5)
dev.off()


t0_comb |> filter(below_line == TRUE)
t0_comb |> filter(id == 'TEAD4')
t0_comb |> mutate(diff = Score.dmso - Score.sndx) |> arrange(desc(diff))



#interesting genes based on semi arbitrary numbers
t0_comb <- t0_comb |>
  filter((Score.sndx < -0.5 & Score.dmso > -0.5) |(below_line==TRUE & Score.sndx < 0 )) |> 
  mutate(diff = Score.dmso - Score.sndx) |> 
  arrange(desc(diff)) |> 
  mutate(id = factor(id, levels = id))
  

t0_comb_hm_dat <- t0_comb |> 
    dplyr::select(id, Score.dmso, Score.sndx, diff) |> 
    pivot_longer(cols = -id, names_to = 'assay', values_to = 'scores')

hm <- grid.grabExpr(tidyHeatmap::heatmap(t0_comb_hm_dat |> head(20*3),
                           .row = id, 
                           .column = assay, 
                           .value = scores, 
                           scale = 'none', 
                           cluster_rows =T, 
                           cluster_columns = F,
                           palette_value = c('blue', 'white', 'red')))

pdf('plots/sndx_screen_hm.pdf', height = 11, width = 8.5)
pageCreate(showGuides = F)
plotGG(grid::grid.grabExpr(hm), x= 0.5, y = 0.5, height = 3.5, width = 3)
dev.off()
hm

pageCreate()
gg_heatmap  <- plotGG(
  plot = hm,
  x = 3, y = 0.5,
  width = 2.5, height = 3
)

  
counts <- read_tsv('data/source_data/crispr/counts/combined_counts.txt')
counts |> 
  pivot_longer(cols = c(-sgRNA, -Gene), names_to = 'samples', values_to = 'counts') |> 
  filter(Gene == 'KDM6A')  |> 
  ggplot(aes(x = samples, y = log2(counts))) + geom_point()
meta <- data.frame(samples = colnames(counts[,3:ncol(counts)]))
meta <- meta[grepl( pattern = 'T0|DMSO|SNDX', meta$sample),]
all_crispr <- read_csv('data/source_data/crispr/all_crispr_data_sets.csv') |> janitor::clean_names()
all_crispr$sample <- paste(all_crispr$cell_line, all_crispr$treatment, all_crispr$replicate,"_" )


# Plots of counts
counts <- read_tsv('data/source_data/crispr/counts/combined_counts.txt')
colSums(counts[,3:ncol(counts)])
counts_norm <- (counts[,3:ncol(counts)]/colSums(counts[,3:ncol(counts)] )* 1e6)
counts_norm <- bind_cols(counts[,1:2], counts_norm)
sample_info <- read_csv('data/source_data/crispr/all_crispr_data_sets.csv') |> janitor::clean_names()
sample_info$sample_id <- with(sample_info, paste(cell_line, treatment, replicate, sep = '_'))
sample_info <- sample_info |> mutate(format = paste0(format, 'D'))
sample_info
counts_long <- counts_norm|> 
  pivot_longer(cols = c(-sgRNA, -Gene), names_to = 'sample_id', values_to = 'counts' ) |> 
  left_join(sample_info, by = 'sample_id')
counts_long


tmp <- t0_comb |> filter(above_line==TRUE) |> pull(id) 
hlf_plot <- counts_long |> 
  filter(cell_line == 'HLF') |> 
  filter(Gene %in%  tmp) |> 
  filter(grepl(x = treatment, pattern = 'DMSO|SNDX|Drug')) |> 
  mutate(sample_grouping = paste(cell_line, treat, format, timepoint, sep = '_') ) |> 
  mutate(timepoint = factor(timepoint, levels = c('Initial', 'Final'))) |> 
  group_by(Gene, format, treat, timepoint, cell_line, sgRNA) |> 
  summarise(mean_count = mean(counts)) |>
  ggplot(aes(x = treat, y = log2(mean_count+1), fill = timepoint)) +
  geom_boxplot(alpha = 0.2, width = 0.5, outlier.color = NA) + 
  geom_point(aes(color = timepoint, group = timepoint), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5) ) +
  facet_wrap( ~ Gene)  +
  theme_bw()  +theme(axis.text = element_text(size = 7, color = 'grey10'), 
                     axis.title = element_text(size = 8, color = 'grey10'), 
                     strip.text = element_text(size = 9, color = 'grey10'),
                     strip.background = element_rect(fill = 'grey90')) + 
  xlab('') + ylab('Log2 Counts + 1') +
  scale_fill_manual(values = c('grey30','steelblue')) +
  scale_color_manual(values = c('grey30','steelblue')) 
hlf_plot


counts_long |> 
  filter(cell_line == 'HLF') |> 
  filter(Gene %in% c('TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'NFYB', 'NFYC', 'MEN1', 'ASH2L', 'NSD1', 'KMT2A', 'PRMT1', 'KDM6A')) |> 
  filter(grepl(x = treatment, pattern = 'DMSO|SNDX|Drug')) |> 
  mutate(sample_grouping = paste(cell_line, treat, format, timepoint, sep = '_') ) |> 
  mutate(timepoint = factor(timepoint, levels = c('Initial', 'Final'))) |> 
  group_by(Gene, format, treat, timepoint, cell_line, sgRNA) |> 
  summarise(mean_count = mean(counts)) |>
  ggplot(aes(x = treat, y = log2(mean_count+1), fill = timepoint)) +
  geom_boxplot(alpha = 0.2, width = 0.5, outlier.color = NA) + 
  geom_point(aes(color = timepoint, group = timepoint), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5) ) +
  facet_wrap( ~ Gene)  +
  theme_bw()  +theme(axis.text = element_text(size = 7, color = 'grey10'), 
                     axis.title = element_text(size = 8, color = 'grey10'), 
                     strip.text = element_text(size = 9, color = 'grey10'),
                     strip.background = element_rect(fill = 'grey90')) + 
  xlab('') + ylab('Log2 Counts + 1') +
  scale_fill_manual(values = c('grey30','steelblue')) +
  scale_color_manual(values = c('grey30','steelblue')) 
hlf_plot


counts_long |> 
  filter(cell_line == 'HLF') |> 
  filter(Gene %in% c('NFYB', 'NFYC', 'LEO1')) |> 
  filter(grepl(x = treatment, pattern = 'DMSO|SNDX|Drug')) |> 
  mutate(sample_grouping = paste(cell_line, treat, format, timepoint, sep = '_') ) |> 
  mutate(timepoint = factor(timepoint, levels = c('Initial', 'Final'))) |> 
  group_by(Gene, format, treat, timepoint, cell_line, sgRNA) |> 
  summarise(mean_count = mean(counts)) |>
  ggplot(aes(x = treat, y = log2(mean_count+1), fill = timepoint)) +
  geom_boxplot(alpha = 0.2, width = 0.5, outlier.color = NA) + 
  geom_point(aes(color = timepoint, group = timepoint), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5) ) +
  facet_wrap( ~ Gene)  +
  theme_bw()  +theme(axis.text = element_text(size = 7, color = 'grey10'), 
                     axis.title = element_text(size = 8, color = 'grey10'), 
                     strip.text = element_text(size = 9, color = 'grey10'),
                     strip.background = element_rect(fill = 'grey90')) + 
  xlab('') + ylab('Log2 Counts + 1') +
  scale_fill_manual(values = c('grey30','steelblue')) +
  scale_color_manual(values = c('grey30','steelblue')) 


t0_comb |> 
  filter(FDR.sndx < 0.05) |> 
  filter(below_line == 'TRUE') |> 
  mutate()