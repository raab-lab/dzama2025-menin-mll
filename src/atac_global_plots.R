# Plot  ATAC data to show no changes

suppressPackageStartupMessages({ 
  library(tidyverse)
  library(plyranges)
  library(DESeq2)
  })


load('data/derived_data/atac/atac_de.Rda')
output_dir <- 'plots/atac/'
if(!dir.exists(output_dir)) { dir.create(output_dir) }
threshold <- log2(1.5) # just ot make it easier to change

res <- lfcShrink(des, coef = 2, type = 'apeglm', format = 'GRanges')
res <- keepStandardChromosomes(res, pruning.mode = 'coarse')

ma_plot <- res |> as_tibble() |> 
  ggplot(aes(y =log2FoldChange, x = log(baseMean)) )+ 
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > threshold) )+ 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(linetype = 'dashed', color = 'grey70'),
        axis.text = element_text(size = 14), 
        axis.title = element_text(size= 16)) + 
  xlab('Log Expression') + ylab ("Log2 Fold Change")  + 
  scale_color_manual(values = c('grey', 'steelblue')) + 
  geom_hline(yintercept = c(-threshold, threshold), linetype = 'dashed', color = 'red3')
  
ma_plot  

ggsave(filename = file.path(output_dir, 'atac_ma_plot.pdf'), plot = ma_plot, width = 7, height = 5)

summary(res$baseMean)
res |> as_tibble() |> arrange(log2FoldChange)


res |> as_tibble() |> filter(baseMean > 50) |> filter(abs(log2FoldChange) > 0.58) |> arrange(pvalue) |> View()
# Create a most differential atac list

up_most   <- res |> as_tibble() |> filter(baseMean > 50) |> filter(abs(log2FoldChange) > threshold) |> arrange(desc(log2FoldChange)) |> filter(log2FoldChange > 0) 
down_most <- res |> as_tibble() |> filter(baseMean > 50) |> filter(abs(log2FoldChange) > threshold) |> arrange(desc(log2FoldChange)) |> filter(log2FoldChange < 0) 

#suprising we see more upregualted than downregulated

write_bed(up_most, file = 'data/derived_data/atac/gain_access_most.bed')
write_bed(down_most, file = 'data/derived_data/atac/lost_access_most.bed')

up_most
down_most
