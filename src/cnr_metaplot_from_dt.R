# Metaplots from deeptools
library(profileplyr)
library(plotgardener)
library(tidyverse)

men_mat <- import_deepToolsMat('data/processed_data/cnr/men1_cov_mat.gz')
mll_mat <- import_deepToolsMat('data/processed_data/cnr/mll_cov_mat.gz')
k4me3_mat <- import_deepToolsMat('data/processed_data/cnr/h3k4me3_cov_mat.gz')
ash2l_mat <- import_deepToolsMat('data/processed_data/cnr/ash2l_cov_mat.gz')
igg_mat  <- import_deepToolsMat('data/processed_data/cnr/igg_cov_mat.gz')
k4_cols <-  '#31a354'
men_cols <-  '#3182bd'
mll_cols <-  '#756bb1'
ash2l_cols <- '#fdc086'
igg_cols  <- 'black'

plot_meta <- function(dt_mat, col) {
  g1_means <- colMeans(assays(dt_mat)[[1]])
  g2_means <- colMeans(assays(dt_mat)[[2]])
  
  df1 <- data.frame(pos = seq(-2000, 2000, length.out = length(g1_means)), vals = g1_means) |> mutate(group = 'DMSO')
  df2 <- data.frame(pos = seq(-2000, 2000, length.out = length(g2_means)), vals = g2_means) |> mutate(group = 'SNDX')
  df <- bind_rows(df1, df2)
  plt <- ggplot(df, aes(x = pos, y = vals, color = group, linetype = group)) + geom_line(linewidth = 0.5) + 
    scale_linetype_manual(values = c('solid', 'dashed')) + 
    scale_color_manual(values = c(col, col) ) + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5), 
          axis.title.x = element_text(size = 6), 
          axis.title.y = element_text(size = 7),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = 'grey90', linetype = 'dashed'))  + 
    xlab('Position Relative to Peak Center') + 
    ylab('Mean Signal')  
  return(plt) 
} 

# plots for figure 2
men_plot <- plot_meta(men_mat, men_cols)
mll_plot <- plot_meta(mll_mat, mll_cols)
k4_plot <-  plot_meta(k4me3_mat, k4_cols)
ash2l_plot <- plot_meta(ash2l_mat, ash2l_cols)
igg_plot  <- plot_meta(igg_mat, igg_cols)+ scale_y_continuous(limits = c(0,4))

pdf('plots/cnr/all_meta.pdf', height = 11, width = 8.5) 
pageCreate(width = 8.5, height = 11, default.units = 'inches', showGuides = F) 
plotGG(men_plot + theme(legend.position = 'none'),  
       width = 1.5, 
       height = 1.5,
       x = 0.5, y = 0.5) 
plotGG(mll_plot + theme(legend.position = 'none') + ylab(''),  
       width = 1.5, 
       height = 1.5,
       x = 2.1, y = 0.5) 
plotGG(k4_plot + theme(legend.position = 'none') + ylab(''),  
       width = 1.5, 
       height = 1.5,
       x = 3.7, y = 0.5) 
#plotGG(ash2l_plot + theme(legend.position = 'none') + ylab(''),
#  width = 1.0, 
#  height = 1.0,
#  x = 5.3, y = 0.5) 

plotGG(igg_plot + theme(legend.position = 'none') + ylab(''),
       width = 1.5, 
       height = 1.5,
       x = 5.3, y = 0.5) 
  

plotText(label = 'Menin', x = 1.4, y = 0.4, just = c('top', 'center') )
plotText(label = 'MLL1', x = 1.4 + 1.6, y = 0.4, just = c('top', 'center') )
plotText(label = 'H3K4me3', x = 1.4+ 3.2, y = 0.4, just = c('top', 'center') ) 
plotText(label = 'IgG', x = 1.4 + 4.8, y = 0.4, just = c('top', 'center'))
dev.off()


# Same plots but for figure 5
men_mat <- import_deepToolsMat('data/processed_data/cnr/men1_cov_nfyb_peaks_mat.gz')
mll_mat <- import_deepToolsMat('data/processed_data/cnr/mll_cov_mat.gz')
k4me3_mat <- import_deepToolsMat('data/processed_data/cnr/h3k4me3_cov_mat.gz')
ash2l_mat <- import_deepToolsMat('data/processed_data/cnr/ash2l_cov_mat.gz')
igg_mat  <- import_deepToolsMat('data/processed_data/cnr/igg_cov_mat.gz')
k4_cols <-  '#31a354'
men_cols <-  '#3182bd'
mll_cols <-  '#756bb1'
ash2l_cols <- '#fdc086'
igg_cols  <- 'black'

plot_meta_group <- function(dt_mat, col) {
  g1_means <- colMeans(assays(dt_mat)[[1]])
  g1_groups <- rowData(dt_mat)$dpGroup
  g2_means <- colMeans(assays(dt_mat)[[2]])
  g2_groups <- rowData(dt_mat)$dpGroup
  df1 <- data.frame(pos = seq(-2000, 2000, length.out = length(g1_means)), vals = g1_means) |> mutate(group = 'DMSO') 
  df2 <- data.frame(pos = seq(-2000, 2000, length.out = length(g2_means)), vals = g2_means) |> mutate(group = 'SNDX') 
  print (df1)
  df <- bind_rows(df1, df2)
  plt <- ggplot(df, aes(x = pos, y = vals, color = group, linetype = group)) + geom_line(linewidth = 0.5) + 
    scale_linetype_manual(values = c('solid', 'dashed')) + 
    scale_color_manual(values = c(col, col) ) + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5), 
          axis.title.x = element_text(size = 6), 
          axis.title.y = element_text(size = 7),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = 'grey90', linetype = 'dashed'))  + 
    xlab('Position Relative to Peak Center') + 
    ylab('Mean Signal')  + 
    facet_wrap(dp)
  return(plt) 
} 

plot_meta_group(men_mat, men_cols)
men_mat
