# test putting deeptools plots on plotgardener
# This code more or less works, but need some finessing to fix lengends and text

library(plotgardener)
library(plyranges)
library(profileplyr)
library(grid)
library(ComplexHeatmap)
library(EnrichedHeatmap)
ht_opt$message <- FALSE

# function to retrun a heatmap list for a pair
hm_fun <- function(dt_mat, col, name) {
  dt_em <- convertToEnrichedHeatmapMat(dt_mat)
  common_min = min(c(dt_em[[1]], dt_em[[2]]) )
  common_max = tail(quantile(c(dt_em[[1]], dt_em[[2]]),probs = seq(0,1, 0.01)), n = 2)[1]
  col_fun = circlize::colorRamp2(c(common_min, common_max), c("white", col))
  
  dmso_plot <- EnrichedHeatmap(dt_em[[1]],
                               top_annotation = NULL, 
                               show_row_names = F, 
                               col = col_fun)
  sndx_plot <- EnrichedHeatmap(dt_em[[2]],
                               top_annotation = NULL,
                               show_row_names = F, 
                               col = col_fun)
  ret_list <- list(dmso_plot, sndx_plot)    
  return(ret_list)
}



k4_mat  <-  import_deepToolsMat('data/processed_data/cnr/h3k4me3_cov_mat.gz')
menin_mat <- import_deepToolsMat('data/processed_data/cnr/men1_cov_mat.gz')
mll_mat   <- import_deepToolsMat('data/processed_data/cnr/mll_cov_mat.gz')
igg_mat <- import_deepToolsMat('data/processed_data/cnr/igg_cov_mat.gz')

k4_plots <- hm_fun(k4_mat, col = '#31a354', name = 'k4')
mll_plots <- hm_fun(mll_mat, col = '#756bb1', name = 'mll')
menin_plots <- hm_fun(menin_mat, col = '#3182bd', name = 'menin')
igg_plots  <- hm_fun(igg_mat , col = 'grey10', name = 'igg')


# This shoudl be done for all plots together to keep rows right
gg_all <- grid.grabExpr(draw(menin_plots[[1]] + menin_plots[[2]] + 
                              mll_plots[[1]] + mll_plots[[2]] + 
                              k4_plots[[1]] + k4_plots[[2]]+ 
                              igg_plots[[1]] + igg_plots[[2]], 
                             show_heatmap_legend = F) )
                                

pageCreate(width = 8.5, height = 11)
plotGG(plot = gg_all, x = 0.5, y = 0.5, just = c('top', 'left'), width = 1.75*4, height = 3.5)
plotText(label = 'DMSO', x = c(1, 1+1.7, 1+1.7*2, 1+1.7*3), y = 0.45)
plotText(label = 'SNDX', x = c(1.8, 1.8+1.7, 1.8+1.7*2, 1.8+1.7*3), y = 0.45)



