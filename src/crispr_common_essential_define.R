# filter for md's genes - and plot heatmap
# 
suppressPackageStartupMessages({
  library(tidyverse)}) 

all_combined <- read_tsv('data/derived_data/crispr/tables/all_scores_combined.rra.tsv')
all_wide     <- read_tsv('data/derived_data/crispr/tables/all_scores_combined_wide.rra.tsv')

all_wide <- all_wide |> 
  mutate(score_filter = ifelse(
    Score_HLF_2D_vs_T0_818 < 0 &  
      Score_HLF_3D_vs_T0_818 < 0 &  
      Score_PLCPRF5_2D_vs_T0_844 < 0 &  
      Score_PLCPRF5_3D_vs_T0_0042 < 0 ,
    'keep', 'drop' )  ) |> 
  mutate(fdr_filter  = ifelse(
    FDR_HLF_2D_vs_T0_818 < 0.05 & 
      FDR_HLF_3D_vs_T0_818 < 0.05 &
      FDR_PLCPRF5_2D_vs_T0_844 < 0.05 & 
      FDR_PLCPRF5_3D_vs_T0_0042 < 0.05, 
    'keep', 'drop') ) 

depmap_genes <- depmap::gene_summary_22Q1()
common_essential <- depmap_genes |> 
  filter(dataset == 'Chronos_Combined') |> 
  filter(common_essential == TRUE)
common_essential

table(all_wide$score_filter == 'keep', all_wide$fdr_filter == 'keep')

all_wide_common <- all_wide |> filter(score_filter == 'keep' & fdr_filter == 'keep') |> 
  filter(!id %in% common_essential$gene_name) 
