# Take the up and downregualted consistant genes (defined in rna_upset) 
# plot a z-scored heatmap
suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(DESeq2)
  library(plotgardener)
  
})

# load the objects with count data
load('data/processed_data/rna/des_obj_jr.Rda')

#  load genes of interest from rna_upset.R 
up_genes <- read_tsv('data/processed_data/rna/upregulated_overlap_menko_sndx.tsv')
down_genes <- read_tsv('data/processed_data/rna/downregulated_overlap_menko_sndx.tsv')

# sanity check things are in right order
table(rownames(des_hlf_ash) == rownames(des_hlf_men))
table(rownames(des_hlf_ash) == rownames(des_plc))

# Make a separate count table for each object
hlf_ash_counts <- counts(des_hlf_ash, norm = T)
hlf_men_counts <- counts(des_hlf_men, norm = T)
plc_men_counts <- counts(des_plc, norm = T)

# get colData for each and combine - for annotations
#hlf_ash_cd <- colData(des_hlf_ash) |> as_tibble()
hlf_men_cd <- colData(des_hlf_men) |> as_tibble()
plc_men_cd <- colData(des_plc) |> as_tibble()

#all_cd <- bind_rows(hlf_ash_cd, hlf_men_cd, plc_men_cd)
# remove the ash2l stuff
all_cd <- bind_rows(hlf_men_cd, plc_men_cd)

# combined all data sets
#all_data <- bind_cols(hlf_ash_counts, hlf_men_counts,plc_men_counts) |> as.matrix()
#exclude ash2l for now
all_data <- bind_cols( hlf_men_counts,plc_men_counts) |> as.matrix()
rownames(all_data) <- rownames(hlf_ash_counts)


# filter down to genes of interest (rows)
all_data_goi <- all_data[rownames(all_data) %in% c(up_genes$rowname, down_genes$rowname),]


# scale data
all_data_goi <- t(scale(t(all_data_goi), center = T, scale = T) )
hm_cols <- circlize::colorRamp2(breaks = c(-2,0,2), colors = c('blue', 'black', 'gold'))
hm_anno <- HeatmapAnnotation( 
                  `Cell Line` = all_cd$cell_line,
                  Treatment = all_cd$treatment,
                  col = list(
                    `Cell Line` = c('PLCPRF5'= 'grey30', 'HLF' = 'grey70'),
                    Treatment = c('SNDX' = '#92c5de', 'DMSO' = '#ca0020', 
                                  'MEN1-KO' = '#0571b0', 'NTC' = '#f4a582')
                  ) )

hm_plot <- Heatmap(all_data_goi, 
        show_row_names = F, 
        col = hm_cols, 
        top_annotation = hm_anno, 
        show_column_names = F, 
        heatmap_legend_param =list( title = "Expression"))

pdf('plots/rna/men_hm.pdf', height = 11, width = 8.5)
pageCreate(showGuides = F)
plotGG(hm_plot, x = 0.5, y= 0.5, height = 4, width = 5) 
dev.off()
