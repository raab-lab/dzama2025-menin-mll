suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(DESeq2)
  library(UpSetR) 
  library(plotgardener)
})
# compariosn of plc and hlf for menin and sndx
hlf_men <- read_tsv('data/processed_data/rna/hlf_menko_v_wt.tsv')
hlf_sndx <- read_tsv('data/processed_data/rna/hlf_sndx_v_dmso.tsv')

plc_men <- read_tsv('data/processed_data/rna/plc_ko_v_wt.tsv')
plc_sndx <- read_tsv('data/processed_data/rna/plc_sndx_v_dmso.tsv')

# Define upregualted and down for each group based on adj p-value only (don't want over filter here)

# up
hlf_men_up <- hlf_men |> filter(padj < 0.05) |> filter(log2FoldChange > 0) 
hlf_sndx_up <- hlf_sndx |> filter(padj < 0.05) |> filter(log2FoldChange > 0) 
plc_men_up <- plc_men |> filter(padj < 0.05) |> filter(log2FoldChange > 0) 
plc_sndx_up <- plc_sndx |> filter(padj <0.05) |> filter(log2FoldChange > 0) 

# down 
hlf_men_down <- hlf_men |> filter(padj < 0.05) |> filter(log2FoldChange < 0) 
hlf_sndx_down <- hlf_sndx |> filter(padj < 0.05) |> filter(log2FoldChange < 0) 
plc_men_down <- plc_men |> filter(padj < 0.05) |> filter(log2FoldChange < 0) 
plc_sndx_down <- plc_sndx |> filter(padj <0.05) |> filter(log2FoldChange < 0) 


# upset of upregulated
up_list <- list('HLF MEN1KO vs WT' = hlf_men_up$rowname, 
                'HLF SNDX vs Veh' = hlf_sndx_up$rowname, 
                'PLC MEN1KO vs WT' = plc_men_up$rowname, 
                'PLC SNDX vs Veh' = plc_sndx_up$rowname)
up_upset <- upset(fromList(up_list), order.by = 'degree',point.size = 1,line.size = 0.6,text.scale = 0.8  )

# about 240 upregulated in all 4 conditions


# upset of downregulated
down_list <- list('HLF MEN1KO vs WT' = hlf_men_down$rowname, 
                  'HLF SNDX vs Veh' = hlf_sndx_down$rowname, 
                  'PLC MEN1KO vs WT' = plc_men_down$rowname, 
                  'PLC SNDX vs Veh' = plc_sndx_down$rowname)

down_upset <- upset(fromList(down_list), order.by = 'degree', point.size = 1, line.size = 0.6, text.scale = 0.8)


# Look at MENKO and ASH2LKO


pdf('plots/rna/upset_up_and_down.pdf', height = 11, width = 8.5)
pageCreate(showGuide = F)
plotGG(up_upset, x = -0.5, y = 0.5, height = 2.25, width = 4, just = c('left', 'top'))
plotGG(down_upset, x = -0.5, y = 2.8, height = 2.25, width = 4, just = c('left', 'top')) 
plotText(label= 'Upregulated', fontsize = 10,rot =  270, x = 3.5, y= 1.15, just = c('left', 'center'))
plotText(label= 'Downregulated', fontsize = 10,rot =  270, x = 3.5, y= 3.25, just = c('left', 'center'))
dev.off()




# Write out the upregulated and downregualted consistent
up_genes <- intersect(hlf_men_up$rowname, intersect(hlf_sndx_up$rowname, intersect(plc_men_up$rowname, plc_sndx_up$rowname)))
up_genes_df <- hlf_men_up |> filter(rowname %in% up_genes) |> select(rowname, gene_name)
down_genes <- intersect(hlf_men_down$rowname, intersect(hlf_sndx_down$rowname, intersect(plc_men_down$rowname, plc_sndx_down$rowname)))
down_genes_df <- hlf_men_down |> filter(rowname %in% down_genes) |> select(rowname, gene_name)
write_tsv(up_genes_df, 'data/processed_data/rna/upregulated_overlap_menko_sndx.tsv')
write_tsv(down_genes_df, 'data/processed_data/rna/downregulated_overlap_menko_sndx.tsv')

