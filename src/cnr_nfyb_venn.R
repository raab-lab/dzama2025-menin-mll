# NFYB Venn Diagrams
#Author: Jesse Raab
# Date: 2025-01-28
# Figure 5A-B

# Menin-MLL-ASH2L overlap DMSO vs SNDX
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(plyranges)
  library(ChIPseeker)
  library(plotgardener)
  library(grid)
  library(ChIPpeakAnno)
  library(clusterProfiler)
  library(msigdbr)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene) 
  library(rGREAT)
  
})
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# Need to clean up their annotations and combine promoter, intron, and exon categories oso there is only one of each. 
fix_category <- function(tib) {
  tib |> 
    mutate(new_cat = case_when(grepl ('Promoter', annotation) ~ 'Promoter', 
                               grepl ('Intron', annotation) ~ 'Intron', 
                               grepl ('Exon', annotation) ~ 'Exon', 
                               .default = as.character(annotation) ) ) 
}

k4_cols <-  '#31a354'
men_cols <-  '#3182bd'
mll_cols <-  '#756bb1'
nfyb_cols <- '#d95f02'
atac_cols <- '#beaed4'

menin_dmso <- read_bed('data/derived_data/peak_sets/cnr_consensus/menin_dmso_consensus.bed') |> keepStandardChromosomes(pruning.mode = 'coarse')
menin_sndx <- read_bed('data/derived_data/peak_sets/cnr_consensus/menin_sndx_consensus.bed') |> keepStandardChromosomes(pruning.mode = 'coarse')


nfyb_dmso <- read_bed('data/derived_data/peak_sets/cnr_consensus/nfyb_dmso_consensus.bed') |> keepStandardChromosomes(pruning.mode = 'coarse')
nfyb_sndx <- read_bed('data/derived_data/peak_sets/cnr_consensus/nfyb_sndx_consensus.bed') |> keepStandardChromosomes(pruning.mode = 'coarse')

atac_dmso  <- read_bed('data/derived_data/peak_sets/atac_consensus/dmso_consensus_atac.bed') |> keepStandardChromosomes(pruning.mode = 'coarse')
atac_sndx <- read_bed('data/derived_data/peak_sets/atac_consensus/sndx_consensus_atac.bed') |> keepStandardChromosomes(pruning.mode = 'coarse')


fill <- c(nfyb_cols, men_cols)
dmso_list <- list(nfyb_dmso, menin_dmso)
sndx_list <- list(nfyb_sndx, menin_sndx)
dmso_plot <-  grid.grabExpr(makeVennDiagram(dmso_list, 
                                            NameOfPeaks = c('NFYB', 'Menin'), 
                                            connectedPeaks = 'merge', 
                                            main.fontfamily  = 'sans', 
                                            fill = fill) )

sndx_plot <-  grid.grabExpr(makeVennDiagram(sndx_list, 
                                            NameOfPeaks = c('NFYB', 'Menin'), 
                                            connectedPeaks = 'merge', 
                                            main.fontfamily  = 'sans', 
                                            fill = fill) )


# Venn of NFYB DMSO vs SNDX
nfyb_list <- list(nfyb_dmso, nfyb_sndx)
nfyb_plot <- grid.grabExpr(makeVennDiagram(nfyb_list, 
                              NameOfPeaks = c('DMSO', 'SNDX'), 
                              connectedPeaks = 'merge', 
                              main.fontfamily  = 'sans', 
                              fill = c(nfyb_cols, '#ef6548')) )

pdf('plots/cnr/fig5_nfyb_venns.pdf', width = 8.5, height = 11)
pageCreate(width = 8.5, height = 11, showGuides = F)
plotGG(dmso_plot, width = 3.5, height = 3.5, x = 0.5, y = 0.5)
plotText(label = 'DMSO Plot', fontsize = 10, x = 4.5, y = 1)
plotGG(sndx_plot, width = 3.5, height = 3.5, x = 0.5, y = 4)
plotText(label = 'SNDX Plot', fontsize = 10, x = 4.5, y= 4.5)
plotGG(nfyb_plot, width = 3.5, height = 3.5, x = 0.5, y = 7.5)
dev.off()

################################################################################
# Supplemental Venn on ATAC + NFYB Overlap
nfyb_dmso <- read_bed('data/derived_data/peak_sets/cnr_consensus/nfyb_union.bed')|> keepStandardChromosomes(pruning.mode = 'coarse')
nfyb_sndx <- read_bed('data/derived_data/peak_sets/cnr_consensus/nfyb_sndx_consensus.bed') |> keepStandardChromosomes(pruning.mode = 'coarse')

atac_dmso <- read_bed('data/derived_data/peak_sets/atac_consensus/union_atac.bed')|> keepStandardChromosomes(pruning.mode = 'coarse')
atac_sndx <- read_bed('data/derived_data/peak_sets/atac_consensus/sndx_consensus_atac.bed') |> keepStandardChromosomes(pruning.mode = 'coarse')


nfyb_atac_dmso_list <- list(nfyb_dmso, atac_dmso)
nfyb_atac_sndx_list <- list(nfyb_sndx, atac_sndx) 
nfyb_occ_up_venn <- grid.grabExpr(
  makeVennDiagram(nfyb_atac_dmso_list,
                  NameOfPeaks = c('nfyb', 'atac'),
                  connectedPeaks = 'merge', 
                  cat.fontfamily = 'sans', 
                  cat.fontsize = 8,
                  fill = c(nfyb_cols, atac_cols) )
)

nfyb_occ_down_venn <- grid.grabExpr(
  makeVennDiagram(nfyb_atac_sndx_list,
                  NameOfPeaks = c('nfyb', 'atac'),
                  connectedPeaks = 'merge', 
                  cat.fontfamily = 'sans', 
                  cat.fontsize = 8,
                  fill = c(nfyb_cols, atac_cols))
)

pdf('plots/cnr/fig5_supp_nfyb_atac_venns.pdf', width = 8.5, height = 11)
pageCreate(width = 8.5, height = 11, showGuides = F)
plotGG(nfyb_occ_up_venn, width = 3.5, height = 3.5, x = 0.5, y = 0.5)
plotText(label = 'DMSO Plot', fontsize = 10, x = 4.5, y = 1)
plotGG(nfyb_occ_down_venn, width = 3.5, height = 3.5, x = 0.5, y = 4)
plotText(label = 'SNDX Plot', fontsize = 10, x = 4.5, y= 4.5)
dev.off()
