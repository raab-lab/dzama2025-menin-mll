# Menin-MLL-ASH2L overlap DMSO vs SNDX
# included
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


menin_dmso <- read_bed('data/derived_data/peak_sets/cnr_consensus/menin_dmso_consensus.bed') |> keepStandardChromosomes(pruning.mode = 'coarse')
mll_dmso   <- read_bed('data/derived_data/peak_sets/cnr_consensus/mll1_dmso_consensus.bed')|> keepStandardChromosomes(pruning.mode = 'coarse')
k4me3_dmso  <- read_bed('data/derived_data/peak_sets/cnr_consensus/k4me3_dmso_consensus.bed') |> keepStandardChromosomes(pruning.mode = 'coarse')


menin_sndx <- read_bed('data/derived_data/peak_sets/cnr_consensus/menin_sndx_consensus.bed') |> keepStandardChromosomes(pruning.mode = 'coarse')
mll_sndx   <- read_bed('data/derived_data/peak_sets/cnr_consensus/mll1_sndx_consensus.bed')|> keepStandardChromosomes(pruning.mode = 'coarse')
k4me3_sndx  <- read_bed('data/derived_data/peak_sets/cnr_consensus/k4me3_sndx_consensus.bed') |> keepStandardChromosomes(pruning.mode = 'coarse')


fill <- c(men_cols, mll_cols)
dmso_list <- list(menin_dmso, mll_dmso)
sndx_list <- list(menin_sndx, mll_sndx)
dmso_plot <-  grid.grabExpr(makeVennDiagram(dmso_list, 
                                                   NameOfPeaks = c('Menin', 'MLL'), 
                                                   connectedPeaks = 'merge', 
                                                   main.fontfamily  = 'sans', 
                                                  fill = fill) )

sndx_plot <-  grid.grabExpr(makeVennDiagram(sndx_list, 
                                            NameOfPeaks = c('Menin', 'MLL'), 
                                            connectedPeaks = 'merge', 
                                            main.fontfamily  = 'sans', 
                                            fill = fill) )

men_list <- list(menin_dmso, menin_sndx)
mll_list <- list(mll_dmso, mll_sndx)
k4me3_list <- list(k4me3_dmso, k4me3_sndx)

men_drug_cols = c('#3182bd', '#9ecae1')
mll_drug_cols  = c('#756bb1', '#bcbddc')
k4me3_drug_cols = c('#31a354', '#a1d99b')


men_drug_plot <- grid.grabExpr(
  makeVennDiagram(men_list,
                  NameOfPeaks = c('DMSO', 'SNDX'),
                  connectedPeaks = 'merge', 
                  cat.fontfamily = 'sans', 
                  cat.fontsize = 8,
                  fill = men_drug_cols)
)

mll_drug_plot <- grid.grabExpr(
  makeVennDiagram(mll_list,
                  NameOfPeaks = c('DMSO', 'SNDX'),
                  connectedPeaks = 'merge', 
                  cat.fontfamily = 'sans', 
                  cat.fontsize = 8,
                  fill = mll_drug_cols)
)

k4me3_drug_plot <- grid.grabExpr(
  makeVennDiagram(k4me3_list,
                  NameOfPeaks = c('DMSO', 'SNDX'),
                  connectedPeaks = 'merge', 
                  cat.fontfamily = 'sans', 
                  cat.fontsize = 8,
                  fill = k4me3_drug_cols)
)


pdf('plots/figs/cnr/fig2_drug_comp_venn.pdf', width = 8.5, height = 11)
pageCreate(height = 11, width = 8.5, showGuide = F)
plotGG(men_drug_plot, x = 0.5, y =0.5, width =2 , height = 2)
plotGG(mll_drug_plot, x = 4, y = 0.5, width = 2, height = 2 )
plotGG(k4me3_drug_plot, x = 4, y = 4, width = 2, height = 2)

dev.off()


pdf('plots/figs/cnr/fig2_all_men_complex.pdf', width = 8.5, height = 11)
pageCreate(height = 11, width = 8.5, showGuide = F)
plotGG(dmso_plot, x = 0.5, y = 0.5, width = 3, height=3)
plotGG(sndx_plot, x = 0.5, y = 5, width = 3, height = 3)
plotText(x = 6, y = 0.5, label = 'DMSO')
plotText(x = 6, y = 5, label = 'SNDX')
dev.off()

# Do a version with just menin/mll/h3k4me3
men_mll_k4_dmso_peaks <- list(menin_dmso, mll_dmso, k4me3_dmso)
names(men_mll_k4_dmso_peaks) <- c('Menin', 'MLL1', 'H3K4me3')
men_mll_k4_sndx_peaks <- list(menin_sndx, mll_sndx, k4me3_sndx)
names(men_mll_k4_sndx_peaks) <- c('Menin', 'MLL1', 'H3K4me3')
plt_cols <- c(men_cols, mll_cols, k4_cols)
men_mll_k4_dmso <- grid.grabExpr(
  makeVennDiagram(men_mll_k4_dmso_peaks, 
                  NameOfPeaks = c('Menin', 'MLL1', 'H3K4me3'), 
                  cat.fontfamily = 'sans', 
                  connectedPeaks = 'merge',
                  cat.font.size = 8, 
                  fill = plt_cols) 
)


men_mll_k4_sndx <- grid.grabExpr(
  makeVennDiagram(men_mll_k4_sndx_peaks, 
                  NameOfPeaks = c('Menin', 'MLL1', 'H3K4me3'), 
                  cat.fontfamily = 'sans', 
                  connectedPeaks = 'merge',
                  cat.font.size = 8, 
                  fill = plt_cols) 
)
pdf('plots/figs/cnr/fig2_venn_men-mll-k4me3.pdf', width = 8.5, height = 11)
pageCreate(height = 11, width = 8.5, showGuide = F)
plotGG(men_mll_k4_dmso, x = 0.5, y = 0.5, width = 3, height=3)
plotGG(men_mll_k4_sndx, x = 0.5, y = 5, width = 3, height = 3)
plotText(label = 'DMSO', x = 1.75, y = 0.4)
plotText(label = 'SNDX', x = 1.75, y = 4.9)
dev.off()

