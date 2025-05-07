# distribution of ATAC peaks# Annotate consensus regions with genomic features
# limiting this to the most gained/lost (LFC > 0.58)
# Overlap with menin/mll/k4me3 peaks


suppressPackageStartupMessages({
  library(tidyverse)
  library(plyranges)
  library(ChIPpeakAnno)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene) 
  library(ChIPseeker)
  library(plotgardener)
})

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# ATAC regions
gain <- read_bed('data/derived_data/atac/gain_access_most.bed') |> keepStandardChromosomes(pruning = 'coarse')
lost <- read_bed('data/derived_data/atac/lost_access_most.bed') |> keepStandardChromosomes(pruning = 'coarse')
gain_all <- read_bed('data/derived_data/atac/gain_access.bed') |> keepStandardChromosomes(pruning = 'coarse')
lost_all <- read_bed('data/derived_data/atac/lost_access.bed') |> keepStandardChromosomes(pruning = 'coarse')
#import regions of interest - using consensus dmso peaks for each group 
menin <- read_bed('data/derived_data/peak_sets/cnr_consensus/menin_dmso_consensus.bed') |> keepStandardChromosomes(pruning = 'coarse')
mll   <- read_bed('data/derived_data/peak_sets/cnr_consensus/mll1_dmso_consensus.bed') |> keepStandardChromosomes(pruning = 'coarse')
h3k4me3 <- read_bed('data/derived_data/peak_sets/cnr_consensus/k4me3_dmso_consensus.bed') |> keepStandardChromosomes(pruning = 'coarse')

gain_anno <- gain |> annotatePeak(TxDb = txdb, tssRegion = c(-2000, 2000)) |> as_granges()
lost_anno <- lost |> annotatePeak(TxDb = txdb, tssRegion = c(-2000, 2000)) |> as_granges()


gain_all_anno <- gain_all |> annotatePeak(TxDb = txdb, tssRegion = c(-2000, 2000)) |> as_granges()
lost_all_anno <- lost_all |> annotatePeak(TxDb = txdb, tssRegion = c(-2000, 2000)) |> as_granges()


# Need to clean up their annotations and combine promoter, intron, and exon categories oso there is only one of each. 
fix_category <- function(tib) {
  tib |> 
    mutate(new_cat = case_when(grepl ('Promoter', annotation) ~ 'Promoter', 
                               grepl ('Intron', annotation) ~ 'Intron', 
                               grepl ('Exon', annotation) ~ 'Exon', 
                               .default = as.character(annotation) ) )
  
  
}

gain_anno <- gain_anno |> fix_category() |> mutate(peakset = 'Gained')
lost_anno <- lost_anno |> fix_category() |> mutate(peakset = 'Lost')

gain_all_anno <- gain_all_anno |> fix_category() |> mutate(peakset = 'Gained')
lost_all_anno <- lost_all_anno |> fix_category() |> mutate(peakset = 'Lost')

anno_gr <- bind_ranges(gain_anno, lost_anno) 
anno_most_plot <- anno_gr |> 
  as_tibble() |> 
  group_by(peakset, new_cat) |> 
  summarise(num = n())  |> 
  mutate(freq = num/sum(num) ) |> 
  ggplot(aes(x = peakset, y = freq, fill = new_cat)) + 
  geom_col(color=  'grey30')+ 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8, color = 'black'), 
        axis.text.y = element_text(size = 8, color = 'black'), 
        axis.title.y = element_text(size = 10, color = 'black'), 
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, 'cm'), 
        legend.title = element_blank() ) + 
  scale_fill_brewer(palette = 'Set3') + 
  ylab('Frequency') + xlab('')

anno_all_gr <- bind_ranges(gain_all_anno, lost_all_anno) 
anno_all_plot <- anno_all_gr |> 
  as_tibble() |> 
  group_by(peakset, new_cat) |> 
  summarise(num = n())  |> 
  mutate(freq = num/sum(num) ) |> 
  ggplot(aes(x = peakset, y = freq, fill = new_cat)) + 
  geom_col(color=  'grey30')+ 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8, color = 'black'), 
        axis.text.y = element_text(size = 8, color = 'black'), 
        axis.title.y = element_text(size = 10, color = 'black'), 
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, 'cm'), 
        legend.title = element_blank() ) + 
  scale_fill_brewer(palette = 'Set3') + 
  ylab('Frequency') + xlab('')



# Look at gained an lost overlap with menin/mll
# combined menin/mll
complex <- bind_ranges(menin, mll) |> reduce_ranges()
complex

# Look at overlap with gain via Venn Diagram
gain_list <- list(gain, complex, h3k4me3)
lost_list <- list(lost, complex, h3k4me3)
gain_all_list  <- list(gain_all, complex, h3k4me3)
lost_all_list  <- list(lost_all, complex, h3k4me3)

# get a list of lost peaks with menin
lost_men <- lost |> filter_by_overlaps(complex)
lost_nomen <- lost |> filter_by_non_overlaps(complex)


lost_men |> as_tibble() |> sample_n(10)
lost_nomen |> as_tibble() |> sample_n(10)

# MAke Plots
gain_most_plot <-grid.grabExpr(makeVennDiagram(gain_list, 
                                               NameOfPeaks = c('Gained', 'Menin.MLL', 'H3K4me3'),
                                               connectedPeaks = 'merge', 
                                               main.fontfamily = 'sans' , 
                                               fill = c('steelblue', 'grey70', '#66c2a4')) )
lost_most_plot <- grid.grabExpr(makeVennDiagram(lost_list, 
                                                 NameOfPeaks = c('Lost', 'Menin.MLL', 'H3K4me3'), 
                                                 connectedPeaks = 'merge', 
                                                 main.fontfamily  = 'sans', 
                                                 fill = c('#bcbddc', 'grey70', '#66c2a4') ) )
# Supplement of peaks without the FC filter
gain_plot <-grid.grabExpr(makeVennDiagram(gain_all_list, 
                                               NameOfPeaks = c('Gained', 'Menin.MLL', 'H3K4me3'),
                                               connectedPeaks = 'merge', 
                                               main.fontfamily = 'sans' , 
                                               fill = c('steelblue', 'grey70', '#66c2a4')) )
lost_plot <- grid.grabExpr(makeVennDiagram(lost_all_list, 
                                                NameOfPeaks = c('Lost', 'Menin.MLL', 'H3K4me3'), 
                                                connectedPeaks = 'merge', 
                                                main.fontfamily  = 'sans', 
                                                fill = c('#bcbddc', 'grey70', '#66c2a4') ) )
# Supplement of peaks without the FC filter
# This makes the plot a better size
# Plot MOST Changed ATAC
pdf('plots/atac/most_changed_regions_anno_plot.pdf', width = 8.5, height = 11, )
pageCreate(height = 11, width = 8.5, showGuide = F)
legend <- cowplot::get_legend(anno_most_plot)
plotGG(anno_most_plot + theme(legend.position = 'none'),
       x = 0.5, y =0.5, width = 2.5, height = 3)
plotGG(legend, x = 3.25, y = 0.5, width = 1, height = 1, just = c('top', 'left')) 
plotGG(gain_most_plot, x= 5, y = 4, width = 3, height = 3) 
plotGG(lost_most_plot, x = 0.5, y =4, width =3 , height = 3)
dev.off()

#Plot all atac that passed padj filter

pdf('plots/atac/all_changed_regions_anno_plot.pdf', width = 8.5, height = 11, )
pageCreate(height = 11, width = 8.5, showGuide = F)
legend <- cowplot::get_legend(anno_all_plot)
plotGG(anno_all_plot + theme(legend.position = 'none'),
       x = 0.5, y =0.5, width = 2.5, height = 3)
plotGG(legend, x = 3.25, y = 0.5, width = 1, height = 1, just = c('top', 'left')) 
plotGG(gain_plot, x= 5, y = 4, width = 3, height = 3) 
plotGG(lost_plot, x = 0.5, y =4, width =3 , height = 3)
dev.off()