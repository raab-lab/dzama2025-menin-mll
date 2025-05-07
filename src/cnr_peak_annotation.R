# Annotate consensus regions with genomic features


suppressPackageStartupMessages({
  library(tidyverse)
  library(plyranges)
  library(ChIPpeakAnno)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene) 
  library(ChIPseeker)
  library(plotgardener)
})

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#import regions of interest - using consensus dmso peaks for each group 
menin <- read_bed('data/derived_data/peak_sets/cnr_consensus/menin_dmso_consensus.bed') |> keepStandardChromosomes(pruning.mode = 'coarse')
mll   <- read_bed('data/derived_data/peak_sets/cnr_consensus/mll1_dmso_consensus.bed')|> keepStandardChromosomes(pruning.mode = 'coarse')
h3k4me3 <- read_bed('data/derived_data/peak_sets/cnr_consensus/k4me3_dmso_consensus.bed')|> keepStandardChromosomes(pruning.mode = 'coarse')



menin_anno <- menin |> annotatePeak(TxDb = txdb, tssRegion = c(-2000, 200) )  |> as_granges()
mll_anno <- mll|> annotatePeak(TxDb = txdb, tssRegion = c(-2000, 200) )  |> as_granges()
h3k4me3_anno <- h3k4me3 |> annotatePeak(TxDb = txdb, tssRegion = c(-2000, 200)) |> as_granges()

# Need to clean up their annotations and combine promoter, intron, and exon categories oso there is only one of each. 
fix_category <- function(tib) {
  tib |> 
    mutate(new_cat = case_when(grepl ('Promoter', annotation) ~ 'Promoter', 
                               grepl ('Intron', annotation) ~ 'Intron', 
                               grepl ('Exon', annotation) ~ 'Exon', 
                               .default = as.character(annotation) ) )
  
  
}

menin_anno   <- menin_anno |> fix_category() |> mutate(peakset = 'Menin')
mll_anno     <- mll_anno |> fix_category() |> mutate(peakset = "MLL1")
h3k4me3_anno <- h3k4me3_anno |> fix_category() |> mutate(peakset = 'H3K4me3')


anno_gr <- bind_ranges(menin_anno, mll_anno, h3k4me3_anno)
anno_plot <- anno_gr |> 
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


# This makes the plot a better size
pdf('plots/cnr/region_anno_plot_all_dmso.pdf', width = 8.5, height = 11, )
pageCreate(height = 11, width = 8.5, showGuides = F)
legend <- cowplot::get_legend(anno_plot)
plotGG(anno_plot + theme(legend.position = 'none'),
       x = 0.5, y =0.5, width = 2.5, height = 3)
plotGG(legend, x = 3.25, y = 0.5, width = 1, height = 1, just = c('top', 'left')) 
pageGuideHide()
dev.off()


