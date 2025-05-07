# Plotting of differential MENIN/MLL/H3K4 peaks
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(plyranges)
  library(ChIPpeakAnno)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene) 
  library(ChIPseeker)
  library(plotgardener)
  library(grid)
})
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
k4_cols <-  '#31a354'
men_cols <-  '#3182bd'
mll_cols <-  '#756bb1'
ash2l_cols <- '#fdc086'

dir.create('plots/figs/cnr', recursive = T)

# Need to clean up their annotations and combine promoter, intron, and exon categories oso there is only one of each. 
# helper function I wrote in atac_peak_annotations.R
fix_category <- function(tib) {
  tib |> 
    mutate(new_cat = case_when(grepl ('Promoter', annotation) ~ 'Promoter', 
                               grepl ('Intron', annotation) ~ 'Intron', 
                               grepl ('Exon', annotation) ~ 'Exon', 
                               .default = as.character(annotation) ) )
  
  
}

# load in the tables
men <- read_tsv('data/derived_data/cnr/menin_diff_menpeaks_res.tsv') |> as_granges() |> keepStandardChromosomes(pruning.mode = 'coarse')
mll <- read_tsv('data/derived_data/cnr/mll_diff_menpeaks_res.tsv') |> as_granges()|> keepStandardChromosomes(pruning.mode = 'coarse')
k4me3  <- read_tsv('data/derived_data/cnr/k4me3_diff_menpeaks_res.tsv') |> as_granges()|> keepStandardChromosomes(pruning.mode = 'coarse')
# no ash2l peaks diffa

table(men$padj < 0.05)
table(mll$padj < 0.05)
table(k4me3$padj < 0.05)

men_anno <- men |>  annotatePeak(TxDb = txdb, tssRegion = c(-2000, 2000)) |> as_granges()
mll_anno <- mll |>  annotatePeak(TxDb = txdb, tssRegion = c(-2000, 2000)) |> as_granges()
k4me3_anno <- k4me3 |>  annotatePeak(TxDb = txdb, tssRegion = c(-2000, 2000)) |> as_granges() 

men_anno

men_anno <- men_anno |> fix_category() |> mutate(peakset = 'Menin')
mll_anno <- mll_anno |> fix_category() |> mutate(peakset = 'MLL')
k4me3_anno <- k4me3_anno |> fix_category() |> mutate(peakset = 'H3K4me3')


anno_gr <- bind_ranges(men_anno, mll_anno, k4me3_anno)
anno_most_plot <- anno_gr |> 
  filter(padj < 0.05) |> 
  filter(log2FoldChange < 0) |> 
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
        legend.title = element_blank(), 
        plot.title = element_text(size = 10)) + 
  scale_fill_brewer(palette = 'Set3') + 
  ylab('Frequency') + xlab('') + 
  ggtitle("Signficantly Changed Peaks")

anno_gr
# Look at overlap between three peak sets of diffeerntial binding

# Look at overlap with gain via Venn Diagram
men_lost <- men_anno |> filter(padj < 0.1, log2FoldChange < 0)
mll_lost <- mll_anno |> filter(padj < 0.1, log2FoldChange < 0)
k4me3_lost <- k4me3_anno |> filter(padj < 0.1, log2FoldChange < 0) 

lost_list <- list(men_lost, mll_lost, k4me3_lost)
fill = c(men_cols, mll_cols, k4_cols)
lost_most_plot <- grid.grabExpr(makeVennDiagram(lost_list, 
                                                NameOfPeaks = c('Menin', 'MLL', 'H3K4me3'), 
                                                connectedPeaks = 'merge', 
                                                main.fontfamily  = 'sans', 
                                                fill = fill) )
                                
                                
pdf('plots/figs/cnr/changed_peaks_anno.pdf', width = 8.5, height = 11)
pageCreate(height = 11, width = 8.5, showGuide = F)
plotText(x = 0.5, y = 0.1, label = 'Figure 2 Bar Graphs')
legend <- cowplot::get_legend(anno_most_plot)
plotGG(anno_most_plot + theme(legend.position = 'none'),
       x = 0.5, y =0.5, width = 2.5, height = 3)
plotText(x = 4, y = 0.5, label = 'Genomic Distribution All Peaks')
plotGG(legend, x = 3.25, y = 0.5, width = 1, height = 1, just = c('top', 'left')) 
plotGG(lost_most_plot, x = 0.5, y =4, width =3 , height = 3)
plotText(x = 4, y=  4, label = 'Overlap of Lost Peaks - supplement?')
dev.off()
                                                
#anno_gr |> 
#  as_tibble() |> 
#  ggplot(aes(x = log2FoldChange, y = -log10(padj))) + 
#  geom_point() +
#  facet_wrap(~peakset)+ 
#  geom_hline(yintercept = -log10(0.05), color = 'red', linetype = 'dashed')
#
