# Plot gardner plots of promtoers of genes of interest
# Takes a GRange and antibody to plot
# if range is specified, will use that
# otherwise will calculate




suppressPackageStartupMessages({
  library(tidyverse)
  library(plotgardener)
  library(plyranges)
})


# leave these here for now
peak_colors <- c('#1b9e77', '#7570b3') # dmso / sndx
plot_out <- 'plots/atac/region_plots'
if (!dir.exists(plot_out)) { dir.create(plot_out, recursive = T)}
bw_path <- file.path('data/source_data/atac/bw/') 
#bw_path_peak <- file.path('data/derived_data/group_norm_peaks/bw/')



# page size width = 4.5, height = 5, default.units = "inches"
#w_path <- bw_path_peak
plot_gene_track <- function(region, 
                            range_max=NULL, 
                            up = 20000, 
                            down = 20000, 
                            genome = 'hg38', 
                            plot_gene_track = TRUE, 
                            peak_colors = c('black', 'grey50'), 
                            offset=0, 
                            tracklabel = NULL) { 
  dmso_signal <- file.path(bw_path, 'groupDMSO_mean.bw') 
  sndx_signal <- file.path(bw_path, 'groupSNDX_mean.bw') 
  
  # 
  chr_str   <- as.character(seqnames(region)) 
  start_str <- as.numeric(start(region))   - up
  end_str   <- as.numeric(end(region))  + down 
  if(offset == 0) { pageCreate(width = 4.5, height = 5, default.units = 'inches', showGuides = F) }
  hk_scale_y <- -0.5 
  ka_scale_y <- 0
  if (is.null(range_max)) { 
    range_max <- 1.1* calcSignalRange(data = list(dmso_signal, sndx_signal), 
                                 chrom = chr_str, 
                                 chromstart = start_str, 
                                 chromend = end_str, 
                                 assembly = genome, 
                                 negData = FALSE) 
  } 
  else { range_max <- c(0, range_max)}
  
  range_label <- paste0( '[', 
                         as.character(round(range_max[1], 1) ), 
                         '-', 
                         as.character(round(range_max[2], 1) ),
                         ']')
  
  
  # create a parameter to hold that range  
  range_d <- pgParams(range = range_max, assembly = genome)
  
  
  params_d <- pgParams(chrom = chr_str, chromstart = start_str, chromend = end_str, 
                       assembly = genome,
                       x = 4, width = 3, default.units = "inches")
  
  ## Plot  signal
  dmso_p <- plotSignal(data = dmso_signal, 
                       params = c(params_d, range_d),
                       fill = peak_colors[1], linecolor = peak_colors[1],
                       y = 1+hk_scale_y+offset, height = 0.25, just = c("right", "bottom"))
  sndx_p <- plotSignal(data = sndx_signal, 
                       params = c(params_d, range_d),
                       fill = peak_colors[2], linecolor = peak_colors[2],
                       y = 1.25+hk_scale_y+offset, height = 0.25, just = c("right", "bottom"))
  
  ## Plot Track Labels
  plotText(label = 'DMSO', fontcolor = "grey30", fontsize = 8,
           x = 1, y = 0.3+offset, just = c("left","top"), default.units = "inches")
  plotText(label = "SNDX", fontcolor = "grey30", fontsize = 8,
           x = 1, y = 0.6+offset, just = c("left","top"), default.units = "inches")
  
  # Plot Track ranges
  plotText(label = range_label, fontcolor = "grey30", fontsize = 8,
           x = 3.9, y = 0.3+offset, just = c("left","top"), default.units = "inches")
  plotText(label = range_label, fontcolor = "grey30", fontsize = 8,
           x = 3.9, y = 0.6+offset, just = c("left","top"), default.units = "inches")
  
  # Plot track labels
  if (!is.null(tracklabel)) {
  plotText(label = tracklabel, fontcolor = "grey30", fontsize = 8,
           x = 0.5 , y = 0.5 + offset, just = c("center","center"), default.units = "inches",
           rot = 90)
  plotSegments(x0 = 0.8, x1 = 0.8, y0 = 0.2 +offset, y1 = 0.6+offset, linecolor = 'grey30', lwd = 3)
  }
  
  ## Plot gene track
  if (plot_gene_track) { 
    
    genes_imr <- plotGenes(params = params_d, stroke = 1, fontsize = 10,
                           strandLabels = FALSE,
                           y = 0.8 +offset, height = 0.4, just = c("right", "top"), 
                           fontcolor = 'grey30', fill = c('grey30', 'grey50') ) 
    ## Annotate genome label
    annoGenomeLabel(plot = genes_imr, params = params_d, 
                    scale = "Kb", fontsize = 10, digits = 0,
                    y = 1.3+offset, just = c("right", "top"))
  } 
  ## Hide page guides
  if (offset == 0) {  pageGuideHide() } 
}





taok3 <- as_granges(data.frame(seqnames = 'chr12', start = 118350000, end =118385000 ))
suds3 <- as_granges(data.frame(seqnames = 'chr12', start = 118360000, end = 118380000))  
fgf5  <- as_granges(data.frame(seqnames = 'chr4',start = 80200000, end =  80320000))
lama2 <- as_granges(data.frame(seqnames = 'chr6', start = 128800000, end = 129000000))
epha5  <- as_granges(data.frame(seqnames = 'chr4', start = 65500000, end = 65700000))
bin3   <- as_granges(data.frame(seqnames = 'chr8', start= 22660000, end = 22680000))
hoxa   <- as_granges(data.frame(seqnames = 'chr7', start= 27080000, end = 27260000))


# Pick somer anges do tplot
up_ranges <- read_tsv('data/derived_data/atac/gain_access.tsv') |> as_granges() |> keepStandardChromosomes(pruning = 'coarse')
down_ranges <- read_tsv('data/derived_data/atac/lost_access.tsv') |> as_granges() |> keepStandardChromosomes(pruning = 'coarse')

up_ranges_best <- up_ranges[order(up_ranges$padj)]
down_ranges_best <- down_ranges[order(down_ranges$padj)]

for (i in 1:10) {
  name <- paste(seqnames(up_ranges_best[i]), start(up_ranges_best[i]), end(up_ranges_best[i]), sep = '_')
  pdf(file.path(plot_out, paste0(name, '_up.pdf')))
  plot_gene_track(up_ranges_best[i], up = 20000, down = 20000) 
  dev.off()
}



for (i in 1:10) {
  name <- paste(seqnames(down_ranges_best[i]), start(down_ranges_best[i]), end(down_ranges_best[i]), sep = '_')
  
  pdf(file.path(plot_out, paste0(name, '_down.pdf')))
    plot_gene_track(down_ranges_best[i], up = 20000, down = 20000) 
  dev.off()
}

plot_gene_track(down_ranges_best[5], up = 300000, down = 300000)
