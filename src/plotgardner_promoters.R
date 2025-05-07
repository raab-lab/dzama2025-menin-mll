# Plot gardner plots of promtoers of genes of interest
# Takes a GRange and antibody to plot
# if range is specified, will use that
# otherwise will calculate


# valid antibodies
# MLL1
# MENIN
# H3K4me3
# H3K27me3
# IgG
# NFYB
# ATAC

suppressPackageStartupMessages({
  library(tidyverse)
  library(plotgardener)
  library(plyranges)
})


# leave these here for now
peak_colors <- c('#1b9e77', '#7570b3') # dmso / sndx
atac_colors <- c('grey30', 'grey60')
plot_out <- 'plots/figs/cnr/region_plots'
if (!dir.exists(plot_out)) { dir.create(plot_out, recursive = T)}

bw_path <- file.path('data/source_data/cnr/bw/')
atac_path <- file.path('data/source_data/atac/bw/')


# page size width = 4.5, height = 5, default.units = "inches"

plot_gene_track <- function(region, 
                            antibody, 
                            range_max=NULL, 
                            up = 20000, 
                            down = 20000, 
                            genome = 'hg38', 
                            plot_gene_track = TRUE, 
                            peak_colors = c('black', 'grey50'), 
                            offset=0, 
                            tracklabel = NULL) { 
     dmso_signal <- file.path(bw_path, paste0('groupHLF_DMSO_', antibody, '_mean.bw') )
     sndx_signal <- file.path(bw_path, paste0('groupHLF_SNDX_', antibody, '_mean.bw') )
     
     # 
     chr_str   <- as.character(seqnames(region)) 
     start_str <- as.numeric(start(region))  
     end_str   <- as.numeric(end(region))  
     if(offset == 0) { pageCreate(width = 4.5, height = 5, default.units = 'inches', showGuides  = F) }
     hk_scale_y <- -0.5 
     ka_scale_y <- 0
     if (is.null(range_max)) { 
       range_max <- calcSignalRange(data = list(dmso_signal, sndx_signal), 
                                         chrom = chr_str, 
                                         chromstart = start_str, 
                                         chromend = end_str, 
                                         assembly = genome, 
                                         negData = FALSE) 
       
       range_top <- max(range_max[2], 5)
       range_max <- c(0, range_top)
     } 
     else { range_max <- c(0,range_max) }
     
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
    plotText(label = tracklabel, fontcolor = "grey30", fontsize = 8,
             x = 0.85 , y = 0.5 + offset, just = c("center","center"), default.units = "inches",
             rot = 90)
    plotSegments(x0 = 0.95, x1 = 0.95, y0 = 0.3 +offset, y1 = 0.8+offset, linecolor = 'grey30', lwd = 3)

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
  
}

plot_atac_track <- function(region, 
                            range_max=NULL, 
                            up = 20000, 
                            down = 20000, 
                            genome = 'hg38', 
                            plot_gene_track = TRUE, 
                            peak_colors = c('black', 'grey50'), 
                            offset=0, 
                            tracklabel = NULL) { 
  dmso_signal <- file.path(atac_path, 'groupATAC_DMSO_mean.bw') 
  sndx_signal <- file.path(atac_path, 'groupATAC_SNDX_mean.bw') 
  
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
                       fill = atac_colors[1], linecolor = peak_colors[1],
                       y = 1+hk_scale_y+offset, height = 0.25, just = c("right", "bottom"))
  sndx_p <- plotSignal(data = sndx_signal, 
                       params = c(params_d, range_d),
                       fill = atac_colors[2], linecolor = peak_colors[2],
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

  
plot_multiple_tracks <- function (region) { 
  plot_atac_track(region = region, 
                  plot_gene_track = F, 
                  peak_colors = c('grey10', 'grey50')) 
  plot_gene_track(region = region,
                  antibody = 'H3K4me3',
                  plot_gene_track = F, 
                  peak_colors = c('#31a354', '#a1d99b'), 
                  offset = 0.75, 
                  tracklabel  = 'H3K4me3')
  plot_gene_track(region = region,
                  antibody = 'Menin',
                  plot_gene_track = F,
                  peak_colors = c('#3182bd','#9ecae1'),
                  offset = 1.5,
                  tracklabel = 'Menin')
  plot_gene_track(region = region, 
                  antibody = 'MLL1', 
                  plot_gene_track = F, 
                  peak_colors = c('#756bb1', '#bcbddc'), 
                  offset = 2.25, tracklabel = 'MLL1' )
  plot_gene_track(region = region, 
                  antibody = 'NFYB', 
                  plot_gene_track = F, 
                  peak_colors = c('#ec7014', '#fec44f'), 
                  offset = 3, tracklabel = 'NF-YB' )
  plot_gene_track(region = region, 
                  antibody = 'IgG', 
                  plot_gene_track = T, 
                  range_max = 4.1, 
                  peak_colors = c('grey10', 'grey50'), 
                  offset = 3.75, 
                  tracklabel = 'IgG') 
}

taok3 <- as_granges(data.frame(seqnames = 'chr12', start = 118350000, end =118385000 ))
suds3 <- as_granges(data.frame(seqnames = 'chr12', start = 118360000, end = 118380000))  
fgf5  <- as_granges(data.frame(seqnames = 'chr4',start = 80200000, end =  80320000))
lama2 <- as_granges(data.frame(seqnames = 'chr6', start = 128800000, end = 129000000))
epha5  <- as_granges(data.frame(seqnames = 'chr4', start = 65500000, end = 65700000))
bin3   <- as_granges(data.frame(seqnames = 'chr8', start= 22660000, end = 22680000))
hoxa   <- as_granges(data.frame(seqnames = 'chr7', start= 27080000, end = 27260000))
myc    <- as_granges(data.frame(seqnames = 'chr8', start = 127730000, end = 127750000))
p53   <- as_granges(data.frame(seqnames = 'chr17', start = 7660000, end = 7690000))
pak3 <- as_granges(data.frame(seqnames = 'chrX', start = 110800000, end = 111230000))
ptp4a1 <- as_granges(data.frame(seqnames = 'chr6', start = 63480000, end = 63620000))
ptp4a3 <- as_granges(data.frame(seqnames = 'chr8', start = 141385000, end = 141435000))
plot_multiple_tracks(suds3)
plot_multiple_tracks(epha5)
plot_multiple_tracks(myc)

pdf(file.path(plot_out, 'taok3_tracks.pdf'))
plot_multiple_tracks(taok3)
plotText(label = 'A', fontsize = 16, x = 0.25, y = 0) 
dev.off()

pdf(file.path(plot_out, 'fgf5_tracks.pdf'))
plot_multiple_tracks(fgf5)
plotText(label = 'B', fontsize = 16, x = 0.25, y = 0) 
dev.off()

pdf(file.path(plot_out, 'epha5_tracks.pdf'))
plot_multiple_tracks(epha5)
plotText(label = 'B', fontsize = 16, x = 0.25, y = 0) 
dev.off()

pdf(file.path(plot_out, 'lama2_tracks.pdf'))
plot_multiple_tracks(lama2)
plotText(label = 'B', fontsize = 16, x = 0.25, y = 0) 
dev.off()

pdf(file.path(plot_out, 'myc_tracks.pdf'))
plot_multiple_tracks(myc)
dev.off()

pdf(file.path(plot_out, 'pak3_tracks.pdf'))
plot_multiple_tracks(pak3)
dev.off() 

pdf(file.path(plot_out, 'p53_tracks.pdf'))
plot_multiple_tracks(p53)
dev.off()

pdf(file.path(plot_out, 'ptp4a1_tracks.pdf'))
plot_multiple_tracks(ptp4a1)
dev.off()

pdf(file.path(plot_out, 'ptp4a3_tracks.pdf'))
plot_multiple_tracks(ptp4a3)
dev.off()
