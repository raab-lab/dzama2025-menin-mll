#  differential accessibility  of NFYB at changed ATAC peaks
# Author: Jesse RAab
# Date: 2025-01-20

# something off with this still - re-running some background stuff -2025-01-21
suppressPackageStartupMessages({
  library(csaw)
  library(plyranges)
  library(tidyverse)
  library(DESeq2)
})
outdir <- 'data/derived_data/cnr/nfyb_analysis/'
dir.create(outdir, recursive = T)

# load union peaks - these are regions of interested for Differential Occupancy
peaks <- read_bed(file.path('data/derived_data/peak_sets/atac_consensus/union_atac.bed'))
nfyb_peaks <- read_bed(file.path('data/derived_data/peak_sets/cnr_consensus/nfyb_union.bed'))
tead_peaks <- read_bed(file.path('data/derived_data/peak_sets/cnr_consensus/tead4_union.bed'))



run_de_cnr <- function(peaks, ab) {
  # get sample_sheet details
  samples <- read_csv('full_cnr_samples.csv') |> janitor::clean_names()
  samples <- samples |> mutate(sample_number = as.character(sample_number))
  samples <- samples |> mutate(experiment_i_id = as.character(experiment_i_id))
  samples <- samples |> dplyr::select(sample_number, experiment_i_id)
  
  # get bam files
  bams <- list.files(path = 'data/source_data/cnr/bams/filtered/', pattern = paste0(ab, '_*.*.bam$'), full.names = T)
  print(bams)
  
  # Sample Sheet
  ss <- data.frame(filename = bams, names = basename(bams))
  
  ss <- ss |> 
    mutate(short_name = str_replace(names, '.aligned.bam', '') ) |> 
    separate(short_name, into = c('airtable_id', 'NGSID', 'cell_line', 'antibody', 'condition', 'rep'), sep = '_' ) 
  ss$condition <- toupper(ss$condition)
  ss <- ss |> left_join(samples, by = c('airtable_id' = 'sample_number'))
  n_exp <- length(unique(ss$experiment_i_id))
  print(n_exp)
  print(ss)
  # CSAW Setup
  param = readParam(pe = 'both', dedup = F, minq = 10) # set this to 10 which matches my bam parameters, higher and we lose some numbers
  peak_counts <- regionCounts(bams, peaks, param = param)
  window_counts <- windowCounts(bams, bin = T, width = 10000) 
  # normalize on the regions we've counted
  
  # Need to set up some sort of sample sheet for comparison
  colData(peak_counts) <- DataFrame(left_join(as.data.frame(colData(peak_counts)), 
                                              ss, by = c('bam.files' = 'filename') ) )
  colData(window_counts) <- DataFrame(left_join(as.data.frame(colData(window_counts)), 
                                                ss, by= c('bam.files' = 'filename') ) )
  
  # 
  if (n_exp > 1) {
    dds <- DESeqDataSet(peak_counts, design = ~  experiment_i_id + condition) 
    dds$condition <- relevel(dds$condition, ref = "DMSO")
    norm_dds <- DESeqDataSet(window_counts, design = ~ experiment_i_id + condition)
  }
  else { 
    dds <- DESeqDataSet(peak_counts, design = ~   condition) 
    dds$condition <- relevel(dds$condition, ref = "DMSO")
    norm_dds <- DESeqDataSet(window_counts, design = ~  condition)
  }
  norm_dds <- estimateSizeFactors(norm_dds)
  sf <- sizeFactors(norm_dds)
  sizeFactors(dds) <- sf
  # comment this next line out if you want bins - TODO: rewrite function to make this a flag
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  return(dds)
}

# run at open chromatin regions
nfyb_de  <- run_de_cnr(peaks, ab = 'NFYB')

nfyb_res <- lfcShrink(nfyb_de, coef = 3, type = 'apeglm', format = 'GRanges')
mll

# Run at NFYB peaks
nfyb_de2 <- run_de_cnr(nfyb_peaks, ab = 'NFYB')
nfyb_res2 <- lfcShrink(nfyb_de2, coef = 3, type = 'apeglm', format = 'GRanges')

# need to fix the function so it normalizes on windows if you run this 
#mll_de  <- run_de_cnr(nfyb_peaks, ab = 'MLL1')
#mll_res <- lfcShrink(mll_de, coef = 2, type = 'apeglm', format = 'GRanges')
#mll_res |> as_tibble() |> ggplot(aes(x = log2FoldChange, y = -log10(padj))) + geom_point()
#mll_res |> as_tibble() |> ggplot(aes(y = log2FoldChange, x = log2(baseMean))) + geom_point()
#
#mll_res |> as_tibble() |> 
  #filter(padj < 0.1) |> 
  #filter(log2FoldChange > 0)
#
# Write out results
write_tsv(nfyb_res |> as.data.frame(), file = file.path(outdir, 'nfyb_diff_at_atac.tsv'))
write_tsv(nfyb_res2 |> as.data.frame(), file = file.path(outdir, 'nfyb_diff_at_nfyb_peaks.tsv'))
save(nfyb_de, nfyb_res, nfyb_de2, nfyb_res2, file = file.path(outdir, 'nfyb_analysis_obj.Rda')) 

