# cnr differential peaks
# Author: Jesse Raab
# Differential occupancy in cut-n-run at each factor's own peaks. 
suppressPackageStartupMessages({
  library(csaw)
  library(plyranges)
  library(tidyverse)
  library(DESeq2)
  library(sva)
})

# load union peaks - these are regions of interested for Differential Occupancy
# Look for differential peaks in MEN/MLL/H3K4me3
men <- read_bed('data/derived_data/peak_sets/cnr_consensus/menin_union.bed')
mll <- read_bed('data/derived_data/peak_sets/cnr_consensus/mll1_union.bed')
h3k4me3 <- read_bed('data/derived_data/peak_sets/cnr_consensus/k4me3_union.bed')


#Write a function that takes a peak set and antibody and returns the DE object
run_de <- function(peaks, ab) {
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
# CSAW Setup
  param = readParam(pe = 'both', dedup = F, minq = 10) # set this to 10 which matches my bam parameters, higher and we lose some numbers
  peak_counts <- regionCounts(bams, peaks, param = param)
  window_counts <- windowCounts(bams, bin = T, width = 10000) 
  # normalize on the regions we've counted

  # Need to set up some sort of sample sheet for comparison
  colData(peak_counts) <- DataFrame(left_join(as.data.frame(colData(peak_counts)), 
                                              ss, by = c('bam.files' ='filename') ) )
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
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  return(dds)
}

# current method controls for batch using airtable experiment ID - this works well for all but ASH2L which does not show signficant differences
men_des <- run_de(men, 'Menin')
print(resultsNames(men_des))
men_res <- lfcShrink(men_des, coef = 3, type = 'apeglm', format  = 'GRanges')

mll_des <- run_de(mll, 'MLL1')
resultsNames(mll_des)
mll_res <- lfcShrink(mll_des, coef = 2, type = 'apeglm', format = 'GRanges')

k4me3_des <- run_de(h3k4me3, 'H3K4me3')
resultsNames(k4me3_des)
k4me3_res <- lfcShrink(k4me3_des, coef = 4, type = 'apeglm', format = 'GRanges')  

# write out the objects
save(mll_des, mll_res, men_des, men_res, k4me3_des, k4me3_res, file = 'data/processed_data/cnr/menin_peaks_diff_all.Rda')
write_tsv(men_res |> as_tibble(), file = 'data/derived_data/cnr/menin_diff_menpeaks_res.tsv')
write_tsv(mll_res |> as_tibble(), file = 'data/derived_data/cnr/mll_diff_mllpeaks_res.tsv')
write_tsv(k4me3_res |> as_tibble(), file = 'data/derived_data/cnr/k4me3_diff_k4me3peaks_res.tsv')
