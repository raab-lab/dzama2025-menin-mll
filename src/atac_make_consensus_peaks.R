# Make consensus atac peaks using mspc
# must have dotnet available
# from console
# module load dotnet 
#make sure you are using bioconductor 3.18 or higher and rsmpc 1.8
suppressPackageStartupMessages({ 
  library(rmspc)
  library(tidyverse)
  library(plyranges)
})

# Output location
output_dir <- 'data/derived_data/peak_sets/atac_consensus/'
if(!dir.exists(output_dir) ) { dir.create(output_dir, recursive = T)} 
# Include an exlcuded region list
cnames <- c('seqnames', 'start', 'end', 'name', 'score', 'strand') 
# I'm excluding these anyway
exclude_regions <- read_tsv('data/external_data/hg38.Nordin.CandRblacklist_hg38.bed', col_names = cnames) |> as_granges()



# Takes pattern for identifying a group of peaks 
# Returns a named list of peaks for that antibody
process_peak_group <- function(pattern) {
  files <- list.files(path = 'data/source_data/atac/peaks/', pattern = pattern, full.names = T)  
  names <- basename(files)
  names <- lapply(names, function(x)  (unlist(str_split(x, pattern = '_') )  ) ) 
  names <- sapply(names, function(x) paste(x[3:6], collapse = '_') )
  peaks <- lapply(files, function(x) read_narrowpeaks(x) |> filter_by_non_overlaps(exclude_regions, maxgap=500) |> mutate(score = pValue) )  
  names(peaks) <- names
  return(peaks) 
}

# Get list of peak calls for each condition
dmso_peaks <- process_peak_group('.*DMSO.*.narrowPeak')
sndx_peaks <- process_peak_group('.*SNDX.*.narrowPeak')

# Create consensus set for each antibody_sample ##############################
# need at least 2 replicates for this to make sense 

create_consensus <- function(peak_list, c) {
  peak_consensus <- mspc(input = GRangesList(peak_list), 
                         replicateType = 'Bio', 
                         stringencyThreshold = 1e-12,
                         weakThreshold = 1e-5, 
                         keep = FALSE,GRanges = TRUE,
                         multipleIntersections = "Lowest",
                         c = c, alpha = 0.05)
  return(peak_consensus$GRangesObjects$ConsensusPeaks)
}

dmso_consensus <- create_consensus(dmso_peaks, c = 3) |> keepStandardChromosomes(pruning.mode = 'coarse')
sndx_consensus <- create_consensus(sndx_peaks, c = 3) |> keepStandardChromosomes(pruning.mode = 'coarse')

# make a union of the consensus peaks
all_peaks <- bind_ranges(dmso_consensus, sndx_consensus) |> reduce_ranges()

# write out peak sets
write_bed(dmso_consensus, file.path(output_dir, 'dmso_consensus_atac.bed'))
write_bed(sndx_consensus, file.path(output_dir, 'sndx_consensus_atac.bed'))
write_bed(all_peaks, file.path(output_dir, 'union_atac.bed'))
        