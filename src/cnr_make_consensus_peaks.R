# MSPC peak calls
# Generate consensus for each group
# Generate Union of all calls

# Try MSPC for defining best peak set

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
output_dir <- 'data/derived_data/peak_sets/cnr_consensus/'
dir.create(output_dir, recursive = T)
# Include an exlcuded region list
cnames <- c('seqnames', 'start', 'end', 'name', 'score', 'strand') 
exclude_regions <- read_tsv('data/external_data/hg38.Nordin.CandRblacklist_hg38.bed', col_names = cnames) |> as_granges()
encode_deny     <- read_tsv('data/external_data/encode_hg38_exclusion_list_encsr636FF.gz', col_names = c('seqnames', 'start', 'end')) |> as_granges()
exclude_regions <- bind_ranges(exclude_regions, encode_deny)

# Load peaks for each set as a list ###########################################

# Takes pattern for identifying a group of peaks 
# Returns a named list of peaks for that antibody
process_peak_group <- function(pattern) {
  files <- list.files(path = 'data/source_data/cnr/peaks/', pattern = pattern, full.names = T, ignore.case = T)  
  names <- basename(files)
  names <- lapply(names, function(x)  (unlist(str_split(x, pattern = '_') )  ) ) 
  names <- sapply(names, function(x) paste(x[3:6], collapse = '_') )
  peaks <- lapply(files, function(x) read_narrowpeaks(x) |> filter_by_non_overlaps(exclude_regions, maxgap=500) |> mutate(score = pValue) )  
  names(peaks) <- names
  return(peaks) 
}

# Get list of peak calls for each antibody/condition
k4_dmso_peaks <- process_peak_group('.*H3K4me3_DMSO.*.narrowPeak')
k4_sndx_peaks <- process_peak_group('.*H3K4me3_SNDX.*.narrowPeak')

menin_dmso_peaks <- process_peak_group('.*Menin_DMSO.*.narrowPeak')
menin_sndx_peaks <- process_peak_group('.*Menin_SNDX.*.narrowPeak')

mll_dmso_peaks <- process_peak_group('.*MLL1_DMSO.*.narrowPeak')
mll_sndx_peaks <- process_peak_group('.*MLL1_SNDX.*.narrowPeak')

tead_dmso_peaks <- process_peak_group('.*TEAD4_DMSO.*.narrowPeak')
tead_sndx_peaks <- process_peak_group('.*TEAD4_SNDX.*.narrowPeak')

nfyb_dmso_peaks <- process_peak_group('.*NFYB_DMSO.*.narrowPeak')
nfyb_sndx_peaks <- process_peak_group('.*NFYB_SNDX.*.narrowPeak')


# Create consensus set for each antibody_sample ##############################
# need at least 2 replicates for this to make sense 

create_consensus <- function(peak_list, c) {
  peak_consensus <- mspc(input = GRangesList(peak_list), 
                         replicateType = 'Bio', 
                         stringencyThreshold = 1e-8,
                         weakThreshold = 1e-4, 
                         keep = FALSE,GRanges = TRUE,
                         multipleIntersections = "Lowest",
                         c = c, alpha = 0.05)
  return(peak_consensus$GRangesObjects$ConsensusPeaks)
}

k4_dmso_consensus <- create_consensus(k4_dmso_peaks, c = 2) |> keepStandardChromosomes(pruning.mode = 'coarse')
k4_sndx_consensus <- create_consensus(k4_sndx_peaks, c = 2) |> keepStandardChromosomes(pruning.mode = 'coarse')

menin_dmso_consensus <- create_consensus(menin_dmso_peaks, c = 1) |> keepStandardChromosomes(pruning.mode = 'coarse')
menin_sndx_consensus <- create_consensus(menin_sndx_peaks, c = 1) |> keepStandardChromosomes(pruning.mode = 'coarse')

mll_dmso_consensus <- create_consensus(mll_dmso_peaks, c = 1) |> keepStandardChromosomes(pruning.mode = 'coarse')
mll_sndx_consensus <- create_consensus(mll_sndx_peaks, c = 1) |> keepStandardChromosomes(pruning.mode = 'coarse')

tead_dmso_consensus <- create_consensus(tead_dmso_peaks, c = 1) |> keepStandardChromosomes(pruning.mode = 'coarse')
tead_sndx_consensus <- create_consensus(tead_sndx_peaks, c = 1) |> keepStandardChromosomes(pruning.mode = 'coarse')

nfyb_dmso_consensus <- create_consensus(nfyb_dmso_peaks, c = 1) |> keepStandardChromosomes(pruning.mode = 'coarse')
nfyb_sndx_consensus <- create_consensus(nfyb_sndx_peaks, c = 1) |> keepStandardChromosomes(pruning.mode = 'coarse')

# Create a union set of these three data sets for statistical testing
union_k4_peaks <- bind_ranges(k4_dmso_consensus, k4_sndx_consensus) |> reduce_ranges() 
union_menin_peaks <- bind_ranges(menin_dmso_consensus, menin_sndx_consensus) |> reduce_ranges() 
union_mll_peaks <- bind_ranges(mll_dmso_consensus, mll_sndx_consensus) |> reduce_ranges()
union_tead_peaks <- bind_ranges(tead_dmso_consensus, tead_sndx_consensus) |> reduce_ranges()
union_nfyb_peaks <- bind_ranges(nfyb_dmso_consensus, nfyb_sndx_consensus) |> reduce_ranges()

# write these out for further analysis
write_bed(k4_dmso_consensus, file = file.path(output_dir, 'k4me3_dmso_consensus.bed'))
write_bed(k4_sndx_consensus, file = file.path(output_dir, 'k4me3_sndx_consensus.bed'))
write_bed(menin_dmso_consensus, file = file.path(output_dir, 'menin_dmso_consensus.bed'))
write_bed(menin_sndx_consensus, file = file.path(output_dir, 'menin_sndx_consensus.bed'))
write_bed(mll_dmso_consensus, file = file.path(output_dir, 'mll1_dmso_consensus.bed'))
write_bed(mll_sndx_consensus, file = file.path(output_dir, 'mll1_sndx_consensus.bed'))
write_bed(tead_dmso_consensus, file = file.path(output_dir, 'tead4_dmso_consensus.bed'))
write_bed(tead_sndx_consensus, file = file.path(output_dir, 'tead4_sndx_consensus.bed'))
write_bed(nfyb_dmso_consensus, file = file.path(output_dir, 'nfyb_dmso_consensus.bed'))
write_bed(nfyb_sndx_consensus, file = file.path(output_dir, 'nfyb_sndx_consensus.bed'))
# write the union sets for each antibody

write_bed(union_k4_peaks, file = file.path(output_dir, 'k4me3_union.bed'))
write_bed(union_menin_peaks, file = file.path(output_dir, 'menin_union.bed'))
write_bed(union_mll_peaks, file = file.path(output_dir, 'mll1_union.bed'))
write_bed(union_tead_peaks, file = file.path(output_dir, 'tead4_union.bed'))
write_bed(union_nfyb_peaks, file = file.path(output_dir, 'nfyb_union.bed'))


