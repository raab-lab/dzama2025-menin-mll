# supplemental figure of fragment length for atac vs non-atac  nfyb peaks
# Author: Jesse Raab 
# Date: 2025-02-25

# Excluded from paper - not useful

library(tidyverse)
library(plyranges)
library(VplotR)
library(plotgardener)
library(memes)
library(universalmotif)
library(grid)
hg38_genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
meme_path <- '/nas/longleaf/rhel8/apps/meme/5.5.2/bin/' # make sure to module load this

nfyb_dar <- read_tsv('data/derived_data/cnr/nfyb_analysis/nfyb_diff_at_nfyb_peaks.tsv') |> as_granges()

atac_up <- read_tsv('data/derived_data/atac/gain_access.tsv') |> as_granges()
atac_down <- read_tsv('data/derived_data/atac/lost_access.tsv') |> as_granges()
atac_unchanged <- read_tsv('data/derived_data/atac/unchanged.tsv') |> as_granges()
atac_all <- read_tsv('data/derived_data/atac/all_data.tsv') |> as_granges()

#nfyb_plus_up
nfyb_atac_up <- nfyb_dar |> 
  filter_by_overlaps(atac_up)|> 
  mutate(group = 'atac_up')

nfyb_atac_down <- nfyb_dar |> 
  filter_by_overlaps(atac_down)|> 
  mutate(group = 'atac_down')

nfyb_atac_unchanged <- nfyb_dar |> 
  filter_by_overlaps(atac_unchanged) |> 
  mutate(group = 'atac_unchanged')


nfyb_atac_all <- c(nfyb_atac_up, nfyb_atac_down, nfyb_atac_unchanged)
# supplemental Figure of nfyb changes by ATAC category
nfyb_volcano <- nfyb_atac_all |> 
  as_tibble() |> 
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = group)) +
  geom_point() + 
  scale_color_manual(values = c('steelblue','grey50', 'red3')) + 
  geom_hline(aes(yintercept = -log10(0.05)), linetype = 'dashed', color = 'grey30')+
  theme_bw() + 
  facet_wrap(~group)+
  theme(panel.grid = element_blank())



# Get a list of nfyb Down Plus ATAC Down
both_up <- nfyb_atac_up |> filter(padj < 0.05) |> filter(log2FoldChange > 0)
# Get a list of nfyb Up Plus ATAC Up
both_down <- nfyb_atac_down |> filter(padj < 0.05) |> filter(log2FoldChange < 0)

nfyb_diff_atac_unchanged <- nfyb_dar |> filter(padj < 0.05) |> filter_by_overlaps(atac_all |> filter(padj > 0.05) ) 

write_bed(both_up, file = 'data/derived_data/cnr/nfyb_analysis/nfyb_and_atac_up.bed')
write_bed(both_down, file = 'data/derived_data/cnr/nfyb_analysis/nfyb_and_atac_down.bed')
write_bed(nfyb_diff_atac_unchanged, file = 'data/derived_data/cnr/nfyb_analysis/nfyb_diff_atac_unchanged.bed')

# relatively small lists

# which DARs do NOT overlap ATAC
nfyb_dar_no_atac <- nfyb_dar |> filter(padj < 0.05) |> filter_by_non_overlaps(atac_all)
nfyb_not_dar_no_atac <- nfyb_dar |> filter(padj > 0.05) |> filter_by_non_overlaps(atac_all)
nfyb_dar_atac    <- nfyb_dar |> filter(padj < 0.05) |> filter_by_overlaps(atac_all)
nfyb_not_dar_atac <- nfyb_dar |> filter(padj > 0.05) |> filter_by_overlaps(atac_all)

  

#lets not look at differential
nfyb_plus_atac <- nfyb_dar |> filter_by_overlaps(atac_all)
nfyb_no_atac   <- nfyb_dar |> filter_by_non_overlaps(atac_all)

nfyb_plus_atac <- nfyb_plus_atac |> anchor_center() |> resize(500)
nfyb_no_atac <-nfyb_no_atac |> anchor_center() |> resize(500)
nfyb_no_atac |> as_tibble() |> sample_n(10)

# supplemental figure
# look at read count distribution at nfyb +/- atac
nfyb_plus_atac <- nfyb_plus_atac |> mutate(group = 'with_atac')
nfyb_no_atac   <- nfyb_no_atac |> mutate(group = 'no_atac')
nfyb_atac_combined <- c(nfyb_plus_atac |> unanchor(), nfyb_no_atac |> unanchor())
plt_coverage_dist <- nfyb_atac_combined |> 
  as_tibble() |> 
  ggplot(aes(x = log2(baseMean), fill = group))  + geom_density(color = 'grey20', alpha = 0.5)  +
  theme_bw() + 
  theme(axis.title = element_text(size = 7, family = 'sans') , 
        axis.text = element_text(size = 6, family = 'sans') )  + 
  xlab('Log2 NFYB Coverage') + ylab ('Density') + 
  scale_fill_manual (values = c('grey40', 'steelblue')) 

#  - ignore this analysis for now -----------------------------------------
nfyb_dar_sig <- nfyb_dar |> filter(padj < 0.05)  
#nfyb_dar_sig_atac <- nfyb_dar_sig |> filter(baseMean > 40) |> filter_by_overlaps(atac_all) |> achor_center() |> resize(500)
#nfyb_dar_sig_no_atac <- nfyb_dar_sig |> filter(baseMean > 40) |> filter_by_non_overlaps(atac_all) |> anchor_center() |> resize(500)
#nfyb_dmso <- list.files('data/source_data/cnr/bams/filtered/', pattern = ".*NFYB_DMSO.*.bam$", full.names = T)
#nfyb_sndx <- list.files('data/source_data/cnr/bams/filtered/', pattern = ".NFYB_SNDX.*.bam$", full.names = T)
#
#
#
#
#dmso_frags <- importPEBamFiles(nfyb_dmso, shift_ATAC_fragments = F, where = c(nfyb_dar_sig_atac, nfyb_dar_sig_no_atac))
#sdnx_frags <- importPEBamFiles(nfyb_sndx, shift_ATAC_fragments = F, where = c(nfyb_dar_sig_atac, nfyb_dar_sig_no_atac))
#
#
#dmso_frags <- do.call('c', dmso_frags)
#sdnx_frags <- do.call('c', sdnx_frags)
#
#
#
#
#list_params <- list(
  #'DMSO nfyb @ nfyb+ATAC' = list(dmso_frags, nfyb_dar_sig_atac),
  #'DMSO nfyb @ nfyb+No ATAC' = list(dmso_frags, nfyb_dar_sig_no_atac),
  #'SNDX nfyb @ nfyb + ATAC'= list(sdnx_frags, nfyb_dar_sig_atac),
  #'SNDX nfyb @ nfyb+No ATAC' = list(sdnx_frags, nfyb_dar_sig_no_atac)  
  #)
#
#p <- plotVmat(list_params, cores = 4, ncol= 2, nrow = 2, normFun = 'zscore')
#p
#--------------------------------------------------------------------------------
## Should find the NFYB motif in each region, center on that. 
nfy_motif <- MotifDb::MotifDb |> MotifDb::query('NFYB') |> convert_motifs() 
nfy_motif <- nfy_motif[[3]]
nfy_motif

nfyb_res <- nfyb_dar_sig |> get_sequence(hg38_genome) |> runFimo(motifs = nfy_motif, thresh = 1e-3, meme_path = meme_path)
           
nfyb_res <- nfyb_res |> 
  mutate(nfyb_range_id = paste(seqnames, start, end, sep = ":")) |>  
  as_granges()

# only a small number of peaks did not have the motif - so I excluded them from further analysis 
nfyb_dar_sig <- nfyb_dar_sig |> 
  mutate(group = ifelse(log2FoldChange > 0, 'Up', 'Down') )|> 
  as_granges() |> 
  join_overlap_left(nfyb_res, suffix= c('.nfyb', '.motif')) |> as.data.frame() |> 
  mutate(range_id = paste(seqnames, start, end, sep = ":")) |> 
  group_by(range_id) |> 
  filter(pvalue.motif == min(pvalue.motif)) |> 
  # If you want just one row even with ties:
  filter(row_number() == 1) 

nfyb_dar_sig <- nfyb_dar_sig |> 
  as_tibble() |> 
  separate(nfyb_range_id, into= c('seqnames_m', 'start_m', 'end_m'), sep =':', remove = F)

nfyb_dar_sig

new_ranges <- with( nfyb_dar_sig,
                   GRanges(seqnames = seqnames_m, ranges = IRanges(start = as.numeric(start_m), end = as.numeric(end_m)), group = group))


new_ranges_atac <- new_ranges |> filter_by_overlaps(atac_all)
new_ranges_noatac <- new_ranges |> filter_by_non_overlaps(atac_all)
new_ranges_atac

# looking at motif centered regions - these look much mbetter
# split up by atac vs not and up vs down
atac_params_up <- list(
 'DMSO nfyb @ nfyb+ATAC' = list(dmso_frags, new_ranges_atac |> filter(group == 'Up')),
 'SNDX nfyb @ nfyb + ATAC'= list(sdnx_frags, new_ranges_atac |> filter(group == 'Up'))
)

noatac_params_up <- list('DMSO nfyb @ nfyb+No ATAC' = list(dmso_frags,new_ranges_noatac |> filter(group == 'Up')),
                       'SNDX nfyb @ nfyb+No ATAC' = list(sdnx_frags, new_ranges_noatac |> filter(group == 'Up')) )


plt_atac_up <- plotVmat(atac_params_up, cores = 4)
plt_noatac_up <- plotVmat(noatac_params_up, cores = 4)



atac_params_down <- list(
  'DMSO nfyb @ nfyb+ATAC' = list(dmso_frags, new_ranges_atac |> filter(group == 'Down')),
  'SNDX nfyb @ nfyb + ATAC'= list(sdnx_frags, new_ranges_atac |> filter(group == 'Down'))
)

noatac_params_down <- list('DMSO nfyb @ nfyb+No ATAC' = list(dmso_frags,new_ranges_noatac |> filter(group == 'Down')),
                         'SNDX nfyb @ nfyb+No ATAC' = list(sdnx_frags, new_ranges_noatac |> filter(group == 'Down')) )


plt_atac_down <- plotVmat(atac_params_down, cores = 4)
plt_noatac_down <- plotVmat(noatac_params_down, cores = 4)
  


# Put all this in a plot
pdf(file = 'plots/integrated/fig5_supplements_nfyb_fragment_lengths.pdf', height = 11, width = 8.5, family = 'Helvetica')
pageCreate(showGuide = F)
plotGG(plt_coverage_dist, x = 4.5, y = 0.5, height = 2, width = 4)
plotGG(plt_noatac_down , x = 0.25, y = 1, width = 4.5, height = 4, just = c('top', 'left'))
plotGG(plt_atac_down, x = 0.25, y = 3.25, width = 4.5, height = 4)
plotGG(plt_noatac_up , x = 0.25, y = 5.5, width = 4.5, height = 4, just = c('top', 'left'))
plotGG(plt_atac_up, x = 0.25, y = 7.75, width = 4.5, height = 4)
dev.off()

# write out motif centers
write_bed(new_ranges |> filter(group == 'Up'), file = 'data/derived_data/integrated/nfyb_motifs_in_sig_up.bed')
write_bed(new_ranges |> filter(group == 'Down'), file = 'data/derived_data/integrated/nfyb_motifs_in_sig_down.bed')

