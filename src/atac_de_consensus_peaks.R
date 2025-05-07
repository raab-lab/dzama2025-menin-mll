#  differential accessibility in atac
# Author: Jesse Raab
# Date: 2025-01-28 (last run)
suppressPackageStartupMessages({
  library(csaw)
  library(plyranges)
  library(tidyverse)
  library(DESeq2)
})

# load union peaks - these are regions of interested for Differential Occupancy

peaks <- read_bed(file.path('data/derived_data/peak_sets/atac_consensus/union_atac.bed'))

# get bam files
bams <- list.files(path = 'data/source_data/atac/bams/filtered/', pattern = '*sorted.bam$', full.names = T)


# Sample Sheet
ss <- data.frame(filename = bams, names = basename(bams))
ss <-
ss <- ss |> 
  mutate(short_name = str_replace(names, '.aligned.bam', '') ) |> 
  separate(short_name, into = c('airtable_id', 'NGSID', 'cell_line', 'condition', 'rep'), sep = '_' ) 
ss


# CSAW Setup
param = readParam(pe = 'both', dedup = F, minq = 10) # set this to 10 which matches my bam parameters, higher and we lose some numbers
peak_counts <- regionCounts(bams, peaks, param = param)

# normalize on the regions we've counted
peak_norm <- normFactors(peak_counts, se.out = T )

# Need to set up some sort of sample sheet for comparison
colData(peak_norm)$bam.files
colData(peak_norm) <- DataFrame(left_join(as.data.frame(colData(peak_norm)), ss, by = c('bam.files' ='filename') ) )
colData(peak_norm)
dds <- DESeq2::DESeqDataSet(peak_norm, design = ~ condition) 
#sizeFactors(dds) <- colData(dds)$norm.factors
dds$condition <- relevel(dds$condition, ref = "DMSO")
des <- DESeq(dds)

dds2 <-estimateSizeFactors(dds)
colData(dds2)
# check output
resultsNames(des)s
ko_v_wt <- lfcShrink(des, coef = 2, format = "GRangesList", type = 'apeglm')
ko_v_wt
# no differences
ko_v_wt |> as_tibble() |> ggplot(aes(  x= log2(baseMean), y = log2FoldChange)) + 
  geom_point(aes(color = padj < 0.05))

ko_v_wt |> as_tibble() |> ggplot(aes(  x = log2FoldChange, y = -log10(padj))) + 
  geom_point()
ko_v_wt |> as_tibble() |> filter(log2FoldChange < 0) |> arrange(padj)  
# no statisticall signficant differencese
table(ko_v_wt$padj < 0.1,ko_v_wt$log2FoldChange <0)
save(des, file = 'data/derived_data/atac/atac_de.Rda')

#write a tsv of everything
write_tsv(as.data.frame(ko_v_wt), file = 'data/derived_data/atac/all_data.tsv')

# unchanged
write_tsv(x = as.data.frame(ko_v_wt) |> filter(padj > 0.05), file = 'data/derived_data/atac/unchanged.tsv')
# WRite out upregulated and down regualted signficant genes
write_tsv(x = as.data.frame(ko_v_wt) |> filter(padj < 0.05) |> filter(log2FoldChange < 0), file = 'data/derived_data/atac/lost_access.tsv', col_names = T)
write_tsv(x = as.data.frame(ko_v_wt) |> filter(padj < 0.05) |> filter(log2FoldChange > 0), file = 'data/derived_data/atac/gain_access.tsv', col_names = T)


write_tsv(x = as.data.frame(ko_v_wt) |> filter(padj < 0.05) |> filter(log2FoldChange < log2(1.5)), file = 'data/derived_data/atac/lost_most_access.tsv', col_names = T)
write_tsv(x = as.data.frame(ko_v_wt) |> filter(padj < 0.05) |> filter(log2FoldChange > -log2(1.5)), file = 'data/derived_data/atac/gain_most_access.tsv', col_names = T)

