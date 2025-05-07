# Figure 5 NFYB analyssis
##### THIS IS AN EARLY VERSION NOT FINAL FIGURE ########
# not in paper

library(tidyverse)
library(DESeq2)
library(rGREAT)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene) 
library(ChIPseeker)
library(clusterProfiler)
library(msigdbr)
library(plyranges)
library(plotgardener)
library(grid)

################################################################################
fix_category <- function(tib) {
  tib |> 
    mutate(new_cat = case_when(grepl ('Promoter', annotation) ~ 'Promoter', 
                               grepl ('Intron', annotation) ~ 'Intron', 
                               grepl ('Exon', annotation) ~ 'Exon', 
                               .default = as.character(annotation) ) )
  
  

}
outdir <- 'data/derived_data/integrated/'
dir.create(outdir, recursive = T)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
load('data/derived_data/cnr/nfyb_analysis/nfyb_analysis_obj.Rda')
nfyb_dar <- read_tsv('data/derived_data/cnr/nfyb_analysis/nfyb_diff_at_atac.tsv')
nfyb_dor <- read_tsv('data/derived_data/cnr/nfyb_analysis/nfyb_diff_at_nfyb_peaks.tsv')
ex <- read_tsv('data/derived_data/rna/d4_sndx_all_genes.tsv') |> as_granges()
seqlevelsStyle(ex) <- 'UCSC'


nfyb_dor_up <- nfyb_dor |> filter(padj < 0.05 & log2FoldChange > 0)
nfyb_dor_down <- nfyb_dor |> filter(padj < 0.05 & log2FoldChange < 0)


#load accessibility regions
load('data/derived_data/atac/atac_de.Rda')
# re-shrink the atac
resultsNames(des)
dar <- lfcShrink(des, type = 'apeglm', coef = 2, format = 'GRanges')
dar
nfyb_dar
nfyb_dor_anno_all <- nfyb_dor |> as_granges() |> annotatePeak(TxDb = txdb, tssRegion = c(-2000, 200), annoDb='org.Hs.eg.db'  )
nfyb_dor_anno_all <- nfyb_dor_anno_all |> as_granges() |> fix_category()
nfyb_dor_anno_all <- nfyb_dor_anno_all |> mutate(group = 
                              case_when( log2FoldChange < 0 & padj < 0.1 ~ 'Down', 
                                         log2FoldChange > 0 & padj < 0.1 ~ 'Up', 
                                         padj > 0.1 ~ 'Stable') )

nfyb_dor_anno_all <- nfyb_dor_anno_all |> as_tibble() |> left_join(ex |> as_tibble(), by = c('SYMBOL'='gene_name'), relationship = 'many-to-many') 

up_sig_dor_num <- nfyb_dor |> as_tibble() |> filter(padj < 0.05 & log2FoldChange > log2(1.5)) |> pull() |> length()
down_sig_dor_num <- nfyb_dor |> as_tibble() |> filter(padj < 0.05 & log2FoldChange < -log2(1.5)) |> pull() |> length() 

# scatter plot of nfyb and accessibilty fold changes
dar_comb <- dar |> as_tibble() |> 
  left_join(nfyb_dar, by = c('seqnames', 'start', 'end', 'strand'), suffix = c('.atac', '.nfyb')) 

n_up_sig <- dar_comb |> 
  filter(padj.nfyb < 0.1 & padj.atac < 0.1) |> 
  filter(log2FoldChange.nfyb > 0) |> 
  filter(sign(log2FoldChange.nfyb) == sign(log2FoldChange.atac)) |> 
  dplyr::count()

n_down_sig <- dar_comb |> 
  filter(padj.nfyb < 0.1 & padj.atac < 0.1) |> 
  filter(log2FoldChange.nfyb < 0) |> 
  filter(sign(log2FoldChange.nfyb) == sign(log2FoldChange.atac)) |> 
  dplyr::count()



#NFYB
nfyb_dmso <- read_bed('data/derived_data/peak_sets/cnr_consensus/nfyb_union.bed')|> keepStandardChromosomes(pruning.mode = 'coarse')
atac_dmso <- read_bed('data/derived_data/peak_sets/atac_consensus/union_atac.bed')|> keepStandardChromosomes(pruning.mode = 'coarse')

nfyb_sndx <- read_bed('data/derived_data/peak_sets/cnr_consensus/nfyb_sndx_consensus.bed') |> keepStandardChromosomes(pruning.mode = 'coarse')
atac_sndx <- read_bed('data/derived_data/peak_sets/atac_consensus/sndx_consensus_atac.bed') |> keepStandardChromosomes(pruning.mode = 'coarse')


### VENN ###
# Venn of upregulated NFYB and ATAC sites (ATAC from sndx) 
nfyb_dor_up <- nfyb_dor |> 
  filter(baseMean > 20) |> 
  filter(log2FoldChange > 0 & padj < 0.1) 
nfyb_up_list <- list(nfyb_dor_up |> as_granges(), atac_sndx)

# Venn of Downregualted NFYB and ATAC sites (ATAC from dmso)
nfyb_dor_down <- nfyb_dor |> 
  filter(baseMean > 20) |> 
  filter(log2FoldChange < 0 & padj < 0.1)

nfyb_down_list <- list(nfyb_dor_down |> as_granges(), atac_dmso)


#nfyb_dor_gr <- as_granges(nfyb_dor)
#nfyb_dor_gr[overlapsAny(nfyb_dor_gr, union(c(atac_sndx, atac_dmso)))]



# Find NFYB up that have ATAC going with it

# this gives me 411
nfyb_dor_up_atac <- nfyb_dor |> 
  as_granges() |> 
  filter(padj < 0.1) |> 
  filter(log2FoldChange > 0) |> 
  filter_by_overlaps(dar |> filter(padj < 0.1& log2FoldChange> 0) )  |> 
  as.data.frame() |> 
  arrange(padj)
  

write_bed(nfyb_dor_up_atac, file = 'data/derived_data/integrated/nfyb_dor_up_atac_up.bed')



nfyb_dor_up_atac_anno <- nfyb_dor |> 
  as_granges() |> 
  join_nearest(ex, suffix = c('.nfyb', '.rna')) |> 
  join_nearest(dar)

nfyb_dor_up_atac_anno |> 
  mutate(group = case_when (log2FoldChange.nfyb > 0.58 & padj.nfyb < 0.05 ~ 'Up',
                            log2FoldChange.nfyb < -0.58 & padj.nfyb < 0.05 ~ 'Down',
                            .default = 'Stable') ) |> 
  as_tibble() |> 
  ggplot(aes(x = group, y = log2FoldChange.rna)) + 
  geom_boxplot()

nfyb_dor_up_atac_anno |> 
  as_tibble() |> 
  ggplot(aes(x = log2FoldChange.nfyb, y = log2FoldChange.rna)) + 
  geom_point() + 
  geom_smooth(method = 'lm')

nfyb_dor_atac_anno <- nfyb_dor_up_atac_anno |> 
  annotatePeak(TxDb = txdb, tssRegion = c(-2000, 200), annoDb='org.Hs.eg.db'  )

nfyb_dor_atac_anno_small <- nfyb_dor_atac_anno |> as.data.frame() |> 
  dplyr::select(seqnames, start, end, 
                baseMean.nfyb, log2FoldChange.nfyb, padj.nfyb, 
                baseMean.rna, log2FoldChange.rna, padj.rna, 
                gene_name, baseMean, log2FoldChange, padj, 
                annotation, distanceToTSS, SYMBOL) |> 
  fix_category()

colnames(nfyb_dor_atac_anno_small)[11:13] <- c('baseMean.atac', 'log2FoldChange.atac', 'padj.atac')
table(nfyb_dor_atac_anno_small$gene_name == nfyb_dor_atac_anno_small$SYMBOL)


# plot just of promotesr# down with concordant RNA
nfyb_dor_atac_anno_small |> 
  mutate(group = case_when (log2FoldChange.nfyb > 0.58 & padj.nfyb < 0.05 ~ 'Up',
                            log2FoldChange.nfyb < -0.58 & padj.nfyb < 0.05 ~ 'Down',
                            .default = 'Stable') ) |> 
  as_tibble() |> 
  filter(new_cat == 'Promoter') |> 
  ggplot(aes(x = group, y = log2FoldChange.rna)) + 
  geom_boxplot()

nfyb_dor_atac_anno_small <- nfyb_dor_atac_anno_small |> 
  filter(new_cat == 'Promoter') |> 
  mutate(rna_chnage = 
           case_when (padj.rna < 0.05 & log2FoldChange.rna > 0 ~ 'Up', 
                      padj.rna < 0.05 & log2FoldChange.rna < 0 ~ 'Down', 
                      .default = 'Unchnaged'))
up_interest <- nfyb_dor_atac_anno_small |> filter(rna_chnage == 'Up') |> filter(padj.nfyb < 0.1)
down_interest <- nfyb_dor_atac_anno_small |> filter(rna_chnage == 'Down') |> filter(padj.nfyb < 0.1)


up_interest |> arrange(desc(log2FoldChange.rna))


nfyb_exp_anno <- nfyb_dor |> as_granges() |> annotatePeak(TxDb = txdb, tssRegion = c(-2000, 200), annoDb='org.Hs.eg.db'  ) |> 
  as_tibble() |> 
  fix_category() |> 
  left_join(ex |> as.data.frame(), by = c('SYMBOL' = 'gene_name'), suffix = c('.nfyb', '.rna'))

nfyb_exp_anno <- nfyb_exp_anno |> mutate(group = case_when (log2FoldChange.nfyb > 0.58 & padj.nfyb < 0.05 ~ 'Up',
                                        log2FoldChange.nfyb < -0.58 & padj.nfyb < 0.05 ~ 'Down',
                                        .default = 'Stable') ) 
nfyb_exp_anno |> 
  ggplot(aes(x = group, y = log2FoldChange.rna)) + 
  geom_boxplot() + 
  theme_bw()

table(nfyb_exp_anno$group)
