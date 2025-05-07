# Define upregualted targets in ATAC data and write bed list to score for motif usage
# Combine with expression to get those that increase


library(tidyverse)
library(plyranges)
library(ChIPseeker)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(biomaRt)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(memes)
hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
motif_db <- '../motif_databases/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme'

# Functions ################################################################################# 
fix_category <- function(tib) {
  tib |> 
    mutate(new_cat = case_when(grepl ('Promoter', annotation) ~ 'Promoter', 
                               grepl ('Intron', annotation) ~ 'Intron', 
                               grepl ('Exon', annotation) ~ 'Exon', 
                               .default = as.character(annotation) ) )
  
  
}

################################################################################

# load RNA data after drug treatment, and get upregulated genes

de_rna <- read_tsv('data/derived_data/rna/d4_sndx_all_genes.tsv') |> 
  filter(padj <0.05) |> 
  filter(log2FoldChange > 0)
de_rna<- as_granges(de_rna, keep_mcols = T) |> keepStandardChromosomes(pruning.mode = 'coarse')
seqlevelsStyle(de_rna) <- 'UCSC'
de_rna

# load ATAC changes after drug treatment, get upregulated sites

de_atac_up <- read_bed('data/derived_data/atac/gain_access.bed') |> keepStandardChromosomes(pruning.mode = 'coarse')
de_atac_up <-  de_atac_up |> annotatePeak(TxDb = txdb, tssRegion = c(-2000, 200), ) |> as_granges(keep_mcols = T) 
de_atac_up <- de_atac_up |> fix_category()
de_atac_up_promoter <- de_atac_up |> filter(new_cat == 'Promoter')
de_atac_up_intron <- de_atac_up |> filter(new_cat == 'Intron')
de_atac_distal <- de_atac_up |> filter(new_cat == 'Distal Intergenic')

# to look at numbers per cvateogry
de_atac_up |> 
  as_tibble() |> 
  count(new_cat)

# some of these may have K27me3 


ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = 'useast')

de_atac_up <- annotatePeakInBatch(de_atac_up, mart = ensembl ) |> addGeneIDs(orgAnn = org.Hs.eg.db, IDs2Add = 'SYMBOL')

# Intersect with  expression changes
simple_rna <- de_rna <- de_rna |> dplyr::select(gene_name, log2FoldChange)
# These are the upregulated atac sites in promoters associated with increased expression
de_atac_plus_rna_up <- de_atac_up |>
  as_tibble() |> 
  left_join(de_rna |> as_tibble(), by = c('SYMBOL'= 'gene_name'), suffix = c('.atac', '.rna'))
# write to a bed file
de_atac_plus_rna_up_bed <- de_atac_plus_rna_up |> 
  dplyr::select(seqnames.atac, start.atac, end.atac)

write_tsv(de_atac_plus_rna_up_bed, file = 'data/derived_data/atac/upregualted_promters_up_rna.bed')

de_atac_plus_rna_up_bed |> sample_n(10)

colnames(de_atac_plus_rna_up_bed) <- c('seqnames', 'start', 'end')

# now use memes to see which of these have ap1 family motifs
dna_in <- de_atac_plus_rna_up_bed |> 
  as_granges() |> 
  plyranges::anchor_center() |> 
  plyranges::mutate(width = 500)

dna_in_string <- dna_in |> as_granges() |>  
                           get_sequence(genome = hg38)

dna_in_string

ap1 <- MotifDb::MotifDb |> MotifDb::query('AP-1') |> universalmotif::convert_motifs()  
fimo_res <- runFimo(dna_in_string, ap1)
colnames(de_atac_plus_rna_up)[1:3] <- c('seqnames', 'start', 'end')
overlap_counts <- count_overlaps(de_atac_plus_rna_up |> as_granges(), fimo_res)
de_atac_plus_rna_up$overlap_ap1 <- overlap_counts
  
de_atac_plus_rna_up |> 
  filter(!is.na(log2FoldChange)) |> 
  filter(overlap_ap1 > 0) |> View()
