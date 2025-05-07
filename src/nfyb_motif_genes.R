# some expression data on NFYB - changed peak

# NFYB 

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
ex <- read_tsv('data/derived_data/rna/d4_sndx_all_genes.tsv') |> as_granges()
seqlevelsStyle(ex) <- 'UCSC'


nfyb_up <- read_bed('data')
