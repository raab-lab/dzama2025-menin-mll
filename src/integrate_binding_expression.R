## Integrate binding and expression
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(csaw)
library(ChIPpeakAnno)
library(plyranges)
library(ggplot2)
library(tidyverse)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
chrs <- paste0("chr", c(1:22, "X", "Y"))
fdr_cut <- 0.1
lfc_cut <- 0

#diffEx <- readRDS("data/derived_data/expression/shrunk_d4SNDX_vs_d4DMSO_diffEx.rds")
#jesse's files
diffEx <- read_tsv('data/processed_data/rna/hlf_sndx_v_dmso.tsv') |> as_granges()
seqlevelsStyle(diffEx) <- "UCSC"
diffEx <- filter(diffEx,
                 seqnames %in% chrs,
                 padj < 0.05) %>%
  promoters()
seqlevels(diffEx) <- seqlevelsInUse(diffEx)

menin <- read_tsv('data/derived_data/cnr/menin_diff_menpeaks_res.tsv') |> as_granges()
k4me3 <- read_tsv('data/derived_data/cnr/k4me3_diff_menpeaks_res.tsv') |> as_granges()

#jesse's overlap
menin_k4_sites <- menin |> filter_by_overlaps(k4me3, maxgap = 50) 
# peytons'
#menin_k4_sites <-
#  read.table(
#    "data/derived_data/binding/menin_H3K4me3_overlap.bed",
#    col.names = c("seqnames", "start", "end", "width", "strand", "FDR", "log2FC")
#  ) %>%
#  as_granges()


## What genes are near sites where k4 and menin are lost?
menin_k4_genes <- join_overlap_inner(diffEx %>% select(gene_name, log2FoldChange, padj),
                                     menin_k4_sites %>% filter(padj < fdr_cut, log2FoldChange < lfc_cut)) %>%
  split(x = ., f = .$gene_name) %>%
  purrr::map(.x = ., .f = ~reduce_ranges(.x,gene_name = unique(gene_name),
                                         geneFC = max(log2FoldChange.x),
                                         padj = max(padj.x),
                                         bindFC = mean(log2FoldChange.y),
                                         FDR = mean(padj.y))) %>%
  bind_ranges() %>%
  mutate(gene_name = as.character(gene_name))

menin_k4_genes

menin_k4_genes |> as_tibble() |> ggplot(aes(x = bindFC, y = geneFC)) + 
   geom_point() + 
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red2')  +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red2')

write.csv(as.data.frame(menin_k4_genes),
          "data/derived_data/integrate/menin_k4_loss_genes.csv",
          row.names = F,
          quote = F)
