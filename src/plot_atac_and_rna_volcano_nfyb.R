#Supplemetnal figure on NFYB change at ATAC grouping
library(tidyverse)
library(plyranges)
library(plotgardener)

############################################################################### #
fix_category <- function(tib) {
  tib |> 
    mutate(new_cat = case_when(grepl ('Promoter', annotation) ~ 'Promoter', 
                               grepl ('Intron', annotation) ~ 'Intron', 
                               grepl ('Exon', annotation) ~ 'Exon', 
                               .default = as.character(annotation) ) )
  
  
  
}
################################################################################
nfyb_dar <- read_tsv('data/derived_data/cnr/nfyb_analysis/nfyb_diff_at_nfyb_peaks.tsv') |> as_granges()
nfyb_dor <- read_tsv('data/derived_data/cnr/nfyb_analysis/nfyb_diff_at_nfyb_peaks.tsv') |> as_granges()

atac_up <- read_tsv('data/derived_data/atac/gain_access.tsv') |> as_granges()
atac_down <- read_tsv('data/derived_data/atac/lost_access.tsv') |> as_granges()
atac_unchanged <- read_tsv('data/derived_data/atac/unchanged.tsv') |> as_granges()
atac_all <- read_tsv('data/derived_data/atac/all_data.tsv') |> as_granges()
################################################################################

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

################################################################################
nfyb_atac_all <- c(nfyb_atac_up, nfyb_atac_down, nfyb_atac_unchanged)
# supplemental Figure of nfyb changes by ATAC category
plt_nfyb_atac_volcano <- nfyb_atac_all |> 
  as_tibble() |> 
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = group)) +
  geom_point() + 
  scale_color_manual(values = c('steelblue','grey50', 'red3')) + 
  geom_hline(aes(yintercept = -log10(0.05)), linetype = 'dashed', color = 'grey30')+
  theme_bw() + 
  facet_wrap(~group)+
  theme(panel.grid = element_blank())
################################################################################
# can I make a similar version of this for RNA
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
ex <- read_tsv('data/derived_data/rna/d4_sndx_all_genes.tsv') |> as_granges()
seqlevelsStyle(ex) <- 'UCSC'

rna_up <- ex |> filter(log2FoldChange > 0 & padj < 0.05)
rna_down <- ex |> filter(log2FoldChange < 0 & padj < 0.05)
rna_stable <- ex |> filter(padj > 0.05)
#Annotate NFYB
nfyb_dor_anno_all <- nfyb_dor |> as_granges() |> 
  annotatePeak(TxDb = txdb, tssRegion = c(-2000, 200), annoDb='org.Hs.eg.db'  ) |> 
  as_tibble() |> 
  fix_category()
nfyb_rna_up <- nfyb_dor_anno_all |> filter(SYMBOL %in% rna_up$gene_name) |> mutate(group = 'RNA Up')
nfyb_rna_down <- nfyb_dor_anno_all |> filter(SYMBOL %in% rna_down$gene_name) |> mutate(group = 'RNA Down')
nfyb_rna_stable <- nfyb_dor_anno_all |> filter(SYMBOL %in% rna_stable$gene_name) |> mutate(group = 'RNA Stable')
nfyb_rna_all <-  bind_rows(nfyb_rna_up, nfyb_rna_down, nfyb_rna_stable)

plt_nfyb_rna_volcano <- nfyb_rna_all |> 
  as_tibble() |> 
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = group)) +
  geom_point() + 
  scale_color_manual(values = c('steelblue','grey50', 'red3')) + 
  geom_hline(aes(yintercept = -log10(0.05)), linetype = 'dashed', color = 'grey30')+
  theme_bw() + 
  facet_wrap(~group)+
  theme(panel.grid = element_blank())

################################################################################
pdf('plots/integrated/fig5_supp_nfyb_sig_at_atacgroup_volcano.pdf', height = 11, width = 8.5)
pageCreate(showGuides = F)
plotGG(plt_nfyb_atac_volcano + theme(legend.position = 'none'), x = 0.5, y = 0.5, width = 4.5, height = 2)
plotGG(plt_nfyb_rna_volcano + theme(legend.position = 'none'),x = 0.5, y = 3, width = 4.5, height = 2)
dev.off()


################################################################################
table(nfyb_rna_all$group, nfyb_rna_all$padj < 0.1)

