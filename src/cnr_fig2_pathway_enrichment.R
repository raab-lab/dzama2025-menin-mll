suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(plyranges)
  library(ChIPseeker)
  library(plotgardener)
  library(grid)
  library(ChIPpeakAnno)
  library(clusterProfiler)
  library(msigdbr)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene) 
  library(rGREAT)
  
})
# Author: Jesse Raab
# updated: 2025-01-20

# use the peak sets from differential expression at menin peaks (we tested each ab specifcally)
menin_diff  <- read_tsv('data/derived_data/cnr/menin_diff_menpeaks_res.tsv') |> as_granges() 
mll_diff  <- read_tsv('data/derived_data/cnr/mll_diff_menpeaks_res.tsv')|> as_granges() 
h3k4me3_diff <- read_tsv('data/derived_data/cnr/k4me3_diff_menpeaks_res.tsv')|> as_granges() 

menin_sig <- menin_diff |> filter(padj < 0.05) # 17002/18586
mll_sig   <- mll_diff |> filter(padj < 0.05) # 78209/18586
h3k4me3_sig <- h3k4me3_diff |> filter(padj < 0.05) # 11106/18586

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# look at GO terms for genes near lost Menin Peaks that lose all 3
lost_sndx <- menin_sig |> filter_by_overlaps(h3k4me3_sig) |> filter_by_overlaps(mll_sig)
# 5989 that fit that criteria

fix_category <- function(tib) {
  tib |> 
    mutate(new_cat = case_when(grepl ('Promoter', annotation) ~ 'Promoter', 
                               grepl ('Intron', annotation) ~ 'Intron', 
                               grepl ('Exon', annotation) ~ 'Exon', 
                               .default = as.character(annotation) ) ) 
}

lost_sndx <- lost_sndx |>  annotatePeak(TxDb = txdb, tssRegion = c(-2000, 200), annoDb = 'org.Hs.eg.db' )  |> as_tibble() |> fix_category() 

# function to run the GREAT analysis and generate table and plot
run_great <- function(peak_gr, T2G){
  great_res <- great(peak_gr, T2G, 'TxDb.Hsapiens.UCSC.hg38.knownGene')
  great_res_tab <- getEnrichmentTable(great_res) 
  plt <- great_res_tab |> 
    ggplot(aes(x = fold_enrichment_hyper, y = -log10(p_adjust_hyper))) + 
    geom_point() + 
    ggrepel::geom_text_repel(aes(label = id |> str_replace('HALLMARK_', '') ), data = great_res_tab |> filter(p_adjust_hyper < 0.05), size = 2, family = 'sans'  )
  return(list(plot=plt, res=great_res_tab))
}

lost_c2cp_out     <- run_great(lost_sndx |> as_granges(), 'msigdb:C2:CP')
lost_hall_out     <- run_great(lost_sndx |> as_granges(), 'MSigDB:H')
lost_c5_out     <- run_great(lost_sndx |> as_granges(), 'msigdb:C5')
lost_c6_out     <- run_great(lost_sndx |> as_granges(), 'MSigDB:C6')
lost_c7_out     <- run_great(lost_sndx |> as_granges(), 'msigdb:C7:IMMUNESIGDB')
lost_c8_out     <- run_great(lost_sndx |> as_granges(), 'msigdb:C8')

lost_c2cp_out$res |> arrange(p_adjust_hyper) |> head(10)
lost_hall_out$res |> arrange(p_adjust_hyper) |> head(10)
lost_c5_out$res |> arrange(p_adjust_hyper) |> head(10)
lost_c6_out$res |> arrange(p_adjust_hyper) |> head(10)
lost_c7_out$res |> arrange(p_adjust_hyper) |> head(10)
lost_c8_out$res |> arrange(p_adjust_hyper) |> head(10)


write_tsv(lost_c2cp_out$res, 'data/derived_data/cnr/lost_c2cp_great.tsv')
write_tsv(lost_hall_out$res, 'data/derived_data/cnr/lost_hall_great.tsv')
write_tsv(lost_c5_out$res, 'data/derived_data/cnr/lost_c5_great.tsv')

# enrichment plot figure 2M
out_plot <- lost_hall_out[[1]]
out_plot <- out_plot + theme_bw() + 
  geom_point(aes(color = p_adjust < 0.05)) + scale_color_manual(values = c('grey20', 'red2')) +
  xlab('Fold Enrichment') + ylab ('-log10 Adjust P-value') + 
  theme(text = element_text(family = 'sans'), 
        axis.title = element_text(size = 9), 
        axis.text = element_text(size = 7), 
        legend.position = 'none'
  ) 


# Make the plot of FC at menin peaks
#merge the original data sets
m1 <- menin_diff |> join_overlap_left(mll_diff, suffix = c('.menin', '.mll'))
m2 <- m1 |> join_overlap_left(h3k4me3_diff, suffix = c('', '.h3k4me3'))
m2 <- m2 |> as_tibble() |> rowwise() |> mutate(minp = min(padj.menin, padj.mll, padj, na.rm = T))
