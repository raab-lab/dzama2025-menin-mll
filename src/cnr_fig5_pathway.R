# Figure 5
# Need and expression and a pathway 

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


df <- read_tsv('data/processed_data/integrated/nfyb_dor_anno_all.tsv')

run_great <- function(peak_gr, T2G){
  great_res <- great(peak_gr, T2G, 'TxDb.Hsapiens.UCSC.hg38.knownGene')
  great_res_tab <- getEnrichmentTable(great_res) 
  plt <- great_res_tab |> 
    ggplot(aes(x = fold_enrichment_hyper, y = -log10(p_adjust_hyper))) + 
    geom_point() + 
    ggrepel::geom_text_repel(aes(label = id |> str_replace('HALLMARK_', '') ), data = great_res_tab |> filter(p_adjust_hyper < 0.05), size = 2, family = 'sans'  )
  return(list(plot=plt, res=great_res_tab))
}

df_up <- df |> filter(group == "Up") |> as_granges()
df_down <- df |> filter(group == "Down") |> as_granges()

                        
lost_c2cp_out     <- run_great(df_down, 'msigdb:C2:CP')
lost_hall_out     <- run_great(df_down, 'MSigDB:H')
lost_c5_out     <- run_great(df_down, 'msigdb:C5')
lost_c6_out     <- run_great(df_down, 'MSigDB:C6')
lost_c7_out     <- run_great(df_down, 'msigdb:C7:IMMUNESIGDB')
lost_c8_out     <- run_great(df_down, 'msigdb:C8')
lost_c2_out    <- run_great(df_down, 'msigdb:C2')


gain_c2cp_out     <- run_great(df_up |> filter(padj.rna >0.1), 'msigdb:C2:CP')
gain_hall_out     <- run_great(df_up |> filter(padj.rna >0.1), 'MSigDB:H')
gain_hall_out     <- run_great(df_up, 'MSigDB:H')
gain_hall_out <- run_great(df_up |> filter(log2FoldChange.rna > 0 & padj.rna < 0.1), 'MSigDB:H')
gain_c5_out     <- run_great(df_up, 'msigdb:C5')
gain_c6_out     <- run_great(df_up, 'MSigDB:C6')
gain_c7_out     <- run_great(df_up, 'msigdb:C7:IMMUNESIGDB')
gain_c8_out     <- run_great(df_up, 'msigdb:C8')
gain_c2_out    <- run_great(df_up|> filter(padj.rna >0.1), 'msigdb:C2')
gain_hall_out <- run_great(df_up |> filter(log2FoldChange.rna > 0 & padj.rna < 0.1), 'msigdb:C2')

df

