  # split data into menin-nfyb-rna-category
# in paper = figure 5K
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
  
  nfyb_dor <- nfyb_dor |> as_granges() |> annotatePeak(TxDb = txdb, tssRegion = c(-2000, 200), annoDb='org.Hs.eg.db'  ) |> as_tibble() |> fix_category()
  nfyb_dor <- nfyb_dor |> as_granges() |> join_nearest(ex, suffix = c('.nfyb', '.rna'))
  nfyb_dor_up <- nfyb_dor |> filter(padj.nfyb < 0.1 & log2FoldChange.nfyb > 0)
  nfyb_dor_down <- nfyb_dor |> filter(padj.nfyb < 0.1 & log2FoldChange.nfyb < 0)
  nfyb_dor_stable <- nfyb_dor |> filter(padj.nfyb > 0.1)
  
  menin_df <- read_tsv('data/derived_data/cnr/menin_diff_menpeaks_res.tsv')
  menin_df <- menin_df |> as_granges() |> annotatePeak(TxDb = txdb, tssRegion = c(-2000, 200), annoDb='org.Hs.eg.db'  ) |> as_tibble() |> fix_category()
  menin_df <- menin_df |> as_granges() |> join_nearest(ex, suffix = c('.menin', '.rna'))
  menin_down <- menin_df |> filter(padj.menin < 0.05 & log2FoldChange.menin < 0) 
  menin_up   <- menin_df |> filter(padj.menin < 0.05 & log2FoldChange.menin > 0)
  menin_stable <- menin_df |> filter(padj.menin > 0.1)
  
  # category 1 
  #menin down , nfyb down, rna down
  both_lost <- menin_down |> filter_by_overlaps(nfyb_dor_down, maxgap = 500) 
  both_lost <- both_lost|> mutate(group = 'alldown')
  
  # cat 2 
  # menin down - nfyb not down
  men_down_nfyb_non_overlap <- menin_down |> filter_by_non_overlaps(nfyb_dor_down, maxgap = 500) |>  filter(padj.rna > 0.05)
  men_down_nfyb_stable <- men_down_nfyb_non_overlap |> filter_by_overlaps(nfyb_dor_stable) # this is about 5765 - menin lost but nfyb remains
  men_down_nfyb_non_overlap <- men_down_nfyb_non_overlap |> filter_by_non_overlaps(men_down_nfyb_stable) # this is about 5444  - so menin indpendent of nfyb-
  men_down_nfyb_stable <- men_down_nfyb_stable |> mutate(group = 'mendown_nfybstable')
  men_down_nfyb_non_overlap <- men_down_nfyb_non_overlap |> mutate(group = 'mendown_nfyb_nonoverlap')
  # NFYB down - not overlap menin down
  nfyb_dor_down_no_menin <- nfyb_dor_down |> filter_by_non_overlaps(menin_down, maxgap=500) 
  nfyb_dor_down_no_menin |> 
    as_tibble() |> ggplot(aes(x= 1, y = log2FoldChange.rna)) + geom_boxplot()
  nfyb_dor_down_no_menin <- nfyb_dor_down_no_menin |> mutate(group = 'nfybdown_nomen')
  
  #NFYB up no menin
  nfyb_dor_up_no_menin <- nfyb_dor_up |> filter_by_non_overlaps(c(menin_stable, menin_down), maxgap = 500) |> mutate(group = 'nfybup_nomen')
  
  all_comb <- c(men_down_nfyb_stable, both_lost, nfyb_dor_up_no_menin)
class_plot <-  all_comb |> 
    as_tibble() |> 
    ggplot(aes(x = group, y = log2FoldChange.rna)) +
    geom_point(position = position_jitter(width = 0.1), alpha = 0.3)+
  geom_boxplot(outliers = F, width = 0.3, fill = NA) + 
  theme_bw()
  
  pairwise.wilcox.test(x = all_comb$log2FoldChange.rna, g = all_comb$group )
  
t.test(both_lost |> as_tibble()  |> pull(log2FoldChange.rna), mu = 0)
t.test(men_down_nfyb_stable|> as_tibble()  |> pull(log2FoldChange.rna), mu = 0)
t.test(nfyb_dor_up_no_menin|> as_tibble()  |> pull(log2FoldChange.rna), mu = 0)

all_comb |> as_tibble() |> ggplot(aes(x = log2FoldChange.rna, -log10(padj.rna))) + geom_point()  
  
  run_great <- function(peak_gr, T2G){
    great_res <- great(peak_gr, T2G, 'TxDb.Hsapiens.UCSC.hg38.knownGene')
    great_res_tab <- getEnrichmentTable(great_res) 
    plt <- great_res_tab |> 
      ggplot(aes(x = fold_enrichment_hyper, y = -log10(p_adjust_hyper))) + 
      geom_point() + 
      ggrepel::geom_text_repel(aes(label = id |> str_replace('HALLMARK_', '') ), data = great_res_tab |> filter(p_adjust_hyper < 0.05), size = 2, family = 'sans'  )
    return(list(plot=plt, res=great_res_tab))
  }
  
  cat1_h   <- run_great(both_lost, 'MSigDB:H')
  cat1_c2     <- run_great(both_lost, 'MSigDB:C2:CGP')
  cat1_pid <- run_great(both_lost, 'MSigDB:C2:CP:PID')
  cat1_c3     <- run_great(both_lost, 'MSigDB:C3:TFT')
  cat1_c6     <- run_great(both_lost, 'MSigDB:C6')
  
  cat2a_h <- run_great(men_down_nfyb_stable, 'MSigDB:H') 
  cat2a_c2   <- run_great(men_down_nfyb_stable, 'MSigDB:C2:CGP') 
  cat2a_pid <- run_great(men_down_nfyb_stable, 'MSigDB:C2:CP:PID')
  cat2a_c3   <- run_great(men_down_nfyb_stable, 'MSigDB:C3:TFT')
  cat2a_c6   <- run_great(men_down_nfyb_stable, 'MSigDB:C6')
  
  cat2b_h <- run_great(men_down_nfyb_non_overlap, 'MSigDB:H')
  cat2b_c2 <- run_great(men_down_nfyb_non_overlap, 'MSigDB:C2:CGP')
  cat2b_pid <- run_great(men_down_nfyb_non_overlap, 'MSigDB:C2:CP:PID')
  cat2b_c3   <- run_great(men_down_nfyb_non_overlap, 'MSigDB:C3:TFT')
  cat2b_c6 <- run_great(men_down_nfyb_non_overlap, 'MSigDB:C6')
  
  cat3_h <- run_great(nfyb_dor_down_no_menin, 'MSigDB:H')
  cat3_c2 <- run_great(nfyb_dor_down_no_menin, 'MSigDB:C2:CGP')
  cat3_pid <- run_great(nfyb_dor_down_no_menin, 'MSigDB:C2:CP:PID')
  cat3_c3   <- run_great(nfyb_dor_down_no_menin, 'MSigDB:C3:TFT')
  cat3_c6 <- run_great(nfyb_dor_down_no_menin, 'MSigDB:C6')
  
  cat4_h <- run_great(nfyb_dor_up_no_menin, 'MSigDB:H')
  cat4_c2 <- run_great(nfyb_dor_up_no_menin, 'MSigDB:C2:CGP')
  cat4_pid <- run_great(nfyb_dor_up_no_menin, 'MSigDB:C2:CP:PID')
  cat4_c3 <- run_great(nfyb_dor_up_no_menin, 'MSigDB:C3:TFT')
  cat4_c6 <- run_great(nfyb_dor_up_no_menin, 'MSigDB:C6')
  
  write_bed(men_down_nfyb_stable, file = 'data/processed_data/cnr/menin_down_nfyb_stable.bed')
  write_bed(nfyb_dor_up_no_menin, file = 'data/processed_data/cnr/nfyb_up_no_menin.bed')

# Dot plot of category
cat1_all <- bind_rows(cat1_h$res, cat1_c2$res, cat1_c6$res)
cat1_hall <- bind_rows(cat1_h$res, cat1_pid$res)
cat1_hall$group <- 'Double-Loss'

cat2a_all <- bind_rows(cat2a_h$res, cat2a_c2$res, cat2a_c6$res)
cat2a_hall <- bind_rows(cat2a_h$res, cat2a_pid$res)
cat2a_hall$group <- 'Menin-loss_NYFB-stable'

cat2b_all <- bind_rows(cat2b_h$res, cat2b_c2$res, cat2b_c6$res)
cat2b_hall <- bind_rows(cat2b_h$res, cat2b_pid$res)
cat2b_hall$group <- 'Menin-loss_No-NFYB'

cat4_all <- bind_rows(cat4_h$res, cat4_c2$res, cat4_c6$res)
cat4_hall <- bind_rows(cat4_h$res, cat4_pid$res)
cat4_hall$group <- 'NFYB-Up'

all_groups <- bind_rows(cat1_hall,cat2a_hall,cat4_hall)
interest <- all_groups |> 
  group_by(group) |> 
  arrange(p_adjust, .by_group = TRUE) |>
  slice_head(n=20) |> 
  ungroup() 
  
great_plot <- all_groups |> 
  filter(id %in% interest$id) |> 
  mutate(id = str_replace(string = id, pattern ='HALLMARK_', replacement = '' )) |> 
  mutate(id = str_replace(string = id, pattern ='PID_', replacement = '' )) |> 
  ggplot(aes(x = group, y = id, fill = fold_enrichment_hyper, size = -log10(p_adjust+1e-14))) + 
  geom_point(color = 'grey30', pch = 21) +
  scale_fill_gradient2(low = 'white', high = 'red3') +
  theme_minimal() + 
  theme(axis.text = element_text(family = 'sans', size = 7),
        axis.title = element_text(family = 'sans', size= 8))



cols = circlize::colorRamp2(breaks = c(0,4), colors = c('white', 'red'))
cols
all_groups |> 
  filter(id %in% interest$id) |> 
  mutate(id = str_replace(string = id, pattern ='HALLMARK_', replacement = '' )) |> 
  mutate(id = str_replace(string = id, pattern ='PID_', replacement = '' )) |> 
  as_tibble() |> 
  tidyHeatmap::heatmap(.row = id, .colum = group, .value = fold_enrichment,
                       scale = 'none', palette_value = cols, na_col = 'white')


# get gene sets 
pdf('plots/integrated/nfyb_class_plots.pdf', width = 8.5, height = 11) 
pageCreate(showGuides = F)
plotGG(great_plot, x = 0.5, y = 0.5, height = 6, width = 7.5)
plotGG(class_plot, x = 0.5, y = 6.5, height = 4, width = 4)
dev.off()

nfyb_dor_up_no_menin |> as_tibble() |> write_tsv(file = 'tmp.nfyb_up_no_menin.tsv')
nfyb_dor_up_no_menin |> as_tibble() |> filter(log2FoldChange.rna> 0.5) |> View()

write_tsv(men_down_nfyb_stable |> as_tibble(), file = 'tmp_men_down_nfyb_stable.tsv')
