# Ontology and comparison of targets

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(msigdbr)
  library(UpSetR)
  library(ggVennDiagram)
  library(ComplexHeatmap)
  library(plotgardener)
})
# label genes as target or nontareet
target <- read_tsv('data/processed_data/cnr/menin_misregualted_genes.tsv')
nontarget <- read_tsv('data/processed_data/cnr/menin_notchanged_genes.tsv')

hlf_sndx <- read_tsv('data/processed_data/rna/hlf_sndx_v_dmso.tsv')
hlf_ko   <- read_tsv('data/processed_data/rna/hlf_menko_v_wt.tsv')

hlf_sndx_mini <- hlf_sndx |> dplyr::select(rowname, gene_name, log2FoldChange, baseMean, padj) 
hlf_ko_mini   <- hlf_ko |> dplyr::select(rowname, gene_name, log2FoldChange, baseMean, padj)

hlf_comb <- left_join(hlf_sndx_mini, hlf_ko_mini, by = c('rowname', 'gene_name'), suffix = c('.hlf.sndx', '.hlf.ko') ) 


plc_sndx <- read_tsv('data/processed_data/rna/plc_sndx_v_dmso.tsv')
plc_ko   <- read_tsv('data/processed_data/rna/plc_ko_v_wt.tsv')

plc_sndx_mini <- plc_sndx |> dplyr::select(rowname, gene_name, log2FoldChange, baseMean, padj)
plc_ko_mini   <- plc_ko |> dplyr::select(rowname, gene_name, log2FoldChange, baseMean, padj)

plc_comb <- left_join(plc_sndx_mini, plc_ko_mini, by = c('rowname', 'gene_name'), suffix = c('.plc.sndx', '.plc.ko') ) 

hlf_sndx_up <- hlf_sndx |> filter(padj < 0.05 & log2FoldChange > 0) 
hlf_ko_up   <- hlf_ko |> filter(padj < 0.05 & log2FoldChange > 0) 

# combine the FC From these
convert_to_mini <- function(x) {
  x <- x |> dplyr::select(gene_name, log2FoldChange, padj)
  return(x)
}

# Look at expression of genes in other conditions
hlf_sndx_targets <- hlf_sndx |> filter(gene_name %in% target$gene_name) |> convert_to_mini() |> mutate(set = 'HLF_SNDX')
hlf_ko_targets <- hlf_ko |> filter(gene_name %in% target$gene_name) |> convert_to_mini() |> mutate(set = 'HLF_MEN1-KO') 
plc_sndx_targets <- plc_sndx |> filter(gene_name %in% target$gene_name) |> convert_to_mini() |> mutate(set = 'PLC_SNDX')
plc_ko_targets   <- plc_ko |> filter(gene_name %in% target$gene_name) |> convert_to_mini() |> mutate(set = 'PLC_MEN1-KO')

# combine into 1  
all_lfc <- bind_rows(hlf_sndx_targets, hlf_ko_targets, plc_sndx_targets, plc_ko_targets)
all_lfc

all_lfc_wide <- all_lfc |> 
  dplyr::select( -padj ) |> 
  pivot_wider(id_cols= gene_name, names_from = set, values_from = log2FoldChange) |> arrange(desc(HLF_SNDX))
all_lfc_wide_rn <- all_lfc_wide$gene_name
all_lfc_wide_mat <- all_lfc_wide[,2:5]
all_lfc_wide_mat <- as.matrix(all_lfc_wide_mat)
rownames(all_lfc_wide_mat) <- all_lfc_wide_rn
hm_cols <- circlize::colorRamp2(breaks = c(-2, 0, 2), colors = c('blue', 'white', 'red'))
Heatmap(all_lfc_wide_mat,col= hm_cols )
pairs(all_lfc_wide_mat)

all_lfc_wide |> 
  janitor::clean_names() |> 
  ggplot(aes( x= plc_men1_ko, y = hlf_men1_ko))+
  geom_point()

# FIlter just for genes bound by menin and whether they are targets
hlf_comb <- hlf_comb |> 
  filter(gene_name %in% c(target$gene_name, nontarget$gene_name)) |> 
  mutate(class = case_when(gene_name %in% target$gene_name  ~ 'Target', 
                           gene_name %in% nontarget$gene_name ~ 'Non-target',
                           .default = 'what') )
plc_comb <- plc_comb |> 
  filter(gene_name %in% c(target$gene_name, nontarget$gene_name)) |> 
  mutate(class = case_when(gene_name %in% target$gene_name  ~ 'Target', 
                           gene_name %in% nontarget$gene_name ~ 'Non-target',
                           .default = 'what') )
upset(fromList(list(hlf_sndx = hlf_sndx_targets |> filter(padj < 0.05 & log2FoldChange < 0) |> pull(gene_name), 
               hlf_men.ko = hlf_ko_targets |> filter(padj < 0.05 & log2FoldChange < 0) |> pull(gene_name), 
               plc_sndx =  plc_sndx_targets |> filter(padj < 0.05 & log2FoldChange < 0) |> pull(gene_name), 
               plc_men.ko = plc_ko_targets |> filter(padj < 0.05 & log2FoldChange < 0 ) |> pull(gene_name)) ) )

upset(fromList(list(hlf_sndx = hlf_sndx_targets |> filter(padj < 0.05 & log2FoldChange < 0) |> pull(gene_name), 
                    hlf_men.ko = hlf_ko_targets |> filter(padj < 0.05 & log2FoldChange < 0) |> pull(gene_name) ) ) )

upset(fromList(list(plc_sndx =  plc_sndx_targets |> filter(padj < 0.05 & log2FoldChange < 0) |> pull(gene_name), 
                    plc_men.ko = plc_ko_targets |> filter(padj < 0.05 & log2FoldChange < 0 ) |> pull(gene_name)) ) )

a = list(hlf_sndx = hlf_sndx_targets |> filter(padj < 0.05 & log2FoldChange < 0) |> as.data.frame() |> pull(gene_name), 
     hlf_men.ko = hlf_ko_targets |> filter(padj < 0.05 & log2FoldChange < 0) |> as.data.frame() |> pull(gene_name)) 

# plot expression of target vs non-target
hlf_sndx_target_plot <- ggplot(hlf_comb, aes(x = class, y = log2FoldChange.hlf.sndx)) + 
  geom_boxplot(width = 0.4) + 
  theme_bw() + 
  xlab('') + 
  ylab('Log2 Fold Change SNDX/DMSO')
wilcox.test(log2FoldChange.hlf.sndx ~ class, data = hlf_comb)$p.value

hlf_ko_target_plot <- ggplot(hlf_comb, aes(x = class, y = log2FoldChange.hlf.ko)) + 
  geom_boxplot(width = 0.4) + 
  theme_bw() +
  xlab('') + 
  ylab('Log2 Fold Change MEN1 KO/WT')
wilcox.test(log2FoldChange.hlf.ko ~ class, data = hlf_comb)$p.value

plc_ko_target_plot <- ggplot(plc_comb, aes(x = class, y = log2FoldChange.plc.ko)) + 
  geom_boxplot(width = 0.4) + 
  theme_bw()+
  xlab('') + 
  ylab('Log2 Fold Change SNDX/DMSO')
wilcox.test(log2FoldChange.plc.ko ~ class, data = plc_comb)$p.value

plc_sndx_target_plot <- ggplot(plc_comb, aes(x = class, y = log2FoldChange.plc.sndx)) + 
  geom_boxplot(width = 0.4) + 
  theme_bw()+
  xlab('') + 
  ylab('Log2 Fold Change MEN1 KO/WT')
wilcox.test(log2FoldChange.plc.sndx ~ class, data = plc_comb)$p.value
# save boxplots of these
pdf('plots/cnr/fig3_target_exp_boxplot.pdf', width = 8.5, height = 11)
pageCreate(height = 11, width = 8.5, showGuide = F)
plotGG(hlf_sndx_target_plot, x = 0.5, y =0.5, width =3 , height = 3)
plotText(label = 'HLF SNDX', x = 2.0, y = 0.4)
plotGG(hlf_ko_target_plot, x = 4, y = 0.5, width = 3, height = 3 )
plotText(label = 'HLF MEN1 KO', x = 5.5, y = 0.4)
plotGG(plc_sndx_target_plot, x = 0.5, y = 4, width = 3, height = 3)
plotText(label = 'PLCPRF5 SNDX', x = 2.0, y = 3.9)
plotGG(plc_ko_target_plot, x = 4, y = 4, width = 3, height = 3)
plotText(label = 'PLCPRF5 MEN1 KO', x = 5.5, y = 3.9)
dev.off()

# enrichment
### DEfine gene sets of interest#######
all <- read_tsv('data/processed_data/rna/hlf_sndx_v_dmso.tsv')
target

# Full data
msig <- msigdbr::msigdbr(species = 'Homo sapiens')

# Hallmark
hallmark_db <- msig |> filter(gs_cat == 'H') |> dplyr::select(gs_name, human_gene_symbol)

# Canonical pathways
cp_db       <- msig |> filter(gs_cat == 'C5') |> dplyr::select(gs_name, human_gene_symbol)

# Cell type
cell_db     <- msig |> filter(gs_cat == 'C8') |> dplyr::select(gs_name, human_gene_symbol)

# C6
cancer_db <- msig |> filter(gs_cat == 'C6') |> dplyr::select(gs_name, human_gene_symbol)

#all

# look at the targets
target_h_enr <- enricher(unique(target$gene_name), universe = all$gene_name, TERM2GENE = hallmark_db)
target_cp_enr <- enricher(target$gene_name, universe= all$gene_name, TERM2GENE = cp_db)
target_cell_enr <- enricher(target$gene_name, universe = all$gene_name, TERM2GENE = cell_db)
target_cancer_enr <- enricher(target$gene_name, universe = all$gene_name, TERM2GENE = cancer_db)
target_all_enr <- enricher(target$gene_name, universe = all|> filter(!is.na(padj)) |> pull(gene_name), TERM2GENE = msig |> dplyr::select(gs_name, human_gene_symbol))
target_h_enr |> as_tibble() |> View()
target_cp_enr |> as_tibble() |> View()
target_cell_enr
target_cancer_enr
target_cancer_enr |> as_tibble() |> View()

################################################################################
# GSEA for each

run_gsea <- function(df, term, pvalcut = 1) {
  df <- df |> filter(!is.na(padj))
  print(nrow(df))
  df_vals <- df |> pull(log2FoldChange)
  names(df_vals) <- df |> pull(gene_name)
  gsout <- GSEA(geneList = rev(sort(df_vals)), TERM2GENE = term, by = 'fgsea', pvalueCutoff = pvalcut) 
  return(gsout)
}

gsea_hlf_ko_h <- run_gsea(hlf_ko_mini, hallmark_db)
gsea_hlf_sndx_h <- run_gsea(hlf_sndx_mini, hallmark_db)
gsea_plc_ko_h <- run_gsea(plc_ko_mini, hallmark_db)
gsea_plc_sndx_h <- run_gsea(plc_sndx_mini, hallmark_db)

gsea_plc_sndx_h_tib <- gsea_plc_sndx_h |> as_tibble() |> mutate(set = 'PLC_SNDX')
gsea_plc_ko_h_tib   <- gsea_plc_ko_h |> as_tibble() |> mutate(set = 'PLC_KO')
gsea_hlf_sndx_h_tib <- gsea_hlf_sndx_h |> as_tibble() |> mutate(set = 'HLF_SNDX')
gsea_hlf_ko_h_tib   <- gsea_hlf_ko_h |> as_tibble() |> mutate(set = 'HLF_KO')

all_h <- bind_rows(gsea_plc_sndx_h_tib, gsea_plc_ko_h_tib, gsea_hlf_sndx_h_tib, gsea_hlf_ko_h_tib)
all_h <- all_h |> mutate(short_id = str_replace(ID, 'HALLMARK_', ''))
pick_groups <- all_h |> 
  group_by(ID) |> 
  summarise(min_q = min(qvalue))
pick_groups
picked <- pick_groups |> filter(min_q < 0.05) |> pull(ID)
h_interest <- all_h |>  
  filter(ID %in% picked)
h_interest
h_interest_df <- h_interest |> pivot_wider(id_cols = short_id, names_from = set, values_from = NES)
h_interest_mat <- as.matrix(h_interest_df[,c(2:5)])
rownames(h_interest_mat) <- h_interest_df$short_id

hm_cols <- circlize::colorRamp2(breaks = c(-2,0,2), colors = c('blue', 'black', 'gold'))
hm <- Heatmap(h_interest_mat)
hm_plot <- Heatmap(h_interest_mat, 
                   show_row_names = T,  
                   col = hm_cols, 
                   row_names_gp = gpar(fontsize = 6),
                   column_names_gp = gpar(fontsize = 8), 
                   heatmap_legend_param =list( title = "NES"))

p <- grid::grid.grabExpr(draw(hm_plot))
pdf('plots/tmp_heatmap.pdf', height = 11, width = 8.5)
pageCreate(width = 8.5, height = 11, showGuides = F)
plotGG(p, height = 4.5, width = 5, x = 0.5, y =0.5)
dev.off()


