# Figures for comparison of RNAseq - 
# See Figure 3

suppressPackageStartupMessages({

  library(DESeq2)
  library(ComplexHeatmap)
  library(clusterProfiler)
  library(msigdbr)
  library(tidyverse)
  library(plotgardener)
  library(grid)
  
})


hlf_sndx <- read_tsv('data/processed_data/rna/hlf_sndx_v_dmso.tsv')
hlf_ko   <- read_tsv('data/processed_data/rna/hlf_menko_v_wt.tsv')

hlf_sndx_mini <- hlf_sndx |> select(rowname, gene_name, log2FoldChange, baseMean, padj) 
hlf_ko_mini   <- hlf_ko |> select(rowname, gene_name, log2FoldChange, baseMean, padj)

hlf_comb <- left_join(hlf_sndx_mini, hlf_ko_mini, by = c('rowname', 'gene_name'), suffix = c('.hlf.sndx', '.hlf.ko') ) 


plc_sndx <- read_tsv('data/processed_data/rna/plc_sndx_v_dmso.tsv')
plc_ko   <- read_tsv('data/processed_data/rna/plc_ko_v_wt.tsv')

plc_sndx_mini <- plc_sndx |> select(rowname, gene_name, log2FoldChange, baseMean, padj)
plc_ko_mini   <- plc_ko |> select(rowname, gene_name, log2FoldChange, baseMean, padj)

plc_comb <- left_join(plc_sndx_mini, plc_ko_mini, by = c('rowname', 'gene_name'), suffix = c('.plc.sndx', '.plc.ko') ) 


# enrichment
### DEfine gene sets of interest#######

# Full data
msig <- msigdbr::msigdbr(species = 'Homo sapiens')

# Hallmark
hallmark_db <- msig |> filter(gs_cat == 'H') |> select(gs_name, human_gene_symbol)

# Canonical pathways
#cp_db       <- msig |> filter(gs_cat == 'C5') |> select(gs_name, human_gene_symbol)

# Cell type
#cell_db     <- msig |> filter(gs_cat == 'C8') |> select(gs_name, human_gene_symbol)

# C6
cancer_db <- msig |> filter(gs_cat == 'C6') |> select(gs_name, human_gene_symbol)

# Load data to check
# HLF
hlf_wt_vs_ko <- hlf_ko_mini |> filter(!is.na(padj))
hlf_wt_vs_ko_vals <- hlf_wt_vs_ko$log2FoldChange
names(hlf_wt_vs_ko_vals) <- hlf_wt_vs_ko$gene_name

hlf_drug_vs_veh <- hlf_sndx_mini |> filter(!is.na(padj))
hlf_drug_vs_veh_vals <- hlf_drug_vs_veh$log2FoldChange
names(hlf_drug_vs_veh_vals) <- hlf_drug_vs_veh$gene_name 

# PLC
plc_wt_vs_ko <- plc_ko_mini |> filter(!is.na(padj))
plc_wt_vs_ko_vals <- plc_wt_vs_ko$log2FoldChange
names(plc_wt_vs_ko_vals) <- plc_wt_vs_ko$gene_name

plc_drug_vs_veh <- plc_sndx_mini |> filter(!is.na(padj))
plc_drug_vs_veh_vals <- plc_drug_vs_veh$log2FoldChange
names(plc_drug_vs_veh_vals) <- plc_drug_vs_veh$gene_name 


# RUN GSEA 
hlf_wvk_h <- GSEA(geneList = rev(sort(hlf_wt_vs_ko_vals)), TERM2GENE = hallmark_db, by = 'fgsea', pvalueCutoff = 1)
hlf_wvk_cancer     <- GSEA(geneList = rev(sort(hlf_wt_vs_ko_vals)), TERM2GENE = cancer_db, by = 'fgsea', pvalueCutoff = 1)

plc_wvk_h <- GSEA(geneList = rev(sort(plc_wt_vs_ko_vals)), TERM2GENE = hallmark_db, by = 'fgsea', pvalueCutoff = 1)
plc_wvk_cancer     <- GSEA(geneList = rev(sort(plc_wt_vs_ko_vals)), TERM2GENE = cancer_db, by = 'fgsea', pvalueCutoff = 1)
#wvk_cp       <- GSEA(geneList = rev(sort(wt_vs_ko_vals)), TERM2GENE = cp_db, by = 'fgsea', pvalueCutoff = 1)
#wvk_cell     <- GSEA(geneList = rev(sort(wt_vs_ko_vals)), TERM2GENE = cell_db, by = 'fgsea', pvalueCutoff = 1)

hlf_dvv_h <- GSEA(geneList = rev(sort(hlf_drug_vs_veh_vals)), TERM2GENE = hallmark_db, by = 'fgsea', pvalueCutoff = 1)
hlf_dvv_cancer     <- GSEA(geneList = rev(sort(hlf_drug_vs_veh_vals)), TERM2GENE = cancer_db, by = 'fgsea', pvalueCutoff = 1)

plc_dvv_h <- GSEA(geneList = rev(sort(plc_drug_vs_veh_vals)), TERM2GENE = hallmark_db, by = 'fgsea', pvalueCutoff = 1)
plc_dvv_cancer     <- GSEA(geneList = rev(sort(plc_drug_vs_veh_vals)), TERM2GENE = cancer_db, by = 'fgsea', pvalueCutoff = 1)

#dvv_cp       <- GSEA(geneList = rev(sort(drug_vs_veh_vals)), TERM2GENE = cp_db, by = 'fgsea', pvalueCutoff = 1)
#dvv_cell     <- GSEA(geneList = rev(sort(drug_vs_veh_vals)), TERM2GENE = cell_db, by = 'fgsea', pvalueCutoff = 1)


# Make a single dataframe for each - drop the c8 and cp
hlf_wvk_all <- bind_rows(hlf_wvk_h |> as_tibble(), hlf_wvk_cancer  |> as_tibble())
hlf_dvv_all <- bind_rows(hlf_dvv_h |> as_tibble(), hlf_dvv_cancer  |> as_tibble() ) 
plc_wvk_all <- bind_rows(plc_wvk_h |> as_tibble(), plc_wvk_cancer |> as_tibble() ) 
plc_dvv_all <- bind_rows(plc_dvv_h |> as_tibble(), plc_dvv_cancer |> as_tibble() ) 

# bring drug and ko together
hlf_all_e <- left_join(hlf_wvk_all, hlf_dvv_all, 
                       by = c('ID', 'Description'), suffix = c('.ko', '.drug') ) 
plc_all_e <- left_join(plc_wvk_all, plc_dvv_all, 
                       by = c('ID', 'Description'), suffix = c('.ko', '.drug') ) 

hlf_all_e |> ggplot(aes(x  = NES.ko, 
                   y = NES.drug) ) + 
  geom_point(data  = hlf_all_e |> filter(p.adjust.ko < 0.1 | p.adjust.drug < 0.1), color = 'steelblue', size= 2) + 
  geom_point(data = hlf_all_e |> filter(p.adjust.ko > 0.1 & p.adjust.drug > 0.1), alpha = 0.2, color = 'grey30')

plc_all_e |> ggplot(aes(x  = NES.ko, 
                   y = NES.drug) ) + 
  geom_point(data  = plc_all_e |> filter(p.adjust.ko < 0.1 | p.adjust.drug < 0.1), color = 'steelblue', size= 2) + 
  geom_point(data = plc_all_e |> filter(p.adjust.ko > 0.1 & p.adjust.drug > 0.1), alpha = 0.2, color = 'grey30')

# convert to long
hlf_all_e_long_nes <- hlf_all_e |> 
  select(ID, NES.ko, NES.drug) |> 
  pivot_longer(cols = c(NES.ko,NES.drug),  
               names_to = 'type', 
               values_to = 'NES') |> 
  separate(type, into = c('nes', 'type'), sep = "\\.")
hlf_all_e_long_padj <- hlf_all_e |> 
  select(ID, p.adjust.ko, p.adjust.drug) |> 
  pivot_longer(cols = c(p.adjust.ko, p.adjust.drug), 
               names_to = 'type',
               values_to = 'p.adjust') |> 
  separate(type, into = c('p','adjust','type'), sep = '\\.') 

hlf_all_e_long <- left_join(hlf_all_e_long_nes, hlf_all_e_long_padj, by = c('ID', 'type') )  |> 
  dplyr::select(ID, type, NES, p.adjust) |> 
  mutate(cell = 'HLF')


plc_all_e_long_nes <- plc_all_e |> 
  select(ID, NES.ko, NES.drug, ) |> 
  pivot_longer(cols = c(NES.ko,NES.drug),  
               names_to = 'type', 
               values_to = 'NES') |> 
  separate(type, into  = c('nes', 'type'), sep  ='\\.') 
plc_all_e_long_padj <- plc_all_e |> 
  pivot_longer(cols = c(p.adjust.ko, p.adjust.drug), 
               names_to = 'type',
               values_to = 'p.adjust') |> 
  separate(type, into = c('p', 'adjust', 'type'), sep = '\\.') 
plc_all_e_long <- left_join(plc_all_e_long_nes, plc_all_e_long_padj, by = c('ID', 'type') ) |> 
  dplyr::select(ID, type, NES, p.adjust) |> 
  mutate(cell = 'PLC')
          



# comibne cell lines
all_e_long <- bind_rows(hlf_all_e_long, plc_all_e_long)
all_e_long <- all_e_long |> mutate(grouping = paste(type, cell, sep = '_') )

# Pick 15 categories based on highest sum or lowest sum of NES
all_e_long_picks_up <- all_e_long |> 
  group_by(ID) |> 
  summarise(sum_nes = sum(NES) ) |> 
  arrange(sum_nes) |> 
  tail(25)

all_e_long_picks_down <- all_e_long |> 
  group_by(ID) |> 
  summarise(sum_nes = sum(NES) ) |> 
  arrange(sum_nes) |> 
  head(25)


all_e_long |> 
  filter(startsWith(ID, 'HALLMARK_')) |> 
 filter(ID %in% c(all_e_long_picks_up$ID, all_e_long_picks_down$ID)) |>
  ggplot(aes(x = grouping, y = ID, color = NES, size = -log10(p.adjust))) + 
  geom_point() +
  scale_color_gradient2(low = 'blue', mid = 'black', high = 'gold') + 
  theme_minimal()
  

best_e_dots <- all_e_long |> 
  filter(ID %in% c(all_e_long_picks_up$ID, all_e_long_picks_down$ID)) |>
  ggplot(aes(x = grouping, y = ID, color = NES, size = -log10(p.adjust))) + 
  geom_point() +
  scale_color_gradient2(low = 'blue', mid = 'black', high = 'gold') + 
  theme_minimal()

best_e_dots

ggsave(filename = 'plots/rna/dots_ontolgies.pdf', best_e_dots)
best_e_dots


# heatmap version - separate hallmark and non-hallmark
sig_e_long_mat_hallmark <- all_e_long |> 
  filter(startsWith(ID, 'HALLMARK_')) |> 
  filter(ID %in% c(all_e_long_picks_up$ID, all_e_long_picks_down$ID)) |>
  dplyr::select(-p.adjust, cell, type) |> 
  pivot_wider(id_cols = ID, names_from = grouping, values_from = NES ) |> 
  mutate(ID = str_replace(ID, 'HALLMARK_', ''))
sig_e_long_mat_hallmark_rn <- sig_e_long_mat_hallmark$ID
sig_e_long_mat_hallmark <- as.matrix(sig_e_long_mat_hallmark[,2:5])
rownames(sig_e_long_mat_hallmark) <- sig_e_long_mat_hallmark_rn
hm_hallmark<- Heatmap(sig_e_long_mat_hallmark,
             col = circlize::colorRamp2(colors = c('blue', 'black', 'gold'), breaks = c(-2,0,2) ),
             column_names_gp = grid::gpar(fontsize = 8),
             row_names_gp = grid::gpar(fontsize = 6),
             show_heatmap_legend = F,
             cluster_columns = F
)


hm_hallmark_g <- grid.grabExpr(draw(hm_hallmark))

# same but not hallmark

sig_e_long_mat_other <- all_e_long |> 
  filter(!startsWith(ID, 'HALLMARK_')) |> 
  filter(ID %in% c(all_e_long_picks_up$ID, all_e_long_picks_down$ID)) |>
  dplyr::select(-p.adjust, cell, type) |> 
  pivot_wider(id_cols = ID, names_from = grouping, values_from = NES ) 
sig_e_long_mat_other_rn <- sig_e_long_mat_other$ID
sig_e_long_mat_other <- as.matrix(sig_e_long_mat_other[,2:5])
rownames(sig_e_long_mat_other) <- sig_e_long_mat_other_rn
hm_other<- Heatmap(sig_e_long_mat_other,
                      col = circlize::colorRamp2(colors = c('blue', 'black', 'gold'), breaks = c(-2,0,2) ),
                      column_names_gp = grid::gpar(fontsize = 8),
                      row_names_gp = grid::gpar(fontsize = 6),
                      show_heatmap_legend = F,
                   cluster_columns = F
)
hm_other_g <- grid.grabExpr(draw(hm_other))
# Do enrichment on just the conserved up and down

up_all <- read_tsv('data/processed_data/rna/upregulated_overlap_menko_sndx.tsv')
down_all <- read_tsv('data/processed_data/rna/downregulated_overlap_menko_sndx.tsv')

# Do enrichment by overlap
up_hall <- enricher(up_all$gene_name, universe = hlf_ko$gene_name, TERM2GENE = hallmark_db)
up_cancer <- enricher(up_all$gene_name, universe = hlf_ko$gene_name, TERM2GENE = cancer_db)
down_hall <- enricher(down_all$gene_name, universe = hlf_ko$gene_name, TERM2GENE = hallmark_db)
down_cancer <- enricher(down_all$gene_name, universe = hlf_ko$gene_name, TERM2GENE = cancer_db)

all_enricher <- bind_rows(up_hall |> as_tibble() |> mutate(category = 'up'),
                          up_cancer |> as_tibble() |> mutate(category = 'up'), 
                          down_hall |> as_tibble() |> mutate(category = 'down'), 
                          down_cancer |> as_tibble() |> mutate(category = 'down'))

all_enricher |> 
  arrange(desc(qvalue)) |> 
  mutate(ID = factor(ID, levels = ID)) |>  
  ggplot(aes(x =ID, y = -log10(qvalue), fill = category)) + geom_col()  +coord_flip()

p1 <- gseaplot(hlf_dvv_h, geneSetID = 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', by = 'runningScore', title = 'DMSO vs SNDX-6318') + 
  theme(axis.title = element_text(size = 8), 
        axis.text.x = element_text(size = 7,color = 'grey10'),
        axis.text.y = element_text(size = 7, color = 'grey10' )) 
p2 <- gseaplot(plc_dvv_h, geneSetID = 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', by = 'runningScore', title = 'DMSO vs SNDX-6318 PLC/PRF/5') + 
  theme(axis.title = element_text(size = 8), 
        axis.text.x = element_text(size = 7,color = 'grey10'),
        axis.text.y = element_text(size = 7, color = 'grey10' )) 
p3 <- gseaplot(hlf_wvk_h, geneSetID = 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',by = 'runningScore', title = 'WT vs KO HLF') + 
  theme(axis.title = element_text(size = 8), 
        axis.text.x = element_text(size = 7,color = 'grey10'),
        axis.text.y = element_text(size = 7, color = 'grey10' )) 
p4 <- gseaplot(plc_wvk_h, geneSetID = 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',by = 'runningScore', title = 'WT vs KO PLC/PRF/5') + 
  theme(axis.title = element_text(size = 8), 
        axis.text.x = element_text(size = 7,color = 'grey10'),
        axis.text.y = element_text(size = 7, color = 'grey10' )) 
t1 <- hlf_dvv_h |> as_tibble() |> filter(ID == 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION')
t2 <- plc_dvv_h |> as_tibble() |> filter(ID == 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION')
t3 <- hlf_wvk_h |> as_tibble() |> filter(ID == 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION')
t4 <- plc_wvk_h |> as_tibble() |> filter(ID == 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION')


pdf('plots/rna/ont_best_cats_hm.pdf', height = 11, width = 8.5)
pageCreate(showGuides = F)
plotGG(plot = hm_hallmark_g, x = 0.25, y = 0.5, width = 3, height = 5)
plotGG(plot = hm_other_g, x = 0.25, y = 6, width = 3, height = 5)
plotGG(p1, x = 4.5, y = 0.5, width = 4, height = 1.5)
plotGG(p2, x = 4.5, y = 2, width = 3.5, height = 1.5)
plotGG(p3, x = 4.5, y = 3.5, width = 3.5, height = 1.5)
plotGG(p4, x = 4.5, y = 5, width = 3.5, height = 1.5) 
plotText(paste0('NES: ', round(t1$NES, 2)), x = 7.75, y = 0.85, just = c('right', 'top'))
plotText(paste0('NES: ', round(t2$NES, 2)), x = 7.75, y = 0.85+1.5, just = c('right', 'top'))
plotText(paste0('NES: ', round(t3$NES, 2)), x = 7.75, y = 0.85+ 3, just = c('right', 'top'))
plotText(paste0('NES: ', round(t4$NES, 2)), x = 7.75, y = 0.85+ 4.5, just = c('right', 'top'))
dev.off()



