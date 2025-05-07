# Figures for comparison of RNAseq - 
# See Figure 2

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(ComplexHeatmap)
  library(clusterProfiler)
  library(msigdbr)
})


hlf_sndx <- read_tsv('data/processed_data/rna/hlf_sndx_v_dmso.tsv')
hlf_ko   <- read_tsv('data/processed_data/rna/hlf_menko_v_wt.tsv')

hlf_sndx_mini <- hlf_sndx |> select(rowname, gene_name, log2FoldChange, baseMean, padj) 
hlf_ko_mini   <- hlf_ko |> select(rowname, gene_name, log2FoldChange, baseMean, padj)

hlf_comb <- left_join(hlf_sndx_mini, hlf_ko_mini, by = c('rowname', 'gene_name'), suffix = c('.sndx', '.ko') ) 


hlf_comb |> 
  filter(!is.na(padj.sndx)) |> 
  filter(!is.na(padj.ko)) |> 
  ggplot(aes(x = log2FoldChange.sndx, y = log2FoldChange.ko) ) + 
  geom_point()


hlf_comb |> 
  filter(!is.na(padj.sndx)) |> 
  filter(!is.na(padj.ko)) |> 
  filter(baseMean.sndx > 250) |> 
  filter(baseMean.ko > 250) |> 
  ggplot(aes(x = log2FoldChange.sndx, y = log2FoldChange.ko) ) + 
  geom_point()

hlf_sndx
# enrichment
### DEfine gene sets of interest#######

# Full data
msig <- msigdbr::msigdbr(species = 'Homo sapiens')

# Hallmark
hallmark_db <- msig |> filter(gs_cat == 'H') |> select(gs_name, human_gene_symbol)

# Canonical pathways
cp_db       <- msig |> filter(gs_cat == 'C5') |> select(gs_name, human_gene_symbol)

# Cell type
cell_db     <- msig |> filter(gs_cat == 'C8') |> select(gs_name, human_gene_symbol)

# C6
cancer_db <- msig |> filter(gs_cat == 'C6') |> select(gs_name, human_gene_symbol)

# Load data to check

wt_vs_ko <- hlf_ko_mini |> filter(!is.na(padj))
wt_vs_ko_vals <- wt_vs_ko$log2FoldChange
names(wt_vs_ko_vals) <- wt_vs_ko$gene_name

drug_vs_veh <- hlf_sndx_mini |> filter(!is.na(padj))
drug_vs_veh_vals <- drug_vs_veh$log2FoldChange
names(drug_vs_veh_vals) <- drug_vs_veh$gene_name 

wvk_h <- GSEA(geneList = rev(sort(wt_vs_ko_vals)), TERM2GENE = hallmark_db, by = 'fgsea', pvalueCutoff = 1)
wvk_cp       <- GSEA(geneList = rev(sort(wt_vs_ko_vals)), TERM2GENE = cp_db, by = 'fgsea', pvalueCutoff = 1)
wvk_cell     <- GSEA(geneList = rev(sort(wt_vs_ko_vals)), TERM2GENE = cell_db, by = 'fgsea', pvalueCutoff = 1)
wvk_cancer     <- GSEA(geneList = rev(sort(wt_vs_ko_vals)), TERM2GENE = cancer_db, by = 'fgsea', pvalueCutoff = 1)

dvv_h <- GSEA(geneList = rev(sort(drug_vs_veh_vals)), TERM2GENE = hallmark_db, by = 'fgsea', pvalueCutoff = 1)
dvv_cp       <- GSEA(geneList = rev(sort(drug_vs_veh_vals)), TERM2GENE = cp_db, by = 'fgsea', pvalueCutoff = 1)
dvv_cell     <- GSEA(geneList = rev(sort(drug_vs_veh_vals)), TERM2GENE = cell_db, by = 'fgsea', pvalueCutoff = 1)
dvv_cancer     <- GSEA(geneList = rev(sort(drug_vs_veh_vals)), TERM2GENE = cancer_db, by = 'fgsea', pvalueCutoff = 1)

# Make a single dataframe for each - drop the c8
wvk_all <- bind_rows(wvk_h |> as_tibble(), wvk_cp |> as_tibble(), wvk_cancer  |> as_tibble())
dvv_all <- bind_rows(dvv_h |> as_tibble(), dvv_cp |> as_tibble(), dvv_cancer  |> as_tibble() ) 

# bring drug and ko together
all_e <- left_join(wvk_all, dvv_all, 
                   by = c('ID', 'Description'), suffix = c('.ko', '.drug') ) 

all_e |> ggplot(aes(x  = NES.ko, 
                   y = NES.drug) ) + 
  geom_point(data  = all_e |> filter(p.adjust.ko < 0.1 | p.adjust.drug < 0.1), color = 'steelblue', size= 2) + 
  geom_point(data = all_e |> filter(p.adjust.ko > 0.1 & p.adjust.drug > 0.1), alpha = 0.2, color = 'grey30')

# convert to long
all_e_long <- all_e |> 
  select(ID, NES.ko, p.adjust.ko, NES.drug, p.adjust.drug) |> 
  pivot_longer(cols = c(NES.ko,NES.drug),  
               names_to = 'type', 
               values_to = 'NES') |> 
  pivot_longer(cols = c(p.adjust.ko, p.adjust.drug), 
               names_to = 'p.type',
               values_to = 'p.adjust')

# this works, but ther eis about 3k to look at
sig_e_long <- all_e |> 
  select(ID, NES.ko, p.adjust.ko, NES.drug, p.adjust.drug) |> 
  filter(p.adjust.ko < 0.1 | p.adjust.drug) |> 
  pivot_longer(cols = c(NES.ko,NES.drug),  
               names_to = 'type', 
               values_to = 'NES') |> 
  pivot_longer(cols = c(p.adjust.ko, p.adjust.drug), 
               names_to = 'p.type',
               values_to = 'p.adjust') |> 
  separate(type, into = c('drop', 'type'), sep = '\\.') |> 
  separate(p.type, into = c('drop', 'drop2', 'p.type'), sep = '\\.') |> 
  mutate(keep  = type == p.type)  |> 
  filter(keep == TRUE) |> 
  select(ID, type, NES, p.adjust)

sig_e_long |> 
  group_by(ID) |> 
  mutate(max_NES = max(abs(NES)) ) |> 
  arrange(desc(max_NES)) |> 
  head(100) |> 
  ggplot(aes(x = type, y = ID, color = NES, size = -log10(p.adjust))) + 
  geom_point() +
  scale_color_gradient2(low = 'blue', mid = 'black', high = 'gold')+
  theme_minimal()
# heatmap version
sig_e_long_mat <- sig_e_long |> 
  select(-p.adjust) |> 
  pivot_wider(id_cols = ID,names_from = type, values_from = NES )
sig_e_long_mat_rn <- sig_e_long_mat$ID
sig_e_long_mat <- as.matrix(sig_e_long_mat[,2:3])
rownames(sig_e_long_mat) <- sig_e_long_mat_rn
Heatmap(sig_e_long_mat)

h_all_long <- h_all |> 
  select(ID, NES.ko, p.adjust.ko, NES.drug, p.adjust.drug) |> 
  pivot_longer(cols = c(NES.ko,NES.drug),  
              names_to = 'type', 
              values_to = 'NES') |> 
  pivot_longer(cols = c(p.adjust.ko, p.adjust.drug), 
               names_to = 'p.type',
               values_to = 'p.adjust')
h_all_long |> 
  ggplot(aes(y = ID, x = type, fill = NES, size = -log10(p.adjust))) + 
  geom_point(pch = 21)  +
  scale_fill_gradient2(low = 'blue', mid = 'black', high = 'gold') + 
  theme_minimal() + 
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size=8))
  