# RNA-seq PCA
library(DESeq2)
library(tidyverse)
library(plotgardener)
load('data/processed_data/rna/des_obj_jr.Rda')

vst_hlf <- varianceStabilizingTransformation(des_hlf_men)
vst_plc <- varianceStabilizingTransformation(des_plc)

hlf_pca <- plotPCA(vst_hlf, intgroup = 'treatment') + scale_color_manual(values= c('grey10', 'grey50', 'darkblue', 'steelblue'))+theme_bw()
plc_pca <- plotPCA(vst_plc, intgroup = 'treatment')+ scale_color_manual(values= c('grey10', 'grey50', 'darkblue', 'steelblue'))+theme_bw()

pdf('plots/rna/rnaseq_pca.pdf', height = 11, width = 8.5)
pageCreate(showGuides = F)
plotGG(hlf_pca, x = 0.5, y = 0.5, height = 3.5, width = 3.5)
plotGG(plc_pca, x = 0.5, y = 4.5, height = 3.5, width = 3.5)
dev.off()



hlf_men_counts <- counts(des_hlf_men, norm = F)
plc_men_counts <- counts(des_plc, norm = F)


write_tsv(x = as.data.frame(hlf_men_counts) |> rownames_to_column(), 'data/derived_data/rna/hlf_rnaseq_rawcounts.tsv')
write_tsv(as.data.frame(plc_men_counts) |> rownames_to_column(), 'data/derived_data/rna/plcprf5_rnaseq_rawcounts.tsv')
