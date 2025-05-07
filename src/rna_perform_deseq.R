# RNAseq for PLC/PRF/5 and HLF WT vs KO and Drug vs Vehicle
# Final RNA seq
# Author: Jesse Raab
# Date Updated: 2025-01-20

suppressPackageStartupMessages({
  library(tidyverse) 
  library(DESeq2) # main package for differential expression
  # Helper functions for reading in count data 
  library(limma)
  library(tximeta)
  library(janitor)
})


# Lazy but do this for the time being
input_dir <- 'data/source_data/rna/quants/'
output_dir <- 'data/processed_data/rna/' # for saving the DES objects
final_dir <- 'data/derived_data/rna/'

for (i in c(output_dir, final_dir)) { 
  if(!dir.exists(i))  {dir.create(i, recursive = T)}
}


# Import 
coldata <- read_csv('full_rna_samples.csv') |> clean_names()
coldata
coldata$names <- paste(coldata$sample_number, coldata$sample_id, coldata$cell_line, coldata$treatment, coldata$replicate, sep = '_') 
coldata$files <- file.path(input_dir, coldata$names, '/quant.sf'  ) 
coldata

# This should be done stratified for ease of interpretiation
cda_plc <- coldata |> filter(sample_number %in% 740:751)
cda_hlf_men <- coldata |> filter(sample_number %in% 540:551)


# use tximeta to import salmon data
txi_plc <- tximeta(coldata = cda_plc, type = 'salmon')
txi_hlf_men <- tximeta(coldata = cda_hlf_men, type = 'salmon')

# sumamrise Tx to Genes
se_plc <- summarizeToGene(txi_plc)
keep_plc <- edgeR::filterByExpr(se_plc, min.count = 10, min.total.count = 100, group = 'treatment')
se_plc <- se_plc[keep_plc]
dds_plc <- DESeqDataSet(se_plc, design =  ~   treatment)
des_plc <- DESeq(dds_plc)

se_hlf_men <- summarizeToGene(txi_hlf_men)
keep_hlf <- edgeR::filterByExpr(se_hlf_men, min.count = 10, min.total.count = 100, group = 'treatment')
dds_hlf_men <- DESeqDataSet(se_hlf_men, design=  ~   treatment)
des_hlf_men <- DESeq(dds_hlf_men)

#shrunken fc
#HLF MEN or SNDX
resultsNames(des_hlf_men)
hlf_sndx_v_dmso<-  lfcShrink(des_hlf_men, coef=4, type = 'apeglm', saveCols = 2, format = 'GRanges')|> keepStandardChromosomes(pruning.mode = 'coarse')
#relevel to look at ko
des_hlf_men$treatment <- relevel(des_hlf_men$treatment, 'NTC')
des_hlf_men <- nbinomWaldTest(des_hlf_men)
resultsNames(des_hlf_men)
hlf_menko_v_wt <- lfcShrink(des_hlf_men, coef = 3, type = 'apeglm', saveCols = 2, format = 'GRanges') |> keepStandardChromosomes(pruning.mode = 'coarse')

# PLCPRF5
resultsNames(des_plc)
plc_sndx_v_dmso <- lfcShrink(des_plc, coef = 4, type = 'apeglm', saveCols = 2, format = 'GRanges') |> keepStandardChromosomes(pruning.mode = 'coarse')
des_plc$treatment <- relevel(des_plc$treatment, 'NTC') 
des_plc <- nbinomWaldTest(des_plc)
resultsNames(des_plc)
plc_ko_v_wt <- lfcShrink(des_plc, coef = 3, type = 'apeglm', saveCols = 2, format = 'GRanges') |> keepStandardChromosomes(pruning.mode = 'coarse')

################################################################################ 

ggplot(plc_ko_v_wt |> as.data.frame(), aes(x = log2(baseMean), y =log2FoldChange, color = padj < 0.05))+ geom_point()
ggplot(plc_sndx_v_dmso |> as.data.frame(), aes(x = log2(baseMean), y =log2FoldChange, color = padj < 0.05))+ geom_point()


plc_all<- left_join(plc_ko_v_wt |> as.data.frame() |> rownames_to_column(), 
                    plc_sndx_v_dmso |> as.data.frame() |> rownames_to_column(), 
                    by= 'rowname', suffix = c('.ko', '.drug'))
plc_all |> ggplot(aes(x = log2FoldChange.ko, y = log2FoldChange.drug)) + 
  geom_point()

cor.test(plc_all$log2FoldChange.ko, plc_all$log2FoldChange.drug)

hlf_men_all<- left_join(hlf_menko_v_wt |> as.data.frame() |> rownames_to_column(), 
                        hlf_sndx_v_dmso |> as.data.frame() |> rownames_to_column(), 
                        by= 'rowname', suffix = c('.men_ko', '.drug'))
hlf_men_all |> ggplot(aes(x = log2FoldChange.men_ko, y = log2FoldChange.drug)) + 
  geom_point(aes(color = (log2FoldChange.men_ko < 0.05)) )


cor.test(hlf_men_all$log2FoldChange.men_ko, hlf_men_all$log2FoldChange.drug)

###############################################################################
# QC plots for supplement
vst_plc <- varianceStabilizingTransformation(dds_plc, blind = T)
vst_hlf_men <- varianceStabilizingTransformation(dds_hlf_men, blind = T)

plotPCA(vst_plc, intgroup = c('treatment', 'replicate'), returnData = T) |> 
  ggplot(aes(x = PC1, y = PC2, color = treatment, shape = as.factor(replicate))) + geom_point(size = 3)

plotPCA(vst_hlf_men, intgroup = c('treatment', 'replicate'), returnData = T) |> 
  ggplot(aes(x = PC1, y = PC2, color = treatment, shape = as.factor(replicate))) + geom_point(size = 3)



# Save the rda from the deseq datasets
save(des_hlf_ash, des_hlf_men, des_plc, file = 'data/processed_data/rna/des_obj_jr.Rda')

write_tsv(plc_ko_v_wt |> as.data.frame() |> rownames_to_column(), file = 'data/processed_data/rna/plc_ko_v_wt.tsv')
write_tsv(plc_sndx_v_dmso |> as.data.frame() |> rownames_to_column(), file = 'data/processed_data/rna/plc_sndx_v_dmso.tsv')
write_tsv(hlf_menko_v_wt |> as.data.frame() |> rownames_to_column(), file = 'data/processed_data/rna/hlf_menko_v_wt.tsv')
write_tsv(hlf_sndx_v_dmso |> as.data.frame() |> rownames_to_column() , file = 'data/processed_data/rna/hlf_sndx_v_dmso.tsv')
write_tsv(hlf_ash2l1_v_wt |> as.data.frame() |> rownames_to_column(), file = 'data/processed_data/rna/hlf_ash2l-1_v_wt.tsv')
write_tsv(hlf_ash2l2_v_wt |> as.data.frame() |> rownames_to_column(), file = 'data/processed_data/rna/hlf_ash2l-2_v_wt.tsv')

