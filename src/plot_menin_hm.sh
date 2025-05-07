#!/bin/bash

#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00



module purge
module load deeptools/3.5.1

OUTDIR=plots/cnr/coverages/
INDIR=data/processed_data/cnr/

plotHeatmap -m ${INDIR}men1_cov_mat.gz \
   --outFileName ${OUTDIR}men1_cov_plot.pdf\
   --sortRegions keep \
   --colorMap Blues\
   --refPointLabel center

plotHeatmap -m ${INDIR}mll_cov_mat.gz \
   --outFileName ${OUTDIR}mll_cov_plot.pdf\
   --sortRegions keep \
   --colorMap Purples\
   --refPointLabel center

plotHeatmap -m ${INDIR}h3k4me3_cov_mat.gz \
   --outFileName ${OUTDIR}h3k4me3_cov_plot.pdf\
   --sortRegions keep \
   --colorMap Greens \
   --refPointLabel center
plotHeatmap -m ${INDIR}igg_cov_mat.gz \
   --outFileName ${OUTDIR}igg_cov_plot.pdf\
   --sortRegions keep \
   --colorMap Greys\
   --refPointLabel center\
   --zMax 3

plotHeatmap -m ${INDIR}nfyb_cov_mat.gz \
   --outFileName ${OUTDIR}nfyb_cov_plot.pdf\
   --sortRegions keep \
   --colorMap Oranges \
   --refPointLabel center

plotHeatmap -m ${INDIR}tead_cov_mat.gz \
   --outFileName ${OUTDIR}tead_cov_plot.pdf\
   --sortRegions keep \
   --colorMap Reds\
   --refPointLabel center 


#plotHeatmap -m ${INDIR}men1_cov_atac_mat.gz \
#   --outFileName ${OUTDIR}men1_cov_atac_plot.pdf\
#   --sortRegions keep \
#   --colorMap Blues\
#   --refPointLabel center
#
#plotHeatmap -m ${INDIR}mll_cov_atac_mat.gz \
#   --outFileName ${OUTDIR}mll_cov_atac_plot.pdf\
#   --sortRegions keep \
#   --colorMap Purples\
#   --refPointLabel center
#
#plotHeatmap -m ${INDIR}h3k4me3_cov_atac_mat.gz \
#   --outFileName ${OUTDIR}h3k4me3_cov_atac_plot.pdf\
#   --sortRegions keep \
#   --colorMap Greens \
#   --refPointLabel center
#
#plotHeatmap -m ${INDIR}igg_cov_atac_mat.gz \
#   --outFileName ${OUTDIR}igg_cov_atac_plot.pdf\
#   --sortRegions keep \
#   --colorMap Greys\
#   --refPointLabel center\
#   --zMax 3

# PLot at most gained/lost
#plotHeatmap -m ${INDIR}men1_cov_most_atac_mat.gz \
#   --outFileName ${OUTDIR}men1_cov_most_atac_plot.pdf\
#   --sortRegions keep \
#   --colorMap Blues\
#   --refPointLabel center
#
#plotHeatmap -m ${INDIR}mll_cov_most_atac_mat.gz \
#   --outFileName ${OUTDIR}mll_cov_most_atac_plot.pdf\
#   --sortRegions keep \
#   --colorMap Purples\
#   --refPointLabel center
#
#plotHeatmap -m ${INDIR}h3k4me3_cov_most_atac_mat.gz \
#   --outFileName ${OUTDIR}h3k4me3_cov_most_atac_plot.pdf\
#   --sortRegions keep \
#   --colorMap Greens \
#   --refPointLabel center
#
#plotHeatmap -m ${INDIR}igg_cov_most_atac_mat.gz \
#   --outFileName ${OUTDIR}igg_cov_most_atac_plot.pdf\
#   --sortRegions keep \
#   --colorMap Greys\
#   --refPointLabel center\
#   --zMax 3
#

