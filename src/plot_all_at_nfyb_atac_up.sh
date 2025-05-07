#!/bin/bash

#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --mem=16G



module purge
module load deeptools/3.5.1

OUTDIR=plots/integrated/coverages/nfyb_peaks/
INDIR=data/processed_data/integrated/
mkdir -p $OUTDIR

plotHeatmap -m ${INDIR}men1_cov_nfyb_peaks_mat.gz \
   --outFileName ${OUTDIR}men1_cov_nfyb_peaksplot.pdf\
   --sortRegions keep \
   --colorMap Blues\
   --refPointLabel center

plotHeatmap -m ${INDIR}mll_cov_nfyb_peaks_mat.gz \
   --outFileName ${OUTDIR}mll_cov_nfyb_peaks_plot.pdf\
   --sortRegions keep \
   --colorMap Purples\
   --refPointLabel center

plotHeatmap -m ${INDIR}h3k4me3_cov_nfyb_peaks_mat.gz \
   --outFileName ${OUTDIR}h3k4me3_cov_nfyb_peaks_plot.pdf\
   --sortRegions keep \
   --colorMap Greens \
   --refPointLabel center
plotHeatmap -m ${INDIR}igg_cov_nfyb_peaks_mat.gz \
   --outFileName ${OUTDIR}igg_cov_nfyb_peaks_plot.pdf\
   --sortRegions keep \
   --colorMap Greys\
   --refPointLabel center\
   --zMax 3

plotHeatmap -m ${INDIR}nfyb_cov_nfyb_peaks_mat.gz \
   --outFileName ${OUTDIR}nfyb_cov_nfyb_peaks_plot.pdf\
   --sortRegions keep \
   --colorMap Oranges \
   --refPointLabel center

plotHeatmap -m ${INDIR}all_cov_nfyb_peaks_mat.gz \
   --outFileName ${OUTDIR}all_cov_nfyb_peaks_plot.pdf\
   --refPointLabel center\
   --colorMap Oranges Oranges Blues Blues Purples Purples Greens Greens PuBu PuBu Reds Reds Greys Greys Greys Greys 




