#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --mail-user=jraab@med.unc.edu
#SBATCH --job-name deeptools-jraab

module purge
module load deeptools/3.5.1

GAIN_peaks='data/derived_data/atac/gain_access.bed'
LOST_peaks='data/derived_data/atac/lost_access.bed'
GAINMOST_peaks='data/derived_data/atac/gain_most_access.tsv'
LOSTMOST_peaks='data/derived_data/atac/lost_most_access.tsv'

INPUTDIR='data/source_data/atac/bw/'
CNRDIR='data/source_data/cnr/bw/'
OUTDIR='data/processed_data/atac/'
OUTPLOT='plots/atac/'

DMSO=${INPUTDIR}groupATAC_DMSO_mean.bw
SNDX=${INPUTDIR}groupATAC_SNDX_mean.bw

MENIN=$(find $CNRDIR -name "groupHLF_*_Menin*mean.bw" |sort)
MLL=$(find $CNRDIR -name "groupHLF_*_MLL*mean.bw" |sort)
H3K4me3=$(find $CNRDIR -name "groupHLF_*_H3K4me3*mean.bw" |sort)
IGG=$(find $CNRDIR -name "groupHLF_*_IgG*mean.bw" |sort)
TEAD=$(find $CNRDIR -name "groupHLF_*_TEAD*mean.bw" |sort)


SORTORDER=${OUTDIR}atac_diff_sort_order.bed
mkdir -p $OUTDIR
mkdir -p $OUTPLOT

computeMatrix reference-point -S $DMSO $SNDX\
      -R $GAIN_peaks $LOST_peaks \
   --outFileName ${OUTDIR}atac_diff_cov_mat.gz\
   --referencePoint center\
   -a 2000\
   -b 2000\
   -bs 10\
   --sortRegions descend\
   -p 8\
   --outFileSortedRegions $OUTDIR/atac_diff_sort_order.bed

computeMatrix reference-point -S $MENIN\
   -R $SORTORDER\
   --outFileName ${OUTDIR}menin_atac_cov_mat.gz\
   --referencePoint center\
   -a 2000\
   -b 2000\
   -bs 10\
   --sortRegions keep\
   -p 8 


computeMatrix reference-point -S $MLL\
   -R $SORTORDER\
   --outFileName ${OUTDIR}mll_atac_cov_mat.gz\
   --referencePoint center\
   -a 2000\
   -b 2000\
   -bs 10\
   --sortRegions keep\
   -p 8 

computeMatrix reference-point -S $H3K4me3\
   -R $SORTORDER\
   --outFileName ${OUTDIR}h3k4me3_atac_cov_mat.gz\
   --referencePoint center\
   -a 2000\
   -b 2000\
   -bs 10\
   --sortRegions keep\
   -p 8 

computeMatrix reference-point -S $IGG\
   -R $SORTORDER\
   --outFileName ${OUTDIR}igg_atac_cov_mat.gz\
   --referencePoint center\
   -a 2000\
   -b 2000\
   -bs 10\
   --sortRegions keep\
   -p 8 


computeMatrix reference-point -S $TEAD\
   -R $SORTORDER\
   --outFileName ${OUTDIR}tead_atac_cov_mat.gz\
   --referencePoint center\
   -a 2000\
   -b 2000\
   -bs 10\
   --sortRegions keep\
   -p 8 

plotHeatmap -m ${OUTDIR}atac_diff_cov_mat.gz\
   --outFileName ${OUTPLOT}/atac_diff.pdf\
   --sortRegions keep\
   --colorMap Greys Greys 
plotHeatmap -m ${OUTDIR}mll_atac_cov_mat.gz\
   --outFileName ${OUTPLOT}/mll_atac_cov.pdf\
   --sortRegions keep\
   --colorMap Purples Purples
plotHeatmap -m ${OUTDIR}menin_atac_cov_mat.gz\
   --outFileName ${OUTPLOT}/menin_atac_cov.pdf\
   --sortRegions keep\
   --colorMap Blues Blues

plotHeatmap -m ${OUTDIR}h3k4me3_atac_cov_mat.gz\
   --outFileName ${OUTPLOT}/h3k4me3_atac_cov.pdf\
   --sortRegions keep\
   --colorMap Greens Greens

plotHeatmap -m ${OUTDIR}igg_atac_cov_mat.gz\
   --outFileName ${OUTPLOT}/igg_atac_cov.pdf\
   --sortRegions keep\
   --colorMap Greys Greys

plotHeatmap -m ${OUTDIR}tead_atac_cov_mat.gz\
   --outFileName ${OUTPLOT}/tead_atac_cov.pdf\
   --sortRegions keep\
   --colorMap Greys Greys





#  computeMatrix reference-point -S $DMSO $SNDX\
#     -R $GAINMOST_peaks $LOSTMOST_peaks \
#     --outFileName ${OUTDIR}atac_most_diff_cov_mat.gz\
#     --referencePoint center\
#     -a 1500\
#     -b 1500\
#     -bs 10\
#     --sortRegions keep\
#     -p 8\
#     --outFileSortedRegions $OUTDIR/atac_most_diff_sort_order.bed
#  
#  plotHeatmap -m ${OUTDIR}atac_most_diff_cov_mat.gz\
#     --outFileName ${OUTPLOT}/atac_most_diff.pdf\
#     --sortRegions keep\
#     --colorMap Greys Greys
#  
#  
