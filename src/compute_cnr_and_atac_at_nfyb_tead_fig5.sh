#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --mail-user=jraab@med.unc.edu
#SBATCH --job-name deeptools-jraab

module purge
module load deeptools/3.5.1

TEAD_PEAKS='data/derived_data/peak_sets/cnr/tead_atac_gained.bed'
NFYB_PEAKS='data/derived_data/peak_sets/cnr/nfyb_atac_gained.bed'

ATAC_MOST=data/processed_data/atac/atac_most_diff_sort_order.bed
TEAD_SORTORDER='data/processed_data/cnr/sort_order_at_tead.bed'
NFYB_SORTORDER='data/processed_data/cnr/sort_order_at_nfyb.bed'
INPUTDIR='data/source_data/cnr/bw/'
OUTDIR='data/processed_data/cnr/'

ATAC_DMSO='data/source_data/atac/bw/groupDMSO_mean.bw'
ATAC_SNDX='data/source_data/atac/bw/groupSNDX_mean.bw'
MEN1=$(find $INPUTDIR -name "groupHLF_*_Menin*mean.bw" | sort)
MLL=$(find $INPUTDIR -name "groupHLF_*_MLL*mean.bw" | sort)
H3K4me3=$(find $INPUTDIR -name "groupHLF_*_H3K4me3*mean.bw"| sort)
NFYB=$(find  $INPUTDIR -name "groupHLF_*_NFYB*mean.bw" |sort)
TEAD=$(find $INPUTDIR -name "groupHLF_*_TEAD4*mean.bw" |sort)
IGG=$(find $INPUTDIR -name "groupHLF_*_IgG*mean.bw" | sort)

PLOTDIR=plots/cnr/coverages/


echo $MEN1
echo $MLL
echo $H3K4me3
echo $IGG

# Get NFYB Gained sort order
computeMatrix reference-point -S $NFYB\
   -R $ATAC_MOST \
   --outFileName ${OUTDIR}nfyb_cov_atac_most_mat.gz\
   --referencePoint center\
   -a 1500\
   -b 1500\
   -bs 10\
   --sortRegions keep\
   -p 8

# Get TEAD Gained sort order 
computeMatrix reference-point -S $TEAD\
   -R $ATAC_MOST \
   --outFileName ${OUTDIR}tead_cov_atac_most_mat.gz\
   --referencePoint center\
   -a 1500\
   -b 1500\
   -bs 10\
   --sortRegions keep\
   -p 8 

#################################################################################
## Get signal for all at NFYB GAINED
## get atac signal at gained vs lost
computeMatrix reference-point -S $ATAC_DMSO $ATAC_SNDX\
   -R $ATAC_MOST \
   --outFileName ${OUTDIR}atac_cov_atac_most_mat.gz\
   --referencePoint center\
   -a 1500\
   -b 1500\
   -bs 10\
   --sortRegions keep\
   -p 8 

computeMatrix reference-point -S $MEN1\
   -R $ATAC_MOST \
   --outFileName ${OUTDIR}men1_cov_atac_most_mat.gz\
   --referencePoint center\
   -a 1500\
   -b 1500\
   -bs 10\
   --sortRegions keep\
   -p 8 

computeMatrix reference-point -S $MLL\
   -R $ATAC_MOST\
   --outFileName ${OUTDIR}mll_cov_atac_most_mat.gz\
   --referencePoint center\
   -a 1500\
   -b 1500\
   -bs 10\
   --sortRegions keep\
   -p 8 

computeMatrix reference-point -S $H3K4me3\
   -R $ATAC_MOST\
   --outFileName ${OUTDIR}h3k4me3_cov_atac_most_mat.gz\
   --referencePoint center\
   -a 1500\
   -b 1500\
   -bs 10\
   --sortRegions keep\
   -p 8 

computeMatrix reference-point -S $IGG\
   -R $ATAC_MOST\
   --outFileName ${OUTDIR}igg_cov_atac_most_mat.gz\
   --referencePoint center\
   -a 1500\
   -b 1500\
   -bs 10\
   --sortRegions keep\
   -p 8 

computeMatrix reference-point -S $TEAD\
   -R $ATAC_MOST\
   --outFileName ${OUTDIR}tead_cov_atac_most_mat.gz\
   --referencePoint center\
   -a 1500\
   -b 1500\
   -bs 10\
   --sortRegions keep\
   -p 8 



##################################################################################
## Get signal for all at TEAD Gained
#computeMatrix reference-point -S $ATAC_DMSO $ATAC_SNDX\
#   -R $TEAD_SORTORDER \
#   --outFileName ${OUTDIR}atac_cov_tead_gained_mat.gz\
#   --referencePoint center\
#   -a 1500\
#   -b 1500\
#   -bs 10\
#   --sortRegions keep\
#   -p 8 
#
#computeMatrix reference-point -S $MEN1\
#   -R $TEAD_SORTORDER \
#   --outFileName ${OUTDIR}men1_cov_tead_gained_mat.gz\
#   --referencePoint center\
#   -a 1500\
#   -b 1500\
#   -bs 10\
#   --sortRegions descend\
#   -p 8 
#
#computeMatrix reference-point -S $MLL\
#   -R $TEAD_SORTORDER\
#   --outFileName ${OUTDIR}mll_cov_tead_gained_mat.gz\
#   --referencePoint center\
#   -a 1500\
#   -b 1500\
#   -bs 10\
#   --sortRegions keep\
#   -p 8 
#
#computeMatrix reference-point -S $H3K4me3\
#   -R $TEAD_SORTORDER\
#   --outFileName ${OUTDIR}h3k4me3_cov_tead_gained_mat.gz\
#   --referencePoint center\
#   -a 1500\
#   -b 1500\
#   -bs 10\
#   --sortRegions keep\
#   -p 8 
#
#computeMatrix reference-point -S $IGG\
#   -R $TEAD_SORTORDER\
#   --outFileName ${OUTDIR}igg_cov_tead_gained_mat.gz\
#   --referencePoint center\
#   -a 1500\
#   -b 1500\
#   -bs 10\
#   --sortRegions keep\
#   -p 8 
#
#computeMatrix reference-point -S $NFYB\
#   -R $TEAD_SORTORDER\
#   --outFileName ${OUTDIR}nfyb_cov_tead_gained_mat.gz\
#   --referencePoint center\
#   -a 1500\
#   -b 1500\
#   -bs 10\
#   --sortRegions keep\
#   -p 8 
#
#
#
################################################################################
#plot heatmaps of these
plotHeatmap -m ${OUTDIR}nfyb_cov_atac_most_mat.gz\
   --outFileName ${PLOTDIR}nfyb_cov_atac_most.pdf\
   --sortRegions keep \
   --colorMap Reds\
   --refPointLabel center\

plotHeatmap -m ${OUTDIR}atac_cov_atac_most_mat.gz\
   --outFileName ${PLOTDIR}atac_cov_atac_most.pdf\
   --sortRegions keep \
   --colorMap Greys\
   --refPointLabel center\

plotHeatmap -m ${OUTDIR}men1_cov_atac_most_mat.gz\
   --outFileName ${PLOTDIR}men1_cov_atac_most.pdf\
   --sortRegions keep \
   --colorMap Blues\
   --refPointLabel center\

plotHeatmap -m ${OUTDIR}mll_cov_atac_most_mat.gz\
   --outFileName ${PLOTDIR}mll_cov_atac_most.pdf\
   --sortRegions keep \
   --colorMap Purples\
   --refPointLabel center\

plotHeatmap -m ${OUTDIR}h3k4me3_cov_atac_most_mat.gz\
   --outFileName ${PLOTDIR}h3k4me3_cov_atac_most.pdf\
   --sortRegions keep \
   --colorMap Greens\
   --refPointLabel center\

plotHeatmap -m ${OUTDIR}igg_cov_atac_most_mat.gz\
   --outFileName ${PLOTDIR}igg_cov_atac_most.pdf\
   --sortRegions keep \
   --colorMap Greys\
   --refPointLabel center\

plotHeatmap -m ${OUTDIR}tead_cov_atac_most_mat.gz\
   --outFileName ${PLOTDIR}tead_cov_atac_most.pdf\
   --sortRegions keep \
   --colorMap Oranges\
   --refPointLabel center\


#plotHeatmap -m ${OUTDIR}tead_cov_tead_gained_mat.gz\
#   --outFileName ${PLOTDIR}tead_cov_tead_gained.pdf\
#   --sortRegions keep \
#   --colorMap Oranges\
#   --refPointLabel center\
#
#plotHeatmap -m ${OUTDIR}atac_cov_tead_gained_mat.gz\
#   --outFileName ${PLOTDIR}atac_cov_tead_gained.pdf\
#   --sortRegions keep \
#   --colorMap Greys\
#   --refPointLabel center\
#
#plotHeatmap -m ${OUTDIR}men1_cov_tead_gained_mat.gz\
#   --outFileName ${PLOTDIR}men1_cov_tead_gained.pdf\
#   --sortRegions keep \
#   --colorMap Blues\
#   --refPointLabel center\
#
#plotHeatmap -m ${OUTDIR}mll_cov_tead_gained_mat.gz\
#   --outFileName ${PLOTDIR}mll_cov_tead_gained.pdf\
#   --sortRegions keep \
#   --colorMap Purples\
#   --refPointLabel center\
#
#plotHeatmap -m ${OUTDIR}h3k4me3_cov_tead_gained_mat.gz\
#   --outFileName ${PLOTDIR}h3k4me3_cov_nfyb_gained.pdf\
#   --sortRegions keep \
#   --colorMap Greens\
#   --refPointLabel center\
#
#plotHeatmap -m ${OUTDIR}igg_cov_tead_gained_mat.gz\
#   --outFileName ${PLOTDIR}igg_cov_nfyb_gained.pdf\
#   --sortRegions keep \
#   --colorMap Greys\
#   --refPointLabel center\
#
#plotHeatmap -m ${OUTDIR}nfyb_cov_tead_gained_mat.gz\
#   --outFileName ${PLOTDIR}nfyb_cov_tead_gained.pdf\
#   --sortRegions keep \
#   --colorMap Reds\
#   --refPointLabel center\
#
#
