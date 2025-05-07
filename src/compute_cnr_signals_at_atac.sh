#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --mail-user=jraab@med.unc.edu
#SBATCH --job-name deeptools-jraab

module purge
module load deeptools/3.5.1

peaks='data/processed_data/atac/atac_diff_sort_order.bed'
most_peaks='data/processed_data/atac/atac_most_diff_sort_order.bed'
SORTORDER='data/processed_data/atac/sort_order_at_atac.bed'
INPUTDIR='data/source_data/cnr/group_norm_bins/bw'
OUTDIR='data/processed_data/atac/'

MEN1=$(find $INPUTDIR -name "groupHLF_*_Menin*mean.bw" | sort)
MLL=$(find $INPUTDIR -name "groupHLF_*_MLL*mean.bw" | sort)
H3K4me3=$(find $INPUTDIR -name "groupHLF_*_H3K4me3*mean.bw"| sort)
IGG=$(find $INPUTDIR -name "groupHLF_*_IgG*mean.bw" | sort)


echo $MEN1
echo $MLL
echo $H3K4me3
echo $IGG

#computeMatrix reference-point -S $MEN1\
#   -R $peaks \
#   --outFileName ${OUTDIR}men1_cov_atac_mat.gz\
#   --referencePoint center\
#   -a 2000\
#   -b 2000\
#   -bs 10\
#   --sortRegions descend\
#   -p 8\
#   --outFileSortedRegions $SORTORDER
#
#
#computeMatrix reference-point -S $MLL\
#   -R $peaks\
#   --outFileName ${OUTDIR}mll_cov_atac_mat.gz\
#   --referencePoint center\
#   -a 2000\
#   -b 2000\
#   -bs 10\
#   --sortRegions keep\
#   -p 8 
#
#computeMatrix reference-point -S $H3K4me3\
#   -R $peaks\
#   --outFileName ${OUTDIR}h3k4me3_cov_atac_mat.gz\
#   --referencePoint center\
#   -a 2000\
#   -b 2000\
#   -bs 10\
#   --sortRegions keep\
#   -p 8 
#
#computeMatrix reference-point -S $IGG\
#   -R $peaks\
#   --outFileName ${OUTDIR}igg_cov_atac_mat.gz\
#   --referencePoint center\
#   -a 2000\
#   -b 2000\
#   -bs 10\
#   --sortRegions keep\
#   -p 8 




# at most gained or lost subset
computeMatrix reference-point -S $MEN1\
   -R $most_peaks \
   --outFileName ${OUTDIR}men1_cov_most_atac_mat.gz\
   --referencePoint center\
   -a 1500\
   -b 1500\
   -bs 10\
   --sortRegions descend\
   -p 8\
   --outFileSortedRegions $SORTORDER


computeMatrix reference-point -S $MLL\
   -R $most_peaks\
   --outFileName ${OUTDIR}mll_cov_most_atac_mat.gz\
   --referencePoint center\
   -a 1500\
   -b 1500\
   -bs 10\
   --sortRegions keep\
   -p 8 

computeMatrix reference-point -S $H3K4me3\
   -R $most_peaks\
   --outFileName ${OUTDIR}h3k4me3_cov_most_atac_mat.gz\
   --referencePoint center\
   -a 1500\
   -b 1500\
   -bs 10\
   --sortRegions keep\
   -p 8 

computeMatrix reference-point -S $IGG\
   -R $most_peaks\
   --outFileName ${OUTDIR}igg_cov_most_atac_mat.gz\
   --referencePoint center\
   -a 1500\
   -b 1500\
   -bs 10\
   --sortRegions keep\
   -p 8 

