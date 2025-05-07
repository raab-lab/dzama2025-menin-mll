#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --mail-user=jraab@med.unc.edu
#SBATCH --job-name deeptools-jraab

module purge
module load deeptools/3.5.1

NFYB_peaks='data/derived_data/integrated/nfyb_dor_up_atac_up.bed'
SORTORDER='data/processed_data/cnr/nfyb_sort_order.bed'
INPUTDIR='data/source_data/cnr/bw'
OUTDIR='data/processed_data/integrated/'
mkdir -p $OUTDIR

MEN1=$(find $INPUTDIR -name "groupHLF_*_Menin*mean.bw" | sort)
MLL=$(find $INPUTDIR -name "groupHLF_*_MLL*mean.bw" | sort)
H3K4me3=$(find $INPUTDIR -name "groupHLF_*_H3K4me3*mean.bw"| sort)
IGG=$(find $INPUTDIR -name "groupHLF_*_IgG*mean.bw" | sort)
NFYB=$(find $INPUTDIR -name "groupHLF_*NFYB*mean.bw" |sort)
TEAD=$(find $INPUTDIR -name "groupHLF_*TEAD*mean.bw" |sort) 
ATAC=$(find data/source_data/atac/bw -name "groupATAC*mean.bw" |sort)
K27AC=$(find $INPUTDIR -name "groupHLF_*H3K27ac*mean.bw" |sort)
K27ME3=$(find $INPUTDIR -name "groupHLF_*H3K27me3*mean.bw"|sort)

echo $MEN1
echo $MLL
echo $H3K4me3
echo $IGG
echo $NFYB
echo $TEAD

computeMatrix reference-point -S $NFYB\
   -R $NFYB_peaks \
   --outFileName ${OUTDIR}nfyb_cov_nfyb_peaks_mat.gz\
   --referencePoint center\
   -a 2000\
   -b 2000\
   -bs 10\
   --sortRegions keep\
   -p 8\
   --outFileSortedRegions $SORTORDER

computeMatrix reference-point -S $MLL\
   -R $SORTORDER\
   --outFileName ${OUTDIR}mll_cov_nfyb_peaks_mat.gz\
   --referencePoint center\
   -a 2000\
   -b 2000\
   -bs 10\
   --sortRegions keep\
   -p 8 

computeMatrix reference-point -S $H3K4me3\
   -R $SORTORDER\
   --outFileName ${OUTDIR}h3k4me3_cov_nfyb_peaks_mat.gz\
   --referencePoint center\
   -a 2000\
   -b 2000\
   -bs 10\
   --sortRegions keep\
   -p 8 

computeMatrix reference-point -S $IGG\
   -R $SORTORDER\
   --outFileName ${OUTDIR}igg_cov_nfyb_peaks_mat.gz\
   --referencePoint center\
   -a 2000\
   -b 2000\
   -bs 10\
   --sortRegions keep\
   -p 8 

computeMatrix reference-point -S $MEN1\
   -R $SORTORDER \
   --outFileName ${OUTDIR}men1_cov_nfyb_peaks_mat.gz\
   --referencePoint center\
   -a 2000\
   -b 2000\
   -bs 10\
   --sortRegions keep\
   -p 8


computeMatrix reference-point -S $NFYB $MEN1 $MLL $H3K4me3 $K27AC $K27ME3 $ATAC $IGG \
   -R $NFYB_peaks \
   --outFileName ${OUTDIR}all_cov_nfyb_peaks_mat.gz\
   --referencePoint center\
   -a 2000\
   -b 2000\
   -bs 10\
   --sortRegions keep\
   -p 8


