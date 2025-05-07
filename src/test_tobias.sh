#!/bin/sh

#SBATCH -c 12
#SBATCH --mem=18G
#SBATCH --time=12:00:00

module load samtools
pixi shell


TEST_CTRL=$(find data/source_data/atac/bams/filtered/ -name "*DMSO*.bam")
TEST_SNDX=$(find data/source_data/atac/bams/filtered/ -name "*SNDX*.bam")

GENOME=/proj/seq/data/hg38_UCSC/Sequence/WholeGenomeFasta/genome.fa
PEAKS=data/derived_data/atac/gain_access.bed
DENY=data/external_data/hg38.Nordin.CandRblacklist_hg38.bed
OUTDIR=data/processed_data/atac/atac_correct/
MOTIFS=data/external_data/hocomov11_core_human_mono_jaspar_format.jaspar

echo $TEST_CTRL
echo $TEST_SNDX
#step 1 ATAC correct

# Merge Bam FILES
#samtools merge -o data/source_data/atac/bams/filtered/DMSO_merged.bam $TEST_CTRL
#samtools merge -o data/source_data/atac/bams/filtered/SNDX_merged.bam $TEST_SNDX


TOBIAS ATACorrect --bam data/source_data/atac/bams/filtered/DMSO_merged.bam\
          --genome $GENOME \
          --peaks $PEAKS \
          --blacklist $DENY \
          --outdir $OUTDIR \
          --cores 12

TOBIAS ATACorrect --bam data/source_data/atac/bams/filtered/SNDX_merged.bam \
             --genome $GENOME \
             --peaks $PEAKS \
             --blacklist $DENY \
             --outdir $OUTDIR \
             --cores 12


#Step 2 scorebw
TOBIAS FootprintScores --signal data/processed_data/atac/atac_correct/DMSO_merged_corrected.bw --regions $PEAKS --output ${OUTDIR}/dmso_atac_footprints.bw --cores 12
TOBIAS FootprintScores --signal data/processed_data/atac/atac_correct/SNDX_merged_corrected.bw --regions $PEAKS --output ${OUTDIR}/sndx_atac_footprints.bw --cores 12
#step 3 BINDetect

TOBIAS BINDetect --motifs $MOTIFS \
   --signals  ${OUTDIR}/dmso_atac_footprints.bw ${OUTDIR}/sndx_atac_footprints.bw \
   --genome $GENOME \
   --peaks $PEAKS \
   --outdir $OUTDIR/BINDetect_output \
   --cond_names DMSO SNDX\
   --cores 12


#step4 - visualize specific TF swith plotAggregates

MOTIF_GAIN=data/derived_data/atac/gain_access.bed
MOTIF_LOSS=data/derived_data/atac/lost_access.bed
TEAD1=data/processed_data/atac/atac_correct/BINDetect_output/_TEAD1_HUMAN.H11MO.0.A/beds/_TEAD1_HUMAN.H11MO.0.A_all.bed
TEAD4=data/processed_data/atac/atac_correct/BINDetect_output/_TEAD4_HUMAN.H11MO.0.A/beds/_TEAD4_HUMAN.H11MO.0.A_all.bed
JUN=data/processed_data/atac/atac_correct/BINDetect_output/_JUN_HUMAN.H11MO.0.A/beds/_JUN_HUMAN.H11MO.0.A_all.bed
NFY=data/processed_data/atac/atac_correct/BINDetect_output/_NFYB_HUMAN.H11MO.0.A/beds/_NFYB_HUMAN.H11MO.0.A_all.bed
KLF1=data/processed_data/atac/atac_correct/BINDetect_output/_KLF1_HUMAN.H11MO.0.A/beds/_KLF1_HUMAN.H11MO.0.A_all.bed
IRF8=data/processed_data/atac/atac_correct/BINDetect_output/_IRF8_HUMAN.H11MO.0.B/beds/_IRF8_HUMAN.H11MO.0.B_all.bed

DMSOBW=data/processed_data/atac/atac_correct/DMSO_merged_corrected.bw
SNDXBW=data/processed_data/atac/atac_correct/SNDX_merged_corrected.bw
#TOBIAS PlotAggregate --TFBS $TEAD1  --signals $DMSOBW $SNDXBW --output $OUTDIR/tead1_all_sdnx_vs_dmso.pdf  --share_y both --plot_boundaries --signal-on-x
#TOBIAS PlotAggregate --TFBS $TEAD4  --signals $DMSOBW $SNDXBW --output $OUTDIR/tead4_all_sdnx_vs_dmso.pdf  --share_y both --plot_boundaries --signal-on-x
#TOBIAS PlotAggregate --TFBS $JUN  --signals $DMSOBW $SNDXBW --output $OUTDIR/jun_all_sdnx_vs_dmso.pdf  --share_y both --plot_boundaries --signal-on-x
#TOBIAS PlotAggregate --TFBS $KLF1  --signals $DMSOBW $SNDXBW --output $OUTDIR/klf1_all_sdnx_vs_dmso.pdf  --share_y both --plot_boundaries --signal-on-x
#TOBIAS PlotAggregate --TFBS $NFY  --signals $DMSOBW $SNDXBW --output $OUTDIR/nfy_all_sdnx_vs_dmso.pdf  --share_y both --plot_boundaries --signal-on-x
#TOBIAS PlotAggregate --TFBS $IRF8  --signals $DMSOBW $SNDXBW --output $OUTDIR/irf8_all_sdnx_vs_dmso.pdf  --share_y both --plot_boundaries --signal-on-x




