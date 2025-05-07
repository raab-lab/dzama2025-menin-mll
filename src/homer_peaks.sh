#!/bin/bash 

# set up script to run HOMER on a range of samples

#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --mail-user=jraab@med.unc.edu
#SBATCH --job-name homer_jraab

module load homer

ATAC_GAIN=data/derived_data/atac/gain_access.bed
ATAC_LOST=data/derived_data/atac/lost_access.bed
UNCHANGED=data/derived_data/atac/unchanged_access.bed
STABLE_ATAC_MENMLL_LOST=data/derived_data/peak_sets/atac/stable_atac_lose_menin_and_mll.bed
NFYB=data/derived_data/peak_sets/cnr_consensus/nfyb_union.bed
TEAD=data/derived_data/peak_sets/cnr_consensus/tead4_union.bed
NFYB_NO_ATAC='data/derived_data/peak_sets/cnr/nfyb_atac_overlap.bed'
NFYB_ATAC='data/derived_data/peak_sets/cnr/nfyb_atac_overlap.bed'
MENDOWN_NFYBSTABLE=data/processed_data/cnr/menin_down_nfyb_stable.bed
NFYB_UP_NOMEN=data/processed_data/cnr/nfyb_up_no_menin.bed


#findMotifsGenome.pl $ATAC_GAIN hg38 data/derived_data/homer/sndx_up_vs_unchanged -size 200 -p 4 -bg $UNCHANGED
#findMotifsGenome.pl $ATAC_LOST hg38 data/derived_data/homer/sndx_down_vs_unchanged -size 200 -p 4 -bg $UNCHANGED
#findMotifsGenome.pl $STABLE_ATAC_MENMLL_LOST hg38 data/derived_data/homer/stable_atac_lose_menin -size 200 -p 4 
#findMotifsGenome.pl $NFYB hg38 data/derived_data/homer/nfyb -size 200 -p 4 
#findMotifsGenome.pl $TEAD hg38 data/derived_data/homer/tead4 -size 200 -p 4 

#findMotifsGenome.pl $NFYB_ATAC hg38 data/derived_data/homer/nfyb_atac_overlap -size 200 -p 4 
#findMotifsGenome.pl $NFYB_NO_ATAC hg38 data/derived_data/homer/nfyb_atac_nooverlap -size 200 -p 4 
$findMotifsGenome.pl $MENDOWN_NFYBSTABLE hg38 data/derived_data/homer/men_down_nfyb_stable -size 200 -p 4 
findMotifsGenome.pl $NFYB_UP_NOMEN hg38 data/derived_data/homer/nfyb_up_nomen -size 200 -p 4 



