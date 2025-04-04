#!/bin/bash

#SBATCH -c 2

#SBATCH --mem=10G

#SBATCH -t 24:00:00

#SBATCH -J NF

module add nextflow
module load anaconda
conda activate mageck 

## MODIFY ALL OF THE PATHS TO YOUR SPECIFIC PROJECT
nextflow run raab-lab/crispr-screens \
      -r main \
		--sample_sheet /proj/jraablab/users/jraab/menin-mll/data/all_crispr_data_sets.csv \
		--grna /proj/jraablab/users/jraab/menin-mll/data/source_data/crispr/epi_lib.tsv \
		--control_guides /proj/jraablab/users/jraab/menin-mll/data/source_data/crispr/controls_sgRNAs.txt \
		--contrasts /proj/jraablab/users/jraab/menin-mll/data/crispr_contrasts.tsv \
		--outdir /proj/jraablab/users/jraab/menin-mll/data/source_data/crispr/ \
		-w /work/users/j/r/jraab/menin-mll/crispr \
		-with-report \
		-N jraab@med.unc.edu \
		-latest \
		-ansi-log false 
