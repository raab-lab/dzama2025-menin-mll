#!/bin/bash

#SBATCH -c 2
#SBATCH --mem=10G
#SBATCH -t 24:00:0
#SBATCH -J NF-RNA
#SBATCH -o ./logs/rna_pipeline.%a.out

module purge
module add nextflow

source ~/.secrets/airtable

## MODIFY THE EXPERIMENT ID, WORK, OUTPUT, AND EMAIL FOR YOUR RUN
nextflow run raab-lab/rnaseq -r main\
	--sample_sheet /proj/jraablab/users/jraab/menin-mll/full_rna_samples.csv\
   --single \
	-w /work/users/j/r/jraab/temp_files/menin-mll/rna   \
	--outdir /proj/jraablab/users/jraab/menin-mll/data/source_data/rna/\
	-with-report \
	-N jraab@med.unc.edu \
	-latest \
	-ansi-log false \
   --profile hg38


