#!/bin/bash

#SBATCH -c 2

#SBATCH --mem=10G

#SBATCH -t 24:00:00

#SBATCH -J NF


module purge
module add nextflow

## MODIFY THE SAMPLESHEET AND EMAIL FOR YOUR RUN
nextflow run raab-lab/cut-n-run \
      -r dev\
		--sample_sheet /proj/jraablab/users/jraab/menin-mll/full_cnr_samples.csv \
      --group_normalize \
		-profile hg38 \
		-w /work/users/j/r/jraab/temp/menin-mll/cnr/ \
		--outdir /proj/jraablab/users/jraab/menin-mll/data/source_data/cnr/\
		-with-report \
		-N jraabs@med.unc.edu \
		-ansi-log false \
		-resume
		##-r main 
		##-latest 
