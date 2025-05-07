#!/bin/bash

#SBATCH -c 2
#SBATCH --mem=10G
#SBATCH -t 24:00:00
#SBATCH -J NF
#SBATCH -o ./logs/output.%a.out

module purge
module add nextflow

# can't get relative paths to work here
SS=/proj/jraablab/users/jraab/menin-mll/atac_samples_norm.csv 

## MODIFY THE SAMPLESHEET AND EMAIL FOR YOUR RUN
nextflow run raab-lab/cut-n-run \
      -r dev\
		--sample_sheet $SS \
		--group_normalize \
		-w /work/users/j/r/jraab/temp_files/2024-01-31/atac \
		--outdir /proj/jraablab/users/share/project_source_data/data/source_data/atac/ \
		-with-report logs/ \
		-N jraab@med.unc.edu \
      --norm_method peaks\
      -profile hg38\
		-latest \
      -resume\
		-ansi-log false 	
