# dzama2025-menin-mll
Code associated with Dzama 2025 - Menin MLL paper

Pipelines for processing data were run using: 
	run_atac.sh
	run_cnr.sh
	run_rna.sh
	run_crispr.sh
sampl sheets for each are located in toplevel directory - cnr_samples.csv, rna_samples.csv, atac_smaples.csv. These require the raw fastq files to be donwloaded from GEO and paths to be updated if desired.

Most downstream analysis was performed using bigWig files normalized within antibody or assay type. Individual bigWig files can be found on GEO under accession numbers

For visualization pruposes, the most useful files are the group*.bw files on GEO which are normalized by group, and averaged across replicates. 

Code for downstream analysis is in src/