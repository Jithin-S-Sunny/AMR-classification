# AMR-classification
This 
Studies used as reference for this framework:
1. Kim, Jee In, et al. "Machine learning for antimicrobial resistance prediction: current practice, limitations, and clinical perspective." Clinical microbiology reviews 35.3 (2022): e00179-21.
2. Santos, João Dourado, et al. "INSaFLU-TELEVIR: an open web-based bioinformatics suite for viral metagenomic detection and routine genomic surveillance." Genome Medicine 16.1 (2024): 61.

## Workflow
Raw reads in FASTQ format for multiple strains
fastp -i strainX_R1.fastq -I strainX_R2.fastq -o clean_R1.fastq -O clean_R2.fastq
spades.py -1 clean_R1.fastq -2 clean_R2.fastq -o spades_output_strainX
amrfinder -n genome_strainX.fasta -o strainX_amr.tsv
or
rgi main --input_sequence genome_strainX.fasta --output_file strainX_rgi --input_type contig
or
resfinder.py -i genome_strainX.fasta -o output_dir_strainX

Combine results from all strains into a single summary matrix (strain × AMR gene presence/absence or abundance)
| Strain | blaTEM-1	| tetA |	mecA	| vanA |
| ------ | -------- | ---- |  ----  | ---- |
| 0001 |
