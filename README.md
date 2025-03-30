# AMR-classification
This workflow outlines a collection of tools utilised to first predict, presence and absence of genes followed by finding mutations in them and then geography 

Studies used as reference for this framework:
1. Kim, Jee In, et al. "Machine learning for antimicrobial resistance prediction: current practice, limitations, and clinical perspective." Clinical microbiology reviews 35.3 (2022): e00179-21.
2. Santos, João Dourado, et al. "INSaFLU-TELEVIR: an open web-based bioinformatics suite for viral metagenomic detection and routine genomic surveillance." Genome Medicine 16.1 (2024): 61.

## Workflow
### Creating gene presence/absence matrix

Raw reads in FASTQ format for multiple strains
fastp -i strain_R1.fastq -I strain_R2.fastq -o clean_R1.fastq -O clean_R2.fastq
spades.py -1 clean_R1.fastq -2 clean_R2.fastq -o spades_output_strain
amrfinder -n genome_strain.fasta -o strain_amr.tsv
or
rgi main --input_sequence genome_strainX.fasta --output_file strainX_rgi --input_type contig
or
resfinder.py -i genome_strainX.fasta -o output_dir_strainX

Combine results from all strains into a single summary matrix (strain × AMR gene presence/absence or abundance)
| Strain | blaTEM-1	| tetA |	mecA	| vanA |
| ------ | -------- | ---- |  ----  | ---- |
| 0001 | 1 |	1 |	0 |	0 |
| 0002 | 1 |	0 |	0 |	0 |
| ............ 

### Generating mutations for each strain

bwa mem REF strain_clean_R1.fastq strain_clean_R2.fastq > strain_aln.sam
samtools view -bS strain_aln.sam | samtools sort -o strain_sorted.bam
samtools index strain_sorted.bam
bcftools mpileup -Ou -f REF strain_sorted.bam 
bcftools call -mv -Ov -o strain_variants.vcf
snpEff ann bacteria_db strain_variants.vcf > strain_annotated.vcf

Now we have: 
1. AMR gene presence/absence matrix — a table of strains × AMR genes
2. strain_annotated.vcf — SNPs (including AMR-associated ones)

| Strain	| blaTEM-1	| tetA	| mecA | vanA | blaTEM-1_mutation |
| ------ | ------ | ------ | ----- | ----- | ------ |
| strain001 |	1 |	1 |	0 |	1 | S83L |
| strain002	| 1	| 0	| 1	| 1	|	E81 |










