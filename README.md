# AMR-surveillance
This workflow outlines a collection of tools utilised to first predict, presence and absence of genes followed by finding mutations in them and then geography 

## Workflow
### Creating gene presence/absence matrix

Raw reads in FASTQ format for multiple strains
1. fastp -i strain_R1.fastq -I strain_R2.fastq -o clean_R1.fastq -O clean_R2.fastq
2. spades.py -1 clean_R1.fastq -2 clean_R2.fastq -o spades_output_strain
3. amrfinder -n genome_strain.fasta -o strain_amr.tsv
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

1. bwa mem REF strain_clean_R1.fastq strain_clean_R2.fastq > strain_aln.sam
2. samtools view -bS strain_aln.sam | samtools sort -o strain_sorted.bam
3. samtools index strain_sorted.bam
4. bcftools mpileup -Ou -f REF strain_sorted.bam 
5. bcftools call -mv -Ov -o strain_variants.vcf
6. snpEff ann bacteria_db strain_variants.vcf > strain_annotated.vcf

Now we have: 
1. AMR gene presence/absence matrix — a table of strains × AMR genes
2. strain_annotated.vcf — SNPs (including AMR-associated ones)

| Strain	| blaTEM-1	| tetA	| mecA | vanA | blaTEM-1_mutation |
| ------ | ------ | ------ | ----- | ----- | ------ |
| 0001 |	1 |	1 |	0 |	1 | S83L |
| 0002	| 1	| 0	| 1	| 1	|	E81 |


### Demographic Information

| Strain | Location | Host | Antibiotic_usage |
| ------ | -------- | ----- | ----- |
| 0001 | UK | Human_01 | xxx |

By combining the three critical layers of information—AMR gene presence/absence, SNPs within resistance-related genes, and demographic metadata (such as host, geographic location, and sample date)—we create a robust multi-dimensional analysis framework. This combination not only supports machine learning-based classification but also allows for advanced statistical methods like Redundancy Analysis (RDA) and Principal Coordinates Analysis (PCoA), which are effective for uncovering patterns and connections within microbial populations.

For instance:
1. RDA (Redundancy Analysis) can be employed to assess how much of the variation in AMR gene profiles or SNP variations is accounted for by metadata variables like host species, geographical location, or antibiotic usage patterns. This aids in pinpointing significant environmental or epidemiological factors contributing to resistance.
2. PCoA (Principal Coordinates Analysis) applied to presence/absence matrices or SNP-generated distance matrices enables us to group strains based on genetic or resistance similarity, while also incorporating metadata to identify regional or host-specific trends.

When these techniques are combined with permutation-based significance tests (such as PERMANOVA), they can statistically confirm associations like:
1. Certain SNPs being notably more prevalent in isolates from particular locations or types of hosts.
2. Specific sets of AMR genes co-occurring in strains from hospital-acquired infections compared to environmental strains. 




