mkdir analysis
cd analysis
mkdir 00_source_data
cd 00_source_data
ln -s ../../bioinformatics_on_the_command_line_files/raw_yeast_rnaseq_data.fastq 
cd ..
tree
mkdir 01_star_index
cd 01_star_index
ln -s ../../bioinformatics_on_the_command_line_files/yeast_genome.fasta
nice STAR --runThreadN 5 --runMode genomeGenerate --genomeDir . --genomeFastaFiles ./yeast_genome.fasta --genomeSAindexNbases 10
cd ..
tree
mkdir 02_aligned_reads
cd 02_aligned_reads
nice STAR --genomeDir ../01_star_index/ --readFilesIn ../00_source_data/raw_yeast_rnaseq_data.fastq --runThreadN 5 --outFileNamePrefix raw_yeast_rnaseq_data. --outSAMtype BAM SortedByCoordinate
cd ..
tree
mkdir 03_coverage
cd 03_coverage
bedtools genomecov -ibam ../02_aligned_reads/raw_yeast_rnaseq_data.Aligned.sortedByCoord.out.bam > raw_yeast_rnaseq_data.genomecov.bg
cd ..
tree
mkdir 04_gene_overlap_counts
cd 04_gene_overlap_counts
ln -s ../../bioinformatics_on_the_command_line_files/yeast_genes.fixed.bed
bedtools intersect -c -a yeast_genes.fixed.bed -b ../02_aligned_reads/raw_yeast_rnaseq_data.Aligned.sortedByCoord.out.bam | awk '$7>0' | sort -k7,7nr > raw_yeast_overlap_data.gene_overlap_counts.bed
cd ..
tree
head -10 04_gene_overlap_counts/raw_yeast_overlap_data.gene_overlap_counts.bed
cd ..


