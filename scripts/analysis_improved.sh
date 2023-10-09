#!/bin/bash
# This is simple RNA-Seq analysis workflow that does the following:
# - Creates a STAR index for the genome fasta file specified in the variable GENOME_FASTA
# - Aligns the raw reads specified in the fastq file referenced by the variable RAW_FASTQ
# - Computes the genome coverage of the aligned reads, producing a bedGraph file
# - Counts the overlaps over the genes specified in the file referenced by the variable GENES_BED

set -euo pipefail


# Global variables ----

RAW_FASTQ=`realpath 'bioinformatics_on_the_command_line_files/raw_yeast_rnaseq_data.fastq'`
GENOME_FASTA=`realpath 'bioinformatics_on_the_command_line_files/yeast_genome.fasta'`
GENES_BED=`realpath 'bioinformatics_on_the_command_line_files/yeast_genes.fixed.bed'`

ls "$RAW_FASTQ" "$GENOME_FASTA" "$GENES_BED" > /dev/null

RAW_FASTQ_BASENAME_PREFIX=`basename $RAW_FASTQ .fastq`
GENOME_FASTA_BASENAME=`basename $GENOME_FASTA`

THREADS=5


# Pipeline commands ----

# Create a directory for the results of the analysis

mkdir analysis
cd analysis

# Link to the fastq containing the raw sequences in '00_source_data'

mkdir 00_source_data
cd 00_source_data
ln -s "$RAW_FASTQ"
cd ..

# Create a STAR index for the genome fasta file, storing all of the results in '01_star_index'

echo 'Creating STAR index...'

mkdir 01_star_index
cd 01_star_index
ln -s "$GENOME_FASTA"
nice STAR --runThreadN "$THREADS" --runMode genomeGenerate --genomeDir . --genomeFastaFiles "$GENOME_FASTA_BASENAME" --genomeSAindexNbases 10
cd ..

# Align the raw sequences using STAR, storing all of the results in '02_aligned_reads'

echo 'Aligning raw reads...'

mkdir 02_aligned_reads
cd 02_aligned_reads
nice STAR --genomeDir ../01_star_index/ --readFilesIn ../00_source_data/"$RAW_FASTQ_BASENAME_PREFIX".fastq --runThreadN "$THREADS" --outFileNamePrefix "$RAW_FASTQ_BASENAME_PREFIX". --outSAMtype BAM SortedByCoordinate
cd ..

# Create a genome coverage file in bedGraph format, and store it in '03_coverage'

echo 'Creating genome coverage file...'

mkdir 03_coverage
cd 03_coverage
bedtools genomecov -ibam ../02_aligned_reads/"$RAW_FASTQ_BASENAME_PREFIX".Aligned.sortedByCoord.out.bam > "$RAW_FASTQ_BASENAME_PREFIX".genomecov.bg
cd ..

# Compute the gene overlap counts, and store them as a bed file in '04_gene_overlap_counts'

echo 'Computing gene overlap counts...'

mkdir 04_gene_overlap_counts
cd 04_gene_overlap_counts
ln -s ../../bioinformatics_on_the_command_line_files/yeast_genes.fixed.bed
bedtools intersect -c -a "$GENES_BED" -b ../02_aligned_reads/"$RAW_FASTQ_BASENAME_PREFIX".Aligned.sortedByCoord.out.bam | awk '$7>0' | sort -k7,7nr > "$RAW_FASTQ_BASENAME_PREFIX".gene_overlap_counts.bed
cd ..

# Create a run report

echo `realpath $0`" run completed successfully." > run_report.txt
date >> run_report.txt
echo 'Software versions:' >> run_report.txt
echo 'STAR' >> run_report.txt
STAR --version >> run_report.txt
echo 'bedtools' >> run_report.txt
bedtools --version | sed 's/^bedtools //' >> run_report.txt
echo 'sort' >> run_report.txt
sort --version | sed 's/^sort //' | head -1 >> run_report.txt
echo 'awk' >> run_report.txt
awk --version | head -1 >> run_report.txt

# If we've got here the pipeline has completed successfully

echo 'Pipeline completed successfully.'
