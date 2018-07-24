#!/usr/bin/env bash

# extract PE reads from BAM file following instructions from https://gist.github.com/darencard/72ddd9e6c08aaff5ff64ca512a04a6dd

module load anaconda/2

# get reads that are properly paired
samtools view -u -f 1 -F 12 --threads 36 my.bam > my_map_map.bam

# R1 unmapped, R2 mapped
samtools view -u -f 4 -F 264 --threads 36 my.bam > my_unmap_map.bam

# R1 mapped, R2 unmapped
samtools view -u -f 8 -F 260 --threads 36 my.bam > my_map_unmap.bam

# R1 & R2 unmapped
samtools view -u -f 12 -F 256 --threads 36 my.bam > my_unmap_unmap.bam

samtools merge --threads 36 -u my_unmapped.bam my_unmap_map.bam my_map_unmap.bam my_unmap_unmap.bam

samtools sort --threads 36 -n my_map_map.bam -o my_mapped.sort

samtools sort --threads 36 -n my_unmapped.bam -o my_unmapped.sort

# convert bam files to fastq
samtools bam2fq --threads 36 my_mapped.sort > my_mapped.fastq
samtools bam2fq --threads 36 my_unmapped.sort > my_unmapped.fastq

# separate fastq files into fwd and rev
cat my_mapped.fastq | grep '^@.*/1$' -A 3 --no-group-separator > my_mapped.1.fastq
cat my_mapped.fastq | grep '^@.*/2$' -A 3 --no-group-separator > my_mapped.2.fastq

cat my_unmapped.fastq | grep '^@.*/1$' -A 3 --no-group-separator > my_unmapped.1.fastq
cat my_unmapped.fastq | grep '^@.*/2$' -A 3 --no-group-separator > my_unmapped.2.fastq

cat my_mapped.1.fastq my_unmapped.1.fastq > my_from_BAM.1.fastq
cat my_mapped.2.fastq my_unmapped.2.fastq > my_from_BAM.2.fastq
