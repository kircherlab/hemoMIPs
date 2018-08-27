#!/bin/bash
bwa mem -L 80 -M -C /fast/projects/cubit/current/static_data/precomputed/BWA/0.7.12/GRCh37/g1k_phase2/hs37d5.fa <( samtools view -F 513 by_sample/${1}.bam | awk 'BEGIN{ OFS="\n"; FS="\t" }{ print "@"$1"\t"$12"\t"$13"\t"$14,$10,"+",$11 }' ) | samtools view -u - | samtools sort - | /fast/groups/ag_kircher/hemophilia/scripts/pipeline2.0/TrimMIParms.py -d hemomips_design.txt -p | samtools reheader new_header.sam - | samtools sort -O bam -o aligned/${1}.bam -
samtools index aligned/${1}.bam
