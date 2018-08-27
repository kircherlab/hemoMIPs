#!/bin/bash
samtools view -u -F 513 -r ${1} sample.bam | /fast/groups/ag_kircher/hemophilia/scripts/pipeline2.0/FilterBAM.py -q --qual_number 5 --qual_cutoff=15 -p > by_sample/${1}.bam
