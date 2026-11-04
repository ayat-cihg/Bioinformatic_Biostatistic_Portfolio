#!/bin/bash

#Perform featureCounts the gene count matrix for differential expression analysis
/data/20_data_ayat/data_jeromine/subread-2.1.1-Linux-x86_64/bin/featureCounts -T 8 -p -a ../gencode.v49.nochr.gtf -o ../RNAseq_fibro/Star_results/counts.txt ../RNAseq_fibro/Star_results/*Aligned.sortedByCoord.out.bam

