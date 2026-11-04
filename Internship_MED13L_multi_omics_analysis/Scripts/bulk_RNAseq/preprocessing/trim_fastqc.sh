#!/bin/bash

# Script to perform trimming, and fastqc on paired-end RNA-seq data
# The script process all *_val_1.fq.gz file and automatically retrives
# The corresponding *_val_2.fq.gz mate for each sample

output_dir="../RNAseq_fibro/Trim_galor"
for files in ../RNAseq_fibro/Trim_galor/RNAseq2-Fibro-*_val_1.fq.gz ;
  do
  ../TrimGalore-0.6.11/trim_galore --paired --gzip --cores 6 --quality 30 --fastqc --output_dir $output_dir $files "${files/_R1_001_val_1.fq.gz/_R2_001_val_2.fq.gz}"
done

