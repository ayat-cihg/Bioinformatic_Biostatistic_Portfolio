#!/bin/bash

#Script to perfom alignment with STAR
#The script process all Trim Galor output file and perform alignment
#using the paired-end files (val_1 and val_2)

mkdir -p ../RNAseq_fibro/Star_results

for files in ../RNAseq_fibro/Trim_galor/*val_1.fq.gz
  do
  base=$(basename "$files" _R1_001_val_1.fq.gz)
  /data/20_data_ayat/data_jeromine/STAR/source/STAR --runThreadN 8 --genomeDir /data/20_data_ayat/data_jeromine/genome --readFilesIn $files "${files/R1_001_val_1.fq.gz/R2_001_val_2.fq.gz}" --readFilesCommand zcat --outFileNamePrefix "../RNAseq_fibro/Star_results/${base}_" --outSAMtype BAM SortedByCoordinate
done
