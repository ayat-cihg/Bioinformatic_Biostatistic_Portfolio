#!/bin/bash

#Alignment quality controle using Qualimap

mkdir -p ../RNAseq_fibro/qualimap

for file in ../RNAseq_fibro/Star_results/*.bam
  do
  base=$(basename "$file" _Aligned.sortedByCoord.out.bam)
  /data/20_data_ayat/data_jeromine/qualimap_v2.3/qualimap rnaseq -pe --bam $file -gtf ../gencode.v49.nochr.gtf --outdir ../RNAseq_fibro/qualimap/$base --java-mem-size=32G
done
