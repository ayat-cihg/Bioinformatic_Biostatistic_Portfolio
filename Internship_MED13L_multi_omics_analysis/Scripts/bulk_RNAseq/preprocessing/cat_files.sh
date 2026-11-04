#!/bin/bash

#Concatenate the paire-end file
#This step was not performed in this project project

for files in ../RNAseq_1_hiPSC_copie/*_R1_001.fastq.gz ;
  do cat $files "${files/_R1_001.fastq.gz/_R2_001.fastq.gz}" > "${files/_R1_001.fastq.gz/_cat.fastq.gz}";
done


