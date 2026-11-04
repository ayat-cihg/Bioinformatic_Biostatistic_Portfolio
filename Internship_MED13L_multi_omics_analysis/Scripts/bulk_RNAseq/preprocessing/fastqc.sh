#!/bin/bash

#Perform fastqc for all file in RNAseq2-Fibro

for files in ../RNAseq_fibro/RNAseq2-Fibro-*;
  do ../FastQC/fastqc $files --outputdir ../RNAseq_fibro/Trim_galor
done
#Effectue le controle qualite et mets le document de sortie dans RNAseq_fibro/Trim_galor

