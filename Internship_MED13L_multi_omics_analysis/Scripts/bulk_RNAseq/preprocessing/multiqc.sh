#!/bin/bash

#Script to perform Multiqc

mkdir -p /data/20_data_ayat/data_jeromine/RNAseq_hiPSc/Multiqc

input_dir="/data/20_data_ayat/data_jeromine/RNAseq_hiPSc/"
output_dir="/data/20_data_ayat/data_jeromine/RNAseq_hiPSc/Multiqc"
multiqc "$input_dir"  --outdir "$output_dir" .
