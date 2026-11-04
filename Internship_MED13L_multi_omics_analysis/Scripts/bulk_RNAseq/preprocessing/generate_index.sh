#!/bin/bash

#Important step befor performing alignment
#Index genome must be generated using STAR
/data/20_data_ayat/data_jeromine/STAR/source/STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /data/20_data_ayat/data_jeromine/genome --genomeFastaFiles /data/10_indices/star/sources/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fasta --sjdbGTFfile /data/20_data_ayat/data_jeromine/gencode.v49.nochr.gtf
