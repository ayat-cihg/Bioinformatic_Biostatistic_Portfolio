#!/bin/bash

#Initialisation des variables
REGION_BED="/data/20_data_ayat/scenicplus_workflow/KO/pycistopic/consensus_peak_calling/consensus_regions.bed"
GENOME_FASTA="/data/20_data_ayat/scenicplus_workflow/cistarget/hg38.fa"
CHROMSIZES="/data/20_data_ayat/scenicplus_workflow/cistarget/hg38.chrom.sizes"
DATABASE_PREFIX="KO"
SCRIPT_DIR="/home/ayat/mamba_scenicplus/ctx_db/create_cisTarget_databases"
CBDIR="/data/20_data_ayat/scenicplus_workflow/cistarget//aertslab_motif_colleciton/v10nr_clust_public/singletons"
OUT_DIR="/data/20_data_ayat/scenicplus_workflow/KO/cistarget/"
mkdir -p $OUT_DIR

#cisTarget va permettre de faire le lien entre les peaks et les TF qui peuvent se fixer. Il va découpé les régions du génome humain en fonction de la position des peaks et va identifié les TF qui peuvent s'y fixer
"${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh" \
        ${GENOME_FASTA} \
        ${CHROMSIZES} \
        ${REGION_BED} \
        ${OUT_DIR}/hg38.WT.with_1kb_bg_padding.fa \
        1000 \
        yes

MOTIF_LIST="${OUT_DIR}/motifs.txt"

ls ${CBDIR} > ${MOTIF_LIST}

FASTA_FILE="/data/20_data_ayat/scenicplus_workflow/cistarget/hg38.CO_AGGR.with_1kb_bg_padding.fa"

"${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
    -f ${FASTA_FILE} \
    -M ${CBDIR} \
    -m ${MOTIF_LIST} \
    -o ${OUT_DIR}/${DATABASE_PREFIX} \
    --bgpadding 1000 \
    -t 20 \
