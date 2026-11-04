#!/bin/bash

bash fastqc.sh
bash trim_fastqc.sh
bash fastqc.sh
bash generate_index.sh
bash alignement.sh
bash qualimap.sh
bash features_counts.sh
bash multiqc.sh
