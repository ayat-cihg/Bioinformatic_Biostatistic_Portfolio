#Chargement des packages

import pycisTopic
import os
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import pickle
from pycisTopic.plotting.qc_plot import plot_sample_stats, plot_barcode_stats
from pycisTopic.qc import get_barcodes_passing_qc_for_sample

#Initialisation des variables
output_dir = "/data/20_data_ayat/scenicplus_workflow/KO/pycistopic"
QC_dir = "/data/20_data_ayat/scenicplus_workflow/KO/qc"
plot_dir = f"{QC_dir}/plots"
os.makedirs(QC_dir, exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)

fragments_dict = { "KO" : "/data/90_data_jerome/donnees_sf/runs/CO_KO/count_CO_KO/outs/atac_fragments.tsv.gz"}

#Lancement du TSS (Transcript Start Site), les régions où la transcription commence, du coup forcément les régions d'ADN sont accessibles 
print("Commencement du TSS")

subprocess.run([
    "pycistopic", "tss", "get_tss",
    "--output", f"{QC_dir}/tss.bed",
    "--name", "hsapiens_gene_ensembl",
    "--to-chrom-source", "ucsc",
    "--ucsc", "hg38"
])

print("TSS fini")

#Lancement du contrôle qualité
print("Début du contrôle qualité")

subprocess.run([
    "pycistopic", "qc",
    "--fragments", fragments_dict["KO"],
    "--regions", f"{output_dir}/consensus_peak_calling/consensus_regions.bed",
    "--tss", f"{QC_dir}/tss.bed",
    "--output", f"{QC_dir}/KO"
])

print("Controle qualité fini")

#Plot
print("Construction des plots")

for sample_id in fragments_dict:
    fig = plot_sample_stats(
        sample_id = sample_id,
        pycistopic_qc_output_dir = QC_dir
    )
    fig.savefig(
        f"{plot_dir}/{sample_id}_sample_stats.png",
        dpi=300,
        bbox_inches="tight"
    )
    plt.close(fig)

sample_id_to_barcodes_passing_filters = {}
sample_id_to_thresholds = {}
for sample_id in fragments_dict:
    (
        sample_id_to_barcodes_passing_filters[sample_id],
        sample_id_to_thresholds[sample_id]
    ) = get_barcodes_passing_qc_for_sample(
            sample_id = sample_id,
            pycistopic_qc_output_dir = QC_dir,
            unique_fragments_threshold = None, # use automatic thresholding
            tss_enrichment_threshold = None, # use automatic thresholding
            frip_threshold = 0,
            use_automatic_thresholds = True,
    )

for sample_id in fragments_dict:
    fig = plot_barcode_stats(
        sample_id = sample_id,
        pycistopic_qc_output_dir = QC_dir,
        bc_passing_filters = sample_id_to_barcodes_passing_filters[sample_id],
        detailed_title = False,
        **sample_id_to_thresholds[sample_id]
    )
    fig.savefig(
        f"{plot_dir}/{sample_id}_barcode_stats.png",
        dpi=300,
        bbox_inches="tight"
    )
    plt.close(fig)

pickle.dump(
    sample_id_to_barcodes_passing_filters,
    open(f"{QC_dir}/barcodes_passing_filters.pkl", "wb")
)
