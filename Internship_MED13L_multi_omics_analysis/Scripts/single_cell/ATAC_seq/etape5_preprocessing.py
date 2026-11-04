#Chargement des packages
import os
import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from pycisTopic.topic_binarization import binarize_topics
from scenicplus.utils import region_names_to_coordinates
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features)

#Initialisation des variables
DATA_dir="/data/20_data_ayat/scenicplus_workflow/WT/pycistopic/"
Plot_dir="/data/20_data_ayat/scenicplus_workflow/WT/qc/plots/"
TEMP_dir="/data/20_data_ayat/tmp"
cell_dir="/data/20_data_ayat/CO_WT/cell_types.csv"

os.makedirs(os.path.join(DATA_dir, "region_sets"), exist_ok=True)

cistopic_obj=pickle.load(
    open(os.path.join(DATA_dir,
        "cistopic_obj_with_model.pkl"), "rb"
        )
)

print(cistopic_obj)
print(cistopic_obj.cell_data.columns)

annot = pd.read_csv("/data/20_data_ayat/CO_WT/cell_types.csv")

annot["Barcode"] = annot["Barcode"].astype(str) + "-WT"
annot = annot.set_index("Barcode")

print(cistopic_obj.cell_data.index[:10])

print(annot.index[:10])

print(len(set(cistopic_obj.cell_data.index) & set(annot.index)))

print(cistopic_obj.cell_data.index[0])
print(annot.index[0])

cistopic_obj.cell_data["cell_types"] = cistopic_obj.cell_data.index.map(
    annot["cell_types"]
).fillna("Unknown")

#Nettoyages des cellules qui sont inconnues

with open(os.path.join(DATA_dir, "cistopic_obj_final.pkl"), "wb") as f:
    pickle.dump(cistopic_obj, f)

print(cistopic_obj.cell_data["cell_types"].value_counts())

cistopic_obj.cell_data = cistopic_obj.cell_data[
    cistopic_obj.cell_data["cell_types"] != "Unknown"
]

region_bin_topics_otsu = binarize_topics(
    cistopic_obj,
    method="otsu",
    num_columns=5,
    plot=True)

for i in plt.get_fignums():
    fig = plt.figure(i)
    fig.savefig(os.path.join(Plot_dir, f"otsu_plot_{i}.pdf"), dpi=300, bbox_inches="tight")
plt.close("all")

region_bin_topics_top3k= binarize_topics(
    cistopic_obj,
    method='ntop',
    ntop=3000,
    num_columns=5,
    plot=True)

for i in plt.get_fignums():
    fig = plt.figure(i)
    fig.savefig(os.path.join(Plot_dir, f"ntop_plot_{i}.pdf"), dpi=300, bbox_inches="tight")
plt.close("all")

os.makedirs(os.path.join(DATA_dir, "region_sets", "topics_otsu"), exist_ok=True)

#Méthode qui permet de prendre les topics les plus significatives selon un seuil, plus de pertinance biologiques
for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index).sort_values(
        ["Chromosome", "Start", "End"]).to_csv(
        os.path.join(DATA_dir, "region_sets", "topics_otsu", f"{topic}.bed"),
        sep="\t", header=False, index=False)

os.makedirs(os.path.join(DATA_dir, "region_sets", "topics_top_3k"), exist_ok=True)

#Méthode qui permet de prendre les 3000 topics les plus élevées, c'est plus comparable
for topic in region_bin_topics_top3k:
    region_names_to_coordinates(
        region_bin_topics_top3k[topic].index).sort_values(
        ["Chromosome", "Start", "End"]).to_csv(
        os.path.join(DATA_dir, "region_sets", "topics_top_3k", f"{topic}.bed"),
        sep="\t", header=False, index=False)

#Rend la matrice moins creuse
imputed_acc_obj = impute_accessibility(
    cistopic_obj,
    selected_cells=cistopic_obj.cell_data.index,
    selected_regions=None,
    scale_factor=10**6
)

mat = imputed_acc_obj.data if hasattr(imputed_acc_obj, "data") else None
print(mat is None)

normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)

#Ici on garde les peak qui varient beaucoup entre cellules, par exemple ouvert dans les astrocytes et fermé dans les photorécepteurs
variable_regions = find_highly_variable_features(
    normalized_imputed_acc_obj,
    min_disp=0.05,
    min_mean=0.0125,
    max_mean=3,
    max_disp = np.inf,
    n_bins=20,
    n_top_features=None,
    plot=True
)

cistopic_obj.cell_data['cell_types'].unique()

#Réalisation des DARs (Differentially Accessible Regions) soit les régions du génomes, dont l'accessibilité varie entre les types cellulaires
markers_dict=find_diff_features(
    cistopic_obj,
    imputed_acc_obj,
    variable='cell_types',
    var_features=variable_regions,
    contrast=None,
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=1,
    _temp_dir=TEMP_dir,
    split_pattern='-'
)
print(len(variable_regions))
print(variable_regions[:20])

os.makedirs(os.path.join(DATA_dir, "region_sets", "DARs_cell_type"), exist_ok=True)
#for cell_type in markers_dict:
for cell_type, df in markers_dict.items():

    if df.empty:
        print(f"{cell_type}: aucun DAR, ignoré")
        continue

    region_names_to_coordinates(df.index).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(
            DATA_dir,
            "region_sets",
            "DARs_cell_type",
            f"{cell_type.replace('/', '')}.bed"
        ),
        sep="\t",
        header=False,
        index=False
    )

print(cistopic_obj.cell_data["cell_types"].value_counts())
print(len(variable_regions))
print(len(cistopic_obj.cell_data["cell_types"].unique()))

print(len(markers_dict))
print(list(markers_dict.keys())[:10])
