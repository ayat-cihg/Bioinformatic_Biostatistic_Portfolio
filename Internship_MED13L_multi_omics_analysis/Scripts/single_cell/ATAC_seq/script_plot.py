#Chargement des packages
import os
import mudata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#Script pour créer le pyramide plots
OUT_dir = "/data/20_data_ayat/scenicplus_workflow/scenicplus_out/"

types_cellulaires = ["Neural_progenitor_cell", "Cortical_neurons", "Pericyte_like_identity", "Oligodendrocytes", "Astrocytes", 
                    "Intermediate_neurons", "Muller_cells", "Non_neural_epithelial_cell", "Mesenchyme", "Unknown", "Photoreceptor"]

correction_labels = {"Mesenchyme_cells" : "Mesenchyme",
                     "neural progenitor cells" : "Neural_progenitor_cells"}
scplus_mdata_WT=mudata.read(os.path.join(OUT_dir,"WT","scplusmdata.h5mu"))
print(scplus_mdata_WT.obs["scATAC_counts:cell_types"].unique())

scplus_mdata_WT.obs["scATAC_counts:cell_types"]=(scplus_mdata_WT.obs["scATAC_counts:cell_types"].replace(correction_labels))
print(scplus_mdata_WT.obs["scATAC_counts:cell_types"].unique())

scplus_mdata_KO=mudata.read(os.path.join(OUT_dir,"KO", "scplusmdata.h5mu"))
print(scplus_mdata_KO.obs["scATAC_counts:cell_types"].unique())

WT_counts=scplus_mdata_WT.obs["scATAC_counts:cell_types"].value_counts()
KO_counts=scplus_mdata_KO.obs["scATAC_counts:cell_types"].value_counts()

print(WT_counts)
print(KO_counts)

cell_counts=pd.DataFrame({
    "WT":WT_counts,
    "KO":KO_counts
})

cell_counts=cell_counts.fillna(0)
print(cell_counts)

cell_counts = cell_counts.sort_values("KO")

# Préparer les données pour la pyramide
cell_types = cell_counts.index

WT_values = -cell_counts["WT"]   # WT à gauche
KO_values = cell_counts["KO"]    # KO à droite

y = range(len(cell_types))

# Création du plot
fig, ax = plt.subplots(figsize=(8, 6))

# Barres horizontales
ax.barh(y, KO_values, label="KO")
ax.barh(y, WT_values, label="WT")

# Labels des cellules
ax.set_yticks(y)
ax.set_yticklabels(cell_types)

# Axe X avec valeurs positives
max_value = max(cell_counts.max())
ax.set_xlim(-max_value*1.1, max_value*1.1)

ticks = ax.get_xticks()
ax.set_xticklabels([abs(int(t)) for t in ticks])

# Ligne centrale
ax.axvline(0)

# Titres
ax.set_xlabel("Nombre de cellules")
ax.set_title("Composition cellulaire WT vs KO")

ax.legend()

plt.tight_layout()
plt.savefig(os.path.join(OUT_dir, "histogram_types_cellulaire.pdf"), bbox_inches="tight")
