#Chargement des packages
import os
import mudata
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
import plotnine as p9
from scenicplus.RSS import (regulon_specificity_scores, plot_rss)
from scenicplus.plotting.dotplot import heatmap_dotplot

#Initialisation des variables
OUT_dir = "/data/20_data_ayat/scenicplus_workflow/scenicplus_out/WT/"
Fig_dir = "/data/20_data_ayat/scenicplus_workflow/scenicplus_out/WT/Figure"

os.makedirs(Fig_dir, exist_ok=True)

scplus_mdata=mudata.read("/data/20_data_ayat/scenicplus_workflow/scenicplus_out/WT/scplusmdata.h5mu")
print(scplus_mdata.uns["direct_e_regulon_metadata"])
print(scplus_mdata.uns["extended_e_regulon_metadata"])
print(scplus_mdata.obs.head())
print(scplus_mdata.obs.columns)

#Création des plots
eRegulon_gene_AUC = anndata.concat(
    [scplus_mdata["direct_gene_based_AUC"], scplus_mdata["extended_gene_based_AUC"]],
    axis = 1,
)

eRegulon_gene_AUC.obs = scplus_mdata.obs.loc[eRegulon_gene_AUC.obs_names]
sc.pp.neighbors(eRegulon_gene_AUC, use_rep = "X")
sc.tl.umap(eRegulon_gene_AUC)
sc.settings.figdir = Fig_dir

#direct_regulons = set(
#    scplus_mdata.uns["direct_e_regulon_metadata"]["eRegulon_name"]
#)

genes = ['SOX2', 'FOXG1', 'TBR2', 'TBR1', 'BCL11B', 'SATB2', 'CUX1', 'MBP','PRDM1',
    'RAX', 'VSX2', 'OTX2', 'CRX','NEUROD1', 'RCVRN', 'RORB', 'THRB', 'USH2A','NRL', 'NR2E3',
    'NR2F1', 'PAX6', 'FOXP2','GAD1','MKI67','FABP7', 'NES', 'MED13L']

#UMAP de l'activité des eRegulons
for i in genes :
    match = [r for r in eRegulon_gene_AUC.var_names if i in r and "_direct_+/-" in r]
    for valide in match:
        sc.pl.umap(
        eRegulon_gene_AUC,
        color=valide,
        save=f"_{i}_eRegulon_gene.pdf"
        )

#Plot avec les 10 eRégulons les plus élevée en fonction des types cellulaires
rss = regulon_specificity_scores(
    scplus_mudata = scplus_mdata,
    variable = "scATAC_counts:cell_types",
    modalities = ["direct_gene_based_AUC"] #, "extended_gene_based_AUC"]
)

plot_rss(
    data_matrix = rss,
    top_n = 10,
    num_columns = 5
)
plt.savefig(os.path.join(OUT_dir, "rss_plot.pdf")) #bbox_inches="tight")
plt.close()

top10_rss = {}

for celltype in rss.index:
    top10_rss[celltype] = (
        rss.loc[celltype]
        .sort_values(ascending=False)
        .head(10)
    )

top10_rss

# Garder uniquement les eRégulons +/+ et -/-
top10_rss_filtered = {}

for celltype, regulons in top10_rss.items():
    filtered = regulons[
        regulons.index.str.contains(r"\+/\+|-/-", regex=True)
    ]
    top10_rss_filtered[celltype] = filtered

top10_rss_filtered

selected_eregulons = list(set(
    er
    for x in top10_rss_filtered
    for er in top10_rss_filtered[x].index
))

import re

selected_direct_clean = [
    re.sub(r"_\(\d+g\)$", "", x)
    for x in selected_eregulons
]
print(selected_direct_clean[:5])
print(scplus_mdata.uns["direct_e_regulon_metadata"]["eRegulon_name"].unique()[:5])

#Heatmap avec les 5 meilleurs eRégulons qui ont un rôle activateur ou inhibiteur
p = heatmap_dotplot(
    scplus_mudata = scplus_mdata,
    color_modality = "direct_gene_based_AUC",
    size_modality = "direct_region_based_AUC",
    group_variable = "scATAC_counts:cell_types",
    eRegulon_metadata_key = "direct_e_regulon_metadata",
    color_feature_key = "Gene_signature_name",
    size_feature_key = "Region_signature_name",
    feature_name_key = "eRegulon_name",
    sort_data_by = "direct_gene_based_AUC",
    subset_feature_names = selected_direct_clean,
    orientation = "horizontal",
    figsize = (25,13),
)

print(type(p))
p = p + p9.theme(
    axis_text_x = p9.element_text(size=15, ha="right"),
    axis_text_y = p9.element_text(size=15),
  #  axis_title_x = p9.element_text(size=16),
  #  axis_title_y = p9.element_text(size=16)
)

p.save(
    os.path.join(Fig_dir, "heatmap_top_RSS.pdf"),
    width=25,
    height=13
)
sc.pl.umap(eRegulon_gene_AUC, color = list(set([x for xs in [rss.loc[ct].sort_values()[0:5].index for ct in rss.index] for x in xs ])), save="_enrichissement_score.pdf")
plt.close()
import plotnine as p9

#p9.theme_set(p9.theme_classic())
#p9.options.figure_size = (35, 15)
#p9.options.limitsize = False

#Heatmap globale
heatmap_dotplot(
    scplus_mudata = scplus_mdata,
    color_modality = "direct_gene_based_AUC",
    size_modality = "direct_region_based_AUC",
    group_variable = "scATAC_counts:cell_types",
    eRegulon_metadata_key = "direct_e_regulon_metadata",
    color_feature_key = "Gene_signature_name",
    size_feature_key = "Region_signature_name",
    feature_name_key = "eRegulon_name",
    sort_data_by = "direct_gene_based_AUC",
    orientation = "horizontal",
    figsize = (25, 15),
    save=os.path.join(Fig_dir, "heatmap.pdf")
)

