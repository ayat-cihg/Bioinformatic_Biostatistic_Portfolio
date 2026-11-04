import scanpy as sc
import os

OUT_DIR = "/data/20_data_ayat/scenicplus_workflow/scenicplus_out/WT"
os.makedirs(OUT_DIR, exist_ok=True)
DATA_DIR = "/data/90_data_jerome/donnees_sf/runs/CO_WT/count_CO_WT/outs/filtered_feature_bc_matrix.h5"

adata = sc.read_10x_h5(DATA_DIR)

#Garder le adata brute
adata.raw = adata.copy()

#Fixer le nom des gènes
adata.var_names = adata.var["gene_ids"]  # ou gene_symbols si dispo

#Retirer les doublons
adata.var_names_make_unique()
adata = adata[:, ~adata.var_names.duplicated()].copy()

#QC
print("dup genes:", adata.var_names.duplicated().sum())

#Normalisation
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

#Ajouter le suffixe si l'ATAC-seq en a
adata.obs_names = adata.obs_names.astype(str) + "-WT"

assert adata.var_names.is_unique
assert adata.obs_names.is_unique

adata.write(os.path.join(OUT_DIR, "adata_clean.h5ad"))
