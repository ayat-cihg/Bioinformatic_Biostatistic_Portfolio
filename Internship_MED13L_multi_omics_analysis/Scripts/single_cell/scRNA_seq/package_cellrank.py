import cellrank as cr
import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import sys
import pathlib
import argparse
import joblib
import matplotlib.pyplot as plt
import math

plt.rcParams["figure.facecolor"]="white"

def fusion_velocyto(adata, loom, data_dir, annotation):
    """
    Fonction qui permet la fusion des données épissées et non épissées dans nos adata, mais aussi l'ajout
    des types cellulaires pour chaque cellule.

    :Paramètres :
    adata : fichier adata des single cell RNA seq
    loom : fichier qui contient les données épissées et non épissées
    data_dir : chemin de sortie
    annotation = fichier csv avec les types cellulaires pour chaque cellule

    :return : le fichier adata complété avec les données d'épissages et d'annotation
    """
    #Distinguer si c'est un fichie KO ou WT
    if "WT" in loom:
      name = "^count_CO_WT:"
      save_name="adata_WT_with_velocyto.h5ad"
      type="cell_types"
    else:
      name = "^count_CO_KO:"
      save_name="adata_KO_with_velocyto.h5ad"
      type="cell_types"

    #Lecture des fichiers loom
    adata_loom=sc.read(loom)
    adata = adata

    #Les noms des cellules ont des suffixes et prefixe, il faut les enlever pour rajouter les données d'épissages dans le adata
    adata_loom.obs_names = adata_loom.obs_names.str.replace(name, "", regex=True)
    adata_loom.obs_names = adata_loom.obs_names.str.replace("x$", "", regex=True)
    adata.obs_names = adata.obs_names.str.replace("-1$", "", regex=True)
    adata.var_names_make_unique()

    #On trie par ordre les noms des cellules
    adata_loom=adata_loom[adata.obs_names,:]

    #On rajoute  les données d'épissages dans le adata
    adata.layers['spliced']=adata_loom.layers['spliced'].copy()
    adata.layers['unspliced']=adata_loom.layers['unspliced'].copy()
    adata.layers['ambiguous']=adata_loom.layers['ambiguous'].copy()

    #Lecture du fichier csv qui contient les types cellulaires pour chaque cellule
    csv_file=pd.read_csv(annotation)
    csv_file = csv_file.set_index("Barcode")
    csv_file.index =csv_file.index.str.replace("-1$", "", regex=True)
    adata.obs["cell_types"]=csv_file[type].reindex(adata.obs_names)
    for ad in [adata]:
      ad.obs["cell_types"] = ad.obs["cell_types"].fillna("Unknown").astype("category")

    print("la fusion des données épissé est un succès")
    print(adata)
    adata.write(data_dir /save_name)
    return adata

def aggregation(adata_1, loom_1, adata_2, loom_2, data_dir, data_name, annotation_1, annotation_2):
    """
    Fusion des données KO et WT en un seul et unique fichier
    
    :Paramètres :
    adata_1, adata_2 : fichiers adata
    loom_1, loom_2 : fichier loom
    data_dir : chemin de sortie
    annotation_1, annotation_2 : fichiers d'annotation

    :return : Fichier adata aggrégé KO et WT
    """
    #Attributions des vairables KO et WT
    if "KO" in str(adata_1):
      KO, WT =adata_1, adata_2
      loom_KO,loom_WT=loom_1, loom_2
    else :
      WT, KO =adata_1, adata_2
      loom_WT, loom_KO=loom_1, loom_2
    #Verifie si le fichier de velocité a déjà été calculé, si c'est le cas, nous n'avons pas besoin de repassé par la fonction fusion, sinon il faut y repassé
    if pathlib.Path("/data/20_data_ayat/CO_KO/adata_KO_with_velocyto_2.h5ad").exists() and pathlib.Path("/data/20_data_ayat/CO_WT/adata_WT_with_velocyto_2.h5ad").exists():
      adata_KO=sc.read("/data/20_data_ayat/CO_KO/adata_KO_with_velocyto.h5ad")
      adata_WT=sc.read("/data/20_data_ayat/CO_WT/adata_WT_with_velocyto.h5ad")
      print(adata_KO.obs_names)
    else :
      read_KO=sc.read_10x_h5(KO)
      read_WT=sc.read_10x_h5(WT)
      if "ko" in annotation_1 :
        adata_KO=fusion_velocyto(read_KO, loom_KO, data_dir, annotation_1)
        adata_WT=fusion_velocyto(read_WT, loom_WT, data_dir, annotation_2)
      else:

        adata_KO=fusion_velocyto(read_KO, loom_KO, data_dir, annotation_2)
        adata_WT=fusion_velocyto(read_WT, loom_WT, data_dir, annotation_1)
    adata_KO.var_names_make_unique()
    adata_WT.var_names_make_unique()

    #Mise en place d'une colonne condition, pour savoir qui sont les KO des WT
    adata_KO.obs['Condition']="KO"
    adata_WT.obs['Condition']="WT"

    adata_KO.obs_names = adata_KO.obs_names + "_KO"
    adata_WT.obs_names = adata_WT.obs_names + "_WT"

    #Concaténation de nos données, et sauvegarde
    adata_CO_AGGR = sc.concat([adata_KO, adata_WT])
    adata_CO_AGGR.write(data_dir / "adata_CO_AGRR_with_velocyto.h5ad")
    return adata_CO_AGGR

def chose_parameters(adata, percentage):
    """
    Fonction qui permet de trier selon un seuil de pourcentage

    :Paramètres :
    adata : fichier adata
    percentage : pourcentage de filtration

    :return : la valeur minimal et maximal de gène à garder
    """
    #Filtration où le nombre de gène est égale à 0
    adata.obs['n_gene']=(adata.X > 0).sum(axis=1)
    #Filtration qui permet de déterminer une fenêtre de gène  valides à garder
    min_gene= int(adata.obs['n_gene'].quantile(percentage/2))
    max_gene= int(adata.obs['n_gene'].quantile(1-(percentage/2)))
    return min_gene, max_gene

def filtration(adata, min_gene, max_gene, data_dir, data_name):
    """
    Filtration, normalisation des données

    :Paramètres
    adata : fichier adata
    min_gene: valeur minimal du nombre de gène à garder
    max_gene: valeur maximal du nombre de gène à garder
    data_dir: chemin où sauvegarder les données
    data_name: nom du dossier, fichier à sauvegarder

    :return fichier adata filtrer est prêt pour cellrank
    """
    #Filtration de nos données
    adata.var_names_make_unique()
    sc.pp.filter_cells(adata, min_genes=min_gene)
    sc.pp.filter_cells(adata, max_genes=max_gene)
    #Normalisation des données
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    #PCA
    sc.tl.pca(adata)
    #KNN
    sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30, random_state=0)
    #Creation de cluster en fonction de la proximité des cellules, s'il n'y a pas de types cellulaire
    sc.tl.leiden(adata, resolution=0.8)
    #Calcule des gènes les plus exprimés en fonction des types cellulaires
    sc.tl.rank_genes_groups(adata, groupby="cell_types", method="wilcoxon")
    #Plot de ces gènes
    sc.pl.rank_genes_groups(adata)
    df = sc.get.rank_genes_groups_df(adata, group=None)
    gene=df.groupby("group").head(20)['names'].tolist()
    #Dotplot des 20 gènes les plus exprimés pour chaque types cellulaires
    sc.pl.dotplot(adata, var_names=gene, groupby="cell_types", save= "plot_rank_gene.pdf")
    
    scv.pp.moments(adata, n_pcs=None, n_neighbors=None)
    sc.tl.umap(adata)
    tab = pd.crosstab(adata.obs["leiden"],adata.obs["cell_types"])
    cluster_to_cell_type=tab.idxmax(axis=1)
    adata.obs["cell_type_leiden"] = adata.obs["leiden"].map(cluster_to_cell_type)
    if "KO" in data_name:
      name="KO"
    elif "WT" in data_name:
      name="WT"
    else :
      name="CO_WT_CO_KO"
    adata.write(data_dir /f"adata_filtred_{name}.h5ad")
    df.to_csv(data_dir/ "top_gene.csv")
    print(f"df file save in  {data_dir}/top_gene.csv")
    print(f"Filtration finish, file are save in {data_dir}/adata_filtered_{name}.h5ad")
    return adata

def recover_chunk(sub_adata):
    """
    Fonction qui va permettre de lancer la parallélisation, si le mode dynamique est utilisé, car très couteux en temps
    
    :Paramètre
    sub_adata : les différents sous_data obtenus via la fonction parallele_recover_dynamics

    :return : les sous-data avec les calculs de la velocité.
    """
    scv.tl.recover_dynamics(sub_adata, n_jobs=1, show_progress_bar=False)
    scv.tl.velocity(sub_adata, mode="dynamical", show_progress_bar=False)
    return sub_adata

def parallel_recover_dynamics(adata, n_splits=4):
    """
    Fonction qui permet de diviser notre adata en 4 sous-data pour réduire le temps de calcule de la velocité avec le mode dynamique

    :Paramètres
    adata
    n_splits : Division de l'adata en ce nombre

    :return : le adata complet et fusionné
    """
    indices = np.array_split(np.arange(adata.n_obs), n_splits)
    adatas = [adata[idx].copy() for idx in indices]

    results = joblib.Parallel(n_jobs=1)(
        joblib.delayed(recover_chunk)(sub_adata)
        for i, sub_adata in enumerate(adatas)
    )
    
    adata_full = results[0].concatenate(results[1:], batch_key="batch").copy()
    print(adata_full)

    for j in adata_full.layers:
        if j not in adata.layers:
            adata.layers[j]= adata_full.layers[j]
    
    for i in adata_full.var.columns:
        if "-" in i :
            base = i.split("-")[0]
            colonnes = [c for c in adata_full.var.columns if c.startswith(base+"-")]
            adata_full.var[base] = adata_full.var[colonnes].mean(axis=1)
            adata_full.var.drop(columns=colonnes, inplace=True)

    for e in adata_full.var: 
        if e not in adata.var:
            adata.var[e]= adata_full.var[e]

    return adata

def pseudotime(adata, kernel, figure_dir, data_dir):
    """
    Fonction qui permet de calculer les pseudotimes

    :Paramètres
    adata
    kernel : Le mode utilisé, la méthode utilisé, ici c'est connectivité
    figure_dire : chemin où les figures seront stockés
    data_dire : chemin où les données sont stockées

    return : Les figures de gene_trends, le adata avec les calculs du pseudotimes
    """
    plt.style.use('_mpl-gallery')
    print("en cours pour l'estimation'")
    estimator = cr.estimators.GPCCA(kernel)
    print("compute macrostate")
    estimator.compute_macrostates(n_states=5)
    print("prediction des états terminal")
    estimator.predict_terminal_states()
    print("calcule de la probabilité des destins cellulaires")
    estimator.compute_fate_probabilities(n_jobs=2)
    adata.write(data_dir/ "adata_with_five_terminal_states.h5ad")
    scv.tl.velocity_pseudotime(adata)
    scv.pl.umap(adata, color="velocity_pseudotime")
    #scv.pl.scatter(adata, color="macrostates_fwd", save=figure_dir /"macrostate.pdf")
    scv.pl.scatter(
        adata,
        color="macrostates_fwd",
        show=False
    )

    plt.savefig(
        figure_dir / "macrostate.pdf",
        dpi=300,
        bbox_inches="tight",
        facecolor="white"
    )

    plt.close()
    model=cr.models.GAM(adata)
    liste_gene = ['PRDM1','RAX', 'VSX2', 'OTX2', 'CRX','NEUROD1', 'RCVRN', 'RORB', 'THRB', 'USH2A',
                  'NRL', 'NR2E3', 'SOX2', 'FOXG1', 'TBR2', 'TBR1', 'BCL11B', 'SATB2', 'CUX1', 'MBP',
                  'NR2F1', 'PAX6', 'FOXP2','GAD1','MKI67','FABP7', 'MED13L']
    valide_gene = [j for j in liste_gene if j in adata.var_names]
    for i in valide_gene:
      cr.pl.gene_trends(
	  adata,
          model=model,
          genes=i,
          same_plot=True,
          ncols=2,
          time_key="velocity_pseudotime",
          hide_cells=True,
          n_jobs=1
      )
      plt.savefig(
        figure_dir/ f"gene_trend_cellrank_{i}.pdf",
        facecolor="white",
        bbox_inches="tight",
        dpi=300
      )
      plt.close()

    return adata



