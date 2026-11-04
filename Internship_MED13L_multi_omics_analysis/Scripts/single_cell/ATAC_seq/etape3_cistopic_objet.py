#Chargement des packages
import os
import pandas as pd
import pycisTopic
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
import polars as pl
import pickle
import scrublet as scr
from pycisTopic.lda_models import run_cgs_models_mallet
import matplotlib.pyplot as plt

#Initialisation des chemins
out_dir = "/data/20_data_ayat/scenicplus_workflow/WT/pycistopic"
path_to_regions = os.path.join(out_dir, "consensus_peak_calling/consensus_regions.bed")
path_to_blacklist = "/home/ayat/uv/pycisTopic/hg38-blacklist.v2.bed"
pycistopic_qc_output_dir = "/data/20_data_ayat/scenicplus_workflow/WT/qc"
plot_dir = os.path.join(pycistopic_qc_output_dir, "plots")
mallet_dir = os.path.join(out_dir, "mallet")
tmp_path="/data/20_data_ayat/tmp"
os.makedirs(mallet_dir, exist_ok=True)

fragments_dict = {"WT" : "/data/90_data_jerome/donnees_sf/runs/CO_WT/count_CO_WT/outs/atac_fragments.tsv.gz"}

sample_id_to_barcodes_passing_filters = pickle.load(open(os.path.join(pycistopic_qc_output_dir, "barcodes_passing_filters.pkl"), "rb"))

#Création de l'object cistopic objet avec les peack trouvé et ceux qui ont passé le contrôle qualité. Une matrice est généré ainsi, pour chaque cellule, il y a l'ensemble des peaks identifier
cistopic_obj_list = []
for sample_id in fragments_dict:
    sample_metrics = pl.read_parquet(
        os.path.join(pycistopic_qc_output_dir, f'{sample_id}.fragments_stats_per_cb.parquet')
    ).to_pandas().set_index("CB").loc[ sample_id_to_barcodes_passing_filters[sample_id] ]
    cistopic_obj = create_cistopic_object_from_fragments(
        path_to_fragments = fragments_dict[sample_id],
        path_to_regions = path_to_regions,
        path_to_blacklist = path_to_blacklist,
        metrics = sample_metrics,
        valid_bc = sample_id_to_barcodes_passing_filters[sample_id],
        n_cpu = 1,
        project = sample_id,
        split_pattern = '-'
    )
    cistopic_obj_list.append(cistopic_obj)

cistopic_obj = cistopic_obj_list[0]
print(cistopic_obj)

pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj.pkl"), "wb")
)

scrub = scr.Scrublet(cistopic_obj.fragment_matrix.T, expected_doublet_rate=0.1)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.savefig(os.path.join(plot_dir, "scrublet_histogram.png"))

scrub.call_doublets(threshold=0.22)
scrub.plot_histogram()
plt.savefig(os.path.join(plot_dir, "scrublet_histogram_nodb.png"))
scrublet = pd.DataFrame([scrub.doublet_scores_obs_, scrub.predicted_doublets_], columns=cistopic_obj.cell_names, index=['Doublet_scores_fragments', 'Predicted_doublets_fragments']).T

cistopic_obj.add_cell_data(scrublet, split_pattern = '-')
sum(cistopic_obj.cell_data.Predicted_doublets_fragments == True)

pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj.pkl"), "wb")
)

# Possibilité de retirer les doublons
singlets = cistopic_obj.cell_data[cistopic_obj.cell_data.Predicted_doublets_fragments == False].index.tolist()

# Créer un objet cistopic sans doublons pour voir si c'est mieux ou non
cistopic_obj_noDBL = cistopic_obj.subset(singlets, copy=True, split_pattern='-')
print(cistopic_obj_noDBL)

pickle.dump(
    cistopic_obj_noDBL,
    open(os.path.join(out_dir, "cistopic_obj_noDBL.pkl"), "wb")
)

os.environ['MALLET_MEMORY'] = '200G'

# Configuration du chemin de l'outil MALLET
MALLET_path="/home/ayat/uv/Mallet-202108/bin/mallet"

#Faire le modèle LDA, qui permet d'identifié des groupes de régions co-accessibles, appelé topics, afin de voir si des programmes chromatiniens sont en cours
models=run_cgs_models_mallet(
    cistopic_obj,
    n_topics=[2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50],
    n_cpu=12,
    n_iter=500,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path= tmp_path,
    save_path=mallet_dir,
    mallet_path=MALLET_path,
)

pickle.dump(
    models,
    open(os.path.join(out_dir, "models.pkl"), "wb")
)

#Prendre le meilleur modèle et l'ajouter dans l'objet, l'idéal est de prendre un modèle ni trop petit ni trop grand
best_model = models[5]

cistopic_obj.add_LDA_model(best_model)

pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj_with_model.pkl"), "wb")
)

#Effecctué avec l'objet sans les doublons, mais le résultas est le même
#models=run_cgs_models_mallet(
#    cistopic_obj_noDBL,
#    n_topics=[2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50],
#    n_cpu=12,
#    n_iter=500,
#    random_state=555,
#    alpha=50,
#    alpha_by_topic=True,
#    eta=0.1,
#    eta_by_topic=False,
#    tmp_path= tmp_path,
#    save_path=mallet_dir,
#    mallet_path=MALLET_path,
#)

#pickle.dump(
#    models,
#    open(os.path.join(out_dir, "models_noDBL.pkl"), "wb")
#)
