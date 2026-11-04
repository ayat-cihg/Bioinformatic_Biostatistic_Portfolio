#Chargement des packages
import os
import pycisTopic
import pandas as pd
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk
from pycisTopic.pseudobulk_peak_calling import peak_calling
from pycisTopic.iterative_peak_calling import get_consensus_peaks
import pyranges as pr

#Initialisation des packages, il est nécessaire d'avoir le fichier fragments.tsv et celui des types cellulaires
data_dir_fragment="/data/90_data_jerome/donnees_sf/runs/CO_KO/count_CO_KO/outs/atac_fragments.tsv.gz"
data_dir = "/data/20_data_ayat/scenicplus_workflow/KO"
output_dir=os.path.join(data_dir, "pycistopic")
cell_data = "/data/20_data_ayat/CO_KO/cell_types.csv"
temp_dir = "/data/20_data_ayat/tmp"
os.makedirs(data_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

#Lecture de nos fichiers
cell_type=pd.read_csv(cell_data)
print(cell_type.head)
cell_type = cell_type.rename(columns={"Barcode" : "cell_id"})
cell_type["cell_types"]=cell_type["cell_types"].fillna('Unknown')
cell_type["cell_types"]=cell_type["cell_types"].str.replace(" ", "_")
print(cell_type.head())

cell_type["sample_id"]="KO"
cell_type = cell_type.set_index("cell_id")
cell_type = cell_type.dropna()

fragment_dict = { "KO" : data_dir_fragment }

#Etape de pseudobulk ATAC-seq
#Cette étape va permettre de récupérer un fichier, où la taille des chromosomes est indiquée, car les les outils d'ATAC-seq doivent le savoir pour éviter de générer des régions impossibles
chromsizes = pd.read_table(
           "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
           header = None,
           names = ["Chromosome", "End"]
           )
chromsizes.insert(1, "Start", 0)
print(chromsizes.head())

os.makedirs(os.path.join(output_dir, "consensus_peak_calling"), exist_ok = True)

#Fichier bed utilisé pour MACS2 et le peak calling
os.makedirs(os.path.join(output_dir, "consensus_peak_calling/pseudobulk_bed_files"), exist_ok = True)
#Fichier bw utilisé pour une visualisation IGV, contient les couverture génomique
os.makedirs(os.path.join(output_dir, "consensus_peak_calling/pseudobulk_bw_files"), exist_ok = True)

print("Les répertoires ont été créer")

#Générer les profiles ATAC seq via la fonction export_pseudobulk
#Cette fonction va rescencer toutes les informations des fragments d'ATAC-seq trouvé.
print("Les profilts ATAC seq vont être générés")

bw_paths, bed_paths = export_pseudobulk(
    input_data = cell_type,
    variable = "cell_types",
    sample_id_col = "sample_id",
    chromsizes = chromsizes,
    bed_path = os.path.join(output_dir, "consensus_peak_calling/pseudobulk_bed_files"),
    bigwig_path = os.path.join(output_dir, "consensus_peak_calling/pseudobulk_bw_files"),
    path_to_fragments = fragment_dict,
    n_cpu = 5,
    normalize_bigwig = True,
    temp_dir = temp_dir,
    split_pattern = "-"
    )

print("Les profiles ont été générés")
#Les fichiers bed continet les données des fragments d'ATAC sous forme de coordonnées, chaque lignes contient une région d'ADN accessible et détécté par l'ATAC-seq
#Les fichiers bw indique combien de fragments ATAC couvrent chaque position du génome, pour la visualisation des profils ATAC et leurs niveaux d'accessibilité entre deux régions.
#Sauvegarder des fichiers bed et bw
print("Sauvegarde des fichiers en cours")

with open(os.path.join(output_dir, "consensus_peak_calling/bed_paths.tsv"), "wt") as f:
    for v in bed_paths:
        _ = f.write(f"{v}\t{bed_paths[v]}\n")

with open(os.path.join(output_dir, "consensus_peak_calling/bw_paths.tsv"), "wt") as f:
    for v in bw_paths:
        _ = f.write(f"{v}\t{bw_paths[v]}\n")

print("Les fichiers bed et bw ont été généré")

#Etapes de MACS
#C'est l'étape de peak calling, c'est à dire identifier les régions ouvertes enrichies des chromatines, et élimine le bruit
print("Etapes de MACS")

bw_paths = {}
with open(os.path.join(output_dir, "consensus_peak_calling/bw_paths.tsv")) as f:
    for line in f:
        v, p = line.strip().split("\t")
        bw_paths.update({v: p})

bed_paths = {}
with open(os.path.join(output_dir, "consensus_peak_calling/bed_paths.tsv")) as f:
    for line in f:
        v, p = line.strip().split("\t")
        bed_paths.update({v: p})

print("Recupération des fichiers bed et bw")

macs_path = "macs2"

os.makedirs(os.path.join(output_dir, "consensus_peak_calling/MACS"), exist_ok = True)

narrow_peak_dict = peak_calling(
    macs_path = macs_path,
    bed_paths = bed_paths,
    outdir = os.path.join(os.path.join(output_dir, "consensus_peak_calling/MACS")),
    genome_size = 'hs',
    n_cpu = 5,
    input_format = 'BEDPE',
    shift = 73,
    ext_size = 146,
    keep_dup = 'auto',
    q_value = 0.01,
    _temp_dir = os.path.join(temp_dir, "ray_spill")
)

print("Etapes MACS fini")

#Autre paramètres, on normalise les données et on utilise le fichier blackliste qui est un fichier qui recense les régions de l'ADN qui sont connues pour être du bruit et ainsi supprimer les peaks touvés ici
print("Obtenir le consensus peak")

peak_half_width=150
path_to_blacklist="/home/ayat/uv/pycisTopic/hg38-blacklist.v2.bed"

#Obtenir le consensus peak
consensus_peaks = get_consensus_peaks(
    narrow_peaks_dict = narrow_peak_dict,
    peak_half_width = peak_half_width,
    chromsizes = chromsizes,
    path_to_blacklist = path_to_blacklist)

# sécurisation support multi-samples
# filtre PyRanges-compatible
if "n_samples" in consensus_peaks.columns:
    consensus_peaks = consensus_peaks[
        consensus_peaks.df["n_samples"] >= 2
    ]

# tri on garde seulement les 220 000 peaks avec un score MACS le plus élevés afin de limiter la mémoire et avoir des données plus significatives
df = consensus_peaks.df

score_col = "Score" if "Score" in df.columns else df.columns[4]

df = df.sort_values(by=score_col, ascending=False).head(220000)

consensus_peaks = pr.PyRanges(df)

#Sauvegarde du consensus

consensus_peaks.to_bed(
    path = os.path.join(output_dir, "consensus_peak_calling/consensus_regions.bed"),
    keep =True,
    compression = 'infer',
    chain = False)

print("Sauvegarde du consensus peak")
