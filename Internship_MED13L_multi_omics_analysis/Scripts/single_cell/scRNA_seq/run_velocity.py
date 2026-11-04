#Chargement des packages
from package_cellrank import *
import argparse
import matplotlib.pyplot as plt

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="Run RNA velocity analysis")

    parser.add_argument("adata", type=str, help="path to h5ad file")
    parser.add_argument("loom", type=str, help="path to velocyto file")
    parser.add_argument("--adata_2", required=False, type=str, help="path to the adata file")
    parser.add_argument("--loom_2", required=False, type=str, help="path to the velocyto file")
    parser.add_argument("--percentage", type=float, help="percentage to elinate the outlier", default=0.01)
    parser.add_argument("--mode", choices=["dynamical", "stochastic", "deterministic"], type=str, help="mode to calculate velocity", default="stochastic")
    parser.add_argument("--split", type=int, help="Number of data split")
    parser.add_argument("--annotation_1", type=str, help="file annotation for cell type")
    parser.add_argument("--annotation_2", type=str, help="file annotation for cell type")

    args=parser.parse_args()

    if args.adata_2:
      data_name="CO_AGGR"
    elif "WT" in args.adata:
      data_name="CO_WT"
    elif "KO" in args.adata:
      data_name="CO_KO"

    data_output_dir = pathlib.Path(f"/data/20_data_ayat/{data_name}")
    figure_dir = pathlib.Path(f"/data/20_data_ayat/{data_name}/Figure")
    data_output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    scv.settings.set_figure_params("scvelo", facecolor="white")
    scv.settings.figdir = str(figure_dir)

    #Si il y a deux adata, on va les aggrégés ensembles
    if args.adata_2:
      adata = aggregation(args.adata, args.loom, args.adata_2, args.loom_2, data_output_dir, data_name, args.annotation_1, args.annotation_2)
    else :
      adata = sc.read_10x_h5(args.adata)
      #Ajout des informations épissées et non épissées indispensable pour le calculs de la vélocité
      if "spliced" not in adata.layers or "unspliced" not in adata.layers:
        adata = fusion_velocyto(adata, args.loom, data_output_dir, args.annotation_1)
    min_gene, max_gene = chose_parameters(adata, args.percentage)
    adata = filtration(adata, min_gene, max_gene, data_output_dir, data_name)
    if args.mode =="dynamical":
      print("dynamical mode is running")
      adata = parallel_recover_dynamics(adata, args.split)
    elif args.mode=="stochastic":
      print("stochastic mode is running")
      scv.tl.velocity(adata, mode="stochastic", show_progress_bar=True)
    else :
      print("deterministic mode is running")
      scv.tl.velocity(adata, mode="deterministic", show_progress_bar=True)

    adata.write(data_output_dir / f"adata_{args.mode}.h5ad")

    scv.tl.velocity_graph(adata)
    scv.tl.velocity_embedding(adata, basis='umap')
    plt.style.use("default")

    #Coloration des UMAP en fonction de l'expression des marqueurs
    scv.pl.velocity_embedding_stream(adata, basis='umap', color="cell_types", save=f"velocity_umap_{data_name}.png")
    cell=['SOX2', 'FOXG1', 'TBR2', 'TBR1', 'BCL11B', 'SATB2', 'CUX1', 'MBP','PRDM1',
    'RAX', 'VSX2', 'OTX2', 'CRX','NEUROD1', 'RCVRN', 'RORB', 'THRB', 'USH2A','NRL', 'NR2E3',
    'NR2F1', 'PAX6', 'FOXP2','GAD1','MKI67','FABP7', 'NES', 'MED13L']
    valide_gene=[g for g in cell if g in adata.var_names]
    for i in valide_gene :
      sc.pl.umap(adata, color=i, s=50, frameon=False, show=False, return_fig=True, cmap='YlGn')
      plt.savefig(figure_dir /f"umap_{i}.pdf")
      plt.close()
    vk = cr.kernels.VelocityKernel(adata)
    vk.compute_transition_matrix()
    ck = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
    pk = 0.8*vk+0.2*ck
    pk.plot_projection(color="cell_types",save=str(figure_dir/f"cellrank_velocity_{data_name}.png"))
    print ("start pseudotime")
    pseudotime(adata, ck, figure_dir, data_output_dir)
