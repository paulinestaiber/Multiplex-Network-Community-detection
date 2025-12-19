# Part 2 of pileline
# generate input file for multilayer community detection for infomap
# How to run it: you need to give directory of the network files, result file of LCC, Output file 
# python3 build_infomap_multilayer.py \  --network-dir /Users/paulinestaiber/Documents/Network/Multiplex-Network-Community-detection/Data/raw/multiplex_InfoWalkR \
#  --lcc-csv /Users/paulinestaiber/Documents/Network/Multiplex-Network-Community-detection/results/lcc_sig_layer.csv \
#  --output-net /Users/paulinestaiber/Documents/Network/Multiplex-Network-Community-detection/results/multiplex_clustering.net \
#  --save-intermediate


import argparse # for reading command-line arguments
from pathlib import Path
import pandas as pd

#Function that loads all layers of a multilayer network stored as .tsv files in a directory
def load_layers(network_dir: Path, layers_to_skip=None):
    """
    Load all .tsv edge files from network_dir, skip some by name,
    and assign clean integer layer IDs starting from 1.
    """
    if layers_to_skip is None:
        layers_to_skip = []

    tsv_files = sorted(network_dir.glob("*.tsv"))
    if not tsv_files:
        raise FileNotFoundError(f"No TSV files found in {network_dir}")

    layers = {}       # layer_id -> DataFrame
    layer_names = {}  # layer_id -> original layer name (stem)

    new_id = 1
    for file in tsv_files:
        layer_name = file.stem
        if layer_name in layers_to_skip:
            continue

        df = pd.read_csv(file, comment="#", header=0, sep="\t", dtype=str)
        df.columns = ["source", "target"]
        df["layer_id"] = new_id

        layer_names[new_id] = layer_name
        layers[new_id] = df

        new_id += 1

    return layers, layer_names

# Function that filters layers based on LCC results
def filter_layers_by_lcc(layers, layer_names, lcc_df):
    """
    Keep only layers whose names appear in lcc_df['network_name'].
    """
    lcc_layer_names = lcc_df["network_name"].tolist()

    filtered_layers = {
        lid: df for lid, df in layers.items()
        if layer_names[lid] in lcc_layer_names
    }
    filtered_layer_names = {
        lid: name for lid, name in layer_names.items()
        if name in lcc_layer_names
    }

    return filtered_layers, filtered_layer_names

# Function that builds a global mapping from node names to unique integer IDs
def build_node_ids(layers):
    """
    Build a global mapping from node name -> node_id (int).
    """
    all_nodes = set()
    for df in layers.values():
        all_nodes.update(df["source"].astype(str))
        all_nodes.update(df["target"].astype(str))

    node_id_map = {name: i for i, name in enumerate(all_nodes, start=1)}
    return node_id_map


def build_vertices_df(node_id_map):
    vertices_df = pd.DataFrame(
        [(i, name) for name, i in node_id_map.items()],
        columns=["node_id", "name"]
    )
    vertices_df = vertices_df.sort_values("node_id")
    return vertices_df

# Funktion that builds intra-layer edges
def build_intra_edges(layers, node_id_map):
    """
    Create intra-layer edges: (layer_id, source_id, target_id, weight).
    """
    intra_rows = []

    for layer_id, df in layers.items():
        df_source = df["source"].astype(str)
        df_target = df["target"].astype(str)

        for s, t in zip(df_source, df_target): #pairs the elements
            u = node_id_map[s]
            v = node_id_map[t]
            w = 1.0  # unweighted
            intra_rows.append((layer_id, u, v, w))

    intra_df = pd.DataFrame(
        intra_rows,
        columns=["layer_id", "source_id", "target_id", "weight"]
    )
    return intra_df


# Function that builds a dictionary of nodes present in each layer (needed for inter-layer edges)
def build_layer_nodes(intra_df):
    """
    For each layer, build the set of node_ids present in that layer.
    """
    layer_nodes = {}
    for layer_id, sub in intra_df.groupby("layer_id"):
        nodes_in_layer = set(sub["source_id"]).union(sub["target_id"])
        layer_nodes[layer_id] = nodes_in_layer
    return layer_nodes

# Function that attaches layer IDs to LCC results
def attach_layer_ids_to_lcc(layer_names, lcc_df):
    """
    Map LCC network_name -> layer_id using layer_names dict,
    and attach layer_id column to lcc_df.
    """
    layer_name_to_id = {name: lid for lid, name in layer_names.items()}
    lcc_df = lcc_df.copy()
    lcc_df["layer_id"] = lcc_df["network_name"].map(layer_name_to_id)
    return lcc_df

# Function that builds a dataframe with layer IDs and z-scores
def build_all_layers_df(layer_names, lcc_df):
    """
    Build a DataFrame with layer_id and z_score for all kept layers.
    """
    all_layers_df = pd.DataFrame({"layer_id": list(layer_names.keys())})
    all_layers_df = all_layers_df.merge(
        lcc_df[["layer_id", "z_score"]],
        on="layer_id",
        how="left"
    )
    return all_layers_df

# Function that builds inter-layer edges
def build_inter_edges(layer_nodes, all_layers_df):
    """
    For each ordered pair of layers (L1, L2, L1 != L2), for each shared node:
    create an inter-layer edge with weight = z_score of L2.
    """
    inter_rows = []

    for L1 in sorted(layer_nodes.keys()):
        for L2 in sorted(layer_nodes.keys()):
            if L1 == L2:
                continue

            shared = layer_nodes[L1] & layer_nodes[L2]
            # weight is z-score of layer L2
            z_vals = all_layers_df.loc[all_layers_df["layer_id"] == L2, "z_score"].values
            if len(z_vals) == 0 or pd.isna(z_vals[0]):
                # if no z_score available,  set default
                omega_L2 = 0.0
            else:
                omega_L2 = float(z_vals[0])

            for nid in shared:
                inter_rows.append((L1, nid, L2, omega_L2))

    inter_df = pd.DataFrame(
        inter_rows,
        columns=["layer_id_1", "node_id", "layer_id_2", "weight"]
    )
    return inter_df


def write_infomap_multilayer(
    output_path: Path,
    vertices_df: pd.DataFrame,
    intra_df: pd.DataFrame,
    inter_df: pd.DataFrame
):
    """
    Write a .net file in Infomap multilayer format:
    *Vertices
    *Multilayer
    with intra- and inter-layer edges.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as f:
        # optional comment
        f.write("# Generated from Python\n")

        # -------- *Vertices --------
        f.write(f"*Vertices {len(vertices_df)}\n")
        f.write("# node_id name\n")
        for _, row in vertices_df.iterrows():
            node_id = int(row["node_id"])
            name = str(row["name"])
            f.write(f'{node_id} "{name}"\n')

        # -------- *Multilayer --------
        f.write("*Multilayer\n")
        f.write("# layer_id node_id layer_id node_id weight\n")

        # ----- intra edges -----
        f.write("# intra\n")

        # unify column naming
        intra_tmp = intra_df.rename(columns={
            "source_id": "node_id_u",
            "target_id": "node_id_v"
        })

        for _, row in intra_tmp.iterrows():
            L = int(row["layer_id"])
            u = int(row["node_id_u"])
            v = int(row["node_id_v"])
            w = float(row["weight"])
            f.write(f"{L} {u} {L} {v} {w}\n")

        # ----- inter edges -----
        f.write("# inter\n")

        for _, row in inter_df.iterrows():
            L1 = int(row["layer_id_1"])
            L2 = int(row["layer_id_2"])
            nid = int(row["node_id"])
            w = float(row["weight"])
            f.write(f"{L1} {nid} {L2} {nid} {w}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Build Infomap multilayer .net file from multiplex TSV layers and LCC results."
    )
    parser.add_argument(
        "--network-dir",
        required=True,
        type=Path,
        help="Directory containing layer .tsv files (edges)."
    )
    parser.add_argument(
        "--lcc-csv",
        required=True,
        type=Path,
        help="CSV file with LCC results (at least columns: network_name, z_score)."
    )
    parser.add_argument(
        "--output-net",
        required=True,
        type=Path,
        help="Output .net file for Infomap multilayer."
    )
    parser.add_argument(
        "--skip-layers",
        nargs="*",
        default=["coex_ITI", "coex_LVR", "coex_MSG", "coex_LCL", "coex_UTR", "reactome_copathway"],
        help="Layer names (stem of TSV file) to skip."
    )
    parser.add_argument(
        "--save-intermediate",
        action="store_true",
        help="If set, also save layer_names, vertices, intra_edges, inter_edges as separate files."
    )

    args = parser.parse_args()

    # 1) Load layers
    layers, layer_names = load_layers(args.network_dir, layers_to_skip=args.skip_layers)
    print(f"Loaded {len(layers)} layers after skipping {len(args.skip_layers)}")

    # 2) Load LCC results
    lcc_df = pd.read_csv(args.lcc_csv)
    print(f"LCC file has {lcc_df.shape[0]} rows")

    # 3) Filter layers to those present in LCC results
    layers, layer_names = filter_layers_by_lcc(layers, layer_names, lcc_df)
    print(f"{len(layers)} layers kept after LCC filtering")

    if not layers:
        raise RuntimeError("No layers left after LCC filtering. Check your LCC CSV and layer names.")

    # 4) Build node IDs and vertices
    node_id_map = build_node_ids(layers)
    print(f"Total unique nodes: {len(node_id_map)}")

    vertices_df = build_vertices_df(node_id_map)

    # 5) Build intra-layer edges
    intra_df = build_intra_edges(layers, node_id_map)
    print(f"Number of intra-layer edges: {len(intra_df)}")

    # 6) Build layer_nodes (which nodes are present in each layer)
    layer_nodes = build_layer_nodes(intra_df)

    # 7) Attach layer_id to LCC table and get z_scores per layer
    lcc_df = attach_layer_ids_to_lcc(layer_names, lcc_df)
    all_layers_df = build_all_layers_df(layer_names, lcc_df)

    # 8) Build inter-layer edges
    inter_df = build_inter_edges(layer_nodes, all_layers_df)
    print(f"Number of inter-layer edges: {len(inter_df)}")

    # 9) Write Infomap multilayer file
    write_infomap_multilayer(args.output_net, vertices_df, intra_df, inter_df)
    print("Written Infomap multilayer file to:", args.output_net)

    # 10) Optionally save intermediate files
    if args.save_intermediate:
        base = args.output_net.with_suffix("")  # remove .net
        layer_names_path = base.parent / f"{base.stem}_layer_names.tsv"
        vertices_path = base.parent / f"{base.stem}_vertices.tsv"
        intra_path = base.parent / f"{base.stem}_intra_edges.csv"
        inter_path = base.parent / f"{base.stem}_inter_edges.csv"

        layer_names_df = pd.DataFrame(
            list(layer_names.items()),
            columns=["layer_id", "layer_name"]
        )
        layer_names_df.to_csv(layer_names_path, sep="\t", index=False)
        vertices_df.to_csv(vertices_path, sep="\t", index=False)
        intra_df.to_csv(intra_path, index=False)
        inter_df.to_csv(inter_path, index=False)

        print("Saved intermediate files:")
        print("  layer names  ->", layer_names_path)
        print("  vertices     ->", vertices_path)
        print("  intra edges  ->", intra_path)
        print("  inter edges  ->", inter_path)


if __name__ == "__main__":
    main()

