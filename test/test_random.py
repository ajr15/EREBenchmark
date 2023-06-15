# script to analyze reaction & specie metrics of enumerated graph
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import OPTICS
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import torinanet as tn

def build_specie_df(rxn_graph: tn.core.RxnGraph, properties: dict):
    g = rxn_graph.to_networkx_graph(use_internal_id=True)
    df = pd.DataFrame()
    for sp in rxn_graph.species:
        d = {name: f(sp, rxn_graph, g) for name, f in properties.items()}
        d["id"] = rxn_graph.specie_collection.get_key(sp)
        df = df.append(d, ignore_index=True)
    return df.set_index("id")

def calc_in_degree(sp: tn.core.Specie, rxn_graph: tn.core.RxnGraph, g: nx.DiGraph):
    key = rxn_graph.specie_collection.get_key(sp)
    return len(list(g.predecessors(key)))

def calc_out_degree(sp: tn.core.Specie, rxn_graph: tn.core.RxnGraph, g: nx.DiGraph):
    key = rxn_graph.specie_collection.get_key(sp)
    return len(list(g.successors(key)))

def calc_specie_mass(sp: tn.core.Specie, rxn_graph: tn.core.RxnGraph, g: nx.DiGraph):
    obmol = sp.parse_identifier()
    return obmol.GetMolWt()

def calc_n_heavy_atoms(sp: tn.core.Specie, rxn_graph: tn.core.RxnGraph, g: nx.DiGraph):
    obmol = sp.parse_identifier()
    return obmol.NumHvyAtoms()

def transformed_properties(col, specie_df, s1, s2=None):
    p1 = specie_df.loc[s1, col]
    p2 = specie_df.loc[s2, col] if s2 is not None else 0
    return [
        p1 + p2,
        p1 * p2 / (p1 + p2) if not (p1 == p2 == 0) else 0
    ]


def embed_reaction(rxn: tn.core.Reaction, rxn_graph: tn.core.RxnGraph, specie_df: pd.DataFrame):
    # adding specie properties (reactants + products)
    embedding = []
    for col in specie_df.columns:
        embedding += transformed_properties(col, specie_df, *[rxn_graph.specie_collection.get_key(s) for s in rxn.reactants])
        embedding += transformed_properties(col, specie_df, *[rxn_graph.specie_collection.get_key(s) for s in rxn.products])
    # adding reaction energy
    embedding.append(rxn.energy())
    return embedding

def build_embedding(rxn_graph: tn.core.RxnGraph):
    """Build a reaction embedding for all reactions in graph as a vector with the following attributes
        - f(r.degree), f(p.degree)
        - g(r.degree), g(p.degree)
        - f(r.mass), f(p.mass)
        - g(r.degree), g(p.degree)
        - f(r.hvy), f(p.hvy)
        - g(r.hvy), g(p.hvy)
    where
        f(x.y) = x1.y + x2.y
        g(x.y) = (x1.y * x2.y) / f(x.y)
    
    f and g are made to ensure symmetry between reactant/products"""
    embeddings = []
    properties = {
                    "in_degree": calc_in_degree,
                    "out_degree": calc_out_degree,
                    "heavy_atoms": calc_n_heavy_atoms,
                    "mass": calc_specie_mass
                }
    specie_df = build_specie_df(rxn_graph, properties)
    rids = []
    for rxn in rxn_graph.reactions:
        rids.append(rxn_graph.reaction_collection.get_key(rxn))
        embeddings.append(embed_reaction(rxn, rxn_graph, specie_df))
    embeddings = np.array(embeddings)
    # returning normalized embeddings
    return rids, (embeddings - np.mean(embeddings, axis=0)) / np.std(embeddings, axis=0)

def find_leafs(rxn_graph):
        """Method to find leaf reactions and species"""
        # collecting all leaf reactions and species
        G = rxn_graph.to_networkx_graph(use_internal_id=True)
        leaf_reactions = []
        for sp in rxn_graph.species:
            sid = rxn_graph.specie_collection.get_key(sp)
            generating_rxns = list(G.predecessors(sid))
            if len(generating_rxns) <= 15:
                leaf_reactions += generating_rxns
        return leaf_reactions

def test_leaf_clusters(path):
    print("reading data...")
    rxn_graph = tn.core.RxnGraph.from_file(path)
    print("making embeddings...")
    rids, embeddings = build_embedding(rxn_graph)
    print("clustering...")
    model = OPTICS()
    categories = model.fit_predict(embeddings)
    rxn_clusters = pd.DataFrame({"rid": rids, "cluster": categories, "rxn": rxn_graph.reactions})
    print("getting leaf species & their similar reactions")
    leaf_rids = find_leafs(rxn_graph)
    leaf_clusters = set([rxn_clusters[rxn_clusters["rid"] == rid].iloc[0, 1] for rid in leaf_rids])
    added_reactions = []
    leaf_rxns = []
    for c in list(leaf_clusters):
        if not c == -1:
            rxns = rxn_clusters[rxn_clusters["cluster"] == c]["rxn"]
            crids = rxn_clusters[rxn_clusters["cluster"] == c]["rid"]
            print("========== CLUSTER {} ==========".format(c))
            for rid, rxn in zip(crids, rxns):
                s = rxn.pretty_string()
                if rid in leaf_rids:
                    print("**", s)
                    leaf_rxns.append(rxn)
                else:
                    print("  ", s)
                    added_reactions.append(rxn)
    print("======================")
    print("N SPECIES", rxn_graph.get_n_species())
    print("N REACTIONS", rxn_graph.get_n_reactions())
    print("N LEAF REACTIONS", len(leaf_rxns))
    print("ADDED REACTIONS", added_reactions)
    print("==== LEAF ONLY REDUCTION ====")
    red_graph = rxn_graph.remove_reactions([rxn for rxn in leaf_rxns if rxn.energy() > 0])
    print("N SPECIES", red_graph.get_n_species())
    print("N REACTIONS", red_graph.get_n_reactions())
    print("==== LEAF+SIMILAR ONLY REDUCTION ====")
    red_graph = rxn_graph.remove_reactions([rxn for rxn in leaf_rxns + added_reactions if rxn.energy() > 0])
    print("N SPECIES", red_graph.get_n_species())
    print("N REACTIONS", red_graph.get_n_reactions())


if __name__ == "__main__":
    test_leaf_clusters("../full_enum/nh3+o2/3/energy_reduced_graph.rxn")