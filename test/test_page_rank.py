import numpy as np
import torinanet as tn

def estimate_rate_constant(rxn: tn.core.Reaction, T: float=600):
        h = 4.135e-15 # eV * sec
        kb = 8.617e-5 # eV / K
        return kb * T / h * np.exp(- max(rxn.energy(), 0) / (kb * T))


def calc_reaction_rate(g, rxn):
    return np.prod([g.nodes[s]["p"] for s in g.predecessors(rxn)]) * g.nodes[rxn]["k"]

def calc_total_out_rate(g, sp):
    return sum([g.nodes[rxn]["rate"] for rxn in g.successors(sp)])

def calc_rate_fraction(g, rxn, sp):
    if g.nodes[sp]["total_rate"] == 0:
        return 0
    else:
        return g.nodes[rxn]["rate"] / g.nodes[sp]["total_rate"]

def page_rank_species(rxn_graph: tn.core.RxnGraph, n_iter: int=10, temperature: float=600):
    """Method to calculate the PageRank metric of a specie"""
    # initializing - converting to networkx graph 
    g = rxn_graph.to_networkx_graph(use_internal_id=True)
    # setting initial probabilities - 1 for source species, 0 for rest
    for node in g:
        g.nodes[node]["p"] = 0
    for sp in rxn_graph.source_species:
        key = rxn_graph.specie_collection.get_key(sp)
        g.nodes[key]["p"] = 1
    # estimating rate constants for reactions
    for rxn in rxn_graph.reactions:
        key = rxn_graph.reaction_collection.get_key(rxn)
        g.nodes[key]["k"] = estimate_rate_constant(rxn, temperature)
    # ==== STARTING PAGERANK ====
    for _ in range(n_iter):
        print("PAGE RANK ITER", _ + 1)
        # assigning reaction rates
        for rxn in rxn_graph.reactions:
            node = rxn_graph.reaction_collection.get_key(rxn)
            g.nodes[node]["rate"] = calc_reaction_rate(g, node)
        # assigning total rates
        for sp in rxn_graph.species:
            node = rxn_graph.specie_collection.get_key(sp)
            g.nodes[node]["total_rate"] = calc_total_out_rate(g, node)
        # calculating new probabilities
        new_ps = {}
        for sp in rxn_graph.species:
            node = rxn_graph.specie_collection.get_key(sp)
            val = 0
            for rxn in g.predecessors(node):
                # taking the product of probability transitions
                val += np.prod([calc_rate_fraction(g, rxn, sp) for sp in g.predecessors(rxn)])
            new_ps[node] = val
        # assigning probabilities
        for sp in rxn_graph.source_species:
            key = rxn_graph.specie_collection.get_key(sp)
            new_ps[key] = 1
        for n, p in new_ps.items():
            g.nodes[n]["p"] = p
    # returning dict with PageRank values on each node
    return new_ps

def get_reactants_page_rank(rxn_graph: tn.core.RxnGraph):
    ranks = page_rank_species(rxn_graph)
    res = []
    for rxn in rxn_graph.reactions:
        ajr = [ranks[rxn_graph.specie_collection.get_key(s)] for s in rxn.reactants]
        # if the reaction has only one reactant (first order), then add 1 to "second reactant"
        # this is made to ensure same length of results
        if len(ajr) < 2:
            ajr.append(1)
        res.append(ajr)
    return res

def main():
    import pandas as pd
    g = tn.core.RxnGraph.from_file("../full_enum/nh3+o2+ch4/2/energy_reduced_graph.rxn")
    ranks = page_rank_species(g, n_iter=4, temperature=600)
    dist_analyzer = tn.analyze.algorithms.ShortestPathAnalyzer(g, prop_func=lambda rxn: 1)
    energy_analyzer = tn.analyze.algorithms.ShortestPathAnalyzer(g, prop_func=lambda rxn: max(0, rxn.energy()))
    res = pd.DataFrame()
    for sp in g.species:
        key = g.specie_collection.get_key(sp)
        energy_dist = sum([rxn.energy() for rxn in energy_analyzer.get_path_to_source(sp)])
        res = res.append({"smiles": sp.identifier, "rank": ranks[key], "distance": dist_analyzer.get_distance_from_source(sp), "energy": energy_dist}, ignore_index=True)
    print(res)
    res.to_csv("page_rank.csv")

if __name__ == "__main__":
    g = tn.core.RxnGraph.from_file("../full_enum/nh3+o2/2/energy_reduced_graph.rxn")
    reducer = tn.analyze.network_reduction.KineticReduction.MolRankReduction(0.01, "k", False, temperature=600)
    g = tn.analyze.kinetics.utils.assign_maximal_rate_constants(g, reducer.temperature, reducer.energy_conversion_factor / 30, reducer.specie_energy_property_name, reducer.rate_constant_property)
    # df = reducer.rank_species(g)
    # df["smiles"] = [sp.identifier for sp in df["specie"].values]
    # df = df.drop(columns=["specie", "visited"])
    # df = df.set_index("smiles")
    # df.to_csv("page_rank.csv")
    print("BEFORE")
    print("species =", g.get_n_species())
    print("reactions =", g.get_n_reactions())
    red_g = reducer.apply(g)
    print("AFTER")
    print("species =", red_g.get_n_species())
    print("reactions =", red_g.get_n_reactions())
    for s in red_g.species:
        print(s.identifier)