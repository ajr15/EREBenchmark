# script with util functions to analyze structure of a reaction graph
# yeilds the following
#   - reaction energy distribution
#   - distribution of reactions by # reactants and # products
#   - distribution of "energy shortest paths" of species
import torinanet as tn
import pandas as pd
import matplotlib.pyplot as plt
import os
from typing import List
from itertools import chain

def calc_reaction_energy(rxn: tn.core.Reaction):
    return sum([s.properties["energy"] for s in rxn.products]) - sum([s.properties["energy"] for s in rxn.reactants])

def reaction_energy_distribution(g: tn.core.RxnGraph, figure_dir: str):
    energies = []
    for r in g.reactions:
        x = calc_reaction_energy(r)
        energies.append(x)
    plt.figure()
    plt.hist(energies, color="gray")
    plt.xlabel("Reaction Energy (Ha)")
    plt.ylabel("N Reactions")
    plt.savefig(os.path.join(figure_dir, "reaction_energy.png"))
    plt.close()

def reactants_and_product_count_distribution(g: tn.core.RxnGraph, figure_dir: str):
    # counting reaction types
    counter = {}
    for r in g.reactions:
        s = "{} -> {}".format(len(r.reactants), len(r.products))
        if s in counter:
            counter[s] += 1
        else:
            counter[s] = 1
    # plotting
    plt.figure()
    plt.bar(
                x=range(1, len(counter) + 1), 
                height=counter.values(),
                color="gray"
            )
    plt.xticks(range(1, len(counter) + 1), counter.keys())
    plt.savefig(os.path.join(figure_dir, "reaction_product_count.png"))
    plt.close()

def shortest_path_distribution(g: tn.core.RxnGraph, figure_dir: str):
    analyzer = tn.analyze.algorithms.ShortestPathAnalyzer(g, prop_func=lambda r: max(0, calc_reaction_energy(r)))
    analyzer.to_dataframe().to_csv(os.path.join(figure_dir, "shortest_path_data.csv"))
    plt.figure()
    plt.hist(analyzer.shortest_path_table["dist"], color="gray")
    plt.xlabel("Shortest Path Energy (Ha)")
    plt.ylabel("N Species")
    plt.savefig(os.path.join(figure_dir, "shortest_path.png"))
    plt.close()

def get_formula_dict(ac_mats: List[tn.core.AcMatrix]) -> dict:
    ajr = {}
    for atom in chain(*[ac.get_atoms() for ac in ac_mats]):
        if atom not in ajr:
            ajr[atom] = 1
        else:
            ajr[atom] += 1
    return ajr

def formula_dict_to_str(formula_d: dict, idx_d: dict) -> str:
    # making sure all keys (atoms) are indexed
    for k in formula_d.keys():
        if not k in idx_d:
            idx_d[k] = len(idx_d)
    # now making string
    ajr = [0 for _ in range(len(idx_d))]
    for k, v in formula_d.items():
        ajr[idx_d[k]] = v
    return ",".join([str(x) for x in ajr])

def proper_formula_str(formula_str: str, idx_d: dict):
    symbols = {1: "H", 6: "C", 7: "N", 8: "O"}
    reverse_idx = {v: k for k, v in idx_d.items()}
    numbers = formula_str.split(",")
    ajr = ""
    for i, n in enumerate(numbers):
        ajr += symbols[reverse_idx[i]]
        ajr += str(n)
    return ajr

def empirical_formula_fit(g: tn.core.RxnGraph, figure_path: str, plot: bool=True):
    """Shows random graph parameters for different stochiometric formulas"""
    # first collecting all different empirical formulas and counting members
    idx_d = {6: 0}
    empirical_d = {}
    for s in g.species:
        ajr = formula_dict_to_str(get_formula_dict([s.ac_matrix]), idx_d)
        if ajr in empirical_d:
            empirical_d[ajr] += 1
        else:
            empirical_d[ajr] = 1
    if plot:
        labels = [proper_formula_str(k, idx_d) for k in empirical_d.keys()]
        values = empirical_d.values()
        for l, v in zip(labels, values):
            print(l, v)
        # sorting by frequency
        # labels = sorted(labels, key=empirical_d.values())
        # values = sorted(empirical_d.values())
        # plotting
        plt.figure()
        plt.bar(range(1, len(labels) + 1), values)
        plt.xticks(range(1, len(labels) + 1), labels)
        plt.savefig(figure_path)
        plt.close()
    return empirical_d
    
    

def run_analysis(graph_path: str, figure_dir: str, analysis_dict: dict):
    print("Analyzing", graph_path)
    g = tn.core.RxnGraph.from_file(graph_path)
    # creating resutls dir
    if not os.path.isdir(figure_dir):
        os.mkdir(figure_dir)
    # running analysis
    for name, func in analysis_dict.items():
        print("Running", name)
        # try:
        func(g, figure_dir)
        # except:
        #     print("Errors with {} analysis, skipping".format(name))
    print("DONE")

def degree_distribution(g: tn.core.RxnGraph, results_dir: str, plot: bool=True):
    G = g.to_networkx_graph(use_internal_id=True)
    # calculating specie degrees
    res = pd.DataFrame()
    for sp in g.species:
        sid = g.specie_collection.get_key(sp)
        degree = len(list(G.neighbors(sid)))
        res = res.append({"smiles": sp.identifier, "degree": degree}, ignore_index=True)
    # plotting and saving
    if plot:
        plt.figure()
        plt.hist(res["degree"], color="gray", bins=50)
        plt.savefig(os.path.join(results_dir, "degree_distribution.png"))
        plt.close()
    res.to_csv(os.path.join(results_dir, "specie_degree.csv"))

def specie_properties_summary(g: tn.core.RxnGraph, results_dir: str):
    """Summary of specie properties - degree (in and out), distance from source (# reactions), # shortest paths, molar mass, # atoms"""
    analyzer = tn.analyze.algorithms.ShortestPathAnalyzer(g, prop_func=lambda rxn: max(0, calc_reaction_energy(rxn)))
    n_paths_df = analyzer.number_of_paths_per_specie()
    n_paths_df = n_paths_df.set_index("smiles")
    G = analyzer.networkx_graph
    # calculating specie degrees
    res = pd.DataFrame()
    for sp in g.species:
        sid = g.specie_collection.get_key(sp)
        in_degree = len(list(G.predecessors(sid)))
        out_degree = len(list(G.successors(sid)))
        obmol = sp.ac_matrix.to_obmol()
        res = res.append({
            "smiles": sp.identifier, 
            "in_degree": in_degree, 
            "out_degree": out_degree, 
            "n_atoms": obmol.NumAtoms(), 
            "mass": obmol.GetExactMass(),
            "min_energy": sum([calc_reaction_energy(r) for r in analyzer.get_path_to_source(sp)]),
            "n_consumed_paths": n_paths_df.loc[sp.identifier, "consumed"],
            "n_created_paths": n_paths_df.loc[sp.identifier, "created"],
            "n_intermediate_paths": n_paths_df.loc[sp.identifier, "intermediate"],
            }, ignore_index=True)
    res.to_csv(os.path.join(results_dir, "specie_properties.csv"))
    plt.figure()
    plt.hist(res["mass"], color="gray")
    plt.xlabel("Mass Distribution")
    plt.ylabel("N")
    plt.savefig(os.path.join(results_dir, "mass_dist.png"))
    plt.close()
    

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("graph_path", type=str, help="full path for graph to analyze")
    parser.add_argument("fig_dir", type=str, help="full path to directory with results")
    args = parser.parse_args()
    if not os.path.isdir(args.fig_dir):
        os.mkdir(args.fig_dir)
    analysis_dict = {
        "reaciton_energy_distribution": reaction_energy_distribution,
        "reactants_and_product_count_distribution": reactants_and_product_count_distribution,
        "shortest_path_distribution": shortest_path_distribution,
        # "empirical_formula_fit": empirical_formula_fit,
        "degree_distribution": degree_distribution,
        "specie_summary": specie_properties_summary
    }
    run_analysis(args.graph_path, args.fig_dir, analysis_dict)
