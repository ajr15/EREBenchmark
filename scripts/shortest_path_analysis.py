# script to analyze reaction pathways using k-shortest paths
import pandas as pd
import numpy as np
from typing import List
from itertools import chain
from torinanet.core import Specie, RxnGraph, Reaction
from torinanet.analyze.network_reduction.KineticReduction import MolRankReduction, assign_maximal_rate_constants
from torinanet.analyze.algorithms import ShortestPathAnalyzer

def initialize_algorithm(rxn_graph: RxnGraph):
    """Initializes the algorithm from a reaction graph object. Calculates rates, concentrations and all other metrices, converts into networkx graph and returns list of pathways"""
    # calculating the rates and concentrations using the MolRank-specie algorithm
    analyzer = MolRankReduction(1, "k", True)
    rxn_graph = assign_maximal_rate_constants(rxn_graph, analyzer.temperature, analyzer.energy_conversion_factor, analyzer.specie_energy_property_name, analyzer.rate_constant_property)
    g = analyzer.rank_species(rxn_graph, return_network=True)
    # collecting information from the graph to the list of species & reactions
    for sp in rxn_graph.species:
        key = rxn_graph.specie_collection.get_key(sp)
        sp.properties["consumption_rate"] = sum([g.nodes[rxn]["rate"] if "rate" in g.nodes[rxn] else 0 for rxn in g.successors(key)])
        sp.properties["production_rate"] = sum([g.nodes[rxn]["rate"] if "rate" in g.nodes[rxn] else 0 for rxn in g.predecessors(key)])
        sp.properties["lifetime"] = g.nodes[key]["p"] / sp.properties["consumption_rate"]
    # now adiding properties to reactions
    for rxn in rxn_graph.reactions:
        key = rxn_graph.reaction_collection.get_key(rxn)
        rxn.properties["rate"] = g.nodes[key]["rate"]
    return rxn_graph

def equivalent_paths(rxn_graph, p1, p2):
    """Method to quickly decide if 2 reaction paths are equal using the hashed collection of the reaction graph"""
    s1 = set([rxn_graph.reaction_collection.get_key(r) for r in p1])
    s2 = set([rxn_graph.reaction_collection.get_key(r) for r in p2])
    return s1 == s2

def yens_k_shortest_path(rxn_graph, target, k, prop_func):
    # initializing
    paths = []
    # adding the 1st shortest path
    analyzer = ShortestPathAnalyzer(rxn_graph, prop_func=prop_func)
    paths.append(list(reversed(analyzer.get_path_to_source(target))))
    # initializing set to store potential kth shortest path
    candidates = []
    # now, going to the main loop of Yen's algorithm
    # going over 
    for i in range(1, k):
        # for each round, we go over all paths in the shortest paths and select one edge from each. 
        # we then remove all of the selected nodes from the parent graph and find the shortest path in the resulting graph
        for j in range(len(paths[-1])):
            # we note the "root path" of edges to select where to start the cutting
            root_path = paths[-1][:j]
            # now we are collecting edges from all paths sharing the same root path & removing them
            reduced_rxn_graph = rxn_graph
            for path in paths:
                if equivalent_paths(rxn_graph, root_path, path[:j]):
                    rxn = path[j]
                    if reduced_rxn_graph.has_reaction(rxn):
                        reduced_rxn_graph = reduced_rxn_graph.remove_reaction(rxn)
            # finally, we find shortest path in the reduced graph (only asuming that target specie is still in graph)
            if reduced_rxn_graph.has_specie(target):
                analyzer = ShortestPathAnalyzer(reduced_rxn_graph, prop_func=prop_func)
                candidate = analyzer.get_path_to_source(target)
                # if candidate is not in candidates, add it
                if all([not equivalent_paths(rxn_graph, candidate, c) for c in candidates]):
                    candidates.append(list(reversed(candidate)))
    # after we added candidates, we now select the kth shortest path
    candidates = sorted(candidates, key=lambda path: sum(prop_func(rxn) for rxn in path))
    paths.append(candidates[i - 1])
    return paths

def _specie_coeff(rxn_graph: RxnGraph, specie: Specie, rxn: Reaction) -> int:
    """Method to get the rate of a specie (production or consumption) with respect to its molar coefficients in it"""
    sid = rxn_graph.specie_collection.get_key(specie)
    return len([1 for s in rxn.products if rxn_graph.specie_collection.get_key(s) == sid]) - len([1 for s in rxn.reactants if rxn_graph.specie_collection.get_key(s) == sid])

def specie_coeff(rxn_graph: RxnGraph, specie: Specie, path: List[Reaction]) -> int:
    c = 0
    for rxn in path:
        c += _specie_coeff(rxn_graph, specie, rxn)
    return c

def path_reactants(rxn_graph: RxnGraph, path: List[Reaction]):
    res = []
    for s in chain(*[r.reactants for r in path]):
        if specie_coeff(rxn_graph, s, path) < 0 and not s in res:
            res.append(s)
    return res

def valid_pathway(rxn_graph: RxnGraph, specie_df: pd.DataFrame, path: List[Reaction]):
    """Method to guarrenty that the suggested max rate path is valid. this is done by ensuring that all consumed species in the pathway are the source species"""
    if None in path:
        return False
    total_path = find_total_path(rxn_graph, specie_df, path)
    return len(total_path) > 0 and all([s in rxn_graph.source_species for s in path_reactants(rxn_graph, total_path)])


def estimate_total_rate(rxn_graph: RxnGraph, specie_df: pd.DataFrame, rxn: Reaction) -> float:
    """Estimate according to Lehmann's formalism the expeceted pathway rate assuming we join reaction with pathways of its reactants.
    Returns maximal rate and reactant to merge paths on"""
    # first calcualte rates 
    lehmann_rate = lambda r, rxn: specie_df.loc[rxn_graph.specie_collection.get_key(r), "dist"] * rxn.properties["rate"] / max([r.properties["consumption_rate"], r.properties["production_rate"]])
    rates = [lehmann_rate(r, rxn) for r in rxn.reactants]
    # get sorted by rate list of reactants and rates
    species_rates = sorted(list(zip(rxn.reactants, rates)), key=lambda x: x[1], reverse=True)
    # going over reactant opetion, if we find valid path we return it
    for reactant, rate in species_rates:
        total_path = _max_rate_path(rxn_graph, specie_df, reactant)
        total_path.append(rxn)
        if valid_pathway(rxn_graph, specie_df, total_path):
            return rate, reactant if reactant not in rxn_graph.source_species else None
    else:
        return - np.inf, None


def find_max_rate_pathways(rxn_graph: RxnGraph, prop_func=None):
    """Method to use Belmann-Ford algorithm together with the Lehmann reaction pathway formalism to find maximal rate pathways for all species in a reaction graph"""
    # initializing algorithm, setting a dataframe with 
    specie_df = {rxn_graph.specie_collection.get_key(s):
                     {"dist": - np.inf, "rxn": None, "reactant": None, "specie": s} for s in rxn_graph.species}
    source_species = [rxn_graph.specie_collection.get_key(s) for s in rxn_graph.source_species]
    for s in source_species:
        specie_df[s].update({"dist": 1})
    specie_df = pd.DataFrame.from_dict(specie_df, orient="index")
    # running belmann-ford algorithm to get max rate pathway
    # the algorithm goes over all possible pathways in network in an efficient way, by going over all edegs exactly the number of species
    for _ in range(len(specie_df) - 1):
        print("====", _, "====")
        for sp in rxn_graph.species:
            print(sp.identifier, valid_pathway(rxn_graph, specie_df, _max_rate_path(rxn_graph, specie_df, sp)))
            for rxn in max_rate_path(rxn_graph, specie_df, sp):
                print(rxn.pretty_string())
        for rxn in rxn_graph.reactions:
            # make sure both reactants are "visited"
            if any(specie_df.loc[rxn_graph.specie_collection.get_key(s), "dist"] < 0 for s in rxn.reactants):
                continue
            # for each reaction estimate
            for specie in rxn.products:
                if not specie in rxn_graph.source_species:
                    key = rxn_graph.specie_collection.get_key(specie)
                    # according to lehmann's formalism, we calculate what will be the total pathway rate if we join the new reaction to the reactants paths
                    pathway_rate, merge_specie = estimate_total_rate(rxn_graph, specie_df, rxn)                    
                    if specie_df.loc[key, "dist"] < pathway_rate:
                        specie_df.loc[key, "dist"] = pathway_rate
                        specie_df.loc[key, "rxn"] = rxn
                        specie_df.loc[key, "reactant"] = merge_specie
    return specie_df

def _max_rate_path(rxn_graph: RxnGraph, specie_df: pd.DataFrame, target: Specie):
    """Method to get the maximal rate pathway"""
    rxns = []
    current_sp = rxn_graph.specie_collection.get_key(target)
    while current_sp is not None:
        rxn = specie_df.loc[current_sp, "rxn"]
        if rxn is None:
            return []
        rxns.append(rxn)
        current_sp = specie_df.loc[current_sp, "reactant"]
        if current_sp is not None:
            current_sp = rxn_graph.specie_collection.get_key(current_sp) 
        else: 
            current_sp = None
    return rxns

def find_total_path(rxn_graph: RxnGraph, specie_df: pd.DataFrame, path: List[Reaction]):
    ajr = []
    covered_species = []
    total_path = path
    while True:
        for s in ajr:
            total_path += _max_rate_path(rxn_graph, specie_df, s)        
        # check for consumed species that are not source species in total path
        ajr = []
        for r in path_reactants(rxn_graph, path):
            if not r in rxn_graph.source_species and not r in covered_species:
                ajr.append(r)
                covered_species.append(r)
        if len(ajr) == 0:
            break
    if all(s in rxn_graph.source_species for s in path_reactants(rxn_graph, path)):
        return total_path
    else: 
        return []

def max_rate_path(rxn_graph: RxnGraph, specie_df: pd.DataFrame, target: Specie):
    # method to build a complete path from source species to target
    path = _max_rate_path(rxn_graph, specie_df, target)
    return find_total_path(rxn_graph, specie_df, path)

if __name__ == "__main__":
    rxn_graph = RxnGraph.from_file("../experimental_networks/curated/okafor2018_nh3+o2.rxn")
    sp = Specie("[O][N][O]", charge=0)
    if not rxn_graph.has_specie(sp):
        raise RuntimeError()
    sp = rxn_graph.add_specie(sp)
    # initialization
    rxn_graph = initialize_algorithm(rxn_graph)
    # finding k shortest paths
    # prop_func = lambda rxn: rxn.properties["rate"]
    # paths = yens_k_shortest_path(rxn_graph, sp, 3, prop_func=prop_func)
    df = find_max_rate_pathways(rxn_graph)
    rxns = max_rate_path(rxn_graph, df, sp)
    # for i, rxns in enumerate(paths):
    #     print("Path", i + 1, "Value =", sum(prop_func(rxn) for rxn in rxns))
    for rxn in rxns:
        print(rxn.pretty_string())
    