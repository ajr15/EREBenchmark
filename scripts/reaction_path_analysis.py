from itertools import chain
from copy import copy
from collections import deque
import networkx as nx
from typing import List, Iterator
import numpy as np
import pandas as pd
from math import gcd
from torinanet.core import RxnGraph, Specie, Reaction
from torinanet.analyze.network_reduction.KineticReduction import MolRankReduction, assign_maximal_rate_constants

class Pathway:

    def __init__(self, rxn_graph: RxnGraph, reactions: List[Reaction], rate: float):
        self.rxn_graph = rxn_graph
        self.rate = rate
        self._reactions = self._build_reactions_df(reactions)
        self._species = self._build_specie_df(self._reactions)
        self.joining_species = [] # list to keep track on what species where used to join the reactions
    
    def _build_reactions_df(self, reactions: List[Reaction]) -> pd.DataFrame:
        df = pd.DataFrame(columns=["obj", "coeff"])
        for rxn in reactions:
            key = self.rxn_graph.reaction_collection.get_key(rxn)
            if key in df.index:
                df.loc[key, "coeff"] += 1
            else:
                df.loc[key] = [rxn, 1]
        return df
    
    def _add_species_from_list(self, df: pd.DataFrame, species: List[Specie], coeff: int):
        for s in species:
            key = self.rxn_graph.specie_collection.get_key(s)
            if key in df.index:
                df.loc[key, "coeff"] += coeff
            else:
                df.loc[key] = [s, coeff]

    def _build_specie_df(self, reactions_df: pd.DataFrame) -> pd.DataFrame:
        df = pd.DataFrame(columns=["obj", "coeff"])
        for rxn, coeff in reactions_df.values:
            self._add_species_from_list(df, rxn.reactants, - coeff)
            self._add_species_from_list(df, rxn.products, coeff)
        return df

    def specie_rate(self, specie: Specie) -> float:
        """Method to find the production / consumption rate of one specie in the pathway"""
        return self.rate * self.specie_coeff(specie)
    
    def has_specie(self, specie: Specie) -> bool:
        """Check if a specie is exists in this pathway"""
        key = self.rxn_graph.specie_collection.get_key(specie)
        return key in self._species.index

    def has_reaction(self, reaction: Reaction) -> bool:
        """Check if a specie is exists in this pathway"""
        key = self.rxn_graph.reaction_collection.get_key(reaction)
        return key in self._reactions.index

    def specie_coeff(self, specie: Specie) -> int:
        """Method to get coefficient of specie in pathway"""
        key = self.rxn_graph.specie_collection.get_key(specie)
        return self._species.loc[key, "coeff"]
    
    def reduce(self):
        """Method to reduce all coefficients to minimum"""
        devisor = gcd(self._reactions["coeff"])
        self._reactions["coeff"] /= devisor
        self._species["coeff"] /= devisor

    def multiply(self, factor: int):
        """Multiply all reaction coefficients by a factor"""
        self._reactions["coeff"] *= factor
        self._species["coeff"] *= factor

    def reactions(self) -> Iterator[Reaction]:
        """Get all reactions as a list, reaction will appear as its multiplicity"""
        for rxn, coeff in self._reacitons.values:
            for _ in range(coeff):
                yield rxn

    def unique_reacitons(self) -> List[Reaction]:
        """Get all reactions as a list each unique reaction appears once"""
        return self._reactions["obj"].tolist()
    
    def reactions_df(self) -> pd.DataFrame:
        """Return the reacions dataframe"""
        return self._reactions
    
    def __repr__(self) -> str:
        """Print nicely pathway information"""
        for rxn, c in self._reactions.values:
            print("{:<3}| {}".format(c, rxn.pretty_string()))



def join_pathways(p1: Pathway, p2: Pathway, specie: Specie) -> Pathway:
    """Join two pathways on a joining specie"""
    # make sure both pathways have the specie
    if not p1.has_specie(specie) or not p2.has_specie(specie):
        raise ValueError("Joining specie must be in both pathways")
    # calculate the rate of the joined pathway
    new_rate = (p1.rate * p2.rate) / max(specie.properties["consumption_rate"], specie.properties["production_rate"])
    # find the reaction coefficients according to lehmann
    p1.multiply(p2.specie_coeff(specie))
    p2.multiply(p1.specie_coeff(specie))
    # now make the new pathway
    new_p = Pathway(p1.rxn_graph, chain(p1.reactions(), p2.reactions()), new_rate)
    new_p.reduce()
    new_p.joining_species = p1.joining_species + [specie] + p2.joining_species
    return new_p


def substract_pathways(p1: Pathway, p2: Pathway) -> Pathway:
    """Method to remove reactions of pathway 2 from pathway now. returning a new pathway"""
    # collect all unique reactions and species to p1
    reactions = [rxn for rxn in p1.unique_reacitons() if not p2.has_reaction(rxn)]
    joining_species = [sp for sp in p1.joining_species if not sp in p2.joining_species]
    # now, rebuild the pathway using the reactions and joining species
    ajr = Pathway(p1.rxn_graph, [reactions[0]], reactions[0].properties["rate"])
    for rxn, sp in zip(reactions[1:], joining_species):
        ajr = join_pathways(ajr, Pathway(p1.rxn_graph, [rxn], rxn.properties["rate"]), sp)
    # the newly created pathway is the difference
    return ajr

def update_pathways(candidate_pathway: Pathway, specie_df: pd.DataFrame, target_species: List[str], parent_specie: str):
    """Method to recursively update reaction pathway of target species (and their children) given a candidate pathway"""
    for specie in target_species:
        sp = specie_df.loc[specie, "specie"]
        current_pathway = specie_df.loc[specie, "pathway"]
        # only if candidate pathway has larger rate, switch it
        if candidate_pathway.specie_rate(sp) > current_pathway.specie_rate(sp):
            # registrating the switched specie as a child of the parent
            specie_df.loc[parent_specie, "children"].append(specie)
            # updating pathway
            specie_df.loc[specie, "pathway"] = candidate_pathway
            for child_specie in copy(specie_df.loc[specie, "children"]):
                child_pathway = specie_df.loc[child_specie, "pathway"]
                # make sure that the child is indeed a child of the parent
                if not child_pathway.has_specie(sp):
                    # if not, remove the child from the children list and continue
                    specie_df.loc[specie, "children"].remove(child_specie)
                    continue
                # now, we can do the same procedure for the child
                # construct a new candidate pathway: this is = (current child pathway) - (old parent pathway) + (new parent pathway)
                p = join_pathways(candidate_pathway, substract_pathways(child_pathway, current_pathway), sp)
                # use now the same function in recursion 
                update_pathways(p, specie_df, [child_specie], specie)


def find_reaction_pathways(rxn_graph: RxnGraph, search_method: str="dfs"):
    """Method to implement modified dijkstra's algorithm for finding distance between specie to the source of the network (reactants)"""
    search_methods = {"bfs": lambda queue: queue.popleft(), "dfs": lambda queue: queue.popright()}
    # converting reaction graph to networkx
    G = rxn_graph.to_networkx_graph(use_internal_id=True)
    graph_objects = nx.get_node_attributes(G, "obj")
    source_species = [rxn_graph.specie_collection.get_key(s) for s in rxn_graph.source_species]
    # initializing 
    # dictionary for each specie its distance from source and making reaction
    specie_df = {rxn_graph.specie_collection.get_key(s):
                     {"pathway": None, "visited": False, "specie": s, "children": []} for s in rxn_graph.species}
    specie_df = pd.DataFrame.from_dict(specie_df, orient="index")
    # initializing the queue of species to explore
    species_queue = deque(source_species)
    touched_species = set(source_species)
    # running main loop
    while not all(specie_df["visited"].values): # runs until all species are visited
        # find the next specie to visit according to search method
        specie = search_methods[search_method](species_queue)
        # chaning flag to visited
        specie_df.loc[specie, "visited"] = True
        # go over all reactions with visited reactants
        for rxn in G.successors(specie):
            pred_species = [s for s in G.predecessors(rxn)]
            if all([specie_df.loc[s, "visited"] for s in pred_species]):
                # collecting neighbor species
                prod_species = [s for s in G.successors(rxn)]
                # adding the species to queue 
                species_queue += [s for s in prod_species if not s in touched_species]
                touched_species.union(prod_species)
                # constructing candidate pathway for species: old pathway + new reaction
                ajr = Pathway(rxn_graph, [graph_objects[rxn]], graph_objects[rxn].properties["rate"])
                if specie_df.loc[specie, "pathway"] is None:
                    pathway = ajr
                else:
                    pathway = join_pathways(specie_df.loc[specie, "pathway"], ajr, graph_objects[specie])
                # recursively updating pathways of species if rate is larger in this pathway
                update_pathways(pathway, specie_df, prod_species, specie)
    # dropping the "visited column"
    specie_df = specie_df.drop(columns=["visited"])
    return specie_df

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

if __name__ == "__main__":
    rxn_graph = RxnGraph.from_file("../experimental_networks/curated/okafor2018_nh3+o2.rxn")
    rxn_graph = initialize_algorithm(rxn_graph)
    df = find_reaction_pathways(rxn_graph, "bfs")