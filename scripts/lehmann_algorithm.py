# script to check out the Lehmann path finding algorithm for reaction networks
import pandas as pd
from itertools import product
from typing import Dict, List
import networkx as nx
import numpy as np
from math import gcd
from dataclasses import dataclass
from torinanet.core import Reaction, Specie, RxnGraph
from torinanet.core.RxnGraph.HashedCollection import IndepCuckooReactionCollection, CuckooSpecieCollection
from torinanet.analyze.network_reduction.KineticReduction import MolRankReduction, assign_maximal_rate_constants

@dataclass
class Pathway:

    rate: float
    reactions: IndepCuckooReactionCollection
    multiplicities: Dict[str, int]

    def specie_coef(self, specie: Specie) -> int:
        ajr = 0
        for key, rxn in self.reactions.items():
            if specie in rxn.reactants:
                ajr += - self.multiplicities[key]
            if specie in rxn.products:
                ajr += self.multiplicities[key]
        return ajr
    
    def reduce_multiplicities(self):
        """Method to reduce pathway to simpler form by making the reaction multiplicities as low as possible"""
        d = gcd(*self.multiplicities.values())
        for k, v in self.multiplicities.items():
            self.multiplicities[k] = int(v / d)

    def has_reaction(self, rxn: Reaction):
        return self.reactions.has(rxn)
    
    def to_str(self):
        """Method to convert pathway to string in a pretty format"""
        return "\n".join([rxn.pretty_string() for rxn in self.reactions.objects()])
    
    def get_specie_coeffs(self):
        """Method to get all specie coeffs in the pathway"""
        species = CuckooSpecieCollection()
        coeffs = {}
        for rxn in self.reactions.objects():
            for s in rxn.reactants + rxn.products:
                if not species.has(s):
                    c = self.specie_coef(s)
                    if c != 0:
                        species.add(s)
                        coeffs[species.get_key(s)] = c
        return species, coeffs
    
    def to_reaction(self) -> Reaction:
        """Convert pathway to an effective reaction"""
        species, coeffs = self.get_specie_coeffs()
        reactants = []
        products = []
        for sp in species.objects():
            if coeffs[species.get_key(sp)] < 0:
                reactants.append(sp)
            else:
                products.append(sp)
        return Reaction(reactants, products, {"rate": self.rate})

        

def join_pathways(p1: Pathway, p2: Pathway, specie: Specie) -> Pathway:
    """Join two pathways insuring that the net production of a given branching point specie is 0"""
    # the new pathway will have the union over the reactions
    new_pathway_rxns = IndepCuckooReactionCollection()
    new_pathway_rxns.union(p1.reactions)
    new_pathway_rxns.union(p2.reactions)
    # now decide on multiplicities (according to original paper formula, new mult = | p1 specie mult | * (p1 reaction mult) + | p2 specie mult | * (p2 reaction mult))
    p1_specie_mult = abs(p1.specie_coef(specie))
    p2_specie_mult = abs(p2.specie_coef(specie))
    new_mults = dict()
    for rxn in new_pathway_rxns.objects():
        mult = 0
        if p1.has_reaction(rxn):
            mult += p1_specie_mult * p1.multiplicities[p1.reactions.get_key(rxn)]
        if p2.has_reaction(rxn):
            mult += p2_specie_mult * p2.multiplicities[p2.reactions.get_key(rxn)]
        new_mults[new_pathway_rxns.get_key(rxn)] = int(mult)
    # at last, we need to calculate the rate of the new pathway
    new_rate = (p1.rate * p2.rate) / max(specie.properties["consumption_rate"], specie.properties["production_rate"])
    # now return the new pathway (reduced)
    p = Pathway(new_rate, new_pathway_rxns, new_mults)
    p.reduce_multiplicities()
    return p


def get_connected_paths(g: nx.DiGraph, specie: str, pathways: List[Pathway]):
    """Method to get the pathways connected to a specie"""
    rxns = [g.nodes[n]["obj"] for n in g.predecessors(specie)] + [g.nodes[n]["obj"] for n in g.successors(specie)]
    sp = g.nodes[specie]["obj"]
    consuming = []
    producing = []
    for rxn in rxns:
        for pathway in pathways:
            if pathway.has_reaction(rxn):
                coeff = pathway.specie_coef(sp)
                if coeff > 0:
                    producing.append(pathway)
                elif coeff < 0:
                    consuming.append(pathway)
    return consuming, producing


def find_branching_point(g: nx.DiGraph, tmin: float):
    """Find branching point specie (specie with shortest lifetime). make sure lifetime is greater than tmin"""
    sp = None
    lt = np.inf
    for node in g.nodes:
        s = g.nodes[node]["obj"]
        if type(s) == Specie:
            slt = s.properties["lifetime"]
            if slt < lt and slt > tmin:
                lt = slt
                sp = node
    return sp, lt

def lehmanns_iteration(g: nx.DiGraph, pathways: List[Pathway], branching_point: str, min_rate: float) -> List[Pathway]:
    """Single iteration in the Lehmann's algorithm. Returns a list of all pathways in the system"""
    # step 1: select branching point - this is the specie node with shortest lifetime
    branching_point_specie = g.nodes[branching_point]["obj"]
    # step 2: merge all consuming and producing pathways
    consuming_pathways, producing_pathways  = get_connected_paths(g, branching_point, pathways)
    # now join all paths together, filtering out all unimportant paths
    for consuming, producing in product(consuming_pathways, producing_pathways):
        new_path = join_pathways(consuming, producing, branching_point_specie)
        if new_path.rate > min_rate:
            # remove old pathways
            if consuming in pathways:
                pathways.remove(consuming)
            if producing in pathways:
                pathways.remove(producing)
            # add new pathway
            if new_path.rate is None or pd.isna(new_path.rate):
                raise RuntimeError("YOU'RE AND IDIOT")
            pathways.append(new_path)
    return pathways

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
    # now building all pathways
    pathways = []
    for rxn in rxn_graph.reactions:
        key = rxn_graph.reaction_collection.get_key(rxn)
        rxns = IndepCuckooReactionCollection()
        rxns.add(rxn)
        new_key = rxns.get_key(rxn)
        pathways.append(Pathway(g.nodes[key]["rate"], rxns, multiplicities={new_key: 1}))
    return g, pathways

def main(rxn_graph: RxnGraph, tmin: float, min_rate: float, niterations: int) -> List[Pathway]:
    """Main run for the Lehmann's algorithm. takes reaciton graph and returns all significant pathways"""
    print("INITIALIZING")
    g, pathways = initialize_algorithm(rxn_graph)
    for i in range(niterations):
        print("ITERATION", i + 1, "OUT OF", niterations)
        branching_point, tmin = find_branching_point(g, tmin)
        print("branching point:", g.nodes[branching_point]["obj"].identifier)
        pathways = lehmanns_iteration(g, pathways, branching_point, min_rate)
        print("# pathways:", len(pathways))
    print("ALL DONE")
    return pathways


def reactions_to_df(reactions: List[Reaction]) -> pd.DataFrame:
    max_r = max([len(x.reactants) for x in reactions])
    max_p = max([len(x.products) for x in reactions])
    data = []
    for rxn in reactions:
        rs = [x.identifier for x in rxn.reactants]
        rs += [None for _ in range(max_r - len(rs))]
        ps = [x.identifier for x in rxn.products]
        ps += [None for _ in range(max_p - len(ps))]
        data.append(rs + ps + [rxn.properties["rate"]])
    return pd.DataFrame(data=data, columns=["r" + str(i) for i in range(max_r)] + ["p" + str(i) for i in range(max_p)] + ["rate"])

if __name__ == "__main__":
    rxn_graph = RxnGraph.from_file("../experimental_networks/curated/miller_nh3+o2.rxn")
    pathways = main(rxn_graph, 0, 1e9, 3)
    df = reactions_to_df([p.to_reaction() for p in pathways])
    df.to_csv("../analysis_results/pathways.csv")