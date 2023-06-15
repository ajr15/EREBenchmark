from typing import List, Optional
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from dataclasses import dataclass
from copy import copy
import os
import openbabel as ob
import json
import sys
sys.path.append("C:\\Users\\Shachar\\OneDrive\\Documents\\repos\\TorinaNet")
sys.path.append("C:\\Users\\Shachar\\OneDrive\\Documents\\repos\\TorinaX")
import torinanet as tn

@dataclass
class DummySpecie:

    symbol: str
    charge: int
    properties: dict

    def __str__(self):
        return self.symbol

    def __hash__(self) -> int:
        return hash(self.symbol)

    def to_specie(self, symbol_dict: dict) -> tn.core.Specie:
        props = copy(self.properties)
        props["symbol"] = self.symbol
        return tn.core.Specie(symbol_dict[self.symbol], charge=self.charge, properties=props)


@dataclass
class DummyReaction:

    reactants: List[DummySpecie]
    products: List[DummySpecie]
    properties: dict

    def __hash__(self):
        return hash("{}->{}".format("+".join(self.reactants), "+".join(self.products)))


def build_rxn_tree(reactions: List[DummyReaction], source_species: List[str], n_iters: int):
    """Build the total reaction graph with the desired source species"""
    allowed_species = set(source_species)
    res = set()
    for _ in range(n_iters):
        new_species = set()
        for reaction in reactions:
            if all([r in allowed_species for r in reaction.reactants]):
                res.add(reaction)
                for r in reaction.products:
                    new_species.add(r)
        allowed_species = allowed_species.union(new_species)
    return list(res)

def get_species_list(rxn_list: List[DummyReaction]):
    species = set()
    for r in rxn_list:
        for s in r.reactants:
            species.add(s)
        for s in r.products:
            species.add(s)
    return list(species)

def init_species_file(specie_list: List[str], specie_file: str):
    with open(specie_file, "w") as f:
        d = {s: "" for s in specie_list}
        json.dump(d, f, sort_keys=True, indent=1)


def dummy_reactions_to_graph(rxn_list: List[DummyReaction], specie_file: str, source_species: List[DummySpecie]) -> tn.core.RxnGraph:
    source_species = [DummySpecie(s, 0, {}) if type(s) is str else s for s in source_species]
    with open(specie_file, "r") as f:
        specie_d = json.load(f)
        rxn_graph = tn.core.RxnGraph()
        source_species = [rxn_graph.add_specie(s.to_specie(specie_d)) for s in source_species]
        rxn_graph.set_source_species(source_species, force=True)
        for rxn in rxn_list:
            reactants = [rxn_graph.add_specie(DummySpecie(s, 0, {}).to_specie(specie_d)) for s in rxn.reactants]
            products = [rxn_graph.add_specie(DummySpecie(s, 0, {}).to_specie(specie_d)) for s in rxn.products]
            reaction = tn.core.Reaction(reactants, 
                                        products, 
                                        properties=rxn.properties)
            rxn_graph.add_reaction(reaction)
        return rxn_graph
    
def update_specie_energies(db_path, rxn_graph) -> tn.core.RxnGraph:
    """Method to update the species energies in the reaction graph from the computation"""
    # # earasing charge data from experimental graph - better match with enumeration
    # for s in rxn_graph.species:
    #     s.charge = None
    # rxn_graph._has_charge = False
    # now, reading energies
    engine = create_engine("sqlite:///{}".format(os.path.abspath(db_path)))
    session = sessionmaker(bind=engine)()
    smiles_energies = session.execute("SELECT smiles, energy FROM species WHERE good_geometry == 1 AND successful == 1")
    for smiles, energy in smiles_energies:
        s = tn.core.Specie(smiles, charge=0)
        if rxn_graph.has_specie(s):
            # workaround to get the specie from the graph and update it. this method returns the "correct" specie.
            # using other methods such as specie_collection.get require proper reading with AC matrix
            ajr = rxn_graph.add_specie(s)
            ajr.properties["energy"] = float(energy)
    for s in rxn_graph.species:
        if not "energy" in s.properties:
            print(s.identifier, "IS A MOTHER F!*&#@ !")
    return rxn_graph
