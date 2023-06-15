from typing import List
import pandas as pd
from utils import DummyReaction, get_species_list, init_species_file, build_rxn_tree, dummy_reactions_to_graph, update_specie_energies

def h2o_parser(row: dict) -> DummyReaction:
    reactants = []
    reactant_charges = {}
    for c in ["R1", "R2", "R3"]:
        if not pd.isna(row[c]):
            reactants.append(row[c])
            reactant_charges[row[c]] = row["C" + c]
    products = []
    product_charges = {}
    for c in ["P1", "P2", "P3", "P4"]:
        if not pd.isna(row[c]):
            products.append(row[c])
            product_charges[row[c]] = row["C" + c]
    return DummyReaction(reactants, products, None, None, None, reactant_charges, product_charges)

def nh3_parser(row: dict) -> DummyReaction:
    reactants = []
    for c in ["R1", "R2"]:
        if not pd.isna(row[c]):
            reactants.append(row[c])
    products = []
    for c in ["P1", "P2", "P3"]:
        if not pd.isna(row[c]):
            products.append(row[c])
    return DummyReaction(reactants, products, row)

def read_all_reactions(csv_path, row_dict_parser) -> List[DummyReaction]:
    df = pd.read_csv(csv_path)
    rxns = []
    for row in df.to_dict(orient="records"):
        rxns.append(row_dict_parser(row))
    return rxns

def init_specie_files(h2o_file, nh3_file):
    # reading from ffcm
    h2o_rxns = read_all_reactions("../raw/H2O.csv", h2o_parser)
    h2o_species = get_species_list(h2o_rxns)
    print(h2o_species)
    init_species_file(h2o_species, h2o_file)
    nh3_rxns = read_all_reactions("../raw/NH3+O2.csv", nh3_parser)
    nh3_species = get_species_list(nh3_rxns)
    init_species_file(nh3_species, nh3_file)

def read_csv(csv_file, parser, specie_file, source_species, n_iter):
    rxns = read_all_reactions(csv_file, parser)
    print("TOTAL N REACTIONS IN GRAPH:", len(rxns))
    print("TOTAL N SPECIES:", len(get_species_list(rxns)))
    tree = build_rxn_tree(rxns, source_species, n_iter)
    print("N REACTIONS IN TREE:", len(tree))
    print("N SPECIES IN TREE:", len(get_species_list(tree)))
    graph = dummy_reactions_to_graph(tree, specie_file, source_species)
    return graph

def make_networks(macro_iteration):
    # print("PARSING H2O")
    # read_csv("../raw/H2O.csv", h2o_parser, "../raw/specie_files/h2o.json", ["H2O"], 4, "../curated_networks/h2o.rxn")
    print("PARSING NH3+O2")
    graph = read_csv("../raw/NH3+O2.csv", nh3_parser, "../raw/specie_files/nh3+o2.json", ["NH3", "O2"], macro_iteration)
    update_specie_energies("../../db_files/nh3+o2.db", graph).save("../curated/nh3+o2_{}.rxn".format(macro_iteration))
    graph = read_csv("../raw/NH3+O2.csv", nh3_parser, "../raw/specie_files/nh3+o2.json", ["NH3", "O2", "H2"], macro_iteration)
    update_specie_energies("../../db_files/nh3+o2.db", graph).save("../curated/nh3+o2+h2_{}.rxn".format(macro_iteration))

if __name__ == "__main__":
    make_networks(1)
    make_networks(2)
    make_networks(3)
    make_networks(8)

