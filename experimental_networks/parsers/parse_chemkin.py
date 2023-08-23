from typing import List
import numpy as np
from utils import DummyReaction, build_rxn_tree, get_species_list, init_species_file, dummy_reactions_to_graph, DummySpecie, update_specie_energies, valid_rxn_graph

def read_species_from_list(species_list: List[str]):
    res = []
    for s in species_list:
        s = s.strip()
        if not "M" in s and len(s) > 0 and not s in ["AR", "HE"]:
            if s.endswith("("):
                s = s[:-1].strip()
            if s[:1].isdigit():
                n = int(s[0])
                for _ in range(n):
                    res.append(s[1:])
            else:
                res.append(s)
    return res

def read_thermo_data(path: str):
    """Read all thermodynamic data from chemkin format
    RETURNS: (dict) dictionary with symbol -> DummySpecie"""
    with open(path, "r") as f:
        res = {}
        block_counter = 1
        sym = None
        low_t_coeffs = np.zeros(7)
        high_t_coeffs = np.zeros(7)
        for line in f.readlines()[2:]:
            if block_counter == 1:
                sym = line.split()[0]
            elif block_counter == 5:
                sym = None
                low_t_coeffs = np.zeros(7)
                high_t_coeffs = np.zeros(7)
                block_counter = 0                
            else:
                line = line.split("!")
                numbers = []
                for i in range(5):
                    s = line[(15 * i):(15 * (i + 1))].strip()
                    if len(s) > 0:
                        numbers.append(float(s))
                if block_counter == 2:
                    for i in range(5):
                        high_t_coeffs[i] = numbers[i]
                elif block_counter == 3:
                    for i in range(5, 7):
                        high_t_coeffs[i] = numbers[i - 5]
                    for i in range(3):
                        low_t_coeffs[i] = numbers[i + 2]
                elif block_counter == 4:
                    for i in range(4):
                        low_t_coeffs[i + 3] = numbers[i]            
            
            block_counter += 1

def read_all_reactions(path: str):
    """Read all reactions from USC Mech v2"""
    with open(path, "r") as f:
        reactions_block = False
        reactions = set()
        for line in f.readlines():
            # skip comment lines
            if line.startswith("!"):
                continue
            # start reaction block
            if "REACTIONS" in line:
                reactions_block = True
                continue
            # end reaction block
            if "END" in line and reactions_block:
                reactions_block = False
                continue
            # parse reaction line
            # make sure its a reaction line
            if "=" in line:
                if "!" in line:
                    line = line.split("!")[0]
                    if not "=" in line:
                        continue
                # make all reaction notations uniform
                line = line.replace(" + ", "+")
                line = line.replace("+ ", "+")
                line = line.replace(" +", "+")
                splited = line.split()
                activation = float(splited[-1])
                beta = float(splited[-2])
                pre_exp = float(splited[-3])
                forward_rxn_props = {"Ea": activation, "beta": beta, "A": pre_exp}
                reverse_rxn_props = {"rEa": activation, "rbeta": beta, "rA": pre_exp}
                if "=>" in line and not "<=>" in line: # non-reversible reaction
                    ajr = line.split("=>")
                    rs = read_species_from_list(ajr[0].split("+"))
                    ps = read_species_from_list(ajr[-1].split()[0].split("+"))
                    # append forward reaction - non reversible reaction
                    reactions.add(DummyReaction(rs, ps, properties=forward_rxn_props))
                else: # "=" or "<=>" in reaction - reversible
                    if "<=>" in line:
                        ajr = line.split("<=>")
                    elif " = " in line:
                        ajr = line.split(" = ")
                    else:
                        ajr = line.split("=")
                    rs = read_species_from_list(ajr[0].split("+"))
                    ps = read_species_from_list(ajr[-1].split()[0].split("+"))
                    # append forward reaction - with kinetic data
                    reactions.add(DummyReaction(rs, ps, properties=forward_rxn_props))
                    # append backward reaction - without kinetic data
                    reactions.add(DummyReaction(ps, rs, properties=reverse_rxn_props))
    return list(reactions)

def init_specie_files(ffcm_file, uscm_file):
    # reading from ffcm
    ffcm_rxns = read_all_reactions("../raw/ffcm1.txt")
    ffcm_species = get_species_list(ffcm_rxns)
    init_species_file(ffcm_species, ffcm_file)
    uscm_rxns = read_all_reactions("../raw/uscm2.txt")
    uscm_species = get_species_list(uscm_rxns)
    init_species_file(uscm_species, uscm_file)

def old_make_networks(macro_iteration: int):
    print("MAKING NETWORKS WITH {} MAX ITERATION".format(macro_iteration))
    ffcm_rxns = read_all_reactions("../raw/ffcm1.txt")
    uscm_rxns = read_all_reactions("../raw/uscm2.txt")
    print("N REACTIONS FFCM v1:", len(ffcm_rxns))
    print("N REACTIONS UCS Mech v2:", len(uscm_rxns))
    # building reaction tree for CH3CH3 pyrolysis
    source_species = ["C2H6"]
    ffcm_ch3ch3 = build_rxn_tree(ffcm_rxns, source_species, macro_iteration)
    print("N REACTIONS CH3CH3-FFCM:", len(ffcm_ch3ch3))
    print("N SPECIES CH3CH3-FFCM:", len(get_species_list(ffcm_ch3ch3)))
    g = dummy_reactions_to_graph(ffcm_ch3ch3, "../raw/specie_files/ffcm.json", [DummySpecie(s, 0, {}) for s in source_species])
    update_specie_energies("../../db_files/ch3ch3.db", g).save("../curated/ffcm_ch3ch3_{}.rxn".format(macro_iteration))
    uscm_ch3ch3 = build_rxn_tree(uscm_rxns,source_species,  macro_iteration)
    print("N REACTIONS CH3CH3-UCSM:", len(uscm_ch3ch3))
    print("N SPECIES CH3CH3-UCSM:", len(get_species_list(uscm_ch3ch3)))
    g = dummy_reactions_to_graph(uscm_ch3ch3, "../raw/specie_files/uscm.json", [DummySpecie(s, 0, {}) for s in source_species])
    update_specie_energies("../../db_files/ch3ch3.db", g).save("../curated/uscm_ch3ch3_{}.rxn".format(macro_iteration))
    # building reaction tree for CH4+H2O
    source_species = ["CH4", "H2O"]
    ffcm_msr = build_rxn_tree(ffcm_rxns, source_species, macro_iteration)
    print("N REACTIONS CH4+H2O-FFCM:", len(ffcm_msr))
    print("N SPECIES CH4+H2O-FFCM:", len(get_species_list(ffcm_msr)))
    g = dummy_reactions_to_graph(ffcm_msr, "../raw/specie_files/ffcm.json", [DummySpecie(s, 0, {}) for s in source_species])
    update_specie_energies("../../db_files/ch4+h2o.db", g).save("../curated/ffcm_ch4+h2o_{}.rxn".format(macro_iteration))
    uscm_msr = build_rxn_tree(uscm_rxns, source_species, macro_iteration)
    print("N REACTIONS CH4+H2O-UCSM:", len(uscm_msr))
    print("N SPECIES CH4+H2O-UCSM:", len(get_species_list(uscm_msr)))
    g = dummy_reactions_to_graph(uscm_msr, "../raw/specie_files/uscm.json", [DummySpecie(s, 0, {}) for s in source_species])
    update_specie_energies("../../db_files/ch4+h2o.db", g).save("../curated/uscm_ch4+h2o_{}.rxn".format(macro_iteration))

def make_network(source: str, source_species: List[str]):
    print("reading {} reacitons...".format(source))
    rxns = read_all_reactions("../raw/{}.txt".format(source))
    print("total number of reactions:", len(rxns))
    print("building reaction graph with {} source species...".format(" + ".join(source_species)))
    g = build_rxn_tree(rxns, source_species, 8)
    g = dummy_reactions_to_graph(g, "../raw/specie_files/{}.json".format(source), source_species)
    print("n reactions:", g.get_n_reactions())
    print("n species:", g.get_n_species())
    print("validating graph...")
    if valid_rxn_graph(g):
        print("adding calculated specie energies when available...")
        g = update_specie_energies("../../db_files/joined.db", g)
        print("saving graph...")
        g.save("../curated/{}_{}.rxn".format(source, "+".join(source_species).lower()))
        print("done!")
    else:
        print("graph has invalid reactions! fix and try againd")

if __name__ == "__main__":
    make_network("konnov2009", ["NH3", "O2"])
    make_network("konnov2009", ["NH3", "O2", "H2"])
    make_network("konnov2009", ["NH3", "O2", "CH4"])
    make_network("okafor2018", ["NH3", "O2"])
    make_network("okafor2018", ["NH3", "O2", "H2"])
    make_network("okafor2018", ["NH3", "O2", "CH4"])