from typing import Optional
import openbabel as ob
import sqlite3
import os
import torinanet as tn

source_parent_dirs = {
    "exp": "/home/shaharpit/Personal/EREBenchmark/experimental_networks/curated",
    "cls": "/home/shaharpit/Personal/EREBenchmark/full_enum",
    "mvc": "/home/shaharpit/Personal/EREBenchmark/mvc_enum"
}

last_reaction_iteration = {
    # "ch3ch3": 2,
    "nh3+o2": 8,
    "nh3+o2+h2": 8,
    "ch4+h2o": 8,
}

mvc_energies = {
    # "ch3ch3": [-1.95, -1.8, -1.65, -1.30],
    "ch4+h2o": [-1.95, -1.8, -1.65, -1.30],
    "nh3+o2": [-0.8, -0.65, -0.5, -0.3],
    "nh3+o2+h2": [-0.8, -0.65, -0.5, -0.3]
}

def get_network(source: str, reaction: str, iteration: int=3, graph_type: str="energy_reduced_graph",  min_electron_energy: float=-8., timestamp: Optional[str]=None):
    if source == "cls":
        internal_path = os.path.join(reaction, str(iteration), graph_type + ".rxn")
    elif source == "mvc":
        internal_path = os.path.join("energy_{}".format(min_electron_energy), reaction, str(iteration), graph_type + ".rxn")
    elif "exp" in source:
        if "_" in source:
            prefix = source.split("_")[-1] + "_"
        else:
            prefix = ""
        return tn.core.RxnGraph.from_file(os.path.join(
            source_parent_dirs["exp"],
            prefix + "_".join([reaction, str(iteration)]) + ".rxn"
        ))
    else:
        raise ValueError("Uknown source {}".format(source))
    if timestamp is not None:
        parent = os.path.join(source_parent_dirs[source], "archive", timestamp)
    else:
        parent = source_parent_dirs[source]
    print("READING FROM", os.path.join(parent, internal_path))
    return tn.core.RxnGraph.from_file(os.path.join(parent, internal_path))

def reaction_energy(rxn: tn.core.Reaction):
    return sum([s.properties["energy"] for s in rxn.products]) - sum([s.properties["energy"] for s in rxn.reactants])

def smiles_to_obmol(smiles: str):
    conv = ob.OBConversion()
    conv.SetInFormat("smi")
    mol = ob.OBMol()
    conv.ReadString(mol, smiles)
    return mol

def get_atom_energies_from_db(db_path: str):
    with sqlite3.connect(db_path) as conn:
        smiles_energies = conn.execute("SELECT smiles, energy FROM species").fetchall()
        res = {}
        for smile, energy in smiles_energies:
            mol = smiles_to_obmol(smile)
            if mol.NumAtoms() == 1:
                atomic_num = mol.GetAtom(1).GetAtomicNum()
                res[atomic_num] = energy
        return res

def get_atom_energies_from_graph(rxn_graph: tn.core.RxnGraph):
    res = {}
    for specie in rxn_graph.species:
        if len(specie.ac_matrix) == 1:
            num = specie.ac_matrix.matrix[0][0]
            res[num] = specie.properties["energy"]
    return res

def atomization_energy(smiles: str, energy: float, atom_energies: dict):
    mol = smiles_to_obmol(smiles)
    atom_energy = 0
    for atom in ob.OBMolAtomIter(mol):
        num = atom.GetAtomicNum()
        atom_energy += atom_energies[num]
    return atom_energy - energy