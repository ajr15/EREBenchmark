import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
import os
import argparse
import enumerator
import torinanet as tn
import torinax as tx

kw_dict = {
    "ch4+h2o": {
              "ac_filters": [
                    tn.iterate.ac_matrix_filters.MaxAtomsOfElement({6: 2, 8: 2}), 
                    tn.iterate.ac_matrix_filters.MaxBondsPerAtom(), 
                    tn.iterate.ac_matrix_filters.MaxComponents(2),
                    tn.iterate.ac_matrix_filters.MaxRingNumber(0)
                ],
              "reaction_energy_th": 0.15,
              "sp_energy_th": 0.25,
        },
    "nh3+o2": {
              "ac_filters": [
                    tn.iterate.ac_matrix_filters.MaxAtomsOfElement({7: 2, 8: 2}), 
                    tn.iterate.ac_matrix_filters.MaxBondsPerAtom(), 
                    tn.iterate.ac_matrix_filters.MaxComponents(2),
                    tn.iterate.ac_matrix_filters.MaxRingNumber(0),
                ],
             "reaction_energy_th": 0.15,
              "sp_energy_th": 0.25,
        },
    "nh3+o2+h2": {
              "ac_filters": [
                    tn.iterate.ac_matrix_filters.MaxAtomsOfElement({7: 2, 8: 2}), 
                    tn.iterate.ac_matrix_filters.MaxBondsPerAtom(), 
                    tn.iterate.ac_matrix_filters.MaxComponents(2),
                    tn.iterate.ac_matrix_filters.MaxRingNumber(0),
                ],
             "reaction_energy_th": 0.15,
              "sp_energy_th": 0.25,
        },
    "nh3+o2+ch4": {
              "ac_filters": [
                    tn.iterate.ac_matrix_filters.MaxAtomsOfElement({6: 2, 7: 2, 8: 2}), 
                    tn.iterate.ac_matrix_filters.MaxBondsPerAtom(), 
                    tn.iterate.ac_matrix_filters.MaxComponents(2),
                    tn.iterate.ac_matrix_filters.MaxRingNumber(0),
                    tn.iterate.ac_matrix_filters.MaxHeavyAtoms(3),
                ],
             "reaction_energy_th": 0.15,
              "sp_energy_th": 0.25
        },
    "ch3ch3": {
              "ac_filters": [
                    tn.iterate.ac_matrix_filters.MaxAtomsOfElement({6: 4}), 
                    tn.iterate.ac_matrix_filters.MaxBondsPerAtom(), 
                    tn.iterate.ac_matrix_filters.MaxComponents(2),
                    tn.iterate.ac_matrix_filters.MaxRingNumber(1),
                    tn.iterate.ac_matrix_filters.BondOrderRingFilter({2: 3, 3: 4})
                ],
              "reaction_energy_th": 0.15,
              "sp_energy_th": 0.25,
        },
}


if __name__ == "__main__":
    # argument parser for different netnworks
    parser = argparse.ArgumentParser()
    parser.add_argument("reaction", type=str, help="reaction to enumerate. options: ch4+h2o, nh3+o2 and ch3ch")
    args = parser.parse_args()
    # init
    os.chdir("./{}".format(args.reaction))
    slurm_client = tx.clients.SlurmClient(4, "8GB", args.reaction)
    dask_client = None
    rxn_graph = tn.core.RxnGraph()
    # setting appropriate source reactants
    if args.reaction == "ch4+h2o":
        reactants = [tn.core.Specie("O"), tn.core.Specie("C")]
    elif args.reaction == "nh3+o2":
        reactants = [tn.core.Specie("O=O"), tn.core.Specie("N")]
    elif args.reaction == "ch3ch3":
        reactants = [tn.core.Specie("CC")]
    elif args.reaction == "nh3+o2+h2":
        reactants = [tn.core.Specie("O=O"), tn.core.Specie("N"), tn.core.Specie("[H][H]")]
    elif args.reaction == "nh3+o2+ch4":
        reactants = [tn.core.Specie("O=O"), tn.core.Specie("N"), tn.core.Specie("C")]
    else:
        raise ValueError("Unknown reaction {}".format(args.reaction))
    # setting source species in graph
    rxn_graph.set_source_species(reactants, force=True)
    # getting threshold
    print("RUNNING FOR {}".format(args.reaction))
    enumer = enumerator.SimpleEnumerator(rxn_graph, 8, ".", slurm_client,
                                         reflect=True,
                                         molrank_temperature=600,
                                         **kw_dict[args.reaction])
    enumer.enumerate()
