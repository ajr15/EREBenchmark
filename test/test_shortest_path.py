import torinanet as tn
import matplotlib.pyplot as plt
import numpy as np

def calc_reaction_energy(rxn: tn.core.Reaction):
    return sum([s.properties["energy"] for s in rxn.products]) - sum([s.properties["energy"] for s in rxn.reactants])

def calc_reaction_energy_weight(rxn: tn.core.Reaction, temp: float):
    R = 8.314 * 2.294e17 / 6.022e23 # in Ha/(mol K)
    return max(calc_reaction_energy(rxn), 0)

def calc_reaction_rate_weight(rxn: tn.core.Reaction, temp: float):
    R = 8.314 * 2.294e17 / 6.022e23 # in Ha/(mol K)
    return np.exp(max(calc_reaction_energy(rxn), 0) / (R * temp))

def analyze_graph(g: tn.core.RxnGraph, prop_func, title):
    print("Running", title)
    analyzer = tn.analyze.algorithms.ShortestPathAnalyzer(g, prop_func=prop_func)
    energies = []
    for sp in g.species:
        rxns = analyzer.get_path_to_source(sp)
        energy = sum([calc_reaction_energy(rxn) for rxn in rxns])
        energies.append(energy)
    plt.title(title)
    plt.hist(energies)
    plt.savefig("./{}.png".format(title))
    plt.close()

if __name__ == "__main__":
    g = tn.core.RxnGraph.from_file("../full_enum/nh3+o2/3/energy_reduced_graph.rxn")
    analyze_graph(g, prop_func=lambda rxn: calc_reaction_energy(rxn), title="ClassicWeight")
    analyze_graph(g, prop_func=lambda rxn: calc_reaction_energy_weight(rxn, 500), title="EnergyWeight")
    analyze_graph(g, prop_func=lambda rxn: calc_reaction_rate_weight(rxn, 500), title="KineticWeight")
