import torinanet as tn
import networkx as nx
import pandas as pd

def calc_reaction_energy(rxn: tn.core.Reaction):
    return sum([s.properties["energy"] for s in rxn.products]) - sum([s.properties["energy"] for s in rxn.reactants])

def reduce_graph(g: tn.core.RxnGraph, name: str):
    print("**** {} ****".format(name))
    print("BEFORE")
    print("N species:", g.get_n_species())
    print("N reactions:", g.get_n_reactions())
    G = g.to_networkx_graph(use_internal_id=True)
    # collecting reactions and species with < 5 creation reactions
    reactions = []
    species = []
    for sp in g.species:
        sid = g.specie_collection.get_key(sp)
        generating_rxns = list(G.predecessors(sid))
        if len(generating_rxns) < 5:
            print(sp.identifier, len(generating_rxns))
            reactions += [G.node[node]["obj"] for node in generating_rxns]
            species.append(sp)
    print("going over {} leaf reactions".format(len(reactions)))
    # removing reactions with high energy
    ajr = []
    for rxn in reactions:
        if calc_reaction_energy(rxn) > 0:
            ajr.append(rxn)
    print("removing {} reactions".format(len(ajr)))
    if len(ajr) > 0:
        g = g.remove_reactions(ajr)
    print("AFTER REACTION REDUCTION")
    print("N species:", g.get_n_species())
    print("N reactions:", g.get_n_reactions())
    analyzer = tn.analyze.algorithms.ShortestPathAnalyzer(g, prop_func=lambda r: max(0, calc_reaction_energy(r)))
    ajr = []
    for sp in species:
        if g.has_specie(sp):
            min_energy = sum([calc_reaction_energy(r) for r in analyzer.get_path_to_source(sp)])
            if min_energy > 0.05:
                ajr.append(sp)
    print("removing {} species".format(len(ajr)))
    if len(ajr) > 0:
        g = g.remove_species(ajr)
    print("AFTER SP REDUCTION")
    print("N species:", g.get_n_species())
    print("N reactions:", g.get_n_reactions())
    g.save("{}.rxn".format(name))

def get_exp_graph_disections(exp_net: tn.core.RxnGraph, cls_net: tn.core.RxnGraph, focus="species"):
    if focus == "species":
        exp = exp_net.specie_collection
        cls_ = cls_net.specie_collection
    else:
        exp = exp_net.reaction_collection
        cls_ = cls_net.reaction_collection
    all_shared = 0
    unique_exp = 0
    for obj in exp.objects():
        if cls_.has(obj):
            all_shared += 1
        else:
            unique_exp += 1
    return {"all_shared": all_shared, 
            "unique_exp": unique_exp, 
            "unique_cls": len(cls_) - all_shared}

def make_graphs():
    g = tn.core.RxnGraph.from_file("../full_enum/ch4+h2o/3/energy_reduced_graph.rxn")
    reduce_graph(g, "CH4+H2O")
    g = tn.core.RxnGraph.from_file("../full_enum/nh3+o2/3/energy_reduced_graph.rxn")
    reduce_graph(g, "NH3+O2")

def compare_graphs():
    pairs = [
        ("NH3+O2", "", "./NH3+O2.rxn", "../experimental_networks/curated/nh3+o2_3.rxn"),
        ("CH4+H2O", "ffcm", "./CH4+H2O.rxn", "../experimental_networks/curated/ffcm_ch4+h2o_3.rxn"),
        ("CH4+H2O", "uscm", "./CH4+H2O.rxn", "../experimental_networks/curated/uscm_ch4+h2o_3.rxn")
    ]
    res = pd.DataFrame()
    for name, source, cls_path, exp_path in pairs:
        print(name, source)
        cls_net = tn.core.RxnGraph.from_file(cls_path)
        print("CLS - N species:", cls_net.get_n_species())
        print("CLS - N reactions:", cls_net.get_n_reactions())
        exp_net = tn.core.RxnGraph.from_file(exp_path)
        for sp in exp_net.species:
            sp.charge = None
            print(sp.properties["energy"])
        print("EXP - N species:", exp_net.get_n_species())
        print("EXP - N reactions:", exp_net.get_n_reactions())
        d = get_exp_graph_disections(exp_net, cls_net, focus="species")
        d["focus"] = "species"
        d["source"] = source
        d["reaction"] = name
        res = res.append(d, ignore_index=True)
        d = get_exp_graph_disections(exp_net, cls_net, focus="reactions")
        d["focus"] = "reactions"
        d["source"] = source
        d["reaction"] = name
        res = res.append(d, ignore_index=True)
    res.to_csv("comparision.csv")
    return res

if __name__ == "__main__":
    make_graphs()
    compare_graphs()