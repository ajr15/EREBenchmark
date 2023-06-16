import torinanet as tn
from matplotlib_venn import venn3_unweighted
import pandas as pd
import matplotlib.pyplot as plt
import os
from typing import List
import utils

def get_pair(reaction, source):
    cls_parent = "/home/shaharpit/Personal/EREBenchmark/full_enum"
    mvc_parent = "/home/shaharpit/Personal/EREBenchmark/mvc_enum"
    exp_parent = "/home/shaharpit/Personal/EREBenchmark/experimental_networks/curated"
    if reaction in ["ch3ch3", "ch4+h2o"]:
        exp_net = tn.core.RxnGraph.from_file(os.path.join(exp_parent, "{}_{}_3.rxn".format(source, reaction)))
        mvc_net = tn.core.RxnGraph.from_file(os.path.join(mvc_parent, reaction, "3", "energy_reduced_graph.rxn"))
        cls_net = tn.core.RxnGraph.from_file(os.path.join(cls_parent, reaction, "3", "energy_reduced_graph.rxn"))
    else:
        exp_net = tn.core.RxnGraph.from_file(os.path.join(exp_parent, "nh3+o2_2.rxn"))
        mvc_net = tn.core.RxnGraph.from_file(os.path.join(mvc_parent, reaction, "2", "energy_reduced_graph.rxn"))
        cls_net = tn.core.RxnGraph.from_file(os.path.join(cls_parent, reaction, "2", "energy_reduced_graph.rxn"))
    # adding data to exp graph, making it more comparable
    for s in exp_net.species:
        # wiping charge data from exp networks
        s.charge = None
        # adding energy data from calculated graphs
        s.properties["energy"] = cls_net.specie_collection.get(s).properties["energy"]
    return exp_net, mvc_net, cls_net

def get_networks(reaction, source, energy):
    
    return ()


def make_venn_diagram(exp_net: tn.core.RxnGraph, cls_net: tn.core.RxnGraph, mvc_net: tn.core.RxnGraph, focus="species"):
    if focus == "species":
        AB = list(exp_net.specie_collection.intersect(cls_net.specie_collection))
        AC = list(exp_net.specie_collection.intersect(mvc_net.specie_collection))
        BC = list(mvc_net.specie_collection.intersect(cls_net.specie_collection))
    else:
        AB = list(exp_net.reaction_collection.intersect(cls_net.reaction_collection))
        AC = list(exp_net.reaction_collection.intersect(mvc_net.reaction_collection))
        BC = list(mvc_net.reaction_collection.intersect(cls_net.reaction_collection))
    ABC = 0
    for x in AB:
        if x in AC and x in BC:
            ABC += 1
    ABc = len(AB) - ABC
    AbC = len(AC) - ABC
    aBC = len(BC) - ABC
    if focus == "species":
        Abc = exp_net.get_n_species() - ABc - AbC - ABC
        aBc = cls_net.get_n_species() - ABc - aBC - ABC
        abC = mvc_net.get_n_species() - aBC - AbC - ABC
    else:
        Abc = exp_net.get_n_reactions() - ABc - AbC - ABC
        aBc = cls_net.get_n_reactions() - ABc - aBC - ABC
        abC = mvc_net.get_n_reactions() - aBC - AbC - ABC
    venn3_unweighted(subsets=[Abc, aBc, ABc, abC, AbC, aBC, ABC], set_labels=["exp_net", "cls_net", "mvc_net"])

def get_exp_graph_disections(exp_net: tn.core.RxnGraph, cls_net: tn.core.RxnGraph, mvc_net: tn.core.RxnGraph, focus="species", use_mvc: bool=True):
    if focus == "species":
        exp = exp_net.specie_collection
        cls_ = cls_net.specie_collection
        mvc = mvc_net.specie_collection
        print_func = lambda obj: print(obj.identifier)
    else:
        exp = exp_net.reaction_collection
        cls_ = cls_net.reaction_collection
        mvc = mvc_net.reaction_collection
        print_func = lambda obj: print(obj.pretty_string(), ",", sum([r.properties["energy"] for r in obj.products]) - sum([r.properties["energy"] for r in obj.reactants]))

    all_shared = 0
    unique_exp = 0
    exp_and_cls = 0
    all_shared_rxns = []
    exp_and_cls_rxns = []
    print(focus)
    for obj in exp.objects():
        if cls_.has(obj):
            if use_mvc and mvc.has(obj):
                all_shared += 1
                all_shared_rxns.append(obj)
            else:
                exp_and_cls += 1
                exp_and_cls_rxns.append(obj)
        else:
            # print_func(obj)
            unique_exp += 1
    # print("ALL SHARED")
    # for obj in all_shared_rxns:
    #     if focus == "species":
    #         print(obj.identifier)
    #     else:
    #         print(obj.pretty_string())
    # print("EXP AND CLS")
    # for obj in exp_and_cls_rxns:
    #     if focus == "species":
    #         print(obj.identifier)
    #     else:
    #         print(obj.pretty_string())
    return {"all_shared": all_shared, 
            "exp_and_cls": exp_and_cls, 
            "unique_exp": unique_exp, 
            "unique_cls": len(cls_) - all_shared - exp_and_cls,
            "unique_mvc": len(mvc) - all_shared}

def missing_reactions_categories(exp_net: tn.core.RxnGraph, cls_net: tn.core.RxnGraph):
    counter = {
        "BAD_FORMAT": 0,
        "ENERGY": 0,
        "PRODUCTS": 0,
        "REACTANTS": 0,
        "OTHER": 0
    }
    for rxn in exp_net.reaction_collection.objects():
        if not cls_net.has_reaction(rxn):
            # if any of the reactants appear also in products - BAD FORMAT
            if any([s in rxn.products for s in rxn.reactants]):
                counter["BAD_FORMAT"] += 1
            # if reaction energy is greater than the barrier - ENERGY
            elif sum([s.properties["energy"] for s in rxn.products]) - sum([s.properties["energy"] for s in rxn.reactants]) > 0.15:
                counter["ENERGY"] += 1
            # if products > 2 - PRODUCTS
            elif len(rxn.products) > 2:
                counter["PRODUCTS"] += 1
            # if reactants > 2 - REACTANTS
            elif len(rxn.reactants) > 2:
                counter["REACTANTS"] += 1
            else:
                counter["OTHER"] += 1
                print(rxn.pretty_string(), sum([s.properties["energy"] for s in rxn.products]) - sum([s.properties["energy"] for s in rxn.reactants]))
    print(counter)
    return counter

def get_summary_df(use_mvc=True):
    pairs = [
        # ("ch4+h2o", "ffcm"),
        # ("ch4+h2o", "uscm"),
        # ("ch3ch3", "ffcm"),
        # ("ch3ch3", "uscm"),
        ("nh3+o2", None),
        ("nh3+o2+h2", None),
    ]
    res = pd.DataFrame()
    for reaction, source in pairs:
        print(reaction, source)
        if source is None:
            source_str = "exp"
        else:
            source_str = "exp_" + source
        cls_net = utils.get_network(source="cls", reaction=reaction, iteration=utils.last_reaction_iteration[reaction], graph_type="leaf_energy_reduced_graph")
        print("CLS SPECIES", cls_net.get_n_species())
        print("CLS RXNS", cls_net.get_n_reactions())
        exp_net = utils.get_network(source=source_str, reaction=reaction, iteration=utils.last_reaction_iteration[reaction])
        # wiping charge data from exp networks
        for s in exp_net.species:
            s.charge = None
            # adding energy data from calculated graphs
            # s.properties["energy"] = cls_net.specie_collection.get(s).properties["energy"]
        if use_mvc:
            for energy in utils.mvc_energies[reaction]:
                print(energy)
                mvc_net = utils.get_network(source="mvc", reaction=reaction, iteration=utils.last_reaction_iteration[reaction], min_electron_energy=energy)
                d = get_exp_graph_disections(exp_net, cls_net, mvc_net, focus="species", use_mvc=use_mvc)
                d["focus"] = "species"
                d["reaction"] = reaction
                d["source"] = source
                d["energy"] = energy
                res = res.append(d, ignore_index=True)
                d = get_exp_graph_disections(exp_net, cls_net, mvc_net, focus="reactions", use_mvc=use_mvc)
                d["focus"] = "reactions"
                d["reaction"] = reaction
                d["source"] = source
                d["energy"] = energy
                res = res.append(d, ignore_index=True)
                res.to_csv("../analysis_results/compare_to_experiments.csv")
        else:
            mvc_net = tn.core.RxnGraph()
            d = get_exp_graph_disections(exp_net, cls_net, mvc_net, focus="species", use_mvc=use_mvc)
            d["focus"] = "species"
            d["reaction"] = reaction
            d["source"] = source
            res = res.append(d, ignore_index=True)
            d = get_exp_graph_disections(exp_net, cls_net, mvc_net, focus="reactions", use_mvc=use_mvc)
            d["focus"] = "reactions"
            d["reaction"] = reaction
            d["source"] = source
            res = res.append(d, ignore_index=True)
            res.to_csv("../analysis_results/compare_to_experiments.csv")
            missing_reactions_categories(exp_net, cls_net)

if __name__ == "__main__":
    # exp_net, mvc_net, cls_net = get_pair("ch4+h2o", "ffcm")
    # mvc_net = tn.core.RxnGraph.from_file("../test/multi_effect.rxn")
    # print("SPECIES")
    # d = get_exp_graph_disections(exp_net, cls_net, mvc_net, focus="species")
    # for k, v in d.items():
    #     print(k, v)
    # print("REACTIONS")
    # d = get_exp_graph_disections(exp_net, cls_net, mvc_net, focus="reactions")
    # for k, v in d.items():
    #     print(k, v)
    get_summary_df(use_mvc=False)
