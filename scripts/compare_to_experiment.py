import torinanet as tn
import pandas as pd
import importlib.util
import sys
import utils

def get_exp_graph_disections(exp_net: tn.core.RxnGraph, cls_net: tn.core.RxnGraph, focus="species"):
    if focus == "species":
        exp = exp_net.specie_collection
        cls_ = cls_net.specie_collection
    else:
        exp = exp_net.reaction_collection
        cls_ = cls_net.reaction_collection
    all_shared = 0
    unique_exp = 0
    exp_and_cls = 0
    for obj in exp.objects():
        if cls_.has(obj):
            exp_and_cls += 1
        else:
            unique_exp += 1
    return {"all_shared": all_shared, 
            "exp_and_cls": exp_and_cls, 
            "unique_exp": unique_exp, 
            "unique_cls": len(cls_) - all_shared - exp_and_cls}

def missing_reactions_categories(name: str, exp_net: tn.core.RxnGraph, cls_net: tn.core.RxnGraph):
    counter = {
        "BAD_SPECIE": 0,
        "BAD_FORMAT": 0,
        "ENERGY": 0,
        "PRODUCTS": 0,
        "REACTANTS": 0,
        "OTHER": 0
    }
    spec = importlib.util.spec_from_file_location("run", "../full_enum/run.py")
    run_main = importlib.util.module_from_spec(spec)
    sys.modules["module.name"] = run_main
    spec.loader.exec_module(run_main)    
    run_kwds = run_main.kw_dict[name]
    for rxn in exp_net.reaction_collection.objects():
        if not cls_net.has_reaction(rxn):
            # if any of the ac filters used in the run do not apply to any specie in the reaction - BAD SPECIE
            if any([any([not f.check(s.ac_matrix) for f in run_kwds["ac_filters"]]) for s in rxn.reactants + rxn.products]):
                counter["BAD_SPECIE"] += 1
            # if any of the reactants appear also in products - BAD FORMAT
            elif any([s in rxn.products for s in rxn.reactants]):
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
    return counter

def get_summary_df():
    reactions = {
        "nh3+o2": ["konnov2009", "okafor2018", "miller"],
        "nh3+o2+h2": ["konnov2009", "okafor2018", "miller"],
        "nh3+o2+ch4": ["konnov2009", "okafor2018"],
    }
    res = pd.DataFrame()
    categories = pd.DataFrame()
    for reaction, sources in reactions.items():
        print("reading", reaction)
        cls_net = utils.get_network(source="cls", reaction=reaction, iteration=-1)
        for source in sources:
            print("comparing with", source)
            if source is None:
                source_str = "exp"
            else:
                source_str = "exp_" + source
            exp_net = utils.get_network(source=source_str, reaction=reaction)
            # wiping charge data from exp networks
            for s in exp_net.species:
                s.charge = None
                # adding energy data from calculated graphs
                # s.properties["energy"] = cls_net.specie_collection.get(s).properties["energy"]
            d = get_exp_graph_disections(exp_net, cls_net, focus="species")
            d["focus"] = "species"
            d["reaction"] = reaction
            d["source"] = source
            res = res.append(d, ignore_index=True)
            d = get_exp_graph_disections(exp_net, cls_net, focus="reactions")
            d["focus"] = "reactions"
            d["reaction"] = reaction
            d["source"] = source
            res = res.append(d, ignore_index=True)
            res.to_csv("../analysis_results/compare_to_experiments.csv")
            d = missing_reactions_categories(reaction, exp_net, cls_net)
            d["reaction"] = reaction
            d["source"] = source
            categories = categories.append(d, ignore_index=True)
            categories.to_csv("../analysis_results/compare_to_experiments_categories.csv")

if __name__ == "__main__":
    get_summary_df()
