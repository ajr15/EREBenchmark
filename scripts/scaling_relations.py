import pandas as pd
from matplotlib import pyplot as plt
import torinanet as tn
import utils


def get_exp_graph_disections(exp_net: tn.core.RxnGraph, cls_net: tn.core.RxnGraph):
    exp = exp_net.reaction_collection
    cls_ = cls_net.reaction_collection
    ajr = []
    for obj in exp.objects():
        if cls_.has(obj):
            obj.properties.update(cls_.get(obj).properties)
            ajr.append(obj)
    return ajr

def calc_reaction_energy(rxn: tn.core.Reaction):
    return sum([s.properties["energy"] for s in rxn.products]) - sum([s.properties["energy"] for s in rxn.reactants])

def gather_scaling_data(rxns):
    data = []
    for rxn in rxns:
        if "Ea" in rxn.properties and not pd.isna(rxn.properties["Ea"]):
            data.append([rxn.pretty_string(), rxn.properties["Ea"], rxn.properties["A"], rxn.properties["beta"], calc_reaction_energy(rxn)])
    return pd.DataFrame(data=data, columns=["rxn", "Ea", "A", "beta", "deltaE"])


if __name__ == "__main__":
    source = "okafor2018"
    print("Reading exp graph...")
    exp_net = utils.get_network(source="exp_{}".format(source), reaction="nh3+o2")
    for s in exp_net.species:
        s.charge = None

    print("Reading cls graph...")
    cls_net = utils.get_network(source="cls", reaction="nh3+o2")
    print("Comparing reactions...")
    rxns = get_exp_graph_disections(exp_net, cls_net)
    print("Collecting scaling data...")
    df = gather_scaling_data(rxns)
    print(df)
    df.to_csv("../analysis_results/{}_scaling_data.csv".format(source))
    print("DONE")
    