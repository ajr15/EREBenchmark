import unittest
import torinanet as tn


def test_two_species():
    max_breaking_bonds = 2
    max_forming_bonds = 2
    g = tn.core.RxnGraph()
    s1 = g._read_specie_with_ac_matrix(tn.core.Specie("[O]O"))
    s2 = g._read_specie_with_ac_matrix(tn.core.Specie("[O]O"))
    iterator = tn.iterate.Iterator(tn.core.RxnGraph(), n_workers=2)
    conv_filters = [
                    tn.iterate.conversion_matrix_filters.MaxChangingBonds(max_breaking_bonds + max_forming_bonds),
                    tn.iterate.conversion_matrix_filters.MaxFormingAndBreakingBonds(max_forming_bonds, max_breaking_bonds),
                    tn.iterate.conversion_matrix_filters.OnlySingleBonds()
                ]
    ac_filters = [
                    tn.iterate.ac_matrix_filters.MaxAtomsOfElement({6: 2, 8: 2}), 
                    tn.iterate.ac_matrix_filters.MaxBondsPerAtom(), 
                    tn.iterate.ac_matrix_filters.MaxComponents(2),
                    tn.iterate.ac_matrix_filters.MaxRingNumber(0)
            ]
    res = iterator.iterate_over_species(g, s1, s2, ac_filters, conv_filters)
    for rxn in res.reactions:
        print("{} -> {}".format(" + ".join([s.ac_matrix.to_specie().identifier.strip() for s in rxn.reactants]),
                                " + ".join([s.ac_matrix.to_specie().identifier.strip() for s in rxn.products])))

if __name__ == "__main__":
    test_two_species()