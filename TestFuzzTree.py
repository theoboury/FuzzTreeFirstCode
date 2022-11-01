import os, glob, pickle
from test import main

def test_graph_where_pattern_is_detected(GP, GTlist, BooleanListGPinGT, nb_samples):
    """
    Input: - A graph Pattern GP that we supposed to be exatly what we are looking for
           - A list of RNA Target Graphs GTlist. For each of these GT, we are looking for GP or a fuzzy version of GP in it.
           - For tests purposes; for each GT graph i, position i of BooleanListGTinGP stores if GP is effectively in GT i or not. 
           - number of samples done for each pattern search nb_samples
    Output: returns the list of samples that were correct injective matching if no mistakes were done and raise an error otherwise.
    """
    resu = []
    for k, GT in enumerate(GTlist):
        mapping = main(GP, GT, nb_samples)
        mapping = [map for (map, length) in mapping if length == len(GP.nodes())]
        assert( (len(mapping) > 0) >= BooleanListGPinGT[k]) #We allow that as which fuzzy graph we are more permissive than the original search.
        resu.append(len(mapping))
    return resu


def iterate_over_graphs_where_pattern_is_detected(GPpath = "rin147.pickle", nb_samples=1000):
    with open(GPpath,'rb') as fP:
        GP = pickle.load(fP)
    path = os.path.abspath(os.getcwd()) + "/RNAstorage"
    for filename in glob.glob(os.path.join(path, '*.pickle')):
        with open(os.path.join(os.getcwd(), filename), 'rb') as fT:
            GT = pickle.load(fT)
