import os, glob, pickle
import networkx as nx
from FuzzTree import main, print_mapping_on_target_graph
import time
from multiprocessing import Process

def first_test():
    with open("rin147.pickle",'rb') as fP:
        GP = pickle.load(fP)
    with open("6v9d.pickle",'rb') as fT:
        GT = pickle.load(fT)
    
    mapping = main(GP, GT)
    print(mapping, nb_samples=100)
    #sur 10 samples on retrouve : ([(1, ('B', 7)), (2, ('B', 8)), (3, ('B', 12)), (4, ('B', 13)), (5, ('B', 17)), (6, ('B', 18)), (7, ('B', 23)), (8, ('B', 24))], 8)

    print_mapping_on_target_graph(GP, GT, 
mapping = [(1, ('B', 7)), (2, ('B', 8)), (3, ('B', 12)), 
(4, ('B', 13)), (5, ('B', 17)), (6, ('B', 18)), (7, ('B', 23)), 
(8, ('B', 24))], name_file="rin147into6v9d", output_format="png")

    #En chronométrant à la main le nombre de samples (sachant qu'ils sont print) ne semblent pas beaucoup impacter la complexité (2 min 20 sec pour 10 samples et 3min12 sec pour 1000  samples)
    # En particulier sur 1000 samples fuzzy sur les labels on arrive à retrouver des samples corrects (8 exactement), dont certains qu'on ne trouver pas avant sur la strand E

    print_mapping_on_target_graph(GP, GT, 
name_file="rin147into6v9d", output_format="png")


def test_graph_where_pattern_is_detected(GPpath = "rin147.pickle", GTlistfolder = "RNAstorage", BooleanListGPinGT = [0,0,0,0,0], nb_samples=1000):
    """
    Input: - A graph Pattern GP file name GPpath that we supposed to be exatly what we are looking for
           - A list of RNA Target Graphs GTlist as a folder of files GTlistfolder. For each of these GT, we are looking for GP or a fuzzy version of GP in it.
           - For tests purposes; for each GT graph i, position i of BooleanListGTinGP stores if GP is effectively in GT i or not. 
           - number of samples done for each pattern search nb_samples
    Output: returns the list of samples that were correct injective matching if no mistakes were done and raise an error otherwise.
    """
    with open(GPpath,'rb') as fP:
        GP = pickle.load(fP)
    path = os.path.abspath(os.getcwd()) + "/" + GTlistfolder
    resu = []
    k = 0
    for filename in glob.glob(os.path.join(path, '*.pickle')):
        with open(os.path.join(os.getcwd(), filename), 'rb') as fT:
            GT = pickle.load(fT)
            def pro():
                return main(GP, GT, nb_samples)
            p = Process(target=pro, name='Process')
            p.start()
            p.join(timeout=300)
            p.terminate()
            mapping = []
            timer = -1
            if p.exitcode is not None:
                timer = time.time()
                mapping = main(GP, GT, nb_samples)
                timer = time.time() - timer
                mapping = [map for (map, length) in mapping if length == len(GP.nodes())]
            assert( (len(mapping) > 0) >= BooleanListGPinGT[k]) #We allow that as which fuzzy graph we are more permissive than the original search.
            k+=1
            resu.append((len(mapping), timer))
    return resu

#first_test()
#test_graph_where_pattern_is_detected()
import RIN
def import_rin(rin_number):
    with open("CaRNAval_1_as_dictionnary.nxpickled",'rb') as ff:
        Gdict = pickle.load(ff)
    GG = Gdict[rin_number].graph
    relabel_mapping = {}
    for k,ind in enumerate(GG.nodes()):
        relabel_mapping[ind] = k + 1
    GG = nx.relabel_nodes(GG, relabel_mapping)
    name = "rin" + str(rin_number) + ".pickle"
    with open(name, 'wb') as ff:
        pickle.dump(GG, ff)
    print(GG.nodes(data=True))
    print(GG.edges(data=True))

#import_rin(12)

test_graph_where_pattern_is_detected(GPpath = "rin12.pickle", GTlistfolder = "RNAstorage", BooleanListGPinGT = [0,0,0,0,0], nb_samples=1000)



