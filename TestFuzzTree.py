import os, glob, pickle
import networkx as nx
from FuzzTree import main
from VarnaDrawing import print_mapping_on_target_graph
import time
from multiprocessing import Process, Queue
from RIN import import_rin

def first_test_mapping():
    #rin147into6v9d  rin169into1Y27
    with open("rin163.pickle",'rb') as fP:
        GP = pickle.load(fP)
    with open("1Y27.nxpickle",'rb') as fT:
        GT = pickle.load(fT)
    timer = time.time()
    mapping = main(GP, GT, nb_samples=100)
    timer = time.time() - timer
    print("mapping", mapping)
    print("time", timer)
   
def first_test_varna():
    with open("rin163.pickle",'rb') as fP:
        GP = pickle.load(fP)
    with open("1Y27.nxpickle",'rb') as fT:
        GT = pickle.load(fT)
    mapping = [(1, ('X', 52)), (2, ('X', 46)), (3, ('X', 47)), 
    (4, ('X', 48)), (5, ('X', 49)), (6, ('X', 50)), 
    (7, ('X', 51)), (8, ('X', 20)), (9, ('X', 53)), 
    (10, ('X', 22)), (11, ('X', 24)), (12, ('X', 25)), 
    (13, ('X', 21)), (14, ('X', 54))]
    print_mapping_on_target_graph([], GT, mapping, output_format = "png", name_file = "rin163into1Y27")



def test_graph_where_pattern_is_detected(GPpath = "rin163.pickle", GTlistfolder = "RNAstorage", BooleanListGPinGT = [0,0,0,0,0], nb_samples=1000, remove_near=True):
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
    print(glob.glob(os.path.join(path, '*.nxpickle')))
    for filename in glob.glob(os.path.join(path, '*.nxpickle')):
        with open(os.path.join(os.getcwd(), filename), 'rb') as fT:
            GT = pickle.load(fT)
            if remove_near:
                edges_to_remove = [(i, j) for (i, j, t) in GT.edges.data() if t['near'] == True]
                print("size of near removal", len(edges_to_remove))
                for (i, j) in edges_to_remove:
                    GT.remove_edge(i, j)
            def pro(queue):
                queue.put(main(GP, GT, nb_samples))
            queue = Queue()
            p = Process(target=pro, args=(queue,), name='Process')
            timer = time.time()
            p.start()
            p.join(timeout=800)
            p.terminate()
            timer = time.time() - timer
            mapping = []
            if p.exitcode is not None:
                mapping = queue.get()
                mapping = [map for (map, length) in mapping if length == len(GP.nodes())]
            assert( (len(mapping) > 0) >= BooleanListGPinGT[k]) #We allow that as which fuzzy graph we are more permissive than the original search.
            k+=1
            print("(len(mapping), timer)", (len(mapping), timer))
            resu.append((len(mapping), timer))
    return resu



import_rin(163)

#first_test_mapping()

#first_test_varna()


test_graph_where_pattern_is_detected(GPpath = "rin163.pickle", GTlistfolder = "RNAstorage", BooleanListGPinGT = [0,0,0,0,0,0,0,0,0,0], nb_samples=100)


#([(1, ('X', 52)), (2, ('X', 46)), (3, ('X', 47)), (4, ('X', 48)), (5, ('X', 49)), (6, ('X', 50)), (7, ('X', 51)), (8, ('X', 20)), (9, ('X', 53)), (10, ('X', 22)), (11, ('X', 24)), (12, ('X', 25)), (13, ('X', 21)), (14, ('X', 54))], 14)
#found for rin163 into 1Y27 in 506.2300956249237 secondes

