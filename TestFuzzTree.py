import os, glob, pickle
from FuzzTree import main
from VarnaDrawing import print_mapping_on_target_graph
import time
from multiprocessing import Process, Queue

DEBUG=1


def first_test_mapping():
    with open("rin163.pickle",'rb') as fP:
        GP = pickle.load(fP)
    #with open("RNAstorage/1Y27-1B53misslabeledCWW1falselabelCHSintoCWH.nxpickle",'rb') as fT:
    #with open("RNAstorage/1Y27-4falselabelsCHSintoCWHandTHWintoTHSandTHWintoTHSandcWWintoTHS.nxpickle",'rb') as fT:
    with open("RNAstorage/1Y27.nxpickle",'rb') as fT:
    #with open("RNAstorage/1Y27-2edgesmissing.nxpickle",'rb') as fT:
    #with open("RNAstorage/1Y27-1gapinside.nxpickle",'rb') as fT:
        GT = pickle.load(fT)

    timer = time.time()
    #mapping = main(GP, GT, 80, 2, 4, nb_samples=100)
    mapping = main(GP, GT, 0, 0, 0, nb_samples=100)
    timer = time.time() - timer
    if DEBUG:
        print("\nmapping", mapping)
        print("\ntime", timer)


def first_test_varna_with_mapping():
    with open("rin163.pickle",'rb') as fP:
        GP = pickle.load(fP)
    with open("RNAstorage/1Y27.nxpickle",'rb') as fT:
        GT = pickle.load(fT)
    mapping = [(1, ('X', 52)), (2, ('X', 46)), (3, ('X', 47)), 
    (4, ('X', 48)), (5, ('X', 49)), (6, ('X', 50)), 
    (7, ('X', 51)), (8, ('X', 20)), (9, ('X', 53)), 
    (10, ('X', 22)), (11, ('X', 24)), (12, ('X', 25)), 
    (13, ('X', 21)), (14, ('X', 54))]
    print_mapping_on_target_graph([], GT, mapping, output_format = "png", name_file = "rin163into1Y27")

def first_test_varna_without_mapping():
    with open("rin163.pickle",'rb') as fP:
        GP = pickle.load(fP)
    with open("RNAstorage/2XNZ.nxpickle", 'rb') as fT:
        GT = pickle.load(fT)
    print_mapping_on_target_graph(GP, GT, output_format = "png", name_file = "rin163into2XNZ", E=80, B=2, A=0)


def test_graph_where_pattern_is_detected(GPpath = "rin163.pickle", GTlistfolder = "RNAstorage", E=0 , B=0, A=0, maxGAPdistance = 3, nb_samples=1000, remove_near=True, timeout=800):
    """
    Input: - A graph Pattern GP file named GPpath that we supposed to be exactly the pattern that we are looking for.
           - A list of RNA Target Graphs GTlist as a folder of files GTlistfolder. For each of these GT, we are looking for GP or a fuzzy version of GP in it.
           - The Fuzzy Parameters E, B, A that are respectively threshold on sum of isostericity, number of edges and sum of gap distances.
           - maxGAPdistance, fuzzy parameter about how far we allow to look for gaps.
           - number of samples done for each searched pattern nb_samples
           - remove_near to True remove all edges labelled "near" and that are not as precise as we want.
           - As we have no return when the process fails, we put a timeout
    Output: returns a list of triplets(filename, mapping, time) with mapping, a mapping from the ones sample by the process. Left empty if no mapping where found in time.
    """
    with open(GPpath,'rb') as fP:
        GP = pickle.load(fP)
    path = os.path.abspath(os.getcwd()) + "/" + GTlistfolder
    resu = []
    for filename in glob.glob(os.path.join(path, '*.nxpickle')):
        with open(os.path.join(os.getcwd(), filename), 'rb') as fT:
            GT = pickle.load(fT)
            if remove_near: #We reove the near edges only if requested.
                edges_to_remove = [(i, j) for (i, j, t) in GT.edges.data() if t['near'] == True]
                if DEBUG:
                    print("size of near removal", len(edges_to_remove))
                for (i, j) in edges_to_remove:
                    GT.remove_edge(i, j)

            #We use an auxiliary process to be able to carry on even if we timeout.
            def pro(queue):
                queue.put(main(GP, GT, E, B, A, maxGAPdistance=maxGAPdistance, nb_samples=nb_samples, respect_injectivity=1))
            queue = Queue()
            p = Process(target=pro, args=(queue,), name='Process')

            timer = time.time()
            p.start()
            p.join(timeout=timeout)
            p.terminate()
            timer = time.time() - timer

            mapping = []
            if p.exitcode is not None:
                mapping = queue.get()
                if mapping:
                    mapping = [map for (map, length) in mapping][0]
            if DEBUG:
                print("filename, mapping, time", (filename, mapping, timer))
            resu.append((filename, mapping, timer))
    return resu


