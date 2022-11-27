import os, glob, pickle
#from FuzzTree_with_dist_preprocessing_and_set_weight import main
#from FuzzTree_without_dist_preprocessing import main
from FuzzTree import main
from VarnaDrawing import print_mapping_on_target_graph
import time
import sys
from multiprocessing import Process, Queue
import networkx as nx

DEBUG=1


def first_test_mapping():
    #TODO : to be removed for final version
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
    #TODO : to be removed for final version
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
    #TODO : to be removed for final version
    with open("rin163.pickle",'rb') as fP:
        GP = pickle.load(fP)
    with open("RNAstorage/2XNZ.nxpickle", 'rb') as fT:
        GT = pickle.load(fT)
    print_mapping_on_target_graph(GP, GT, output_format = "png", name_file = "rin163into2XNZ", E=80, B=2, A=0)


def test_graph_where_pattern_is_detected(GPpath, GTlistfolder = "RNAstorage", E=0 , B=0, A=0, maxGAPdistance = 3, nb_samples=1000, remove_near=True, timeout=800, D = 5):
    """
    Input: - A graph Pattern GP file named GPpath that we supposed to be exactly the pattern that we are looking for.
           - A list of RNA Target Graphs GTlist as a folder of files GTlistfolder. For each of these GT, we are looking for GP or a fuzzy version of GP in it.
           - The Fuzzy Parameters E, B, A that are respectively threshold on sum of isostericity, number of edges and sum of gap distances.
           - maxGAPdistance and D, fuzzy parameter about how far we allow to look respectively for gaps and missing edges.
           - number of samples done for each searched pattern nb_samples
           - remove_near to True remove all edges labelled "near" and that are not as precise as we want.
           - As we have no return when the process fails, we put a timeout
    Output: returns a list of triplets(filename, mapping, time) with mapping, a mapping from the ones sample by the process. Left empty if no mapping where found in time.
    """
    with open(GPpath,'rb') as fP:
        GP = pickle.load(fP)
    path = os.path.abspath(os.getcwd()) + "/" + GTlistfolder
    resu = []
    if DEBUG:
        print("File exploration", glob.glob(os.path.join(path, '*.nxpickle')))
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
                loc = main(GP, GT, E, B, A, maxGAPdistance=maxGAPdistance, nb_samples=nb_samples, respect_injectivity=1, D = D)
                queue.put(loc)
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
                if isinstance(mapping, Exception):
                    mapping = []
                elif mapping:
                    mapping = [map for (map, length) in mapping][0]
            filename = (filename.split('/'))[-1]
            if DEBUG:
                print("filename, mapping, time", (filename, mapping, timer))
            resu.append((filename, mapping, timer))
    return resu


def similar_mapping(mapping_ref, mapping, GT):
    """
    Input: - mapping_ref, the perfect mapping that we wan to compare with. 
           - mapping, a mapping that was found in the RNA sampled by the algorithm.
           - GT, the graph target.
    Output: Returns 1, if the mappings are the same and 0 otherwise.
    """
    mapping = mapping[0]
    #print("mapping_ref, mapping",mapping_ref, mapping)
    for (i, j) in mapping:
        ii,t = [(ii, tt['pdb_position']) for ((ii,jj), tt) in GT.nodes.data() if (ii,jj) == j][0]
        #print('i, t', i, t)
        if (i, (ii, t)) not in mapping_ref:
            return 0
    return 1


def test_perfect_mapping(perfect_mapping, GPpath, E=0 , B=0, A=0, maxGAPdistance = 3, nb_samples=10, remove_near=True, timeout=800, D = 5):
    """
    Input: - A graph Pattern GP file named GPpath that we supposed to be exactly the pattern that we are looking for.
           - A perfect_mapping that icontains the list of RNA Target Graphs, but also for each of them the "perfect mapping" that we want to find nad that is indicated in the RNA data.
           - The Fuzzy Parameters E, B, A that are respectively threshold on sum of isostericity, number of edges and sum of gap distances.
           - maxGAPdistance and D, fuzzy parameter about how far we allow to look respectively for gaps and missing edges.
           - number of samples done for each searched pattern nb_samples
           - remove_near to True remove all edges labelled "near" and that are not as precise as we want.
           - As we have no return when the process fails, we put a timeout
    Output: returns a list of triplets(filename, mapping, time) with mapping, a mapping from the ones sample by the process. Left empty if no mapping where found in time.
    """
    with open(GPpath,'rb') as fP:
        GP = pickle.load(fP)
    path = os.path.abspath(os.getcwd()) + "/bigRNAstorage/"
    path_list = []
    for (RNAname, _) in perfect_mapping:
        path_list.append(path+ RNAname + '.nxpickle')
    if DEBUG:
        print("perfect_mapping", perfect_mapping)
        print("list of studied RNA files", path_list)
    resu = []
    for index, filename in enumerate(path_list):
        with open(filename, 'rb') as fT:
            GT = pickle.load(fT)
            if remove_near: #We reove the near edges only if requested.
                edges_to_remove = [(i, j) for (i, j, t) in GT.edges.data() if t['near'] == True]
                if DEBUG:
                    print("size of near removal", len(edges_to_remove))
                for (i, j) in edges_to_remove:
                    GT.remove_edge(i, j)
            chains = []
            for mappinger in perfect_mapping[index][1]:
                (cha, num) = mappinger[1]
                if cha not in chains:
                    chains += cha #We retrieve the letter of the chain as we will look only at the objective chain in order to study a smaller graph
            print("chains", chains)
            if DEBUG:
                print("nb_nodes_GT_before", len(GT.nodes.data()),"nb_edges_GT_before", len(GT.edges.data()))
            Gnew=nx.DiGraph() #Initiate the new GT graph.
            for ((i, ii),t) in GT.nodes.data():
                if i in chains:
                    Gnew.add_node((i, ii), pdb_position = t['pdb_position'], atoms = t['atoms'])
            for ((i, ii),(j, jj),t) in GT.edges.data():
                if i in chains and j in chains:
                    Gnew.add_edge((i, ii),(j, jj), label=t['label'], near=t['near'])
            GT = Gnew
            if DEBUG:
                print("nb_nodes_GT_after", len(GT.nodes.data()),"nb_edges_GT_after", len(GT.edges.data()))
            #We use an auxiliary process to be able to carry on even if we timeout.
            def pro(queue):
                loc = main(GP, GT, E, B, A, maxGAPdistance=maxGAPdistance, nb_samples=nb_samples, respect_injectivity=1, D = D)
                queue.put(loc)
            queue = Queue()
            p = Process(target=pro, args=(queue,), name='Process')
            timer = time.time()
            p.start()
            p.join(timeout=timeout)
            p.terminate()
            timer = time.time() - timer
            proportion = -1
            mapping = []
            if p.exitcode is not None:
                mapping = queue.get()
                if isinstance(mapping, Exception):
                    mapping = []
                elif mapping:
                    #We compute the proportion of mappings that correspond to the "perfect mapping"
                    print('mapping_first', mapping[0])
                    proportion = len([mapp for mapp in mapping if similar_mapping(perfect_mapping[index][1], mapp, GT) ])/len(mapping)
            filename = (filename.split('/'))[-1]
            if DEBUG:
                print("filename, proportion, time", (filename, proportion, timer))
            resu.append((filename, mapping, timer))
    return resu

#TODO: add the tests here for final version