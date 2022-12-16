import os, glob, pickle
from FuzzTree import main
from VarnaDrawing import print_mapping_on_target_graph
from SliceInCubes import slicer
import time
from multiprocessing import Pool, Process, Queue
import networkx as nx
import matplotlib.pyplot as plt
import func_timeout

DEBUG=1

def open_graph(graph_path):
    with open(graph_path,'rb') as f:
        G = pickle.load(f)
    return G

def near_removal(GT):
    edges_to_remove = [(i, j) for (i, j, t) in GT.edges.data() if t['near'] == True]
    if DEBUG:
        print("size of near removal", len(edges_to_remove))
    for (i, j) in edges_to_remove:
        GT.remove_edge(i, j)  
    return GT
 
def outer_chain_removal(GT, chains):
    Gnew=nx.DiGraph() #Initiate the new GT graph.
    for ((i, ii),t) in GT.nodes.data():
        if i in chains:
            Gnew.add_node((i, ii), pdb_position = t['pdb_position'], atoms = t['atoms'])
    for ((i, ii),(j, jj),t) in GT.edges.data():
        if i in chains and j in chains:
            Gnew.add_edge((i, ii),(j, jj), label=t['label'], near=t['near'])
    if DEBUG:
        print("nb_nodes_GT_before", len(GT.nodes.data()),"nb_edges_GT_before", len(GT.edges.data()))
        print("nb_nodes_GT_after", len(Gnew.nodes.data()),"nb_edges_GT_after", len(Gnew.edges.data()))
    return Gnew

def weak_similar_mapping(mapping_ref, mapping, GT, strong_mapping):
    """
    Input: - mapping_ref, the perfect mapping that we wan to compare with. 
           - mapping, a mapping that was found in the RNA sampled by the algorithm.
           - GT, the graph target.
           - strong_mapping, the proportio, of the mapping that we ensure to correspond exactly to say that the mapping are "weakly" similar
    Output: Returns 1, if the mappings are weakly the same given a proportion and 0 otherwise.
    """
    #mapping = mapping[0]
    #print("mapping_ref, mapping",mapping_ref, mapping)
    quantity_nodes_valid = 0
    for (i, j) in mapping:
        ii,t = [(ii, tt['pdb_position']) for ((ii,jj), tt) in GT.nodes.data() if (ii,jj) == j][0]
        #print('i, t', i, t)
        if (i, (ii, t)) in mapping_ref:
            quantity_nodes_valid +=1
    if quantity_nodes_valid >= strong_mapping*len(mapping):
        return 1
    else:
        return 0

def wrapper_main(filename_local_mapping_strong_mapping_timeout_GP_GT_E_B_A_maxGAPdistance_nb_samples_respect_injectivity_D_Distancer_preprocessed):
    (filename, local_mapping, strong_mapping, timeout, GP, GT, E, B, A, maxGAPdistance, nb_samples, respect_injectivity, D, Distancer_preprocessed) = filename_local_mapping_strong_mapping_timeout_GP_GT_E_B_A_maxGAPdistance_nb_samples_respect_injectivity_D_Distancer_preprocessed
    timer = time.time()
    if timeout == -1:
        mapping = main(GP, GT, E, B, A, maxGAPdistance=maxGAPdistance, nb_samples=nb_samples, respect_injectivity=respect_injectivity, D = D, Distancer_preprocessed = Distancer_preprocessed)
        mapping = [mapp for (mapp, _) in mapping]
    else:
        def computation():
            loc = main(GP, GT, E, B, A, maxGAPdistance=maxGAPdistance, nb_samples=nb_samples, respect_injectivity=respect_injectivity, D = D, Distancer_preprocessed = Distancer_preprocessed)
            return [mapp for (mapp, _) in loc]
        try:
            mapping = func_timeout.func_timeout(timeout, computation)
        except: #func_timeout.FunctionTimedOut:
            mapping = []
    timer = time.time() - timer
    if local_mapping:
        proportion = []
        for loc_mapp in local_mapping:
            if mapping:
                proportion.append(len([mapp for mapp in mapping if weak_similar_mapping(loc_mapp, mapp, GT, strong_mapping) ])/len(mapping))
            else:
                proportion.append(0)
    else:
        proportion = [1.01] #It means here that we have no way to verify that the mappings are the ones that we want here.
    if DEBUG:
        print("filename", filename, "timer", timer, "proportion", proportion, "number_mapping", len(mapping))
    return (filename, timer, proportion, mapping)


def fusion_resu_cube(resu): #resu empty never happens
    full_filename = resu[0][0]
    full_timer = max([timer for (_, timer, _, _) in resu])
    full_mapping = []
    full_proportion = [0]*len(resu[0][2])
    for (_, _, proportion, mapping) in resu:
        for k in range(len(proportion)):
            full_proportion[k] += proportion[k] * len(mapping)
        full_mapping+= mapping
    if full_mapping:
        for k in range(len(full_proportion)):
            full_proportion[k] = full_proportion[k]/len(full_mapping)
    return (full_filename, full_timer, full_proportion, full_mapping)

def initialise_perfect_mapping(perfect_mapping, motifs_mapping):
    if motifs_mapping:
        new_motifs_mapping = [-1 for i in range(len(motifs_mapping) + 1)]
        for (k, l) in motifs_mapping:
            new_motifs_mapping[l] = k
    new_perfect_mapping = {}
    for index in range(len(perfect_mapping)):
        chains = []
        for mappinger in perfect_mapping[index][1]:
            (cha, num) = mappinger[1]
            if cha not in chains:
                chains.append(cha) #We retrieve the letter of the chain as we will look only at the objective chain in order to study a smaller graph
        chains.sort()
        chains = tuple(chains)
        if motifs_mapping:
            local_mapping = []
            for (num, stay) in perfect_mapping[index][1]:
                local_mapping.append((new_motifs_mapping[num], stay))
        else:
            local_mapping = perfect_mapping[index][1]
        if (perfect_mapping[index][0], chains) in new_perfect_mapping.keys():
            new_perfect_mapping[(perfect_mapping[index][0], chains)] =new_perfect_mapping[(perfect_mapping[index][0], chains)] + [local_mapping]
        else:
            new_perfect_mapping[(perfect_mapping[index][0], chains)] = [local_mapping]
    return new_perfect_mapping



def test_mapping(GPpath, GTpath, E, B, A, maxGAPdistance=3, nb_samples=1000, respect_injectivity=1, D = 5, Distancer_preprocessed = {}):
    GP = open_graph(GPpath)
    GT = open_graph(GTpath)
    timer = time.time()
    mapping = main(GP, GT, E, B, A, maxGAPdistance=maxGAPdistance, nb_samples=nb_samples, respect_injectivity=respect_injectivity, D = D, Distancer_preprocessed = Distancer_preprocessed)
    timer = time.time() - timer
    if DEBUG:
        print("\nmapping", mapping)
        print("\ntime", timer)
    return mapping

def test_varna(name_file,GPpath, GTpath, show=1, output_format='png', E = 0, B = 0, A = 0, maxGAPdistance=3, nb_samples=1000, respect_injectivity=1, D = 5, Distancer_preprocessed = {}, mapping = []):
    GP = open_graph(GPpath)
    GT = open_graph(GTpath)
    if mapping:
        print_mapping_on_target_graph([], GT, mapping=mapping, output_format = output_format, name_file = name_file, show=show)
    else:
        print_mapping_on_target_graph(GP, GT, mapping = [], output_format = output_format, name_file = name_file, show=show, E=E, B=B, A=A, maxGAPdistance=maxGAPdistance, nb_samples=nb_samples, respect_injectivity=respect_injectivity, D = D)



def test_GP_into_multiples_GT(GPpath, GTlistfolder = "bigRNAstorage", threshold_bigGT = 500, strong_mapping = 1, respect_injectivity=1, E=0 , B=0, A=0, maxGAPdistance = 3, nb_samples=1000, remove_near=True, timeout=-1, D = 5, nb_procs = 32, perfect_mapping = [], motifs_mapping = []):
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
    #TODO : update doc string
    GP = open_graph(GPpath)
    path = os.path.abspath(os.getcwd()) + "/" + GTlistfolder
    smallGT = []
    bigGT = []
    if perfect_mapping:
        perfect_mapping = initialise_perfect_mapping(perfect_mapping, motifs_mapping)
        path_list = []
        chains_list = []
        RNA_list = []
        for (RNAname, chains) in perfect_mapping.keys():
            path_list.append(path+ "/" + RNAname + '.nxpickle')
            RNA_list.append(RNAname)
            chains_list.append(list(chains))
    else:
        path_list = glob.glob(os.path.join(path, '*.nxpickle'))
    if DEBUG:
        if perfect_mapping:
            print("Perfect_mapping", perfect_mapping)
        print("File exploration", path_list)
    for num_list, filename in enumerate(path_list):
        GT = open_graph(os.path.join(os.getcwd(), filename))
        compact_filename = (filename.split('/'))[-1][:-9]
        local_mapping = []
        if remove_near: #We remove the near edges only if requested.
            GT = near_removal(GT)
             
        if perfect_mapping:
            chains = chains_list[num_list]
            GT = outer_chain_removal(GT, chains)
            local_mapping = perfect_mapping[(RNA_list[num_list], tuple(chains))]
        if len(list(GT.nodes())) > threshold_bigGT:
            if DEBUG:
                print("Big GT", compact_filename, "size", len(GT.nodes()), "\n")
            bigGT.append((compact_filename, GT.copy(), local_mapping))
        else:
            if DEBUG:
                print("Small GT", compact_filename, "size", len(GT.nodes()), "\n")
            smallGT.append((compact_filename, GT.copy(), local_mapping))
    entry = []
    for (filename, GT, local_mapping) in smallGT:
        entry.append((filename, local_mapping, strong_mapping, timeout, GP, GT, E, B, A, maxGAPdistance, nb_samples, respect_injectivity, D, {}))
    with Pool(nb_procs) as pool:
        resu = list(pool.imap_unordered(wrapper_main, entry))
    for (filename, GT, local_mapping) in bigGT:
        entry = []
        graph_grid, Distancer = slicer(GP, GT,  size_cube_versus_radius=1, filename=filename) #instead of 0.5 for now to have less cubes
        for GTsmall in graph_grid:
            local_nb_samples = int(nb_samples/len(graph_grid)) + 1
            entry.append((filename, local_mapping, strong_mapping, timeout, GP, GTsmall, E, B, A, maxGAPdistance, local_nb_samples, respect_injectivity, D, Distancer))
        with Pool(nb_procs) as pool:
            resu_big = list(pool.imap_unordered(wrapper_main, entry))
        resu.append(fusion_resu_cube(resu_big))
    compact_resu = [(i,j,k) for (i,j,k,_) in resu]
    print("compact_resu", compact_resu)
    print("\n resu", resu)
    return resu










def newmain4(GP_GTsmall_ind_debug_E_B_A_maxGAPdistance_nb_samples_D_Distancer):
    (GP, GTsmall, ind_debug, E, B, A, maxGAPdistance, nb_samples, D, Distancer) = GP_GTsmall_ind_debug_E_B_A_maxGAPdistance_nb_samples_D_Distancer
    try:
        mapping = main(GP, GTsmall, E, B, A, maxGAPdistance=maxGAPdistance, nb_samples=nb_samples, respect_injectivity=1, D = D, Distancer_preprocessed=Distancer)
        print("I got a cube done", ind_debug, "len:", len(mapping), "\n")
    except:
        mapping = []
        print("I got a void cube done",ind_debug, "\n")
    if isinstance(mapping, Exception):
        mapping = []
    return mapping

def test_perfect_mapping_multiprocess_oneRNA_sliced(perfect_mapping_one_RNA, GPpath, pattern_name, E=0 , B=0, A=0, maxGAPdistance = 3, nb_samples=10, remove_near=True, timeout=800, D = 5, motifs_mapping = [], strong_mapping = 1):
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
    chains = []
    for mappinger in perfect_mapping_one_RNA[1]:
        (cha, num) = mappinger[1]
        if cha not in chains:
            chains.append(cha) #We retrieve the letter of the chain as we will look only at the objective chain in order to study a smaller graph
    chains.sort()
    RNAname = perfect_mapping_one_RNA[0]
    with open(GPpath,'rb') as fP:
        GP = pickle.load(fP)
    path = os.path.abspath(os.getcwd()) + "/bigRNAstorage/" + RNAname + '.nxpickle'
    perfect_mapping_one_RNA = perfect_mapping_one_RNA[1]
    if DEBUG:
        print("perfect_mapping", perfect_mapping_one_RNA)
        print("Studied RNA file", path)
    resu = []
    entry = []
    with open(path, 'rb') as fT:
        GT = pickle.load(fT)
        if remove_near: #We reove the near edges only if requested.
            edges_to_remove = [(i, j) for (i, j, t) in GT.edges.data() if t['near'] == True]
            if DEBUG:
                print("size of near removal", len(edges_to_remove))
            for (i, j) in edges_to_remove:
                GT.remove_edge(i, j)
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
        if motifs_mapping:
            local_mapping = []
            new_motifs_mapping = [-1 for i in range(len(motifs_mapping) + 1)]
            for (k, l) in motifs_mapping:
                new_motifs_mapping[l] = k
            pre_local_mapping = perfect_mapping_one_RNA
            for (num, stay) in pre_local_mapping:
                local_mapping.append((new_motifs_mapping[num], stay))
        else:
            local_mapping = perfect_mapping_one_RNA
        if DEBUG:
            print("local_mapping", local_mapping)
        graph_grid, Distancer = slicer(GP, GT,  size_cube_versus_radius=1, filename=RNAname) #instead of 0.5 for now to have less cubes
        for ind_debug, GTsmall in enumerate(graph_grid):
            entry.append((GP, GTsmall, ind_debug, E, B, A, maxGAPdistance, nb_samples, D, Distancer))
    timer = time.time()
    with Pool(32) as pool:
        resu = list(pool.imap_unordered(newmain4, entry))
    if DEBUG:
        print("resu", resu)
    resubis = []
    for ii in resu:
        resubis+= ii
    resu = resubis
    proportion = len([mapp for mapp in resu if weak_similar_mapping(local_mapping, mapp, GT, strong_mapping) ])/len(resu)
    timer = time.time() - timer
    if DEBUG:
        print("RNAname", RNAname, "proportion", proportion, "timer", timer, "resu", resu)
    return RNAname, proportion, timer, resu
    
def bar_graph_3proportions_1time_by_filename(resu, title, bar_length = 0.3):
    """
    Input: - resu, the list of triplets (filename, proportion, time) obtained during the test.
           - title, a title for the plot.
           - bar_length, size of bar drawn.
    Output: returns a bar plot for the proportion and the time for each RNA studied during the tests.
    """
    plt.rc('xtick', labelsize=7)
    y11 = []
    y12 = []
    y13 = []
    y2 = []
    graph_names = []
    max_proportion = 1
    max_time = 1
    for (name, time, proportion, _) in resu:
        max_proportion = max(max_proportion, len(proportion))
        max_time = max(max_time, time)
    if max_proportion > 3:
        print("TOO MUCH PROPORTION TO HANDLE LIMITED TO 3 FOR NOW")
        return 
    for (name, time, porportion, _) in resu:
        if len(proportion) >= 1:
            y11.append(proportion[0])
        else:
            y11.append(0)
        if len(proportion) >= 2:
            y12.append(proportion[1])
        else:
            y12.append(0)
        if len(proportion) >= 3:
            y13.append(proportion[2])
        else:
            y13.append(0)
        y2.append(time)
        graph_names.append(name)

    x1 = range(len(y11)) 

    x2 = [i + bar_length for i in x1]


    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.set_ylabel('Proportion of "perfect" mapping', color='orange')
    ax1.set(ylim=(0, 1))
    ax2.set(ylim=(1, max_time + 1))
    ax2.set_ylabel('Time', color='blue')
    ax1.bar(x1, y11, width = bar_length, color = 'orange', edgecolor = 'black', linewidth = 2)
    ax1.bar(x1, y12, width = bar_length, color = 'green', bottom=y11, edgecolor = 'black', linewidth = 2)
    ax1.bar(x1, y13, width = bar_length, color = 'red', bottom=[y11[i] + y12[i] for i in range(len(y11))], edgecolor = 'black', linewidth = 2)
    ax2.bar(x2, y2, width = bar_length, color = 'blue', edgecolor = 'black', linewidth = 2, log = True)
    ax1.set_facecolor(color='white')
    ax2.set_facecolor(color='white')
        # fontsize of the x and y labels
    fig.set_facecolor(color='white')
    plt.xticks([r + bar_length / 2 for r in range(len(y11))], graph_names)
    plt.title(title)
    plt.show()

