import os, glob, pickle
from FuzzTree import main
from VarnaDrawing import print_mapping_on_target_graph
from SliceInCubes import slicer
import time
from multiprocessing import Pool
import networkx as nx
import matplotlib.pyplot as plt
import func_timeout

DEBUG=1

def open_graph(graph_path):
    """
    Input: A graph path graph_path.
    Output: Returns graph G extracted from graphpath in a networkx format.
    """
    with open(graph_path,'rb') as f:
        G = pickle.load(f)
    return G

def near_removal(GT):
    """
    Input: A target graph GT.
    Output: Returns GT where every edges labeled "near" are removed.
    """
    edges_to_remove = [(i, j) for (i, j, t) in GT.edges.data() if t['near'] == True]
    if DEBUG:
        print("size of near removal", len(edges_to_remove))
    for (i, j) in edges_to_remove:
        GT.remove_edge(i, j)  
    return GT
 
def outer_chain_removal(GT, chains):
    """
    Input: - A target graph GT.
           - A list of chains chains subset of the ones present in GT.
    Output: Returns a copy of GT, where every chains not in chains list are removed.
    """
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

def rename_author_position(GT, GTref):
    """
    Input:  A target graph GT.
    Output: Change name of author position in pdb_position
    """
    GTref = near_removal(GTref)
    #print(GT.nodes.data())
    #print(GTref.nodes.data())
    Gnew=nx.DiGraph() #Initiate the new GT graph.
    #print(GT[('A', 20)])
    auth = {}
    for ((i, j), t) in GT.nodes.data():
        if 'author_position' in t.keys():
            auth[i] = t['author_chain']
    for ((i, j), t) in GT.nodes.data():
        if 'author_position' not in t.keys():
            authpos = GTref.nodes.data()[(auth[i], j)]['pdb_position']
        else:
            authpos = t['author_position']
        if 'atoms' not in t.keys():
            at = []
        else:
            #print("\nHOLE", t['atoms'])
            at = []
            for (type, label, blub, pos1, pos2, pos3) in t['atoms']:
                tt = {}
                tt['type'] = type
                tt['label']= label
                tt['position'] = (pos1, pos2, pos3)
                at.append(tt)
        Gnew.add_node((auth[i], j), pdb_position = authpos, atoms = at)
    for ((i, ii), (j, jj), t) in GT.edges.data():
        if t['label'] in ['B53', 'CHH', 'TWH', 'THW', 'CWW', 'THS', 'TSH', 'CWS', 'CSW', 'CSS', 'CWH','CHW', 'CHS','CSH','TWS','TSW','TSS','TWW','THH']:
            Gnew.add_edge((auth[i], ii), (auth[j], jj), label=t['label'], near=t['near'])
    for (i, j, t) in Gnew.edges.data():

        if t['near'] == True:
            if (i,j) not in GTref.edges():
                #print("AHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH", GTref[i][j]['label'], Gnew[i][j]['label'])
                GTref.add_edge(i, j, label = t['label'], near = t['near'])
    return GTref
    #return Gnew
   
def exact_similar_mapping(mapping_ref, mapping, GT):
    """
    Input: - mapping_ref, the perfect mapping that we want to compare with. 
           - mapping, a mapping that was found in the RNA sampled by the algorithm.
           - GT, the graph target.
    Output: Returns 1, if the mappings are exactly the same with same number of nucleotides and 0 otherwise.
    """
    for (i, j) in mapping:
        ii,t = [(ii, tt['pdb_position']) for ((ii,jj), tt) in GT.nodes.data() if (ii,jj) == j][0]
        if (i, (ii, t)) not in mapping_ref:
            return 0
    return 1


def weak_similar_mapping(mapping_ref, mapping, GT, strong_mapping):
    """
    Input: - mapping_ref, the perfect mapping that we want to compare with. 
           - mapping, a mapping that was found in the RNA sampled by the algorithm.
           - GT, the graph target.
           - strong_mapping, the proportio, of the mapping that we ensure to correspond exactly to say that the mapping are "weakly" similar
    Output: Returns 1, if the mappings are weakly the same given a proportion and 0 otherwise.
    """
    quantity_nodes_valid = 0
    ref_len = len(mapping_ref)
    motif_len = len(mapping)
    mapping_unfold = []
    for (i, j) in mapping:
        ii,t = [(ii, tt['pdb_position']) for ((ii,jj), tt) in GT.nodes.data() if (ii,jj) == j][0]
        mapping_unfold.append((ii, t))
    mapping_unfold = list(set(mapping_unfold))
    ref_unfold = [j for (i, j) in mapping_ref]
    for i in mapping_unfold:
        if i in ref_unfold:
            quantity_nodes_valid +=1
    if quantity_nodes_valid >= strong_mapping*min(ref_len, motif_len):
        return 1
    else:
        return 0

def best_effort_similar_mapping(mapping_ref, mapping, GT):
    """
    Input: - mapping_ref, the perfect mapping that we want to compare with. 
           - mapping, a mapping that was found in the RNA sampled by the algorithm.
           - GT, the graph target.
    Output: Returns 1, if the mappings are cover each others depending on which one is the larger one and 0 otherwise.
    """
    return weak_similar_mapping(mapping_ref, mapping, GT, 1)
         

def wrapper_main(filename_local_mapping_strong_mapping_timeout_GP_GT_E_B_A_maxGAPdistance_nb_samples_respect_injectivity_D_Distancer_preprocessed):
    """
    A wrapper around the main of Fuzztree.py to compute it in a multiprocess fashion with two techniques :
    - for small RNA graph, compute them all in parallel
    - for large RNA graphs, proceed eah RNA one after the other. We slice each in small spheres and compute spheres themselves in parallel.
    """
    (filename, local_mapping, strong_mapping, timeout, GP, GT, E, B, A, maxGAPdistance, nb_samples, respect_injectivity, D, Distancer_preprocessed) = filename_local_mapping_strong_mapping_timeout_GP_GT_E_B_A_maxGAPdistance_nb_samples_respect_injectivity_D_Distancer_preprocessed
    timer = time.time()
    if timeout == -1:
        mapping = main(GP, GT, E, B, A, maxGAPdistance=maxGAPdistance, nb_samples=nb_samples, respect_injectivity=respect_injectivity, D = D, Distancer_preprocessed = Distancer_preprocessed, nb_procs = 1)
        mapping = [mapp for (mapp, _) in mapping]
    else:
        def computation():
            loc = main(GP, GT, E, B, A, maxGAPdistance=maxGAPdistance, nb_samples=nb_samples, respect_injectivity=respect_injectivity, D = D, Distancer_preprocessed = Distancer_preprocessed, nb_procs = 1)
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
                exact = len([mapp for mapp in mapping if exact_similar_mapping(loc_mapp, mapp, GT)])/len(mapping)
                best_effort = len([mapp for mapp in mapping if best_effort_similar_mapping(loc_mapp, mapp, GT)])/len(mapping)
                weak = len([mapp for mapp in mapping if weak_similar_mapping(loc_mapp, mapp, GT, strong_mapping)])/len(mapping)
                proportion.append([exact, best_effort, weak])
            else:
                proportion.append([0, 0, 0])
    else:
        proportion = [[1.01, 1.01, 1.01]] #It means here that we have no way to verify that the mappings are the ones that we want here.
    if DEBUG:
        print("filename", filename, "timer", timer, "proportion", proportion, "number_mapping", len(mapping))
    return (filename, timer, proportion, mapping)


def fusion_resu_cube(resu): #resu empty never happens
    """
    Input: A list of results resu of the form (filename, time, proportion, mappings) for a single RNA.
    Output: Returns the fusion of this results in one single result of the form    
            (filename, full_time, full_proportion, full_mappings).
    """
    full_filename = resu[0][0]
    full_timer = max([timer for (_, timer, _, _) in resu])
    full_mapping = []
    full_proportion = [[0, 0, 0]]*len(resu[0][2])
    for (_, _, proportion, mapping) in resu:
        for k in range(len(proportion)):
            for l in range(3):
                full_proportion[k][l] += proportion[k][l] * len(mapping)
        full_mapping+= mapping
    if full_mapping:
        for k in range(len(full_proportion)):
            for l in range(3):
                full_proportion[k][l] = full_proportion[k][l]/len(full_mapping)
    return (full_filename, full_timer, full_proportion, full_mapping)

def initialise_perfect_mapping(perfect_mapping, motifs_mapping):
    """
    Input: - perfect_mapping, for each RNA GT graphs, the list of mappings that correspond to the research mappings.
           - motifs_mapping, the mapping between the motifs themselves to allow search for a motif in another family.
    Output: Return perfect mapping in a more convenient way and ubstitute motif mapping if needed in a more convenient way.
    """
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



def test_mapping(GPpath, GTpath, E, B, A, maxGAPdistance=3, nb_samples=1000, respect_injectivity=1, D = 5, Distancer_preprocessed = {}, nb_procs = 1):
    """
    Input: - A graph Pattern GP file named GPpath that we supposed to be exactly the pattern that we are looking for.
           - A RNA Target Graph GT file named GTpath. We are looking for GP or a fuzzy version of GP in it.
           - The Fuzzy Parameters E, B, A that are respectively threshold on sum of isostericity, number of edges and sum of gap distances.
           - maxGAPdistance and D, fuzzy parameter about how far we allow to look respectively for gaps and missing edges.
           - number of samples done for each searched pattern nb_samples.
           - A boolean respect_injectivity to ask if we want to ensure that the injectivity is respected or if we allow mapping with doublons.
           - strong_mapping indicates the proportion of the mapping that we want to be correct and the one that we allow to be "faulty".
           - If Distance in GT are already preprocessed, we cant ake them in input to avoid further computation in Distancer_preprocessed.
           - nb_procs, the number of processors allowed for precomputations.
    Output: Return the nb_samples mappings obtained for this instance.
    """
    GP = open_graph(GPpath)
    GT = open_graph(GTpath)
    timer = time.time()
    mapping = main(GP, GT, E, B, A, maxGAPdistance=maxGAPdistance, nb_samples=nb_samples, respect_injectivity=respect_injectivity, D = D, Distancer_preprocessed = Distancer_preprocessed, nb_procs = nb_procs)
    timer = time.time() - timer
    if DEBUG:
        print("\nmapping", mapping)
        print("\ntime", timer)
    return mapping

def test_varna(name_file,GPpath, GTpath, show=1, output_format='png', E = 0, B = 0, A = 0, maxGAPdistance=3, nb_samples=1000, respect_injectivity=1, D = 5, mapping = [], Distancer_preprocessed = {}, nb_procs = 1):
    """
    Input: - A filename for the graph name_file.
           - show to show the graph with matplotlib.pyplot.
           - output_format, specification of the type of storage for the graph.
           - mapping, the mapping that we ant to draw, if affected we discard the instance parameter, f not we compute a list of mappings for curent instance.
           - A graph Pattern GP file named GPpath that we supposed to be exactly the pattern that we are looking for.
           - A RNA Target Graph GT file named GTpath. We are looking for GP or a fuzzy version of GP in it.
           - The Fuzzy Parameters E, B, A that are respectively threshold on sum of isostericity, number of edges and sum of gap distances.
           - maxGAPdistance and D, fuzzy parameter about how far we allow to look respectively for gaps and missing edges.
           - number of samples done for each searched pattern nb_samples.
           - A boolean respect_injectivity to ask if we want to ensure that the injectivity is respected or if we allow mapping with doublons.
           - strong_mapping indicates the proportion of the mapping that we want to be correct and the one that we allow to be "faulty".
           - If Distance in GT are already preprocessed, we cant ake them in input to avoid further computation in Distancer_preprocessed.
           - nb_procs, the number of processors allowed for precomputations.
    Output: Return the Varna drawing for the mapping specified or for the first mapping in the list of mappings obtained for this instance.
    """
    GP = open_graph(GPpath)
    GT = open_graph(GTpath)
    if mapping:
        print_mapping_on_target_graph([], GT, mapping=mapping, output_format = output_format, name_file = name_file, show=show)
    else:
        print_mapping_on_target_graph(GP, GT, mapping = [], output_format = output_format, name_file = name_file, show=show, E=E, B=B, A=A, maxGAPdistance=maxGAPdistance, nb_samples=nb_samples, respect_injectivity=respect_injectivity, D = D, Distancer_preprocessed = Distancer_preprocessed, nb_procs = nb_procs)



def test_GP_into_multiples_GT(GPpath, GTlistfolder = "bigRNAstorage", threshold_bigGT = 500, strong_mapping = 1, respect_injectivity=1, E=0 , B=0, A=0, maxGAPdistance = 3, nb_samples=1000, remove_near=True, timeout=-1, D = 5, nb_procs = 32, perfect_mapping = [], motifs_mapping = [], slice = -1):
    """
    Input: - A graph Pattern GP file named GPpath that we supposed to be exactly the pattern that we are looking for.
           - A list of RNA Target Graphs GTlist as a folder of files GTlistfolder. For each of these GT, we are looking for GP or a fuzzy version of GP in it.
           - The threshold on the number of nucleotides threshold_bigGT above which we consider the RNA GT graph as "big".
           - The Fuzzy Parameters E, B, A that are respectively threshold on sum of isostericity, number of edges and sum of gap distances.
           - maxGAPdistance and D, fuzzy parameter about how far we allow to look respectively for gaps and missing edges.
           - number of samples done for each searched pattern nb_samples.
           - remove_near to True remove all edges labelled "near" and that are not as precise as we want.
           - A boolean respect_injectivity to ask if we want to ensure that the injectivity is respected or if we allow mapping with doublons.
           - strong_mapping indicates the proportion of the mapping that we want to be correct and the one that we allow to be "faulty".
           - As we have no return when the process fails, we put a timeout.
           - perfect_mapping, for each RNA GT graphs, the list of mappings that correspond to the research mappings.
           - motifs_mapping, the mapping between the motifs themselves to allow search for a motif in another family.
    Output: returns a list of quadruplets(filename, time, proportion, mappings) with mappings, mappings from the ones sample by the process. Left empty if no mapping where found in time.
    """
    GP = open_graph(GPpath)
    path = os.path.abspath(os.getcwd()) + "/" + GTlistfolder
    pathbis = os.path.abspath(os.getcwd()) + "/bigRNAstorage" 
    smallGT = []
    bigGT = []
    if GTlistfolder == "bigRNAstorage":
        addon = ".nxpickle"
    else:
        addon = ".pickle"
    if perfect_mapping:
        perfect_mapping = initialise_perfect_mapping(perfect_mapping, motifs_mapping)
        path_list = []
        path_listbis = []
        chains_list = []
        RNA_list = []
        for (RNAname, chains) in perfect_mapping.keys():
            path_list.append(path+ "/" + RNAname + addon)
            path_listbis.append(pathbis + "/" + RNAname + ".nxpickle")
            RNA_list.append(RNAname)
            chains_list.append(list(chains))
    else:
        path_list = glob.glob(os.path.join(path, '*' + addon))
        path_listbis= glob.glob(os.path.join(path,  "*.nxpickle"))
    if DEBUG:
        if perfect_mapping:
            print("Perfect_mapping", perfect_mapping)
        print("File exploration", path_list)
    for num_list, filename in enumerate(path_list):
        #print("filnemae", filename)
        #if  filename != "/home/uqamportable/Documents/FuzzTreeFirstCode/bigRNAstoragenear/3SIU.pickle":
        #    continue
        GT = open_graph(os.path.join(os.getcwd(), filename))
        GTref = open_graph(os.path.join(os.getcwd(), path_listbis[num_list]))
        if GTlistfolder != "bigRNAstorage":
            
            GT = rename_author_position(GT, GTref)

            compact_filename = (filename.split('/'))[-1][:-7]
        else:
            compact_filename = (filename.split('/'))[-1][:-9]
        #print("nodes", [i for (i, t) in GTref.nodes().data() if i not in GT.nodes()] +[i for (i, t) in GT.nodes().data() if i not in GTref.nodes()]  )     
        #print("edges", [i for (i, j, t) in GTref.edges().data() if (i, j) not in GT.edges()] +[i for (i, j, t) in GT.edges().data() if (i, j) not in GTref.edges()] )
        #print([t for (i, j, t) in GT.edges.data()])
        #TODO : long range missing et albel bizarre à modifier
        local_mapping = []
        if remove_near: #We remove the near edges only if requested.
            GT = near_removal(GT)
         # #
        if perfect_mapping:
            chains = chains_list[num_list]
            #print("chains", chains)
            GT = outer_chain_removal(GT, chains)
            #GTref = near_removal(GTref)
            #GTref = outer_chain_removal(GTref, chains)
            #print('GTnodes', [(i, t['atoms']) for (i, t) in GT.nodes.data()][:3])
            #print('GTnodesref', [(i, t['atoms']) for (i, t) in GTref.nodes.data()][:3])
            #print('GTB53', [(i, j) for (i, j, t) in GT.edges.data() if t['label'] == 'B53'])
            #print('GTB53ref', [(i, j) for (i, j, t) in GTref.edges.data() if t['label'] == 'B53'])
            #print('GTlink', len([(i, j, t['label']) for (i, j, t) in GT.edges.data() if t['label'] != 'B53']))
            #print('GTlinkref', len([(i, j, t['label']) for (i, j, t) in GTref.edges.data() if t['label'] != 'B53']))
            #print('GTnear', [(i, j) for (i, j, t) in GT.edges.data() if t['near'] == True])
            #print('GTnearref', [(i, j) for (i, j, t) in GTref.edges.data() if t['near'] == True])
            #GT = near_removal(GT)
            local_mapping = perfect_mapping[(RNA_list[num_list], tuple(chains))]
            #print("loc", local_mapping)
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
    #resu = []
    #for ent in entry:
    #    resu.append(wrapper_main(ent))
    for (filename, GT, local_mapping) in bigGT:
        entry = []
        if slice == -1:
            slice = A/4
        graph_grid, Distancer = slicer(GP, GT, nb_procs, filename=filename, D = slice) #instead of 0.5 for now to have less cubes
        for GTsmall in graph_grid:
            local_nb_samples = max(10, int(nb_samples/len(graph_grid)) + 1) 
            #TODO: Is the maximum necessary here to avoid unlucky mismatch ? Furthermore, should we remove multiple indentical samples. In this case, by what should we replace the proportion ?
            entry.append((filename, local_mapping, strong_mapping, timeout, GP, GTsmall, E, B, A, maxGAPdistance, local_nb_samples, respect_injectivity, D, Distancer))
        with Pool(nb_procs, maxtasksperchild=1) as pool:
            resu_big = list(pool.imap_unordered(wrapper_main, entry))
        resu.append(fusion_resu_cube(resu_big))
    compact_resu = [(i,j,k) for (i,j,k,_) in resu]
    print("compact_resu", compact_resu)
    #print("\n resu", resu)
    return resu

    
def bar_graph_3proportions_1time_by_filename(resu, title, proportion_choice='exact', bar_length = 0.3):
    """
    Input: - resu, the list of triplets (filename, time, proportion, mappings) obtained during the test.
           - title, a title for the plot.
           - bar_length, size of bar drawn.
    Output: returns a bar plot for the proportion and the time for each RNA studied during the tests.
    """
    if proportion_choice == 'exact':
        ind_pro = 0
    elif proportion_choice == 'best_effort':
        ind_pro = 1
    elif proportion_choice == 'weak':
        ind_pro = 2
    else:
        print('Not a valid choice of proportion !\n')
        return
    plt.rc('xtick', labelsize=7)
    y11 = []
    y12 = []
    y13 = []
    y14 = []
    y15 = []
    y2 = []
    graph_names = []
    max_proportion = 1
    max_time = 1
    for (name, time, proportion, _) in resu:
        max_proportion = max(max_proportion, len(proportion))
        max_time = max(max_time, time)
    if max_proportion > 5:
        print("TOO MUCH PROPORTION TO HANDLE LIMITED TO 5 FOR NOW\n")
        return 
    for (name, time, proportion, _) in resu:
        if len(proportion) >= 1:
            y11.append(proportion[0][ind_pro])
        else:
            y11.append(0)
        if len(proportion) >= 2:
            y12.append(proportion[1][ind_pro])
        else:
            y12.append(0)
        if len(proportion) >= 3:
            y13.append(proportion[2][ind_pro])
        else:
            y13.append(0)
        if len(proportion) >= 4:
            y14.append(proportion[3][ind_pro])
        else:
            y14.append(0)
        if len(proportion) >= 5:
            y15.append(proportion[4][ind_pro])
        else:
            y15.append(0)
        y2.append(time)
        graph_names.append(name)
    print(y11, y12, y13, y14, y15)
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
    ax1.bar(x1, y14, width = bar_length, color = 'blue', bottom=[y11[i] + y12[i] + y13[i] for i in range(len(y11))], edgecolor = 'black', linewidth = 2)
    ax1.bar(x1, y15, width = bar_length, color = 'pink', bottom=[y11[i] + y12[i] + y13[i] + y14[i] for i in range(len(y11))], edgecolor = 'black', linewidth = 2)
    ax2.bar(x2, y2, width = bar_length, color = 'blue', edgecolor = 'black', linewidth = 2, log = True)
    ax1.set_facecolor(color='white')
    ax2.set_facecolor(color='white')
    fig.set_facecolor(color='white')
    plt.xticks([r + bar_length / 2 for r in range(len(y11))], graph_names)
    plt.title(title)
    plt.show()

