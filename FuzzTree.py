# FuzzTree
# Copyright (C) 2023 THEO BOURY 

import math
import networkx as nx

#!export PYTHONPATH="${PYTHONPATH}:/home/uqamportable/miniconda3/lib/python3.9/site-packages/infrared/"
import infrared as ir
from infrared import def_constraint_class, def_function_class
from multiprocessing import Pool

DEBUG = 0
infinite = 100

def check(GP, GT, Mapping, L, E, G, D, IDI_matrix, Distancer):
    """ 
    Input : - Pattern and target graphs GP and GT
            - Fuzzy thresholds L, E, G, D
            - A mapping found in GT for GP
            - IDI_matrix for label isostericity computation and Distancer that stores distances between edges of GT.
    Output : Returns 1 if the Mapping is admissible and respects all of the fuzzy thresholds and 0 otherwise.
    """
    list_edges_GT = list(GT.edges())
    mapped = {}
    for (i, j) in Mapping:
        mapped[i] = j
    sum_label = 0
    sum_edge_missing = 0
    sum_dist_gap = 0
    for (i, j, t) in GP.edges.data():
        if (mapped[i], mapped[j]) not in list_edges_GT:
            if t['label'] == 'B53':
                return 0
            if Distancer(mapped[i], mapped[j]) > D:
                return 0
            sum_edge_missing+=1
        else:
            (ii, jj, tt) = [(ii, jj, tt) for (ii, jj, tt) in GT.edges.data() if ii == mapped[i] and jj == mapped[j]][0]
            if t['label'] != tt['label'] and (t['label'] == 'B53' or tt['label'] == 'B53'): #To avoid problem with B53 not in IDI matrix
                return 0
            elif t['label'] != tt['label']:
                sum_label += IDI_matrix[t['label']][tt['label']]
            sum_dist_gap += tt['dist'] #Can be done in any case as this label is set to 0 when it is not a "false" edge.
    if (sum_label > L) or (sum_edge_missing > E) or (sum_dist_gap > G):
        return 0
    return 1

def interaction_to_number(interact_char):
    """ 
    Input : A non-canonical interaction/label as a string, interact_char.
    Output : The corresponding index of this interaction in the Distance_matrix/ IDI_matrix 
            and -1 if it is a backbone.
    """
    if interact_char == 'B53':
        return -1
    if interact_char[2] == 'W':
        interact_char = interact_char[0] + interact_char[2] + interact_char[1]
    if interact_char[1] == 'S' and interact_char[2] == 'H':
        interact_char = interact_char[0] + 'HS'
    traduction = ['CHH', 'TWH', 'CWW', 'THS', 'CWS', 'CSS', 'CWH', 'CHS','TWS','TSS','TWW','THH']
    return traduction.index(interact_char)


def augment_graph(GT, maxGapallowed, Distancer):
    """
    Input : - A graph GT that must be augmented with "false" edges between the B53 chains to account for gaps. 
    	    - The length of BM_by_gap number_gaps_allowed which for each consecutive gap account for the cost in big mistakes to use "false" edges.
    Output : - A modified GT graph that account gaps with "false edges" with two additional labels.
            * label "dist", the distance between the two nodes (i,j) in the edges minus the distance between first node and its neighbor link with a B53:
                dist = D(i,j) - D(i, N_{B53}(i))
            * label "correspondant_nodes" that correspond to the list of shorcuted nodes. They will serve in the post procesing as shortcuted nodes must not be taken in mappings. 
    """
    Gnew=nx.DiGraph() #Initiate the new GT graph.
    for (i,t) in GT.nodes.data():
        Gnew.add_node(i, atoms = t['atoms'])
    for (i,j,t) in GT.edges.data():
        Gnew.add_edge(i,j, label=t['label'], near=t['near'], correspondant_nodes=[], dist=0)
    for node in GT.nodes():
        iter_node = node
        correspondant_nodes = []
        B53_neighbors=[n for n in GT.successors(iter_node) if GT[iter_node][n]['label'] == 'B53']
        if len(B53_neighbors) > 1: #It means that two backbones start from iter_node, which is not biologically admissible.
            print("THE WORLD BLOWS UP")
        if len(B53_neighbors) == 0: #It means that we are at the end of the strand already
            continue
        iter_node=B53_neighbors[0]
        correspondant_nodes.append(iter_node)
        first_dist = Distancer[node][iter_node]
        dist = first_dist

        while dist < maxGapallowed: #We do not take into account nodes that are too far from the origin node. 
        #Can be replaced here by  "while dist - first_dist < maxGapallowed:", if we only want to take into account everything except the first gap distance.
            B53_neighbors=[n for n in GT.successors(iter_node) if GT[iter_node][n]['label'] == 'B53']
            if len(B53_neighbors) > 1:
                print("THE WORLD BLOWS UP")
            if len(B53_neighbors) == 0:
                break #no nodes after, on this strand, to carry on.
            iter_node = B53_neighbors[0]
            dist = distance(node, iter_node, GT)
            if dist < maxGapallowed:
                Gnew.add_edge(node, iter_node, label='B53', near='False', correspondant_nodes=correspondant_nodes.copy(), dist = max(0, dist - first_dist))
            correspondant_nodes.append(iter_node)
    return Gnew


def filter(GP, GTaugment, mapping):
    """
    Input: - A graph Pattern GP
           - A Target graph GTaugment augmented with "false" edges to allow gaps
    Output: Return 1 if the mapping have effectively all correspondant_node unmapped which allows this mapping to be feasible. 
    """
    banned = []
    for k1,(p1,t1) in enumerate(mapping):
        for k2, (p2, t2) in enumerate(mapping):
            if k1 != k2:
                if (p1, p2) in GP.edges() and (t1, t2) in GTaugment.edges():
                    banned += GTaugment[t1][t2]['correspondant_nodes'] #If edges used in the mapped graph have correspondants, they must be unmapped.
    banned = list(set(banned))
    GTmapped = [t for (_, t) in mapping]
    for b in banned:
        if b in GTmapped:
            return 0
    return 1

 
def preL2distance(triplet1, triplet2):
    """
        Compute square of L2 distances between the three coordinates triplet1=(x1,y1,z1) and triplet2=(x2,y2,z2).
    """
    (x1,y1,z1) = triplet1
    (x2,y2,z2) = triplet2
    (x1,y1,z1) = float(x1), float(y1), float(z1)
    (x2,y2,z2) = float(x2), float(y2), float(z2)
    return (x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2


def distance(node1, node2, GT):
    """
    Input : - Two nodes node1 and node2 in GT
            - Graph Target GT
    Output :
            - The minimal distance between node1 and node2 in GT by considering the atoms that are the more close to each other.
    """
    # distance here is computed on the base of the two closests atoms and not on "centrality" of atoms.
    li1 = [t['atoms'] for i,t in GT.nodes.data() if i == node1][0]
    li2 = [t['atoms'] for i,t in GT.nodes.data() if i == node2][0]
    dist = math.inf
    for atom1 in li1:
        for atom2 in li2:
            dist = min(dist, preL2distance(atom1['position'], atom2['position']))
    return math.sqrt(dist)


def wrapper_distance(node1_GT_k1_li):
    """
    A wrapper to compute distance with multiprocessing between each pairs of nodes
    """
    (node1, GT, k1, li) = node1_GT_k1_li
    resu = []
    for k2, node2 in enumerate(li):
        if k2 >= k1:
            value = distance(node1, node2, GT)
            resu.append((node1, node2, value))
    return resu

def precompute_distance(GT, nb_procs):
    """
    Input : - The target graph GT on which we preprocess the distances between each node.
            - The number of processors nb_procs to multiprocess the distance precomputing.
    Output : Return the Distancer dictionnary that contains dictionnary such as 
             Distancer[node1][node2] is the distance in Angstrom between the two nodes.
    """
    Distancer = {}
    entry = []
    li = list(GT.nodes())
    for k1, node1 in enumerate(li):
            entry.append((node1, GT, k1, li))
    if nb_procs == 1:
        resu = []
        for elem in entry:
            resu.append(wrapper_distance(elem))
    else:
        with Pool(nb_procs) as pool:
            resu= list(pool.imap_unordered(wrapper_distance, entry))
    for li in resu:
        for (node1, node2, value) in li:
            if node1 not in Distancer.keys():
                Distancer[node1] = {}
            if node2 not in Distancer.keys():
                Distancer[node2] = {}
            Distancer[node1][node2] = value
            Distancer[node2][node1] = value
    if DEBUG:
        print("Distancer precomputing done\n")
    return Distancer

def _EdgeRespect(x, y, label_edge_ij, nodes_target, edges_target, Distancer, E, D):
    """
    Input : - Two nodes x and y from the GT respectively the mappings of i and j from the GP.
            - The list of the nodes in GT nodes_target to map the index of infrared functions with real names of nodes in the graph.
            - The list of all edges in GT, edges_target.
            - Storage of the Target graph GT to compute distance.
            - E, the threshold on number of missing edges.
            - D, the maximum distance above which we do not consider edges that are missing
    Output : A value between 0 and +inf that indicates if the edge is present, returns 0 if it is the case, 1/E if it is absent and +inf if it is eliminatory
    """
    if x == y: #Necessary condition of Injectivity : If neighbors are mapped to the same node in GT, we can already reject  
        return infinite
    if (nodes_target[x], nodes_target[y]) in edges_target: #We check that this edge is present or not in GT
        return 0
    if E == 0 or label_edge_ij == 'B53' or (Distancer[nodes_target[x]][nodes_target[y]] > D):
        return infinite
    return 1


def _LabelRespect(x, y, label_edge_ij, nodes_target, label_edge_target, edges_target, IDI_matrix):
    """
    Input : - Two nodes x and y from the GT respectively the mappings of i and j from the GP.
            - The label of edge between i and j, label_edge_ij.
            - The list of the nodes in GT nodes_target to map the index of infrared functions with real names of nodes in the graph.
            - The list of all label of the edges in GT, label_edge_target.
            - The list of all edges in GT, edges_target.
            - A IDI Matrix giving distances between non-canonical interactions/labels (for instance the IDI matricx for isostericity distance).
    Output : A value between 0 and +inf that indicates  how much labels are close from each other.
    """
    if (nodes_target[x], nodes_target[y]) in edges_target: #We check that this edge is present or not in GT as no label must be checked otherwise.
        (xx, yy) = interaction_to_number(label_edge_ij), interaction_to_number(label_edge_target[edges_target.index((nodes_target[x], nodes_target[y]))])  #Convert canonical labels to numbers.
        if xx == yy: #it includes the case where xx == yy == -1.
            return 0 
        if xx == -1 or yy == -1:
            return infinite #Backbone is replaced by something else, it is eliminatory.
        return IDI_matrix[xx][yy] 
    return 0 #Nothing to add if we do not talk about an edge in GT, eliminatory case for couple not in GT alreay in _EdgeRespect function.


def _GapRespect(x, y, nodes_target, edges_target, GT, G):
    """
    Input : - Two nodes x and y from the GT respectively the mappings of i and j from the GP.
            - The list of the nodes in GT nodes_target to map the index of infrared functions with real names of nodes in the graph.
            - The list of all edges in GT, edges_target.
            - Storage of the Target graph GT to compute distance.
            - A, the maximal sum of distances of gaps allowed in total.
    Output : A value between 0 and +inf that indicates the sum of all "false" edges (gaps) distances. Return +inf for eliminatory condition.
    """
    if (nodes_target[x], nodes_target[y]) in edges_target: 
        if GT[nodes_target[x]][nodes_target[y]]['correspondant_nodes'] != []: #This field is filled only ig it represents a gap.
            if G == 0: #If no gap allowed.
                return infinite
            return GT[nodes_target[x]][nodes_target[y]]['dist']
    return 0  


class EdgeRespect(ir.infrared.WeightedFunction):
    """
    Constrain complementarity mapping of any pair of nodes (i,j) an check if it is an edge in graph_target.
    ```
    EdgeRespect(i, j, label_edge_ij, nodes_target, edges_target, GT, E, D).
    The constraint is satisfied if the edge (x, y) mapped to (i,j) is an edge in GT.
    """
def_function_class('EdgeRespect', lambda i, j, label_edge_ij, nodes_target, edges_target, Distancer, E, D: [i, j],
            lambda x, y, label_edge_ij, nodes_target, edges_target, Distancer, E, D: _EdgeRespect(x, y, label_edge_ij, nodes_target, edges_target, Distancer, E, D),
            module=__name__)


class LabelRespect(ir.infrared.WeightedFunction):
    """
    Constrain complementarity mapping of any pair of nodes (i,j) that are an edge in graph_pattern.
    ```
    LabelRespect(i, j,label_edge_ij, nodes_target, label_edge_target, edges_target, IDI_matrix).
    ```
    The constraint is satisfied if label of the edge (x, y) mapped to (i,j) has label near or equals to the label of (i, j).
    """
def_function_class('LabelRespect', lambda i, j, label_edge_ij, nodes_target, label_edge_target, edges_target, IDI_matrix: [i, j],
            lambda x, y, label_edge_ij, nodes_target, label_edge_target, edges_target, IDI_matrix: _LabelRespect(x, y, label_edge_ij, nodes_target, label_edge_target, edges_target, IDI_matrix),
            module=__name__)


class GapRespect(ir.infrared.WeightedFunction):
    """
    Constrain the length of any pair of nodes (x,y) that are mapped to a "false" edge/ gap in graph_target.
    ```
    GapRespect(i, j,  nodes_target, edges_target, IDI_matrix, GT, G).
    ```
    The constraint report the cost of the use of such "false edge".
    """
def_function_class('GapRespect', lambda i, j, nodes_target, edges_target, GT, G: [i, j],
            lambda x, y, nodes_target, edges_target, GT, G: _GapRespect(x, y, nodes_target, edges_target, GT, G),
            module=__name__)


def main(GP, GT, L, E, G, maxGAPdistance=3, nb_samples=1000, respect_injectivity=1, D = 5, Distancer_preprocessed = {}, nb_procs = 1): 
    """
    Input : - Two graphs, the pattern graph GP and the target graph GT.
            - L, the threshold in term of sum of isostericity allowed
            - E, the threshold in term of number of missing edges allowed
            - G, the threshold in term of sum in Angstrom of gaps allowed
            - maxGAPdistance the maximal distance after which we are not looking for gap anymore
    	    - nb_samples, maximum numbers of samples allowed to look for a pattern in GT.
            - respect_injectivity is a boolean, if set to true we filter patterns that are not injective by postprocessing
    Output : List of (pre)mapping between GP and GT and the number of nodes that are effectively mapped.
        As injectivity is not respected yet and with gaps conditions, the results are preprocessed before being output to return only valid premapping. 
    """

    #The matrix for the isostericity between non-canonical interactions
    IDI_matrix = [[8.9, 12.0, 14.7, 14.0, 13.7, 12.7, 15.1, 14.7, 16.2, 16.6, 16.2, 14.0], 
    [12.0, 2.6, 10.6, 9.7, 14.3, 15.6, 11.2, 15.2, 13.8, 15.4, 11.9, 11.4], 
    [14.7, 10.6, 4.1, 8.2, 9.2, 13.1, 14.5, 16.0, 12.4, 11.3, 11.1, 15.5], 
    [14.0, 9.7, 8.2, 2.1, 7.0, 12.7, 12.0, 12.1, 10.0, 11.9, 13.1, 15.8], 
    [13.7, 14.3, 9.2, 7.0, 3.5, 7.4, 14.9, 12.3, 10.9, 10.8, 14.6, 17.7], 
    [12.7, 15.6, 13.1, 12.7, 7.4, 1.3, 15.8, 12.9, 13.8, 12.0, 17.1, 19.0], 
    [15.1, 11.2, 14.5, 12.0, 14.9, 15.8, 3.2, 8.8, 8.4, 11.5, 10.6, 10.8], 
    [14.7, 15.2, 16.0, 12.1, 12.3, 12.9, 8.8, 2.4, 7.9, 11.2, 14.7, 14.9], 
    [16.2, 13.8, 12.4, 10.0, 10.9, 13.8, 8.4, 7.9, 3.4, 6.4, 9.6, 14.4], 
    [16.6, 15.4, 11.3, 11.9, 10.8, 12.0, 11.5, 11.2, 6.4, 2.2, 9.0, 14.4], 
    [16.2, 11.9, 11.1, 13.1, 14.6, 17.1, 10.6, 14.7, 9.6, 9.0, 3.8, 9.0],
    [14.0, 11.4, 15.5, 15.8, 17.7, 19.0, 10.8, 14.9, 14.4, 14.4, 9.0, 4.0]]

    if Distancer_preprocessed != {}:
        Distancer = Distancer_preprocessed
    else:
        Distancer = precompute_distance(GT, nb_procs)

    #We enrich the target Graph with False Edges that account for gaps
    GT = augment_graph(GT, maxGAPdistance, Distancer)
    #We modify IDI_matrix by the isostericity to avoid taking into account isostericy interned to families.
    if L != 0:
        for i in range(len(IDI_matrix)):
            storage = IDI_matrix[i][i]
            for j in range(len(IDI_matrix[0])):
                IDI_matrix[i][j] = (IDI_matrix[i][j]- storage)
    #Intialisation.
    n_pattern = len(GP.nodes)
    n_target = len(GT.nodes)

    edges_pattern = list(GP.edges())
    edges_target = list(GT.edges())

    nodes_pattern = list(GP.nodes())
    nodes_target = list(GT.nodes())

    label_edge_pattern = [t['label'] for (_, _, t) in GP.edges.data()]
    label_edge_target = [t['label'] for (_, _, t) in GT.edges.data()]
    #No label on vertices for now, can be introduced if needed.
    if DEBUG: #To print graph infos at each call to main
        print("Graphs Infos\n n_pattern:", n_pattern, "\n n_target:", 
    n_target,"\n edges_pattern:", edges_pattern,"\n edges_target:", 
    edges_target,"\n nodes_pattern:", nodes_pattern,"\n nodes_target:", 
    nodes_target,"\n label_edge_pattern:", label_edge_pattern,
    "\n label_edge_target:", label_edge_target)


    model = ir.Model(n_pattern, n_target)

    #We define Edge respect function on each edges in the pattern to check if the mapping of this couple is effectively an edge in the target.
    model.add_functions([EdgeRespect(nodes_pattern.index(i), nodes_pattern.index(j), label_edge_pattern[k], nodes_target,  edges_target, Distancer, E, D) for k, (i,j) in enumerate(edges_pattern)], 'EdgeRespect')
    #For now max distance to consider missing edge is consider as max distance with which we consider gap 

    #We define Label respect function on each edges in the pattern to check if the label in the target correspond to it.
    model.add_functions([LabelRespect(nodes_pattern.index(i), nodes_pattern.index(j), label_edge_pattern[k], nodes_target, label_edge_target,  edges_target, IDI_matrix) for k, (i,j) in enumerate(edges_pattern)], 'LabelRespect')

    #We define Gap respect function on each edges in the pattern to check if we mapped "false edges" that represent gaps and we take into account the distance between nodes that they induced.
    model.add_functions([GapRespect(nodes_pattern.index(i), nodes_pattern.index(j), nodes_target,  edges_target, GT, G) for k, (i,j) in enumerate(edges_pattern)], 'GapRespect')
   
    #We now take some samples of results given the built model.
    sampler = ir.Sampler(model)

    #Next we ensure that edges in the pattern are represented enough
    sampler.set_target(E/2, E/2, 'EdgeRespect') 
    #Next we ensure that the labels are not too far in term of isostericity.
    sampler.set_target(L/2, L/2, 'LabelRespect') 
    #Finally we ensure that the sum of gaps does not exceed a certain Angstrom value.
    sampler.set_target(G/2, G/2, 'GapRespect')
    if DEBUG:
        print("treewidth", sampler.treewidth())
    samples = []
    for _ in range(nb_samples):
        try:
            samples.append(sampler.targeted_sample())
        except:
            if DEBUG:
                print("Inconsistancy error detected ! No solution computable in this case.")
            return []
    resu = [([(nodes_pattern[k], nodes_target[x]) for k,x in enumerate(sample.values())], len(list(set(sample.values())))) for sample in samples]
    
    #We postprocess now the gap procedure as shorcuted nodes must be left unaffected by the mapping.
    resu = [(r, r2) for r, r2 in resu if filter(GP, GT, r) == 1]

    #We postprocess here by checking number of nodes mapped to same nodes in GT in order to obtain the injectivity.
    if respect_injectivity:
        resu = [(r, length) for r, length  in resu if length == n_pattern]

    return resu


