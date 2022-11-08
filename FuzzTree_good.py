import infrared as ir
#import infrared.rna as rna
#import matplotlib.pyplot as plt

#import re
#import collections
#import math
#import sys
#import os

#!export PYTHONPATH="${PYTHONPATH}:/home/uqamportable/miniconda3/lib/python3.9/site-packages/infrared/"
from infrared import def_constraint_class, def_function_class

#import pickle
import networkx as nx
#import varnaapi

from FuzzynessParameters import FuzzyParameters
DEBUG = 1

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
    


def DistanceEdgesLabel(x, y, Distance_matrix, len_edges_pattern, order=0):
    """
    Input : - Two edges labels x and y respectively from GP and GT.
               - A Distance Matrix between non-canonical interactions/labels (for instance the IDI matricx for isostericity distance).
               - len_edges_patterns to scale the impact of errors with the size of the graph pattern.
    Output : A value between 0 and 1 that indicates how much these labels are close from each other.
    """
    (number_BM_allowed, elim_iso_threshold, BM_iso_threshold, B53_missing_elim, allow_iso_nonBM, _, BM_by_gap, scale_BM_with_P_size) = FuzzyParameters()
    if scale_BM_with_P_size:
        number_BM_allowed = number_BM_allowed*len_edges_pattern/20
    (xx, yy) = interaction_to_number(x), interaction_to_number(y)
    if xx == yy: #including the case where xx == -1 (is a backbone)
        #return 1
        if xx != -1 or BM_by_gap == [] or (number_BM_allowed == 0):
            return 1 #not about backbone or no gap allowed
        BM_by_gap = [0] + BM_by_gap
        return 1 - BM_by_gap[order]/number_BM_allowed #equals to 1 for order 0
    if (number_BM_allowed == 0) and (Distance_matrix[xx][yy] >= min(BM_iso_threshold, elim_iso_threshold)):
        return 0
    if xx == -1 or yy == -1: 
        if B53_missing_elim or (number_BM_allowed == 0):
            return 0 #backbone not respected or BM not allowed
        return 1 - 1/number_BM_allowed
    if Distance_matrix[xx][yy] >= elim_iso_threshold: 
        return 0 #Too different 
    if Distance_matrix[xx][yy] >= BM_iso_threshold: 
        return 1 - 1/number_BM_allowed #Just a big mistake
    if allow_iso_nonBM and number_BM_allowed: #If small mistakes (non BM) are allowed
        return 1 - Distance_matrix[xx][yy]/(20*number_BM_allowed)
    return 1

def _EdgeLabelRespect(x, y, label_edge_ij, len_edges_pattern, nodes_target, label_edge_target, order_edge_target, edges_target, IDI_matrix, oriented):
    """
    Input : - Two nodes x and y from the GT respectively the mappings of i and j from the GP.
            - The label of edge between i and j, label_edge_ij.
            - len_edges_patterns to scale the impact of errors with the size of the graph pattern.
            - The list of the nodes in GT nodes_target to map the index of infrared functions with real names of nodes in the graph.
            - The list of all label of the edges in GT, label_edge_target.
            - The list of all edges in GT, edges_target.
            - A IDI Matrix giving distances between non-canonical interactions/labels (for instance the IDI matricx for isostericity distance).
            - oriented specifies if the graph is oriented (1 if it is the case)
    Output : A value between 0 and 1 that indicates if the edge is present or not but also how much these labels are close from each other.
    """
    #print("label_edge_target", label_edge_target)
    #print("order_edge_target", label_edge_target)
    (number_BM_allowed, _, _, _, _, BM_by_missing_edge, _, scale_BM_with_P_size) = FuzzyParameters()
    if scale_BM_with_P_size:
        number_BM_allowed = number_BM_allowed*len_edges_pattern/20
    if x == y: #If neighbors are mapped to the same node in GT, we can already reject 
        return 0 
        #TODO: Etendre ce critère d'injectivité à + ou - 1 le long du backbone ?
    if (not oriented):
        x, y =  min(x, y), max(x, y)
    if (nodes_target[x], nodes_target[y]) in edges_target: #We check that this edge is present or not in GT
        index = edges_target.index((nodes_target[x], nodes_target[y]))
        #print("index", index)
        return DistanceEdgesLabel(label_edge_ij,label_edge_target[edges_target.index((nodes_target[x], nodes_target[y]))], IDI_matrix, len_edges_pattern, order = order_edge_target[index])       
    if BM_by_missing_edge == -1 or number_BM_allowed == 0:
        return 0 #currently 0 here simply eliminates this pattern when the edge is not present
    return 1 - BM_by_missing_edge/number_BM_allowed




class EdgeLabelRespect(ir.infrared.WeightedFunction):
    """
    Constrain complementarity mapping of any pair of nodes (i,j) that are an edge in graph_pattern.
    ```
    Label(i, j,label_edge_ij, len_edges_pattern, nodes_target, label_edge_target, edges_target, IDI_matrix, oriented).
    ```
    The constraint is satisfied if label of the mapped edge (i,j) is an edge and if its label near or equal to the label of (i, j).
    """
def_function_class('EdgeLabelRespect', lambda i, j, label_edge_ij, len_edges_pattern, nodes_target, label_edge_target, order_edge_target, edges_target, IDI_matrix, oriented: [i, j],
            lambda x, y, label_edge_ij, len_edges_pattern, nodes_target, label_edge_target, order_edge_target, edges_target, IDI_matrix, oriented: _EdgeLabelRespect(x, y, label_edge_ij, len_edges_pattern, nodes_target, label_edge_target, order_edge_target, edges_target, IDI_matrix, oriented),
            module=__name__)


def main(GP, GT, nb_samples=1000): 
    """
    Input : - Two graphs, the pattern graph GP and the target graph GT.
    	    - nb_samples, maximum numbers of samples allowed to look for a pattern in GT.
    Output : List of "pre"mapping between GP and GT (as injectivity is not respected yet) with the numbers of nodes of GP that are correctly mapped
            The mapping is correct if this number is equals to n_pattern.
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
    
    (_, _, _, _, _, _, BM_by_gap, _) = FuzzyParameters()
    #TODO : Gap incomping
    #print("lenGTnodes", len(list(GT.nodes())))
    #print("lenGTedges", len(list(GT.edges())))
    GT = augment_graph(GT, len(BM_by_gap))
    n_pattern = len(GP.nodes)
    n_target = len(GT.nodes)
    #print("lenGTnodesnew", len(list(GT.nodes())))
    #print("lenGTedgesnew", len(list(GT.edges())))
    edges_pattern = list(GP.edges())
    edges_target = list(GT.edges())

    nodes_pattern = list(GP.nodes())
    nodes_target = list(GT.nodes())
    label_edge_pattern = [t['label'] for (_, _, t) in GP.edges.data()]
    label_edge_target = [t['label'] for (_, _, t) in GT.edges.data()]
    order_edge_target = [t['order'] for (_, _, t) in GT.edges.data()]
    #No label on vertices for now, can be introduced if needed.
    if DEBUG:
        print("Graphs Infos\n n_pattern:", n_pattern, "\n n_target:", 
    n_target,"\n edges_pattern:", edges_pattern,"\n edges_target:", 
    edges_target,"\n nodes_pattern:", nodes_pattern,"\n nodes_target:", 
    nodes_target,"\n label_edge_pattern:", label_edge_pattern,
    "\n label_edge_target:", label_edge_target, "\norder_edge_target:", order_edge_target)

   
    model = ir.Model(n_pattern, n_target)


    #We define Label respect function on each edges in the pattern to check if the label in the target correspond to it.
    model.add_functions([EdgeLabelRespect(nodes_pattern.index(i), nodes_pattern.index(j), label_edge_pattern[k], len(edges_pattern), nodes_target, label_edge_target, order_edge_target, edges_target, IDI_matrix, oriented=1) for k, (i,j) in enumerate(edges_pattern)], 'EdgeLabelRespect')
    #len(edges_pattern) here will allow us to scale the error depending on the size of the motif without touching the tolerancy indicated in the set_target part here.
    
    #We now take some samples of results given the built model.
    sampler = ir.Sampler(model)

    #Next we ensure that edges in the pattern are represented correctly for now (or enough in the future) in the target graph.
    sampler.set_target(1*len(edges_pattern), 0.999 , 'EdgeLabelRespect') 
    #With the tolerancy that is not scale on length of edges_pattern here, we ensure namely that if we put 0 as a return for EdgeLabelRespect then the pattern is automatically eliminated.
    #Other fuzzy "mistakes" must return values higher than 0 and "close" to 1 to stay in the long run as an admissible matching.
    #print("HELLO1")
    samples = [sampler.targeted_sample() for _ in range(nb_samples)]
    resu = [([(nodes_pattern[k], nodes_target[x]) for k,x in enumerate(sample.values())], len(list(set(sample.values())))) for sample in samples]
    #print("HELLO2")
    #TODO : Gap incomping
    resu = [r for r in resu if filter(GP, GT, r) == 1]
    return resu

def augment_graph(GT, number_gaps_allowed):
    """
    Input : - A graph GT that must be augmented with "false" edges between the B53 chains to account for gaps 
    	    - The length of BM_by_gap number_gaps_allowed which for each consecutive gap account for the cost in big mistakes to use "false" edges
    Output : - A modified GT graph that account gaps with "false edges" with two additionals labels
            * label "order" on edges that indicate the "depth of the gap
            * label "correspondant_nodes" that correspond to the list of shorcuted nodes. They will serve in the post procesing as shortcuted nodes must not be taken in mappings 
    """
    Gnew=nx.DiGraph()
    Gnew.add_nodes_from(GT)
    #print("len in process", len(GT.edges.data()))
    for (i,j,t) in GT.edges.data():
        Gnew.add_edge(i,j, label=t['label'], long_range=t['long_range'], near=t['near'], order=0, correspondant_nodes=[])
    #print("len in process", len(Gnew.edges.data()))
    for node in GT.nodes():
        iter_node = node
        correspondant_nodes = []
        B53_neighbors=[n for n in GT.successors(iter_node) if GT[iter_node][n]['label'] == 'B53']
        if len(B53_neighbors) > 1:
            print("THE WORLD BLOWS UP")
        if len(B53_neighbors) == 0:
            continue
        iter_node=B53_neighbors[0]
        correspondant_nodes.append(iter_node)
        for order in range(1,number_gaps_allowed+1):
            B53_neighbors=[n for n in GT.successors(iter_node) if GT[iter_node][n]['label'] == 'B53']
            if len(B53_neighbors) > 1:
                print("THE WORLD BLOWS UP")
            if len(B53_neighbors) == 0:
                break #no nodes after on this strand to carry on
            Gnew.add_edge(node, B53_neighbors[0], label='B53', long_range='False', near='False', order=order, correspondant_nodes=correspondant_nodes.copy())
            iter_node = B53_neighbors[0]
            correspondant_nodes.append(iter_node)
    return Gnew

def filter(GP, GTaugment, mapping):
    """
    Input: - A graph Pattern
           - A Target graph GTaugment augmented with edges to allow gaps
    Output: Return 1 if the mapping have effectively all correspondant_node unmapped which allows this mapping as feasible 
    """
    banned = []
    for k1,(p1,t1) in mapping:
        for k2, (p2, t2) in mapping:
            if k1 != k2:
                if (p1, p2) in GP.edges():
                    banned += GTaugment[t1][t2]['correspondant_nodes']
    banned = list(set(banned))
    GTmapped = [t for (_, t) in mapping]
    for b in banned:
        if b in GTmapped:
            return 0
    return 1


#import pickle
#with open("1Y27.nxpickle",'rb') as fT:
#    GT = pickle.load(fT)
#    print(GT.nodes.data())
#    print(GT.edges.data())
#GG = augment_graph(GT, 2)
#print(GG.nodes.data())
#print(GG.edges.data())

