import infrared as ir
#import infrared.rna as rna
#import matplotlib.pyplot as plt

#import re
#import collections
import math
#import sys
#import os

#!export PYTHONPATH="${PYTHONPATH}:/home/uqamportable/miniconda3/lib/python3.9/site-packages/infrared/"
from infrared import def_constraint_class, def_function_class

#import pickle
import networkx as nx
#import varnaapi

from FuzzynessParameters import FuzzyParameters
DEBUG = 1
infinite = 1000
#Des qu on autorise des fuzzy edges tout devient plus long on passe de 80 secondes à 250 pour un graphe où il ne manque aucun edge
#TODO: faut il precomputer toutes les distances ?
#TODO: WARNING gros probleme avec l infini.... l'infini est une constante ? (locale ou globale ?) #Autre probleme si l'infini est 1000000 -> inconsistancy error mais 100000 converge par exemple Mais si l'infini est placé à 10 c'est plus long 114 secondes contre 77 environ dans le cas où infinite = 1000 ou 10000 81 secs pour 100
#TODO: WARNING Infrared (ou moi) semble définitivement avoir un problème avec les valeurs grandes j'ai pour le moment quotienté par 100 les thresholds et valeurs retournées en conséquence???
#Ilse peut que ce soit tout simplement un problme de . flotant qui fait qu il retourne une difference de matricce inexac te ou alors il cherche tou les 10 puissance - qqc ce qui lu fait trop grand a explorere avec des entiers assez grand
def _EdgeRespect(x, y, nodes_target, edges_target, GT, B, D):
    """
    Input : - Two nodes x and y from the GT respectively the mappings of i and j from the GP.
            - The list of the nodes in GT nodes_target to map the index of infrared functions with real names of nodes in the graph.
            - The list of all edges in GT, edges_target.
    Output : A value between 0 and +inf that indicates if the edge is present, returns 0 if it is the case, 1 if it is absent and +inf if it is eliminatory
    """
    #if x == y: #If neighbors are mapped to the same node in GT, we can already reject 
    #    return math.inf#math.inf 
    if (nodes_target[x], nodes_target[y]) in edges_target: #We check that this edge is present or not in GT
        return 0 #Change perhaps when it is a backbone ?
    if B == 0 or (distance(nodes_target[x], nodes_target[y], GT) > D):
        return infinite#math.inf #To be modified to allow fuzzy edges
    return 1/B


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
    #if x == y: #If neighbors are mapped to the same node in GT, we can already reject 
    #    return math.inf#math.inf 
    if (nodes_target[x], nodes_target[y]) in edges_target: #We check that this edge is present or not in GT as no label to check if not present
        (xx, yy) = interaction_to_number(label_edge_ij), interaction_to_number(label_edge_target[edges_target.index((nodes_target[x], nodes_target[y]))])
        if xx == yy: 
            return 0 #1 #0
        if xx == -1 or yy == -1:
            return infinite#math.inf#0#math.inf #backbone replaces by something else, it is eliminatory
        #print("\nIDI", IDI_matrix[xx][yy] - IDI_matrix[xx][xx])
        return IDI_matrix[xx][yy]#(IDI_matrix[xx][yy] - IDI_matrix[xx][xx])/100
    return 0#1 #0

def _GapRespect(x, y, nodes_target, edges_target, GT, A):
    """
    Input : - Two nodes x and y from the GT respectively the mappings of i and j from the GP.
            - The list of the nodes in GT nodes_target to map the index of infrared functions with real names of nodes in the graph.
            - The list of all edges in GT, edges_target.
    Output : A value between 0 and +inf that indicates  the sum of all "false" edge use distances.
    """
    #if x == y: #If neighbors are mapped to the same node in GT, we can already reject 
    #    return math.inf#math.inf 
    if (nodes_target[x], nodes_target[y]) in edges_target: #We check that this edge is present or not in GT
        if GT[nodes_target[x]][nodes_target[y]]['dist'] != 0:
            if A == 0:
                return infinite
            return GT[nodes_target[x]][nodes_target[y]]['dist']/A
    return 0
class EdgeRespect(ir.infrared.WeightedFunction):
    """
    Constrain complementarity mapping of any pair of nodes (i,j) that are an edge in graph_pattern.
    ```
    EdgeRespect(i, j,label_edge_ij, nodes_target, edges_target).
    The constraint is satisfied if the edge (x, y) mapped to (i,j) is an edge in GT.
    """
def_function_class('EdgeRespect', lambda i, j, nodes_target, edges_target, GT, B, D: [i, j],
            lambda x, y, nodes_target, edges_target, GT, B, D: _EdgeRespect(x, y, nodes_target, edges_target, GT, B, D),
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
    Constrain the length of any pair of nodes (x,y) that are mapped to a "false" edge in graph_target.
    ```
    GapRespect(i, j,label_edge_ij, nodes_target, label_edge_target, edges_target, IDI_matrix).
    ```
    The constraint report the cost of the use of such "false edge".
    """
def_function_class('GapRespect', lambda i, j, nodes_target, edges_target, GT, A: [i, j],
            lambda x, y, nodes_target, edges_target, GT, A: _GapRespect(x, y, nodes_target, edges_target, GT, A),
            module=__name__)


def main(GP, GT, E, B, A, maxGAPdistance=3, nb_samples=1000): 
    """
    Input : - Two graphs, the pattern graph GP and the target graph GT.
            - E, the threshold in term of sum of isostericity allowed
            - B, the threshold in term of number of missing edges allowed
            - A, the threshold in term of sum in Angstrom of gaps allowed
            - maxGAPdistance the maximal distance after which we are not looking for gap anymore
    	    - nb_samples, maximum numbers of samples allowed to look for a pattern in GT.
    Output : List of mapping between GP and GT 
        As injectivity is not respected yet and with gaps conditions, the results are preprocessed before being output to return only valid mapping. 
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

    #We enrich the target Graph with False Edges that account for gaps
    print("GTnodes", len(GT.nodes()))
    GT = augment_graph(GT, maxGAPdistance)
    print("GTnodesbis", len(GT.nodes()))
    if E != 0:
        for i in range(len(IDI_matrix)):
            storage = IDI_matrix[i][i]
            for j in range(len(IDI_matrix[0])):
                IDI_matrix[i][j] = (IDI_matrix[i][j]- storage)/E
    n_pattern = len(GP.nodes)
    n_target = len(GT.nodes)

    edges_pattern = list(GP.edges())
    edges_target = list(GT.edges())

    nodes_pattern = list(GP.nodes())
    nodes_target = list(GT.nodes())

    label_edge_pattern = [t['label'] for (_, _, t) in GP.edges.data()]
    label_edge_target = [t['label'] for (_, _, t) in GT.edges.data()]

    #No label on vertices for now, can be introduced if needed.
    if DEBUG:
        print("Graphs Infos\n n_pattern:", n_pattern, "\n n_target:", 
    n_target,"\n edges_pattern:", edges_pattern,"\n edges_target:", 
    edges_target,"\n nodes_pattern:", nodes_pattern,"\n nodes_target:", 
    nodes_target,"\n label_edge_pattern:", label_edge_pattern,
    "\n label_edge_target:", label_edge_target)

   
    model = ir.Model(n_pattern, n_target)

    #We define Edge respect function on each edges in the pattern to check if the mapping of this couple is effectively an edge in the target.
    model.add_functions([EdgeRespect(nodes_pattern.index(i), nodes_pattern.index(j), nodes_target,  edges_target, GT, B, maxGAPdistance) for k, (i,j) in enumerate(edges_pattern)], 'EdgeRespect')
    #For now max distance to consider missing edge is consider as max distance with which we consider gap
    #We define Label respect function on each edges in the pattern to check if the label in the target correspond to it.
    model.add_functions([LabelRespect(nodes_pattern.index(i), nodes_pattern.index(j), label_edge_pattern[k], nodes_target, label_edge_target,  edges_target, IDI_matrix) for k, (i,j) in enumerate(edges_pattern)], 'LabelRespect')

    #We define Label respect function on each edges in the pattern to check if the label in the target correspond to it.
    model.add_functions([GapRespect(nodes_pattern.index(i), nodes_pattern.index(j), nodes_target,  edges_target, GT, A) for k, (i,j) in enumerate(edges_pattern)], 'GapRespect')
   
    #We now take some samples of results given the built model.
    sampler = ir.Sampler(model)

    #Next we ensure that edges in the pattern are represented enough
    sampler.set_target(0, 1, 'EdgeRespect') 
    #Next we ensure that the labels are not too far in term of isostericity.
    sampler.set_target(0, 1, 'LabelRespect') 
    #sampler.set_target(1*len(edges_pattern), 0.001, 'LabelRespect') 
    #Finally we ensure that the sum of gaps does not exceed a certain Angstrom value.
    sampler.set_target(0, 1, 'GapRespect') 

    samples = [sampler.targeted_sample() for _ in range(nb_samples)]

    #We postprocess here by checking number of nodes mapped to same nodes in GT in order to obtain the injectivity
    resu = [([(nodes_pattern[k], nodes_target[x]) for k,x in enumerate(sample.values())], len(list(set(sample.values())))) for sample in samples]
    
    #We postprocess now the gap procedure as shorcuted nodes must be left unaffected by the mapping.
    resu = [(r, r2) for r, r2 in resu if filter(GP, GT, r) == 1]

    return resu

def augment_graph(GT, maxGapallowed):
    """
    Input : - A graph GT that must be augmented with "false" edges between the B53 chains to account for gaps 
    	    - The length of BM_by_gap number_gaps_allowed which for each consecutive gap account for the cost in big mistakes to use "false" edges
    Output : - A modified GT graph that account gaps with "false edges" with one additional label
            * label "correspondant_nodes" that correspond to the list of shorcuted nodes. They will serve in the post procesing as shortcuted nodes must not be taken in mappings 
    """
    Gnew=nx.DiGraph()
    for (i,t) in GT.nodes.data():
        Gnew.add_node(i, pdb_position = t['pdb_position'], atoms = t['atoms'])
    for (i,j,t) in GT.edges.data():
        Gnew.add_edge(i,j, label=t['label'], long_range=t['long_range'], near=t['near'], correspondant_nodes=[], dist=0)
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
        first_dist = distance(node, iter_node, GT)
        dist = first_dist
        while dist < maxGapallowed: #can be replaced here by  while dist - first_dist < maxGapallowed: #if we only want to take into account everything but the first gap distance
            B53_neighbors=[n for n in GT.successors(iter_node) if GT[iter_node][n]['label'] == 'B53']
            if len(B53_neighbors) > 1:
                print("THE WORLD BLOWS UP")
            if len(B53_neighbors) == 0:
                break #no nodes after on this strand to carry on
            iter_node = B53_neighbors[0]
            dist = distance(node, iter_node, GT)
            if dist < maxGapallowed: #if dist - first_dist < maxGapallowed:
                Gnew.add_edge(node, iter_node, label='B53', long_range='False', near='False', correspondant_nodes=correspondant_nodes.copy(), dist = max(0, dist - first_dist))
            correspondant_nodes.append(iter_node)
    return Gnew

def filter(GP, GTaugment, mapping):
    """
    Input: - A graph Pattern
           - A Target graph GTaugment augmented with edges to allow gaps
    Output: Return 1 if the mapping have effectively all correspondant_node unmapped which allows this mapping as feasible 
    """
    banned = []
    for k1,(p1,t1) in enumerate(mapping):
        for k2, (p2, t2) in enumerate(mapping):
            if k1 != k2:
                if (p1, p2) in GP.edges() and (t1, t2) in GTaugment.edges():
                    banned += GTaugment[t1][t2]['correspondant_nodes']
    banned = list(set(banned))
    GTmapped = [t for (_, t) in mapping]
    for b in banned:
        if b in GTmapped:
            return 0
    return 1

def preL2distance(x,y,z):
    return x**2 + y**2 + z**2
 
def distance(node1, node2, GT):
    li1 = [t['atoms'] for i,t in GT.nodes.data() if i == node1][0]
    li2 = [t['atoms'] for i,t in GT.nodes.data() if i == node2][0]
    dist = math.inf
    for atom1 in li1:
        for atom2 in li2:
            (x1, y1, z1) = atom1['position']
            (x2, y2, z2) = atom2['position']
            (x1, y1, z1) = (float(x1), float(y1), float(z1))
            (x2, y2, z2) = (float(x2), float(y2), float(z2))
            dist = min(dist, preL2distance(x1 - x2, y1 - y2, z1 - z2))
    return math.sqrt(dist)
#import pickle
#with open("1Y27.nxpickle",'rb') as fT:
#    GT = pickle.load(fT)
#    print(GT.nodes.data())
#    print(GT.edges.data())
#print(distance(('X', 46), ('X', 55), GT))
#GG = augment_graph(GT, 3)
#print(GG.nodes.data())
#print(GG.edges.data())

