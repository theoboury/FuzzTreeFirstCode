import infrared as ir
import infrared.rna as rna
import matplotlib.pyplot as plt

import re
import collections
import math
import sys
import os
import infrared
import pickle
#!export PYTHONPATH="${PYTHONPATH}:/home/uqamportable/miniconda3/python3.9/site-packages/infrared/"
from infrared import def_constraint_class, def_function_class
import networkx as nx
import graphviz as gv

DEBUG = 0

def init_from_graph(graph_pattern, label_pattern, graph_target, label_target, oriented):
    n_pattern = len(graph_pattern)
    n_target = len(graph_target)
    if oriented:
        edges_pattern = [(i, j) for i in range(n_pattern) for j in range(n_pattern) if j in graph_pattern[i]]
        edges_target = [(i, j) for i in range(n_target) for j in range(n_target) if j in graph_target[i]]
    else:
        edges_pattern = [(i, j) for i in range(n_pattern) for j in range(n_pattern) if j in graph_pattern[i] and i < j]
        edges_target = [(i, j) for i in range(n_target) for j in range(n_target) if j in graph_target[i] and i < j]
    return n_pattern, edges_pattern, n_target, edges_target

#WARNING : Should we define another injectivity criteria ? Or is a restiction of the entry couple enough (for now we can focus on consecutive backbone for instance !)

def _Injectivity(x, y):
    if DEBUG:
        print("Positions of mapping tried by Injectivity", x, y)
    if x == y:
        return 0
    else:
        if DEBUG:
            print("that succeded !")
        return 1


class Injectivity(infrared.Constraint):
    """
    Constrain complementarity mapping of any pair of nodes (i,j)

    ```
    Injectivity(i, j)
    ```
    The constraint is satisfied if values of the ma#def DistanceEdgesLabel(x, y, Distance_matrix, len_edges_pattern):
#    """
#    Input : - Two edges labels x and y respectively from GP and GT.
#               - A Distance Matrix between non-canonical interactions/labels (for instance the IDI matricx for isostericity distance).
#               - len_edges_patterns to scale the impact of errors with the size of the graph pattern.
#    Output : A value between 0 and 1 that indicates how much these labels are close from each other.
#    """
#    (number_BM_allowed, elim_iso_threshold, BM_iso_threshold, B53_missing_elim, allow_iso_nonBM, _, BM_by_gap, scale_BM_with_P_size) = FuzzyParameters()
#    if scale_BM_with_P_size:
#        number_BM_allowed = number_BM_allowed*len_edges_pattern/20
#    (xx, yy) = interaction_to_number(x), interaction_to_number(y)
#    if xx == yy: #including the case where xx == -1 (is a backbone)
#        return 1
#    if xx == -1 or yy == -1: 
#        if B53_missing_elim:
#            return 0 #backbone not respected
#        else:
#            if number_BM_allowed:
#                return 1 - 1/number_BM_allowed
#            else:
#                return 0
#    if Distance_matrix[xx][yy] >= elim_iso_threshold: 
#        return 0 #Too different 
#    if Distance_matrix[xx][yy] >= BM_iso_threshold: 
#        if number_BM_allowed:
#            return 1 - 1/number_BM_allowed #Just a big mistake
#        else:
#            return 0
#    if allow_iso_nonBM and number_BM_allowed: #If small mistakes (non BM) are allowed
#        return 1 - Distance_matrix[xx][yy]/(20*number_BM_allowed)
#    return 1pping for i and j are distincts.
    """
def_constraint_class('Injectivity', lambda i, j: [i, j],
                     lambda x, y: _Injectivity(x, y),
                     module=__name__)



#def _EdgeRespect(x, y, edges_target):
#    if DEBUG:
#        print("Positions of mapping tried by EdgeRespect", x, y)
#    if (x, y) in edges_target or (y, x) in edges_target:
#        if DEBUG:
#            print("that succeded !")
#        return 1
#    else:
#        return 0


def _LabelRespect(x, label_i, label_target):
    if DEBUG:
        print("Positions of mapping tried by LabelRespect", x)
    if (label_i == label_target[x]):
        if DEBUG:
            print("that succeded !")
        return 1
    else:
        return 0

def interaction_to_number(interact_char):
    if interact_char == 'B53':
        return -1
    if interact_char[2] == 'W':
        interact_char = interact_char[0] + interact_char[2] + interact_char[1]
    if interact_char[1] == 'S' and interact_char[2] == 'H':
        interact_char = interact_char[0] + 'HS'
    traduction = ['CHH', 'TWH', 'CWW', 'THS', 'CWS', 'CSS', 'CWH', 'CHS','TWS','TSS','TWW','THH']
    return traduction.index(interact_char)
    
def DistanceEdgesLabel(x, y, Distance_matrix, len_edges_pattern):
    (xx, yy) = interaction_to_number(x), interaction_to_number(y)
    if xx == yy:
        return 1
    return 0
    if xx == -1 and yy == -1:
        return 1
    if xx == -1 or yy == -1 or Distance_matrix[xx][yy] > 10: #To change it will not be these "toy" characters forever
        return 0 #Too different or backbone not respected
    else:
        return 1 #- Distance_matrix[xx][yy]/(10*len_edges_pattern)
    #TODO : carry on by distinguishing cases depending on distance on label
    #Les bombres placés devront être de la forme 0.5 /len(edges_pattern) pour garder la logique de pourcentages ?


def _EdgeLabelRespect(x, y, label_edge_ij, len_edges_pattern, nodes_target, label_edge_target, edges_target, IDI_matrix, oriented):
    if DEBUG:
        print("Positions of mapping tried by EdgeLabelRespect", x, y)
    if (not oriented):
        x, y =  min(x, y), max(x, y)
    print("Positions of mapping tried by EdgeLabelRespect", x, y)
    print("Positions of mapping tried by EdgeLabelRespect", nodes_target[x], nodes_target[y])
    if (nodes_target[x], nodes_target[y]) in edges_target:
        return DistanceEdgesLabel(label_edge_ij,label_edge_target[edges_target.index((nodes_target[x], nodes_target[y]))], IDI_matrix, len_edges_pattern)
        #if (label_edge_ij == label_edge_target[edges_target.index((x, y))]):
            #if DEBUG:
            #    print("that succeded !")
            #return 1
            
    return 0 #currently 0 simply eliminate this pattern #-math.inf

def rev_matrix(mat):
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            if i < j:
                mat[j][i] = mat[i][j]




    
#class EdgeRespect(infrared.infrared.WeightedFunction):
#    """
#    Constrain complementarity mapping of any pair of nodes (i,j) that are an edge in graph_pattern
#
#    ```
#    EdgeRespect(i, j, edges_target)
#    ```
#    The constraint is satisfied if values of the mapping (i; j) form an edge in graph_target.
#    """
#def_function_class('EdgeRespect', lambda i, j, edges_target: [i, j],
#            lambda x, y, edges_target: _EdgeRespect(x, y, edges_target),
#            module=__name__)


class LabelRespect(infrared.infrared.WeightedFunction):
    """
    Constrain complementarity mapping of a node i

    ```
    Label(i,label_pattern, label_target)
    ```
    The constraint is satisfied if label of the mapping i is near or equal to the label of i.
    """
def_function_class('LabelRespect', lambda i, label_i, label_target: [i],
            lambda x, label_i, label_target: _LabelRespect(x, label_i, label_target),
            module=__name__)

class EdgeLabelRespect(infrared.infrared.WeightedFunction):
    """
    Constrain complementarity mapping of any pair of nodes (i,j) that are an edge in graph_pattern
    ```
    Label(i, j,label_edge_ij label_edge_target, edges_target)
    ```
    The constraint is satisfied if label of the mapped edge (i,j) is an edge and if its label near or equal to the label of (i, j).
    """
def_function_class('EdgeLabelRespect', lambda i, j, label_edge_ij, len_edges_pattern, nodes_target, label_edge_target, edges_target, IDI_matrix, oriented: [i, j],
            lambda x, y, label_edge_ij, len_edges_pattern, nodes_target, label_edge_target, edges_target, IDI_matrix, oriented: _EdgeLabelRespect(x, y, label_edge_ij, len_edges_pattern, nodes_target, label_edge_target, edges_target, IDI_matrix, oriented),
            module=__name__)


def main_old(graph_pattern, label_pattern, label_edge_pattern , graph_target, label_target, label_edge_target, oriented=0):
    
    IDI_matrix = [[8.9, 12.0, 14.7, 14.0, 13.7, 12.7, 15.1, 14.7, 16.2, 16.6, 16.2, 14.0],
            [0, 2.6, 10.6, 9.7, 14.3, 15.6, 11.2, 15.2, 13.8, 15.4, 11.9, 11.4],
            [0, 0, 4.1, 8.2, 9.2, 13.1, 14.5, 16.0, 12.4, 11.3, 11.1, 15.5],
            [0, 0, 0, 2.1, 7.0, 12.7, 12.0, 12.1, 10.0, 11.9, 13.1, 15.8],
            [0, 0, 0, 0, 3.5, 7.4, 14.9, 12.3, 10.9, 10.8, 14.6, 17.7],
            [0, 0, 0, 0, 0, 1.3, 15.8, 12.9, 13.8, 12.0, 17.1, 19.0],
            [0, 0, 0, 0, 0, 0, 3.2, 8.8, 8.4, 11.5, 10.6, 10.8],
            [0, 0, 0, 0, 0, 0, 0, 2.4, 7.9, 11.2, 14.7, 14.9],
            [0, 0, 0, 0, 0, 0, 0, 0, 3.4, 6.4, 9.6, 14.4],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 2.2, 9.0, 14.4],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.8, 9.0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4.0]]
    rev_matrix(IDI_matrix)
    
    
    
    n_pattern, edges_pattern, n_target, edges_target = init_from_graph(graph_pattern, label_pattern, graph_target, label_target, oriented)
    model = ir.Model(n_pattern, n_target)

    #First we put constraint on Injectivity to have a proper mapping
    #model.add_constraints(Injectivity(i,j) for i in range(n_pattern) for j in range(n_pattern) if i != j)
    
    #We define Edge respect function on each couples of edges in the pattern to check if these edges appeared in the target
    #model.add_functions([EdgeRespect(i, j, edges_target) for (i,j) in edges_pattern], 'EdgeRespect')
    #We could have put the list on edge already but set_weight not fully understood by me ?
    #model.set_feature_weight(1, 'EdgeRespect')

    #We define Label respect function on each nodes in the pattern to check if the label in the target correspond to it.
    model.add_functions([LabelRespect(i, label_pattern[i], label_target) for i in range(n_pattern)], 'LabelRespect')

    #We define Label respect function on each edges in the pattern to check if the label in the target correspond to it.
    model.add_functions([EdgeLabelRespect(i, j, label_edge_pattern[k], len(edges_pattern), label_edge_target, edges_target, IDI_matrix, oriented) for k, (i,j) in enumerate(edges_pattern)], 'EdgeLabelRespect')
    #len(edges_pattern) here will allow us to scale th error depending on the size of the motif wihtout touching the tolerancy indicated in the set_target part
    #We now take some samples of results
    sampler = ir.Sampler(model)

    #Next we ensure that edges in the pattern are represented correctly for now (or enough in the future) n the target
    #sampler.set_target(targetting value; tolerance; name of group of functions)
    #sampler.set_target(1*len(edges_pattern), 0, 'EdgeRespect')
    #Same for labels on nodes for now
    sampler.set_target(1*n_pattern, 0, 'LabelRespect')
    #Same for labels on edges
    sampler.set_target(1*len(edges_pattern), 0.999 , 'EdgeLabelRespect') #All faulty xase returns a 0 which eliminate the corresponding matching
    #Other fuzzy mistake ust returns something high and close to 1 to stay in the long run in the pattern
    samples = [sampler.targeted_sample() for _ in range(10)]
    resu = [sample.values() for sample in samples]
    print('edges_pattern', edges_pattern, "edges_target", edges_target)
    return resu


#graph_pattern = [[1], [0]] #graph in adjacency matrix here
#label_pattern = [0, 0]

#graph_target = [[1, 2], [0], [0]]
#label_target = [0, 0, 0]

#graph_pattern = [[1], [0]] #graph in adjacency matrix here
#label_pattern = [0, 0]

#graph_target = [[1, 2], [0], [0]]
#label_target = [0, 0, 1]



#A first example that are nothing to do with RNA but is a good toy example
graph_pattern = [[1, 2, 3], [0, 2], [0, 1, 3, 4], [0, 2], [2]] 
label_pattern = ["O", "O", "S", "O", "T"]
label_edge_pattern = ["L", "D", "S", "K", "S", "L"]
graph_target = [[], [2, 3, 4], [1, 4, 5], [1, 7],  [1, 2],
                [2, 6, 7, 10, 12], [5, 10, 11], 
                [3, 5, 8, 10], [7, 9], [8, 11], [5, 6, 7, 11],
                [6, 9, 10], [5]]
label_target = ["BLUB", "O", "O", "T", "S", "O", "O",
                "T","O","O","S","O","S"]
label_edge_target = ["L", "L", "D", "K", "L", "L","S", "S", "S", "S", "D", "L", "L", "L", "L", "L", "K"]


#In case of oriented graph, the adjacency matrix should be constructed accordingly and put orieted = 1 in main entry
#main_old(graph_pattern, label_pattern, label_edge_pattern, graph_target, label_target, label_edge_target)
#sans l'injectivité et avec même edge label l'algo peut confondre le noeud 5 et le noeud 11



def main(GP, GT):
    
    IDI_matrix = [[8.9, 12.0, 14.7, 14.0, 13.7, 12.7, 15.1, 14.7, 16.2, 16.6, 16.2, 14.0],
            [0, 2.6, 10.6, 9.7, 14.3, 15.6, 11.2, 15.2, 13.8, 15.4, 11.9, 11.4],
            [0, 0, 4.1, 8.2, 9.2, 13.1, 14.5, 16.0, 12.4, 11.3, 11.1, 15.5],
            [0, 0, 0, 2.1, 7.0, 12.7, 12.0, 12.1, 10.0, 11.9, 13.1, 15.8],
            [0, 0, 0, 0, 3.5, 7.4, 14.9, 12.3, 10.9, 10.8, 14.6, 17.7],
            [0, 0, 0, 0, 0, 1.3, 15.8, 12.9, 13.8, 12.0, 17.1, 19.0],
            [0, 0, 0, 0, 0, 0, 3.2, 8.8, 8.4, 11.5, 10.6, 10.8],
            [0, 0, 0, 0, 0, 0, 0, 2.4, 7.9, 11.2, 14.7, 14.9],
            [0, 0, 0, 0, 0, 0, 0, 0, 3.4, 6.4, 9.6, 14.4],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 2.2, 9.0, 14.4],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.8, 9.0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4.0]]
    rev_matrix(IDI_matrix)
    
    n_pattern = len(GP.nodes)
    n_target = len(GT.nodes)
    edges_pattern = list(GP.edges())#[(x-1, y-1) for (x,y) in list(GP.edges())]
    edges_target = list(GT.edges())#[(x-1, y-1) for (x,y) in list(GT.edges())]
    nodes_pattern = list(GP.nodes())
    nodes_target = list(GT.nodes())
    label_edge_pattern = [t['label'] for (_, _, t) in GP.edges.data()]
    label_edge_target = [t['label'] for (_, _, t) in GT.edges.data()]
    print("HELLO\n", n_pattern, n_target, edges_pattern, edges_target, nodes_pattern, nodes_target, label_edge_pattern, label_edge_target)
    #No label on vertices needed for now    
    model = ir.Model(n_pattern, n_target)

    #First we can put constraint on Injectivity to have a proper mapping and rewrok them around backbone edges
    #model.add_constraints(Injectivity(i,j) for i in range(n_pattern) for j in range(n_pattern) if i != j)

    #We define Label respect function on each edges in the pattern to check if the label in the target correspond to it.
    model.add_functions([EdgeLabelRespect(nodes_pattern.index(i), nodes_pattern.index(j), label_edge_pattern[k], len(edges_pattern), nodes_target, label_edge_target, edges_target, IDI_matrix, oriented=1) for k, (i,j) in enumerate(edges_pattern)], 'EdgeLabelRespect')
    #len(edges_pattern) here will allow us to scale th error depending on the size of the motif wihtout touching the tolerancy indicated in the set_target part
    #We now take some samples of results
    sampler = ir.Sampler(model)
    
    #Next we ensure that edges in the pattern are represented correctly for now (or enough in the future) n the target
    sampler.set_target(1*len(edges_pattern), 0.999 , 'EdgeLabelRespect') #All faulty case returns a 0 which eliminates the corresponding matching
    #Other fuzzy mistakes must return something high and close to 1 to stay in the long run in the pattern
    samples = [sampler.targeted_sample() for _ in range(1)]
    print("blub")
    resu = [([nodes_target[x] for x in sample.values()], len(list(set(sample.values())))) for sample in samples]
    print('edges_pattern', edges_pattern, "edges_target", edges_target)
    return resu


#Graphviz for Python test

#dot = gv.Digraph(comment='The Round Table')
#dot.node('A', 'King Arthur')
#dot.node('B', 'Sir Bedevere the Wise')
#dot.node('L', 'Sir Lancelot the Brave')
#dot.edges(['AB', 'AL'])
#dot.edge('B', 'L', constraint='false')
#print(dot.source)
#dot.render('doctest-output/round-table.gv', view=True) 




#GP = nx.read_gpickle("rin147-GP.pickle")
#GT = nx.read_gpickle("6V9D-GT.pickle")

with open("rin147-GP.pickle",'rb') as fP:
    GP = pickle.load(fP)

with open("1A3M-GT.nxpickle",'rb') as fT:
    GT = pickle.load(fT)


print(GP.edges(data=True))
#print(GP.nodes(data=True))
print(GT.edges(data=True))
#print(GT.nodes(data=True))

main(GP, GT)

#def DistanceEdgesLabel(x, y, Distance_matrix, len_edges_pattern):
#    """
#    Input : - Two edges labels x and y respectively from GP and GT.
#               - A Distance Matrix between non-canonical interactions/labels (for instance the IDI matricx for isostericity distance).
#               - len_edges_patterns to scale the impact of errors with the size of the graph pattern.
#    Output : A value between 0 and 1 that indicates how much these labels are close from each other.
#    """
#    (number_BM_allowed, elim_iso_threshold, BM_iso_threshold, B53_missing_elim, allow_iso_nonBM, _, BM_by_gap, scale_BM_with_P_size) = FuzzyParameters()
#    if scale_BM_with_P_size:
#        number_BM_allowed = number_BM_allowed*len_edges_pattern/20
#    (xx, yy) = interaction_to_number(x), interaction_to_number(y)
#    if xx == yy: #including the case where xx == -1 (is a backbone)
#        return 1
#    if xx == -1 or yy == -1: 
#        if B53_missing_elim:
#            return 0 #backbone not respected
#        else:
#            if number_BM_allowed:
#                return 1 - 1/number_BM_allowed
#            else:
#                return 0
#    if Distance_matrix[xx][yy] >= elim_iso_threshold: 
#        return 0 #Too different 
#    if Distance_matrix[xx][yy] >= BM_iso_threshold: 
#        if number_BM_allowed:
#            return 1 - 1/number_BM_allowed #Just a big mistake
#        else:
#            return 0
#    if allow_iso_nonBM and number_BM_allowed: #If small mistakes (non BM) are allowed
#        return 1 - Distance_matrix[xx][yy]/(20*number_BM_allowed)
#    return 1



def print_mapping_on_target_graph(GP, GT, mapping = [], output_format = "pdf", name_file = "", show=1):
    """
    Input : - Two graphs, the pattern graph GP and the target graph GT 
            - A mapping between GP and GT, if this mapping is not specified, it is calculated here
            - output_format is the type of file that we want to save, specify None for no save.
            - name_file is the name of the output graph if specified.
            - show is set to 1 if we want to observe the output with matplotlib.
    Output : Build the GT graph with a red colored part for each node in GP mapped in GT, can save it as a file or plot it depending on the optons selected.
    """
    if mapping == []:
        mapping = main(GP, GT)
        mapping = [map for (map, length) in mapping if length == len(GP.nodes())]
        mapping = mapping[0] #We take only the first result
    GG = GT.copy()
    pos = nx.nx_agraph.graphviz_layout(GG, prog="sfdp") #The layout can be modified is the graph is too "packed" visually
    #pos = nx.spring_layout(GG, k=0.2)  #An alternative layout
    
    nx.draw(GG, pos, width=1, linewidths=1)
    
    # nodes drawing
    options = {"edgecolors": "tab:black", "node_size": 500, "alpha": 0.9}
    mapped = [j for (_, j) in mapping]
    nx.draw_networkx_nodes(GG, pos, nodelist=mapped, node_color="tab:red")
    nx.draw_networkx_nodes(GG, pos, nodelist=[i for i in list(GG.nodes()) if i not in mapped])

    # edges drawing
    nx.draw_networkx_edges(GG, pos)

    # labels drawing
    labels = {}
    for (i, j) in mapping:
        labels[j] = i
    nx.draw_networkx_labels(GG, pos, labels, font_color="whitesmoke")
    labels_edges = {}
    for (i, j, t) in GT.edges.data():
        labels_edges[(i,j)] = t['label']
    nx.draw_networkx_edge_labels(GG, pos, edge_labels = labels_edges, font_size=5)

    plt.axis("off")
    if output_format:
        if name_file == "":
            name = "GPintoGT" + "." + output_format
        else:
            name = name_file + "." + output_format
        plt.savefig(name) 
    if show:
        plt.show()
