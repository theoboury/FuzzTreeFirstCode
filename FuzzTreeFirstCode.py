import infrared as ir
import infrared.rna as rna
import matplotlib.pyplot as plt

import re
import collections
import math
import sys
import os
import infrared
#!export PYTHONPATH="${PYTHONPATH}:/home/uqamportable/miniconda3/python3.9/site-packages/infrared/"
from infrared import def_constraint_class, def_function_class

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

#WARNING : Should we define another injectvity criteria

#def _Injectivity(x, y):
#    if DEBUG:
#        print("Positions of mapping tried by Injectivity", x, y)
#    if x == y:
#        return 0
#    else:
#        if DEBUG:
#            print("that succeded !")
#        return 1


#class Injectivity(infrared.Constraint):
#    """
#    Constrain complementarity mapping of any pair of nodes (i,j)
#
#    ```
#    Injectivity(i, j)
#    ```
#    The constraint is satisfied if values of the mapping for i and j are distincts.
#    """
#def_constraint_class('Injectivity', lambda i, j: [i, j],
#                     lambda x, y: _Injectivity(x, y),
#                     module=__name__)



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

def DistanceEdgesLabel(x, y, Distance_matrix, DistanceMax = -1):
    if DistanceMax == -1:
        DistanceMax = -math.inf
        for i in range(len(Distance_matrix)):
            DistanceMax = max(DistanceMax, max(Distance_matrix[i]))
    if Distance_matrix[ord(x)%12][ord(y)%12] > 10: #To change it will not be these "toy" characters forever
        return 0
    else:
        return 1
    #TODO : carry on by distinguishing cases depending on distance on label
    #Les bombres placés devront être de la forme 0.5 /len(edges_pattern) pour garder la logique de pourcentages ?


def _EdgeLabelRespect(x, y, label_edge_ij, len_edges_pattern, label_edge_target, edges_target, IDI_matrix, oriented):
    if DEBUG:
        print("Positions of mapping tried by EdgeLabelRespect", x, y)
    if (not oriented):
        x, y =  min(x, y), max(x, y)
    if (x, y) in edges_target:
        return DistanceEdgesLabel(label_edge_ij,label_edge_target[edges_target.index((x, y))], IDI_matrix)
        #if (label_edge_ij == label_edge_target[edges_target.index((x, y))]):
            #if DEBUG:
            #    print("that succeded !")
            #return 1
            
    return 0#-math.inf

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
def_function_class('EdgeLabelRespect', lambda i, j, label_edge_ij, len_edges_pattern, label_edge_target, edges_target, IDI_matrix, oriented: [i, j],
            lambda x, y, label_edge_ij, len_edges_pattern, label_edge_target, edges_target, IDI_matrix, oriented: _EdgeLabelRespect(x, y, label_edge_ij, len_edges_pattern, label_edge_target, edges_target, IDI_matrix, oriented),
            module=__name__)


def main(graph_pattern, label_pattern, label_edge_pattern , graph_target, label_target, label_edge_target, oriented=0):
    
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
main(graph_pattern, label_pattern, label_edge_pattern, graph_target, label_target, label_edge_target)
#sans l'injectivité et avec même edge label l'algo peut confondre le noeud 5 et le noeud 11
