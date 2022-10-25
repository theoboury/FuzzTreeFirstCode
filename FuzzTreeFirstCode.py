import infrared as ir
import infrared.rna as rna
import matplotlib.pyplot as plt

import re
import collections
import math
import infrared
!export PYTHONPATH="${PYTHONPATH}:/home/uqamportable/miniconda3/python3.9/site-packages/infrared/"
from infrared import def_constraint_class, def_function_class

DEBUG = 0

def init_from_graph(graph_pattern, label_pattern, graph_target, label_target):
    n_pattern = len(graph_pattern)
    edges_pattern = [(i, j) for i in range(n_pattern) for j in range(n_pattern) if j in graph_pattern[i] and i < j]
    n_target = len(graph_target)
    edges_target = [(i, j) for i in range(n_target) for j in range(n_target) if j in graph_target[i] and i < j]
    return n_pattern, edges_pattern, n_target, edges_target


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
    The constraint is satisfied if values of the mapping for i and j are distincts.
    """
def_constraint_class('Injectivity', lambda i, j: [i, j],
                     lambda x, y: _Injectivity(x, y),
                     module=__name__)



def _EdgeRespect(x, y, edges_target):
    if DEBUG:
        print("Positions of mapping tried by EdgeRespect", x, y)
    if (x, y) in edges_target or (y, x) in edges_target:
        if DEBUG:
            print("that succeded !")
        return 1
    else:
        return 0


def _LabelRespect(x, label_i, label_target):
    if DEBUG:
        print("Positions of mapping tried by LabelRespect", x)
    if (label_i == label_target[x]):
        if DEBUG:
            print("that succeded !")
        return 1
    else:
        return 0


class EdgeRespect(infrared.infrared.WeightedFunction):
    """
    Constrain complementarity mapping of any pair of nodes (i,j) that are an edge in graph_pattern

    ```
    EdgeRespect(i, j, edges_target)
    ```
    The constraint is satisfied if values of the mapping (i; j) form an edge in graph_target.
    """
def_function_class('EdgeRespect', lambda i, j, edges_target: [i, j],
            lambda x, y, edges_target: _EdgeRespect(x, y, edges_target),
            module=__name__)


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



def main(graph_pattern, label_pattern, graph_target, label_target):
    n_pattern, edges_pattern, n_target, edges_target = init_from_graph(graph_pattern, label_pattern, graph_target, label_target)
    model = ir.Model(n_pattern, n_target)

    #First we put constraint on Injectivity to have a proper mapping
    model.add_constraints(Injectivity(i,j) for i in range(n_pattern) for j in range(n_pattern) if i != j)
    
    #We define Edge respect function on each couples of edges in the pattern to check if these edges appeared in the target
    model.add_functions([EdgeRespect(i, j, edges_target) for (i,j) in edges_pattern], 'EdgeRespect')
    #We could have put the list on edge already but set_weight not fully understood by me ?
    #model.set_feature_weight(1, 'EdgeRespect')

    #We define Label respect function on each nodes in the pattern to check if the label in the target correspond to it.
    model.add_functions([LabelRespect(i, label_pattern[i], label_target) for i in range(n_pattern)], 'LabelRespect')

    #We now take some samples of results
    sampler = ir.Sampler(model)

    #Next we ensure that edges in the pattern are represented correctly for now (or enough in the future) n the target
    #sampler.set_target(targetting value; tolerance; name of group of functions)
    sampler.set_target(1*len(edges_pattern), 0, 'EdgeRespect')
    #Same for labels on nodes for now
    sampler.set_target(1*n_pattern, 0, 'LabelRespect')
    samples = [sampler.targeted_sample() for _ in range(10)]
    resu = [sample.values() for sample in samples]
    print(resu)
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

graph_target = [[], [2, 3, 4], [1, 4, 5], [1, 7],  [1, 2],
                [2, 6, 7, 10, 12], [5, 10, 11], 
                [3, 5, 8, 10], [7, 9], [8, 11], [5, 6, 7, 11],
                [6, 9, 10], [5]]
label_target = ["BLUB", "O", "O", "T", "S", "O", "O",
                "T","O","O","S","O","S"]



main(graph_pattern, label_pattern, graph_target, label_target)
