
import math
import os, glob, pickle
import networkx as nx

#!export PYTHONPATH="${PYTHONPATH}:/home/uqamportable/miniconda3/lib/python3.9/site-packages/infrared/"
import infrared as ir
from infrared import def_constraint_class, def_function_class

DEBUG = 0
#infinite = 1000
infinite = 100

def preL2distance(triplet1, triplet2):
    """
        Compute square of L2 distances between the three coordinates triplet1=(x1,y1,z1) and triplet2=(x2,y2,z2).
    """
    (x1,y1,z1) = triplet1
    (x2,y2,z2) = triplet2
    (x1,y1,z1) = float(x1), float(y1), float(z1)
    (x2,y2,z2) = float(x2), float(y2), float(z2)
    return (x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2

def distance_cube(node, cube,  Distancer):
    """
    Input : - One node and one cube in GT.
            - Preprocessed distance between nodes of GT Distancer.
    Output :
            - The minimal distance between node and the nodes in the cube in GT by considering the atoms that are the more close to each other.
    """
    dist = math.inf
    for node_cube in cube:
            dist = min(dist, Distancer[node][node_cube])
    return dist

def get_diameter(GP):
    diam = 0
    for node1 in GP.nodes.data():
        for node2 in GP.nodes.data():
            atoms1 = node1[1]['atoms']
            atoms2 = node2[1]['atoms']
            for a1 in atoms1:
                for a2 in atoms2:
                    diam = max(diam, preL2distance(a1['position'], a2['position']))
    return math.sqrt(diam)


def full_allocate(GT, cutoff):
    grid = {}
    for node in GT.nodes.data():
        id = node[0]
        atoms = node[1]['atoms']
        positions = [t['position'] for t in atoms]
        for  (x, y, z) in positions:
            (xx, yy, zz) = (float(x), float(y), float(z))
            row = int(xx / cutoff), int(yy / cutoff), int(zz / cutoff)
            if row not in grid.keys():
                grid[row] = []
            if id not in grid[row]:
                grid[row].append(id)
    return grid
import itertools
def full_cube(GT, cutoff, Distancer):
    grid = full_allocate(GT, cutoff)
    new_grid = {}
    for row in grid:
        new_grid[row] = grid[row]
        print("blub")
        for dx, dy, dz in itertools.product([-1, 0, 1], repeat=3):
            drow = (row[0] + dx, row[1] + dy, row[2] + dz)
            if drow in grid:
                preresu = [node for node in grid[drow] if distance_cube(node, grid[row], Distancer) <= diam]
                resu = list(set(new_grid[row] + preresu))
                new_grid[row] = resu
    return new_grid

with open("ALLkinkturntarget/22IL_29549.9into5J7L.nxpickle",'rb') as fP:
    GP = pickle.load(fP)
#print(GP.nodes.data())
with open("/home/uqamportable/Documents/FuzzTreeFirstCode/bigRNAstorage/4WF9.nxpickle", 'rb') as fT:
    GT = pickle.load(fT)
    edges_to_remove = [(i, j) for (i, j, t) in GT.edges.data() if t['near'] == True]
    for (i, j) in edges_to_remove:
        GT.remove_edge(i, j)
    chains = ['X']
    Gnew=nx.DiGraph() #Initiate the new GT graph.
    for ((i, ii),t) in GT.nodes.data():
        if i in chains:
            Gnew.add_node((i, ii), pdb_position = t['pdb_position'], atoms = t['atoms'])
    for ((i, ii),(j, jj),t) in GT.edges.data():
        if i in chains and j in chains:
            Gnew.add_edge((i, ii),(j, jj), label=t['label'], near=t['near'])
    GT = Gnew
from FuzzTree import distance
Distancer = {}
for node1 in GT.nodes():
    loc_distance = {}
    for node2 in GT.nodes():
        if node2 in Distancer.keys():
            if node1 in Distancer[node2].keys():
                loc_distance[node2] = Distancer[node2][node1]
            else:
                loc_distance[node2] = distance(node1, node2, GT)
        else:
            loc_distance[node2] = distance(node1, node2, GT)
    Distancer[node1] = loc_distance.copy()
diam = get_diameter(GP)
print(diam)
grid = full_allocate(GT, diam)
print(len(grid))
print(([len(grid[i]) for i in grid.keys()]))
grid = full_cube(GT, diam, Distancer)
print(len(grid))
print(([len(grid[i]) for i in grid.keys()]))