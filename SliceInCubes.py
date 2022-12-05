
import math
import os, glob, pickle
import networkx as nx

#!export PYTHONPATH="${PYTHONPATH}:/home/uqamportable/miniconda3/lib/python3.9/site-packages/infrared/"
import infrared as ir
from infrared import def_constraint_class, def_function_class
import itertools
from FuzzTree import distance
DEBUG = 1
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

def get_radius(GP):
    diam = math.inf
    for node1 in GP.nodes.data():
        diam_loc = 0
        for node2 in GP.nodes.data():
            atoms1 = node1[1]['atoms']
            atoms2 = node2[1]['atoms']
            for a1 in atoms1:
                for a2 in atoms2:
                    diam_loc = max(diam_loc, preL2distance(a1['position'], a2['position']))
        diam = min(diam, diam_loc)
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


def full_cube_diam_rad(GT, cutoff, radius, Distancer):
    grid = full_allocate(GT, cutoff)
    new_grid = {}
    for row in grid:
        new_grid[row] = grid[row]
        #print("blub")
        for dx, dy, dz in itertools.product([-1, 0, 1], repeat=3):
            drow = (row[0] + dx, row[1] + dy, row[2] + dz)
            if drow in grid:
                preresu = [node for node in grid[drow] if distance_cube(node, grid[row], Distancer) <= radius]
                resu = list(set(new_grid[row] + preresu))
                new_grid[row] = resu
    return new_grid

def extract_graph(GT, list_nodes):
    Gnew=nx.DiGraph()
    for ((i, ii),t) in GT.nodes.data():
        if (i, ii) in list_nodes:
            Gnew.add_node((i, ii), pdb_position = t['pdb_position'], atoms = t['atoms'])
    for ((i, ii),(j, jj),t) in GT.edges.data():
        if (i, ii) in list_nodes and (j, jj) in list_nodes:
            Gnew.add_edge((i, ii),(j, jj), label=t['label'], near=t['near'])
    return Gnew

def slicer(pattern_name, GT,  size_cube_versus_radius=0.5):
    with open("ALLkinkturntarget/" + pattern_name + ".nxpickle",'rb') as fP:
        GP = pickle.load(fP)
    diam = get_diameter(GP)
    rad = get_radius(GP)
    if DEBUG:
        print("Diameter", diam, "Radius", rad)
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
    grid = full_cube_diam_rad(GT, size_cube_versus_radius*rad, rad, Distancer)
    if DEBUG:
        print("Number of cubes", len(grid), "Max size cube", max([len(grid[i]) for i in grid.keys()]), "Size of each cube", ([len(grid[i]) for i in grid.keys()]))
    pre_graph_grid = [grid[i] for i in grid.keys()]
    pre_graph_grid.sort(key=len)
    graph_grid = []
    for list_nodes in pre_graph_grid:
        graph_grid.append(extract_graph(GT, list_nodes))
    return graph_grid, Distancer

#with open("ALLkinkturntarget/22IL_29549.9into5J7L.nxpickle",'rb') as fP:
#    GP = pickle.load(fP)
#print(GP.nodes.data())
#diam = get_diameter(GP)
#rad = get_radius(GP)
#print(diam, rad)
#with open("/home/uqamportable/Documents/FuzzTreeFirstCode/bigRNAstorage/4WF9.nxpickle", 'rb') as fT:
#    GT = pickle.load(fT)
#    edges_to_remove = [(i, j) for (i, j, t) in GT.edges.data() if t['near'] == True]
#    for (i, j) in edges_to_remove:
#        GT.remove_edge(i, j)
#    chains = ['X']
#    Gnew=nx.DiGraph() #Initiate the new GT graph.
#    for ((i, ii),t) in GT.nodes.data():
#        if i in chains:
#            Gnew.add_node((i, ii), pdb_position = t['pdb_position'], atoms = t['atoms'])
#    for ((i, ii),(j, jj),t) in GT.edges.data():
#        if i in chains and j in chains:
#            Gnew.add_edge((i, ii),(j, jj), label=t['label'], near=t['near'])
#    GT = Gnew

#Distancer = {}
#for node1 in GT.nodes():
#    loc_distance = {}
#    for node2 in GT.nodes():
#        if node2 in Distancer.keys():
#            if node1 in Distancer[node2].keys():
#                loc_distance[node2] = Distancer[node2][node1]
#            else:
#                loc_distance[node2] = distance(node1, node2, GT)
#        else:
#            loc_distance[node2] = distance(node1, node2, GT)
#    Distancer[node1] = loc_distance.copy()
#diam = get_diameter(GP)
#rad = get_radius(GP)
#print(diam, rad)
#grid = full_allocate(GT, diam)
#print(len(grid))
#print(([len(grid[i]) for i in grid.keys()]))
#grid = full_cube(GT, diam, Distancer)
#print(len(grid))
#print(([len(grid[i]) for i in grid.keys()]))
#grid = full_cube_diam_rad(GT, rad/2, rad, Distancer)
#print(len(grid))
#print(([len(grid[i]) for i in grid.keys()]))

#116
#[213, 734, 503, 816, 444, 647, 514, 657, 963, 1090, 286, 278, 214, 78, 210, 413, 241, 612, 475, 645, 1009, 753, 1263, 842, 632, 336, 311, 240, 466, 286, 529, 135, 145, 92, 185, 95, 246, 102, 89, 427, 242, 358, 68, 71, 83, 357, 1177, 584, 389, 151, 611, 163, 872, 342, 148, 189, 226, 295, 218, 419, 635, 715, 632, 358, 29, 551, 550, 313, 196, 150, 210, 107, 70, 39, 44, 173, 79, 267, 327, 144, 71, 97, 120, 57, 62, 62, 454, 582, 521, 499, 261, 313, 177, 102, 88, 167, 327, 319, 194, 226, 154, 441, 137, 281, 91, 379, 140, 48, 82, 97, 94, 122, 339, 114, 84, 268]
#with the modification with distances