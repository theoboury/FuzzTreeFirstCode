
import math
import networkx as nx
from FuzzTree import preL2distance, precompute_distance
from multiprocessing import Pool
DEBUG = 1


def get_radius(GP):
    """
    Input : A pattern graph GP.
    Output : Locate the node the more "central" node c in GP. returns the maximal distance between c and the others nucleotides in GP.
    """
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


def wrapper_sphere(GT_Distancer_cube_cutoff_sphere):
    """
    A werapper around the creation of spheres to allow multiprocessing"""
    (GT, Distancer_cube, cutoff_sphere) = GT_Distancer_cube_cutoff_sphere
    preresu = [node for node in GT.nodes() if Distancer_cube[node] <= cutoff_sphere]
    return list(set(preresu))

def allocate_sphere(GT, cutoff_sphere, Distancer, nb_procs):
    """
    Input : - A target graph GT.
            - cutoff_sphere to define the size of the sphere around each nodes.
            - Distancer the precomputed distance between nodes of GT.
            - nb_procs the number of allowed processors to multiprocess this precomputation.
    Output :  We output a grid that contains multiple list of nodes that are all spheres around nodes of GT.
    """
    sphere_grid = []
    entry = []
    for node in GT.nodes.data():
        entry.append((GT, Distancer[node[0]].copy(), cutoff_sphere))
    with Pool(nb_procs) as pool:
        resu= list(pool.imap_unordered(wrapper_sphere, entry))
    for li in resu:
        lili = li
        lili.sort()
        if lili not in sphere_grid:
            sphere_grid.append(lili)
    if DEBUG:
        print("Sphere grid done\n")
    return sphere_grid


def extract_small_sphere_graph(GT, list_nodes):
    """
    Input : - A target graph GT.
            - A list of nodes list_nodes.
    Output :  We extract from GT the minimal subgraph that contains all the nodes in list_nodes. 
    """
    Gnew=nx.DiGraph()
    for ((i, ii),t) in GT.nodes.data():
        if (i, ii) in list_nodes:
            Gnew.add_node((i, ii), pdb_position = t['pdb_position'], atoms = t['atoms'])
    for ((i, ii),(j, jj),t) in GT.edges.data():
        if (i, ii) in list_nodes and (j, jj) in list_nodes:
            Gnew.add_edge((i, ii),(j, jj), label=t['label'], near=t['near'])
    return Gnew

def slicer(GP, GT, nb_procs, filename = "", D = 5, A = 0):
    """
    Input : - A pattern graph GP and a target graph GT.
            - size_cube_versus_radius serves to quantify the size of the buce compare to the radius of sphere as we are free 
              about the size of the cubes but it can have impact on the performances in practise.
            - filename, for debug purposes.
    Output :  We slice in tcube the GT graph depending on diameter and.or radius of the pattern graph. We return the list of 
              all subgraphs of GT obtained by slicing and also the precomputed distance between all nodes of GT as it already serves at this point.
    """
    rad = get_radius(GP)
    if DEBUG:
        print("Radius", rad)
    Distancer = precompute_distance(GT, nb_procs) 
    grid = allocate_sphere(GT, 20, Distancer, nb_procs)#TODO : here we cheated, should be rad + A + D instead of rad or 20 ?
    if DEBUG:
        print("filename", filename, "Number of cubes", len(grid), "Max size cube", max([len(grid[i]) for i in range(len(grid))]), "Size of each cube", ([len(grid[i]) for i in range(len(grid))]))
    pre_graph_grid = [grid[i] for i in range(len(grid))]
    pre_graph_grid.sort(key=len)
    graph_grid = []
    for list_nodes in pre_graph_grid:
        if len(list_nodes) >= len(GP.nodes.data()):
            graph_grid.append(extract_small_sphere_graph(GT, list_nodes))
    if DEBUG:
        print("filename", filename, "Number of cubes after", len(graph_grid), "Max size cube after", max([len(graph_grid[i].nodes.data()) for i in range(len(graph_grid))]), "Size of each cube after", [len(graph_grid[i].nodes.data()) for i in range(len(graph_grid))])
    return graph_grid, Distancer
