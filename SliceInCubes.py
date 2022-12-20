
import math
import networkx as nx
import itertools
from FuzzTree import preL2distance, precompute_distance
from multiprocessing import Pool
DEBUG = 1


def distance_cube(node, cube, Distancer):
    """
    Input : - One node and one cube in GT.
            - Preprocessed distance between nodes of GT, Distancer.
    Output :
            - The minimal distance between input node and the nodes in the cube in GT by considering the atoms that are the more close to each other.
    """
    dist = math.inf
    for node_cube in cube:
            dist = min(dist, Distancer[node][node_cube])
    return dist


def get_diameter(GP):
    """
    Input : A pattern graph GP.
    Output : The maximal distance between two nucleotides in GP.
    """
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


def full_allocate_cube(GT, cutoff_cube):
    """
    Input : - A target graph GT.
            - cutoff_cube to define the size of the cube.
    Output :  We output a grid that contains multiple list of nodes that are all the cubes of size GT that pave GT.
    """
    grid = {}
    for node in GT.nodes.data():
        id = node[0]
        atoms = node[1]['atoms']
        positions = [t['position'] for t in atoms]
        for  (x, y, z) in positions:
            (xx, yy, zz) = (float(x), float(y), float(z))
            row = int(xx / cutoff_cube), int(yy / cutoff_cube), int(zz / cutoff_cube)
            if row not in grid.keys():
                grid[row] = []
            if id not in grid[row]:
                grid[row].append(id)
    return grid

def wrapper_sphere(row_grid_GT_Distancer_cutoff_sphere):
    (row, grid, GT, Distancer, cutoff_sphere) = row_grid_GT_Distancer_cutoff_sphere
    preresu = [node for node in GT.nodes() if distance_cube(node, grid[row], Distancer) <= cutoff_sphere]
    return (row, list(set(grid[row] + preresu)))

def full_allocate_cube_and_sphere(GT, cutoff_cube, cutoff_sphere, Distancer, nb_procs):
    """
    Input : - A target graph GT.
            - cutoff_cube to define the size of the cube.
            - cutoff_sphere to define the size of the sphere around the cuve.
    Output :  We output a grid that contains multiple list of nodes that are all the cubes of size GT that pave GT extended by the sphere around each cube.
    """
    print("Starting grid\n")
    grid = full_allocate_cube(GT, cutoff_cube)
    sphere_grid = {}
    entry = []
    print("Starting sphere\n")
    for row in grid:
        entry.append((row, grid, GT, Distancer, cutoff_sphere))
    print("Entry sphere done\n")
    with Pool(nb_procs) as pool:
        resu= list(pool.imap_unordered(wrapper_sphere, entry))
    print("Resu sphere done\n")
    for (row, li) in resu:
        sphere_grid[row] = li
    print("Sphere done\n")
    return sphere_grid
def full_allocate_cube_and_sphere_old(GT, cutoff_cube, cutoff_sphere, Distancer):
    """
    Input : - A target graph GT.
            - cutoff_cube to define the size of the cube.
            - cutoff_sphere to define the size of the sphere around the cuve.
    Output :  We output a grid that contains multiple list of nodes that are all the cubes of size GT that pave GT extended by the sphere around each cube.
    """
    grid = full_allocate_cube(GT, cutoff_cube)
    sphere_grid = {}
    for row in grid:
        sphere_grid[row] = grid[row]
        for dx, dy, dz in itertools.product([-1, 0, 1], repeat=3):
            drow = (row[0] + dx, row[1] + dy, row[2] + dz)
            if drow in grid:
                preresu = [node for node in grid[drow] if distance_cube(node, grid[row], Distancer) <= cutoff_sphere]
                resu = list(set(sphere_grid[row] + preresu))
                sphere_grid[row] = resu
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

def slicer(GP, GT, nb_procs, size_cube_versus_radius=0.5, filename = "", D = 5, A = 0):
    """
    Input : - A pattern graph GP and a target graph GT.
            - size_cube_versus_radius serves to quantify the size of the buce compare to the radius of sphere as we are free 
              about the size of the cubes but it can have impact on the performances in practise.
            - filename, for debug purposes.
    Output :  We slice in tcube the GT graph depending on diameter and.or radius of the pattern graph. We return the list of 
              all subgraphs of GT obtained by slicing and also the precomputed distance between all nodes of GT as it already serves at this point.
    """
    diam = get_diameter(GP)
    rad = get_radius(GP)
    if DEBUG:
        print("Diameter", diam, "Radius", rad)
    Distancer = precompute_distance(GT, nb_procs) 
    grid = full_allocate_cube_and_sphere(GT, size_cube_versus_radius*rad, 20, Distancer, nb_procs)#rad + A + D
    if DEBUG:
        print("filename", filename, "Number of cubes", len(grid), "Max size cube", max([len(grid[i]) for i in grid.keys()]), "Size of each cube", ([len(grid[i]) for i in grid.keys()]))
    pre_graph_grid = [grid[i] for i in grid.keys()]
    pre_graph_grid.sort(key=len)
    graph_grid = []
    for list_nodes in pre_graph_grid:
        if len(list_nodes) >= len(GP.nodes.data()):
            graph_grid.append(extract_small_sphere_graph(GT, list_nodes))
    if DEBUG:
        print("filename", filename, "Number of cubes after", len(graph_grid), "Max size cube after", max([len(graph_grid[i].nodes.data()) for i in range(len(graph_grid))]), "Size of each cube after", [len(graph_grid[i].nodes.data()) for i in range(len(graph_grid))])
    return graph_grid, Distancer
