import os, glob, pickle
import networkx as nx

DEBUG=1

def extract(G, list_nodes, cutting_edges):
    Gnew=nx.DiGraph() #Initiate the new GT graph.
    Gnewnx=nx.DiGraph() #Initiate the new GT graph.
    for (i, t) in G.nodes.data():
        if i in list_nodes:
            Gnew.add_node(list_nodes.index(i) + 1, pdb_position = t['pdb_position'], atoms = t['atoms'], nt = t['nt'])
            Gnewnx.add_node(i, pdb_position = t['pdb_position'], atoms = t['atoms'], nt = t['nt'])
    for (i,j,t) in G.edges.data():
        if i in list_nodes or j in list_nodes:
            print( i, j, t['label'])
        if i in list_nodes and j in list_nodes:
            print("AH", t['label'])
            Gnew.add_edge(list_nodes.index(i) + 1, list_nodes.index(j) + 1, label=t['label'], near=t['near'])
            Gnewnx.add_edge(i, j, label=t['label'], near=t['near'])
    for i in range(len(list_nodes) - 1):
        if (i + 1, i + 2) not in Gnew.edges() and (list_nodes[i], list_nodes[i + 1]) not in cutting_edges:
            Gnew.add_edge(i + 1, i + 2, label='B53', near=False)
            Gnewnx.add_edge(list_nodes[i], list_nodes[i + 1], label='B53', near=False)
    return Gnew, Gnewnx

def extract_with_pdb_number(G, list_nodes, cutting_edges):
    Gnew=nx.DiGraph() #Initiate the new GT graph.
    Gnewnx=nx.DiGraph() #Initiate the new GT graph.
    list_nodes_pdb = {}
    for ((i, j), t) in G.nodes.data():
        if (i, int(t['pdb_position'])) in list_nodes:
            list_nodes_pdb[(i, j)] = (i, int(t['pdb_position']))
            Gnew.add_node(list_nodes.index((i, int(t['pdb_position']))) + 1, pdb_position = t['pdb_position'], atoms = t['atoms'], nt = t['nt'])
            Gnewnx.add_node((i, int(t['pdb_position'])), pdb_position = t['pdb_position'], atoms = t['atoms'], nt = t['nt'])
    for (i,j,t) in G.edges.data():
        if i in list_nodes_pdb.keys() and j in list_nodes_pdb.keys() :
            Gnew.add_edge(list_nodes.index(list_nodes_pdb[i]) + 1, list_nodes.index(list_nodes_pdb[j]) + 1, label=t['label'], near=t['near'])
            Gnewnx.add_edge(list_nodes_pdb[i], list_nodes_pdb[j], label=t['label'], near=t['near'])
    for i in range(len(list_nodes) - 1):
        if (i + 1, i + 2) not in Gnew.edges() and (list_nodes[i], list_nodes[i + 1]) not in cutting_edges:
            Gnew.add_edge(i + 1, i + 2, label='B53', near=False)
            Gnewnx.add_edge(list_nodes[i], list_nodes[i + 1], label='B53', near=False)
    return Gnew, Gnewnx

def extractor(Gpath, Gnewname, list_nodes, cutting_edges, with_pdb_num=0):
    with open(Gpath,'rb') as fG:
        G = pickle.load(fG)
        print("Gpath", Gpath)
    if with_pdb_num:
        Gnew, Gnewnx = extract_with_pdb_number(G, list_nodes,cutting_edges)
    else:
        Gnew, Gnewnx = extract(G, list_nodes,cutting_edges)
    if DEBUG:
        print([(i) for (i,t) in Gnew.nodes.data()])
        #print([(i,j, t['label']) for (i,j,t) in Gnew.edges.data() if t['label'] != 'B53'])
        print([(i,j, t['label']) for (i,j,t) in Gnew.edges.data() if t['label'] != 'B53'])
        print([(i,j, t['label']) for (i,j,t) in Gnewnx.edges.data() if t['label'] != 'B53'])
        #print([(i) for (i,t) in Gnewnx.nodes.data()])
        #print([(i,j, t['label']) for (i,j,t) in Gnewnx.edges.data()])
    with open("kinkturnpattern/" +Gnewname + ".pickle", 'wb') as ff:
        pickle.dump(Gnew, ff)
    with open("kinkturntarget/" + Gnewname + ".nxpickle", 'wb') as ff2:
        pickle.dump(Gnewnx, ff2)


