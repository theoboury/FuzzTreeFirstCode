import os, glob, pickle
import networkx as nx
import csv

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

def extract_with_pdb_number(G, list_nodes,  list_nodes_clean, cutting_edges):
    Gnew=nx.DiGraph() #Initiate the new GT graph.
    Gnewnx=nx.DiGraph() #Initiate the new GT graph.
    list_nodes_pdb = {}
    for ((i, j), t) in G.nodes.data():
        if (i, t['pdb_position']) in list_nodes:
            list_nodes_pdb[(i, j)] = (i, t['pdb_position'])
            index = list_nodes.index((i, t['pdb_position']))
            Gnew.add_node(index + 1, pdb_position = t['pdb_position'], atoms = t['atoms'], nt = t['nt'])
            Gnewnx.add_node(list_nodes_clean[index], pdb_position = t['pdb_position'], atoms = t['atoms'], nt = t['nt'])
    for (i,j,t) in G.edges.data():
        if i in list_nodes_pdb.keys() and j in list_nodes_pdb.keys() :
            indexi = list_nodes.index(list_nodes_pdb[i])
            indexj = list_nodes.index(list_nodes_pdb[j])
            Gnew.add_edge(indexi + 1, indexj + 1, label=t['label'], near=t['near'])
            Gnewnx.add_edge(list_nodes_clean[indexi], list_nodes_clean[indexj], label=t['label'], near=t['near'])
    for i in range(len(list_nodes) - 1):
        if (i + 1, i + 2) not in Gnew.edges() and (list_nodes[i], list_nodes[i + 1]) not in cutting_edges:
            Gnew.add_edge(i + 1, i + 2, label='B53', near=False)
            Gnewnx.add_edge(list_nodes_clean[i], list_nodes_clean[i + 1], label='B53', near=False)
    return Gnew, Gnewnx

def extractor(Gpath, Gnewname, list_nodes, cutting_edges, with_pdb_num=0, pattern_place="kinkturnpattern/", target_place ="kinkturntarget/", list_nodes_clean = []):
    with open(Gpath,'rb') as fG:
        G = pickle.load(fG)
        print("Gpath", Gpath)
    if with_pdb_num:
        Gnew, Gnewnx = extract_with_pdb_number(G, list_nodes, list_nodes_clean, cutting_edges)
    else:
        Gnew, Gnewnx = extract(G, list_nodes, cutting_edges)
    if DEBUG:
        print([(i) for (i,t) in Gnew.nodes.data()])
        #print([(i,j, t['label']) for (i,j,t) in Gnew.edges.data() if t['label'] != 'B53'])
        print([(i,j, t['label']) for (i,j,t) in Gnew.edges.data() if t['label'] != 'B53'])
        print([(i,j, t['label']) for (i,j,t) in Gnew.edges.data() if t['label'] == 'B53'])
        #print([(i) for (i,t) in Gnewnx.nodes.data()])
        #print([(i,j, t['label']) for (i,j,t) in Gnewnx.edges.data()])
    with open(pattern_place +Gnewname + ".pickle", 'wb') as ff:
        pickle.dump(Gnew, ff)
    with open(target_place + Gnewname + ".nxpickle", 'wb') as ff2:
        pickle.dump(Gnewnx, ff2)

def csv_parse(family_name, break_list):
    #break_list is the couples of nucleotides between which there is a break.
    # We consider that nucleotides are numeroted from 1 to n
    with open(family_name + '.csv', newline='', encoding='utf-8') as f:
        reader = csv.reader(f)
        for k,row in enumerate(reader):
            symmetry = 0
            RNAtarget = row[0].split('|')[0]
            Gpath = "bigRNAstorage/" + RNAtarget + ".nxpickle"
            print(Gpath)
            list_nodes = []
            list_nodes_clean = []
            cutting_edges = []
            for nucleo in row:
                explore_list = nucleo.split('|')
                #print(explore_list)
                if len(explore_list) > 8:
                    symmetry = 1
                if len(explore_list) > 7:
                    pdb_string = explore_list[4] + explore_list[7]
                else:
                    pdb_string = explore_list[4] 
                list_nodes_clean.append((explore_list[2], int(explore_list[4])))
                list_nodes.append((explore_list[2], pdb_string))
            #print(list_nodes)
            #print(list_nodes_clean)
            for (i, j) in break_list:
                cutting_edges.append((list_nodes[i - 1], list_nodes[j - 1]))
            #print(cutting_edges)
            if symmetry == 0:
                extractor(Gpath, str(k) + family_name + 'into' + RNAtarget, list_nodes, cutting_edges, with_pdb_num=1, pattern_place="ALLkinkturnpattern/", target_place ="ALLkinkturntarget/", list_nodes_clean = list_nodes_clean)
            else:
                print("SKIPPED DUE TO SYMMETRY")
