# FuzzTree
# Copyright (C) 2023 THEO BOURY 

import os, glob, pickle
import networkx as nx
import csv

DEBUG=0


def extract_with_pdb_number(G, list_nodes, list_nodes_clean, cutting_edges):
    """
    Input : - G, the graph from which we extract the pattern graph. 
            - list_nodes, the list of the nodes that we want to extract using pdb format
            - list_nodes_clean, the list of the nodes that we want to extract but still as exctrated in the FR3D UNIT format, useful for labelling.
            - cutting_edges, the place where the backbone B53 is broken in the pattern that we want to extract expressed as a list of couples in [|1, n|]**2.
    Output : Twos graphs, the Gnew graph serves directly as pattern and the Gnewnx graph serves as target for research of pattern inside other patterns.
    """
    Gnew=nx.DiGraph() #Initiate the new GT graph with label from 1 to n to serve as pattern.
    Gnewnx=nx.DiGraph() #Initiate the new GT graph with label similar to the one in the RNA to serve as target.
    list_nodes_pdb = {}
    for ((i, j), t) in G.nodes.data():
        if (i, t['pdb_position']) in list_nodes: #Important to use the pdb_position here to check and not j
            list_nodes_pdb[(i, j)] = (i, t['pdb_position']) #We keep track of which pdb position correspond to which j
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


def extract_with_pdb_number_and_gaps(G, list_nodes, list_nodes_clean, cutting_edges):
    """
    Input : - G, the graph from which we extract the pattern graph. 
            - list_nodes, the list of the nodes that we want to extract using pdb format
            - list_nodes_clean, the list of the nodes that we want to extract but still as exctrated in the FR3D UNIT format, useful for labelling.
            - cutting_edges, the place where the backbone B53 is broken in the pattern that we want to extract expressed as a list of couples in [|1, n|]**2.
    Output : Twos graphs, the Gnew graph serves directly as pattern and the Gnewnx graph serves as target for research of pattern inside other patterns.
    """
    Gnew=nx.DiGraph()
    node_list = [-1]
    index = 1
    for (i, j) in list_nodes:
        k, t = [(kk,tt) for ((ii, kk), tt) in G.nodes.data() if tt['pdb_position'] == j and i == ii][0]
        Gnew.add_node(index, pdb_position = t['pdb_position'], atoms = t['atoms'], nt = t['nt'], chain = i)
        node_list.append((i, k))
        if (i, j) == list_nodes[-1]:
            iter_node = None
        elif ((i,j),list_nodes[list_nodes.index((i, j)) + 1]) in cutting_edges:
            iter_node = None
        else:
            B53_neighbors=[n for n in G.successors((i, k)) if G[(i, k)][n]['label'] == 'B53']
            if len(B53_neighbors) > 1: #It means that two backbones start from iter_node, which is not biologically admissible.
                print("THE WORLD BLOWS UP1")
            if len(B53_neighbors) == 0: 
                print("THE WORLD BLOWS UP2")
            iter_node = B53_neighbors[0]
            succ1, succ2 = list_nodes[list_nodes.index((i,j)) + 1]
            succ = (succ1, int(succ2))
        index +=1
        while iter_node:
            c1, blub = iter_node
            t = [tt for (ii, tt) in G.nodes.data() if ii == iter_node][0]
            if (c1, int(t['pdb_position'])) == succ:
                iter_node = None
            else:
                Gnew.add_node(index, pdb_position = t['pdb_position'], atoms = t['atoms'], nt = t['nt'], chain = c1)     
                node_list.append(iter_node)
                index +=1
                B53_neighbors=[n for n in G.successors(iter_node) if G[iter_node][n]['label'] == 'B53']
                if len(B53_neighbors) > 1: #It means that two backbones start from iter_node, which is not biologically admissible.
                    print("THE WORLD BLOWS UP3")
                if len(B53_neighbors) == 0: 
                    print("THE WORLD BLOWS UP4")
                iter_node = B53_neighbors[0]
    for i in node_list:
        for j in node_list:
            potential_edge = [(ii,jj,tt) for (ii,jj,tt) in G.edges.data() if i==ii and j==jj]
            if len(potential_edge) ==  1:
                (i,j,t) = [(ii,jj,tt) for (ii,jj,tt) in G.edges.data() if i==ii and j==jj][0]
                indexi = node_list.index(i)
                indexj = node_list.index(j)
                Gnew.add_edge(indexi, indexj, label=t['label'], near=t['near'])
            elif len(potential_edge) > 1:
                print("THE WORLD BLOWS UP")
    return Gnew, Gnew

def extractor(Gpath, Gnewname, list_nodes, cutting_edges, pattern_place="kinkturnpattern/", target_place ="kinkturntarget/", list_nodes_clean = [], withgaps = 0):
    """
    Input : - Gpath, the graph path from which we extract the pattern graph. 
            - Gnewname, the name for the pattern and target graphs that we are extracted
            - list_nodes, the list of the nodes that we want to extract using pdb format
            - list_nodes_clean, the list of the nodes that we want to extract but still as exctrated in the FR3D UNIT format, useful for labelling.
            - cutting_edges, the place where the backbone B53 is broken in the pattern that we want to extract expressed as a list of couples in [|1, n|]**2.
            - pattern_place and target_place, the places where respectively put the extracted pattern and target graphs.
    Output : Dump twos graphs in a pickle and nx pickle format, the Gnew graph serves directly as pattern and the Gnewnx graph serves as target for research of pattern inside other patterns.
    """
    with open(Gpath,'rb') as fG:
        G = pickle.load(fG)
        if DEBUG:
            print("Gpath", Gpath)
    if withgaps:
        Gnew, Gnewnx = extract_with_pdb_number_and_gaps(G, list_nodes, list_nodes_clean, cutting_edges)
    else:
        Gnew, Gnewnx = extract_with_pdb_number(G, list_nodes, list_nodes_clean, cutting_edges)
    if DEBUG:
        print([(i) for (i,t) in Gnew.nodes.data()])
        print([(i,j, t['label']) for (i,j,t) in Gnew.edges.data() if t['label'] != 'B53'])
        print([(i,j, t['label']) for (i,j,t) in Gnew.edges.data() if t['label'] == 'B53'])

    with open(pattern_place +Gnewname + ".pickle", 'wb') as ff:
        pickle.dump(Gnew, ff)
    with open(target_place + Gnewname + ".nxpickle", 'wb') as ff2:
        pickle.dump(Gnewnx, ff2)

def csv_parse(family_name, break_list_entry, RNAstorage = "bigRNAstorage/", csvlocation = "RNAcsv/", pattern_place="ALLkinkturnpattern/", target_place ="ALLkinkturntarget/", withgaps = 0):
    """
    Input : - family_name, the name of the family of RNA from which we extract the graphs
            - Gnewname, the name for the pattern and target graphs that we are extracted
            - break_list, a precursor for cutting_edges with a more natural numerotation
            - csvlocation, place where csv are stored.
            - pattern_place and target_place, the places where respectively put the extracted pattern and target graphs.
    Output : Extract the pattern graphs from the csv and return also the list of all mapping that we are waiting for for this family.
    """
    #break_list is the couples of nucleotides between which there is a break.
    # We consider that nucleotides are numeroted from 1 to n
    if break_list_entry == -1:
        annotated = 1
    else:
        annotated = 0
    with open(csvlocation + family_name + '.csv', newline='', encoding='utf-8') as f:
        reader = csv.reader(f)
        resu = []
        for k,prerow in enumerate(reader):
            if annotated:
                for i, e in enumerate(prerow):
                    if e != '':
                        prebreak = e.split('|')
                        stoper = i
                break_list = [int(elem) for elem in prebreak]
                row = prerow[:stoper]
            else:
                break_list = [i for (i, j) in break_list_entry]
                row = prerow
            symmetry = 0
            RNAtarget = row[0].split('|')[0]
            Gpath = RNAstorage + RNAtarget + ".nxpickle"
            list_nodes = []
            list_nodes_clean = []
            cutting_edges = []
            perfect_mapping = []
            for kk, nucleo in enumerate(row):
                explore_list = nucleo.split('|')
                if len(explore_list) > 8:
                    symmetry = 1
                if len(explore_list) > 7:
                    pdb_string = explore_list[4] + explore_list[7]
                else:
                    pdb_string = explore_list[4] 
                list_nodes_clean.append((explore_list[2], int(explore_list[4])))
                list_nodes.append((explore_list[2], pdb_string))
                perfect_mapping.append((kk + 1, (explore_list[2], pdb_string)))
            if DEBUG:
                print(perfect_mapping)
            for i in break_list:
                cutting_edges.append((list_nodes[i - 1], list_nodes[i]))
            if symmetry == 0:
                extractor(Gpath, str(k) + family_name + 'into' + RNAtarget, list_nodes, cutting_edges, pattern_place=pattern_place, target_place =target_place, list_nodes_clean = list_nodes_clean, withgaps = withgaps)
                resu.append((RNAtarget, perfect_mapping))
            else:
                if DEBUG:
                    print("SKIPPED DUE TO SYMMETRY")
        return resu

