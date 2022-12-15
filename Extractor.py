import os, glob, pickle
import networkx as nx
import csv

DEBUG=1


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

def extractor(Gpath, Gnewname, list_nodes, cutting_edges, pattern_place="kinkturnpattern/", target_place ="kinkturntarget/", list_nodes_clean = []):
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
    Gnew, Gnewnx = extract_with_pdb_number(G, list_nodes, list_nodes_clean, cutting_edges)
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

def csv_parse(family_name, break_list, RNAstorage = "bigRNAstorage/", csvlocation = "RNAcsv/", pattern_place="ALLkinkturnpattern/", target_place ="ALLkinkturntarget/"):
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
    with open(csvlocation + family_name + '.csv', newline='', encoding='utf-8') as f:
        reader = csv.reader(f)
        resu = []
        for k,row in enumerate(reader):
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
            for (i, j) in break_list:
                cutting_edges.append((list_nodes[i - 1], list_nodes[j - 1]))
            if symmetry == 0:
                extractor(Gpath, str(k) + family_name + 'into' + RNAtarget, list_nodes, cutting_edges, pattern_place=pattern_place, target_place =target_place, list_nodes_clean = list_nodes_clean)
                resu.append((RNAtarget, perfect_mapping))
            else:
                if DEBUG:
                    print("SKIPPED DUE TO SYMMETRY")
        return resu

