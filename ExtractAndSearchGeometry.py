import networkx as nx 
import networkx.algorithms.isomorphism as iso
import pickle
from TestFuzzTree import outer_chain_removal
from collections import defaultdict
        



def edge_match(d_g1, d_g2):
    return d_g1['label'][:3] == d_g2['label'][:3]


def line_graph_to_graph(mapping, g):
    new_g = nx.DiGraph()
    for k, v in mapping.items():
        new_g.add_edge(k[0], k[1], **g.get_edge_data(*k))
    d = {}
    for n in new_g.nodes():
        d[n] = g.nodes(data=True)[n]
        #new_g.node[n] = g.node[n]
    nx.set_node_attributes(new_g, d)

    return new_g


def make_line_graph(g):
    lg = nx.line_graph(g)
    d = {n:g.get_edge_data(*n)['label'] for n in lg.nodes}
    nx.set_node_attributes(lg, values=d, name='label')
    return lg


def has_subgraph(g1, g2):
    """G2 is included in G1"""
    M = nx.isomorphism.DiGraphMatcher(make_line_graph(g1),
                                        make_line_graph(g2), node_match=edge_match)
    return M.subgraph_is_isomorphic()


def find_all_subgraph(g, subgraph):
    found = set()
    new_data = []
    M = nx.isomorphism.DiGraphMatcher(make_line_graph(g),
                                      make_line_graph(subgraph),
                                      node_match=edge_match)
    for m in M.subgraph_isomorphisms_iter():
        iso = line_graph_to_graph(m, g)

        name = tuple(sorted(iso.edges()))
        if name not in found:
            found.add(name)
            new_data.append(iso)

    return new_data

def VF2(GP, GT):
    return find_all_subgraph(GT, GP)

def eliminate_similar_geometry(motifs_to_search):
    """
    Input: - A  list of geometry that can serve as motif to be search exhaustievely in same RNA and in others RNAs. 
    Output: Return same list limited to unique geometry by removing the ones that are indicated multiple times.
    """
    em = iso.categorical_edge_match("label", "cWW")
    resu = []
    for i1 in range(len(motifs_to_search)):
        toadd = 1
        G1 = motifs_to_search[i1]
        for i2 in range(i1 + 1, len(motifs_to_search)):
            if nx.is_isomorphic(G1, motifs_to_search[i2], edge_match=em):
                toadd = 0
                break
        if toadd:
            resu.append(G1)
    return resu

def ff(tuple):
    (a,b) = tuple
    return a

def abstract_in_geometry(GT, mappings, cutting_edges):#, GTlistfolder = "bigRNAstorage/"):
    """
    Input: - The graph target GT.
           - The list of mappings obtained as a request from a pattern graph in a target graph.
           - cutting_edges, list of numbers where the backbone has to be separated.
    Output: A minimal list of geometry that can serve as motif to be search exhaustievely in same RNA and in others RNAs. 
    """
    #GTpath = GTlistfolder + GTpath + ".nxpickle"
    #with open(GTpath,'rb') as f:
    #    GT = pickle.load(f)
    new_mappings = []
    for mapp in mappings:
        mapp.sort()
        new_mapp = []
        for (i,j) in mapp:
            new_mapp.append(j)
        new_mappings.append(new_mapp)
    #TODO: we assume for now that mappings have a length above 2.
    motifs_to_search = []
    for mapp in new_mappings:
        Gnew=nx.DiGraph()
        node_list = [-1]
        index = 1
        for i in mapp:
            c, blub = i
            t = [tt for (ii, tt) in GT.nodes.data() if ii == i][0]
            Gnew.add_node(index, pdb_position = t['pdb_position'], atoms = t['atoms'])
            node_list.append(i)
            if mapp.index(i) + 1 in cutting_edges or i == mapp[-1]: #mapp.index(i) + 1(c, t['pdb_position'])
                iter_node = None
            else:
                B53_neighbors=[n for n in GT.successors(i) if GT[i][n]['label'] == 'B53']
                if len(B53_neighbors) > 1: #It means that two backbones start from iter_node, which is not biologically admissible.
                    print("THE WORLD BLOWS UP")
                if len(B53_neighbors) == 0: 
                    print("THE WORLD BLOWS UP")
                iter_node = B53_neighbors[0]
                succ = mapp[mapp.index(i) + 1]
            index +=1
            while iter_node:
                t = [tt for (ii, tt) in GT.nodes.data() if ii == iter_node][0]
                if iter_node == succ:
                    iter_node = None
                else:
                    Gnew.add_node(index, pdb_position = t['pdb_position'], atoms = t['atoms'])     
                    node_list.append(iter_node)
                    index +=1
                    B53_neighbors=[n for n in GT.successors(iter_node) if GT[iter_node][n]['label'] == 'B53']
                    if len(B53_neighbors) > 1: #It means that two backbones start from iter_node, which is not biologically admissible.
                        print("THE WORLD BLOWS UP")
                    if len(B53_neighbors) == 0: 
                        print("THE WORLD BLOWS UP")
                    iter_node = B53_neighbors[0]
        for i in node_list:
            for j in node_list:
                potential_edge = [(ii,jj,tt) for (ii,jj,tt) in GT.edges.data() if i==ii and j==jj]
                if len(potential_edge) ==  1:
                    (i,j,t) = [(ii,jj,tt) for (ii,jj,tt) in GT.edges.data() if i==ii and j==jj][0]
                    indexi = node_list.index(i)
                    indexj = node_list.index(j)
                    Gnew.add_edge(indexi, indexj, label=t['label'], near=t['near'])
                elif len(potential_edge) > 1:
                    print("THE WORLD BLOWS UP")

        motifs_to_search.append(Gnew)
    print(len(mappings))
    print(len(motifs_to_search))
    motifs_to_search = eliminate_similar_geometry(motifs_to_search)
    print(len(motifs_to_search))
    return motifs_to_search

def look_at_all_occurences(GT, chains, mappings, cutting_edges):
    """
    Input: - The graph target GT.
           - The list of mappings obtained as a request from a pattern graph in a target graph.
           - cutting_edges, list of numbers where the backbone has to be separated.
    Output: Seach with VF2 all occurences of the minimal list of geometry obtained from mappings. 
    """
    #GTpath = GTlistfolder + GTpath + ".nxpickle"

    motifs_to_search = abstract_in_geometry(GT, mappings, cutting_edges)#, GTlistfolder)
    resu = []
    for GP in motifs_to_search:
        print(GP.nodes())
        inst = VF2(GP, GT) 
        for G in inst:
            print(G.nodes())
            local_mapping = []
            falseind = 0
            for (((i, j), t)) in G.nodes.data():
                local_mapping.append((falseind,(i, j)))
                falseind+=1
        resu.append(local_mapping.copy())
        #TODO : VF2 to be implemented
    #TODO ; still have to purge "doublons" ? not sure as geometric form are distinct.
    print("resu", resu)
    return resu



def compute_metrics(ref_mappings, GT, occurences, strong_mapping = 1 ):
    ref_unfold = []
    mapping_unfold = []
    TP = 0
    for mapping_ref in ref_mappings:
        ref_unfold+= [j for (i,j) in mapping_ref]
    ref_unfold = list(set(ref_unfold))
    for mapping in occurences:
        for (i, j) in mapping:
            ii,t = [(ii, tt['pdb_position']) for ((ii,jj), tt) in GT.nodes.data() if (ii,jj) == j][0]
            mapping_unfold.append((ii, t))
    mapping_unfold = list(set(mapping_unfold))
    TP = len([i for i in mapping_unfold if i in ref_unfold])
    specificity = TP / len(mapping_unfold)
    sensitivity = TP / len(ref_unfold)
    F = 2 * TP /(len(mapping_unfold) + len(ref_unfold))
    return(specificity, sensitivity, F)

example = [[(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 4)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 4)), (5, ('B', 7)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 5)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 5)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 4)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 4)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 5)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 4)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 4)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 4)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 4)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 4)), (5, ('B', 7)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 4)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 4)), (5, ('B', 7)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 5)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 5)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 4)), (5, ('B', 7)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 4)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 5)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 4)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 5)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 4)), (5, ('B', 7)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 4)), (5, ('B', 7)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 4)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 4)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 4)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 4)), (5, ('B', 7)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 4)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 5)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 4)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 5)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 4)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 4)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 4)), (5, ('B', 7)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 5)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 4)), (5, ('B', 7)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 4)), (5, ('B', 7)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 5)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 4)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 4)), (5, ('B', 7)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 4)), (5, ('B', 7)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 4)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 4)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 4)), (5, ('B', 7)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 5)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 4)), (5, ('B', 7)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 5)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 5)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 4)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 3)), (7, ('A', 4)), (5, ('B', 6)), (8, ('A', 6)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 11)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 4)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 1))], [(12, ('A', 10)), (11, ('A', 9)), (6, ('A', 2)), (7, ('A', 3)), (5, ('B', 7)), (8, ('A', 4)), (9, ('A', 7)), (10, ('A', 8)), (4, ('B', 5)), (3, ('B', 4)), (2, ('B', 3)), (1, ('B', 2))]]

#abstract_in_geometry('5G4T', example, [5], GTlistfolder = "bigRNAstorage/")

#perfect_mapping = csv_parse("IL_29549.9", [(5,6)])
#perfect_mapping = [perfect_mapping[i] for i in range(len(perfect_mapping)) if perfect_mapping[i][0] in ['5G4T']]

#occ = look_at_all_occurences('5G4T', ['A', 'B'], example, [5], GTlistfolder = "bigRNAstorage/")
#print("\nperfect", perfect_mapping)
#compute_metrics(perfect_mapping, occ, strong_mapping = 1, GTlistfolder = "bigRNAstorage/")
from Extractor import csv_parse
from TestFuzzTree import initialise_perfect_mapping
import os
def full_metrics( dict_mappings, GTlistfolder = "bigRNAstorage"):
    resu = []
    perfect_mapping, full_cut = csv_parse("kink_turn", -1, return_cutting_edges=1)
    perfect_mapping = [perfect_mapping[i] for i in range(len(perfect_mapping)) if perfect_mapping[i][0] in ['5G4T']]
    perfect_mapping = initialise_perfect_mapping(perfect_mapping, [])
    GT_list = []
    chains_list = []
    RNA_list = []
    cutting_list = []
    path = os.path.abspath(os.getcwd()) + "/" + GTlistfolder
    for (RNAname, chains) in perfect_mapping.keys():
        with open(path+ "/" + RNAname + '.nxpickle','rb') as f:
            GT = pickle.load(f)
        GT = outer_chain_removal(GT, chains)
        GT_list.append(GT)
        RNA_list.append(RNAname)
        chains_list.append(list(chains))
        precut = [li for (nam, li) in full_cut if nam == RNAname]
        cut = [ ]
        for li in precut:
            cut+=li
        cutting_list.append(cut)
    print("BLUB", GT_list, chains_list, RNA_list, cutting_list)
    for i in range(len(GT_list)):
        mappings = dict_mappings[RNA_list[i]]
        occ = look_at_all_occurences(GT_list[i], chains_list[i], mappings, cutting_list[i])
        ref_mappings = perfect_mapping[(RNA_list[i], tuple(chains_list[i]))]
        loc = compute_metrics(ref_mappings, GT_list[i], occ, strong_mapping = 1)
        resu.append((RNA_list[i], chains_list[i],loc))
    print("resu", resu)
    return resu

dicto = {}
dicto['5G4T'] = example
full_metrics( dicto, GTlistfolder = "bigRNAstorage")

dicto = {}
dicto['4LFB'] = example2[3]
full_metrics( dicto, GTlistfolder = "bigRNAstorage")
