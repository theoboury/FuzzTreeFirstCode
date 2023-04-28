# FuzzTree
# Copyright (C) 2023 THEO BOURY 

import networkx as nx 
import networkx.algorithms.isomorphism as iso
import pickle
from multiprocessing import Pool
from TestFuzzTree import outer_chain_removal, initialise_perfect_mapping
from Extractor import csv_parse
import os 

DEBUG=0

from TestFuzzTree import rename_author_position

def edge_match(d_g1, d_g2):
    return d_g1['label'][:3] == d_g2['label'][:3]

def line_graph_to_graph(mapping, g):
    new_g = nx.DiGraph()
    for k, v in mapping.items():
        new_g.add_edge(k[0], k[1], **g.get_edge_data(*k))
    d = {}
    for n in new_g.nodes():
        d[n] = g.nodes(data=True)[n]
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


def abstract_in_geometry(GT, mappings, cutting_edges):
    """
    Input: - The graph target GT.
           - The list of mappings obtained as a request from a pattern graph in a target graph.
           - cutting_edges, list of numbers where the backbone has to be separated.
    Output: A minimal list of geometry that can serve as motif to be search exhaustievely in same RNA and in others RNAs. 
    """
    new_mappings = []
    for mapp in mappings:
        mapp.sort()
        new_mapp = []
        for (i,j) in mapp:
            new_mapp.append(j)
        new_mappings.append(new_mapp)
    #We assume here that mappings have a length above 2.
    motifs_to_search = []
    for mapp in new_mappings:
        Gnew=nx.DiGraph()
        node_list = [-1]
        index = 1
        for i in mapp:
            t = [tt for (ii, tt) in GT.nodes.data() if ii == i][0]
            Gnew.add_node(index, pdb_position = t['pdb_position'], atoms = t['atoms'])
            node_list.append(i)
            if mapp.index(i) + 1 in cutting_edges or i == mapp[-1]: 
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
    motifs_to_search = eliminate_similar_geometry(motifs_to_search)
    return motifs_to_search

def look_at_all_occurences(GT, chains, mappings, cutting_edges):
    """
    Input: - The graph target GT.
           - The list of mappings obtained as a request from a pattern graph in a target graph.
           - cutting_edges, list of numbers where the backbone has to be separated.
    Output: Seach with VF2 all occurences of the minimal list of geometry obtained from mappings. 
    """
    motifs_to_search = abstract_in_geometry(GT, mappings, cutting_edges)
    resu = []
    for GP in motifs_to_search:
        inst = VF2(GP, GT) 
        for G in inst:
            local_mapping = []
            falseind = 0
            for (((i, j), t)) in G.nodes.data():
                local_mapping.append((falseind,(i, j)))
                falseind+=1
            resu.append(local_mapping.copy())
    #No need to purge "doublons" as geometric forms are distinct.
    return resu



def compute_metrics(ref_mappings, GT, occurences):
    """
    Input: - ref_mappings, the mapping of references that is known for a given Kink-Turn motif. 
           - The graph target GT in which we are looking this Kink-Turn motif.
           - occurences, the list of all occurences found for this Kink-Turn motif.
    Output: A tuple precision, specificity, sensitivity, F score, number of found motifs, list of found motifs 
    for current instance composed of one Kink-Turn motif and all the occurences that we found for it with our method.
    """
    if len(occurences) == 0:
        return (0, 0, 0, 0, 0, [])
    ref_unfold = []
    mapping_unfold = []
    TP = 0
    for mapping_ref in ref_mappings:
        ref_unfold+= [j for (i,j) in mapping_ref]
    ref_unfold = list(set(ref_unfold))
    occ_clean = []
    for mapping in occurences:
        mapping_bis = mapping.copy()
        mapping_bis.sort()
        if mapping_bis not in occ_clean:
            occ_clean.append(mapping_bis)
    for mapping in occurences:
        for (i, j) in mapping:
            ii,t = [(ii, tt['pdb_position']) for ((ii,jj), tt) in GT.nodes.data() if (ii,jj) == j][0]
            mapping_unfold.append((ii, t))
    mapping_unfold = list(set(mapping_unfold))
    TP = len([i for i in mapping_unfold if i in ref_unfold])
    FP = len(mapping_unfold) - TP
    NotReferenced = len(GT.nodes()) - len(ref_unfold)
    precision = TP / len(mapping_unfold)
    specificity = (NotReferenced - FP)/ NotReferenced
    sensitivity = TP / len(ref_unfold)
    F = 2 * TP /(len(mapping_unfold) + len(ref_unfold))
    return (precision, specificity, sensitivity, F, len(occ_clean), occ_clean)


def wrapper_metrics(mappings_ref_mappings_GT_listi_chains_listi_cutting_listi_RNA_listi):
    """
    A simple wrapper arond the computation of the metrics in order to be able to multiprocess the computation.
    """
    (mappings, ref_mappings, GT_listi, chains_listi, cutting_listi, RNA_listi) = mappings_ref_mappings_GT_listi_chains_listi_cutting_listi_RNA_listi
    occ = look_at_all_occurences(GT_listi, chains_listi, mappings, cutting_listi)
    loc = compute_metrics(ref_mappings, GT_listi, occ)
    print("\nresu_by_RNA", (RNA_listi, chains_listi, loc))
    return (RNA_listi, chains_listi,loc)
        
def full_metrics(dict_mappings, GTlistfolder = "bigRNAstorage", csvtostudy = "kink_turn", nb_procs = 1, cutting_edge = []):
    """
    Input: - dict_mappings, the list of all found mappings for studied Kink-Turn motifs obtained with the FuzzTree method. 
           - GTlistfolder, location of full RNA graphs.
           - csvtostudy, the csv from which the motif of references known for the Kink-Turn are extracted.
           - nb_procs, the number of processors to use to reduce time of computation.
           - cutting_edges, list of numbers where the backbone has to be separated in the pattern graph that wes used to found the motifs in the FuzzTree method.
    Output: Compute for each Kink-Turn motif of reference a tuple precision, specificity, sensitivity, F score, number of found motifs, list of found motifs. 
    """
    resu = [("precision", "specificity", "sensitivity", "F", "nb_found_motifs", "found_motifs")]
    perfect_mapping = csv_parse(csvtostudy, -1)
    perfect_mapping = initialise_perfect_mapping(perfect_mapping, [])
    new_perfect_mapping = {}
    for (RNAname, chains) in perfect_mapping.keys():
        new_mapping = perfect_mapping[(RNAname, chains)]
        if RNAname in new_perfect_mapping.keys():
            (cha, mapper) = new_perfect_mapping[RNAname]
            new_perfect_mapping[RNAname] = (tuple(list(cha) + list(chains)), mapper + new_mapping)
        else:
            new_perfect_mapping[RNAname] = (chains, new_mapping)
    perfect_mapping = {}
    for RNAname in new_perfect_mapping.keys():
        (cha, mapper) = new_perfect_mapping[RNAname]
        perfect_mapping[(RNAname, cha)] = mapper
    GT_list = []
    chains_list = []
    RNA_list = []
    path = os.path.abspath(os.getcwd()) + "/" + GTlistfolder
    pathbis = os.path.abspath(os.getcwd()) + "/bigRNAstorage"
    for (RNAname, chains) in perfect_mapping.keys():
        if GTlistfolder != "bigRNAstorage": 
            with open(path+ "/" + RNAname + '.pickle','rb') as f1:
                GT = pickle.load(f1)
            with open(pathbis + "/" + RNAname + '.nxpickle','rb') as f2:
                GTref = pickle.load(f2)
            GT = rename_author_position(GT, GTref)
        else:
            with open(path+ "/" + RNAname + '.nxpickle','rb') as f:
                GT = pickle.load(f)
        GT = outer_chain_removal(GT, chains)
        GT_list.append(GT)
        RNA_list.append(RNAname)
        chains_list.append(list(chains))
    if DEBUG:
        print("Temporary lists to study", [G.nodes() for G in GT_list], chains_list, RNA_list)
    entry = []
    for i in range(len(GT_list)):
        mappings = dict_mappings[RNA_list[i]]
        ref_mappings = perfect_mapping[(RNA_list[i], tuple(chains_list[i]))]
        entry.append((mappings, ref_mappings, GT_list[i], chains_list[i], cutting_edge, RNA_list[i]))
    with Pool(nb_procs, maxtasksperchild=1) as pool:
        resu += list(pool.imap_unordered(wrapper_metrics, entry))
    print("\nresu", resu)
    return resu




