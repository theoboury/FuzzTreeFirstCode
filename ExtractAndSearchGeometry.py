import networkx as nx 
import networkx.algorithms.isomorphism as iso

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
        new_mapp = []
        for (i,j) in mapp:
            new_mapp.append(j)
        new_mappings.append(new_mapp)
    new_mappings = list(set(new_mappings))
    #TODO: we assume for now that mappings have a length above 2.
    motifs_to_search = []
    for mapp in new_mappings:
        Gnew=nx.DiGraph()
        node_list = [-1]
        index = 1
        for i in mapp:
            t = [tt for (ii, tt) in GT.nodes.data() if ii == i][0]
            Gnew.add_node(index, pdb_position = t['pdb_position'], atoms = t['atoms'], nt = t['nt'])
            node_list.append(i)
            if index in cutting_edges or i == mapp[-1]:
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
                Gnew.add_node(index, pdb_position = t['pdb_position'], atoms = t['atoms'], nt = t['nt'])     
                node_list.append(iter_node)
                index +=1
                if iter_node == succ:
                    iter_node = None
                else:
                    B53_neighbors=[n for n in GT.successors(iter_node) if GT[iter_node][n]['label'] == 'B53']
                    if len(B53_neighbors) > 1: #It means that two backbones start from iter_node, which is not biologically admissible.
                        print("THE WORLD BLOWS UP")
                    if len(B53_neighbors) == 0: 
                        print("THE WORLD BLOWS UP")
                    iter_node = B53_neighbors[0]
        for i in node_list:
            for j in node_list:
                (i,j,t) = [(ii,jj,tt) for (ii,jj,tt) in GT.edges.data() if i==ii and j==jj][0]
                indexi = node_list.index(i)
                indexj = node_list.index(j)
                Gnew.add_edge(indexi, indexj, label=t['label'], near=t['near'])
        motifs_to_search.append(Gnew)
    motifs_to_search = eliminate_similar_geometry(motifs_to_search)
    return motifs_to_search
