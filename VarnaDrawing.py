# FuzzTree
# Copyright (C) 2023 THEO BOURY 

import varnaapi
from varnaapi.param import BaseAnnotation
from FuzzTree import main

DEBUG=0


def traduction(label_char):
    """
    Input : - An interaction label as 3 capital letters label_char
    Output : - An interaction label foamlized as something readable by VARNA.
    """
    resu = []
    if label_char[0] == 'C':
        resu.append('cis')
    else:
        resu.append('trans')
    for i in range(1, 3): 
        if label_char[i] == 'W':
            resu.append('wc')
        elif label_char[i] == 'H':
            resu.append('h')
        else:
            resu.append('s')
    return resu[1], resu[2], resu[0]

def draw_varna(GT, nodes_target, mapping):
    """
        Input :
               - Graph Target GT.
               - The list of all nodes in GT nodes_target, nodes_target should be ordered in lexicographic order.
               - The mapping that we want to represent on GT.
        Output : 
        - Return GT graph with mapping drawn in red using VARNA.
    """
    secondary_structures = [(nodes_target.index(i), nodes_target.index(j)) for (i, j, t) in GT.edges.data() if t['label'] == "CWW"]
    noncanonicaledges = [(nodes_target.index(i), nodes_target.index(j), t['label']) for (i, j, t) in GT.edges.data() if (nodes_target.index(i), nodes_target.index(j)) not in secondary_structures if t['label'] !='B53']
    sequence  = []
    if nodes_target != []:
        strain = nodes_target[0][0]
    for node in nodes_target:
        (i, t) = [(i, t) for (i, t) in GT.nodes.data() if i == node][0]
        if node[0] != strain:
            sequence += "&" #We add this symbol when we change of strain in VARNA
            strain = node[0]
        sequence += t['nt'] #Sequence is otherwise the list of nucleotides in order
    sequence = "".join(sequence)
    v  = varnaapi.Structure(sequence, secondary_structures)
    for (i,j,t) in noncanonicaledges:
        e5,e3,ster = traduction(t)
        if DEBUG:
            print("(i, j) =", (i, j), "t", t, "edge5 = ", e5,"edge3 = ", e3, "stericity = ", ster)
        if i > j:
            v.add_aux_BP(j, i, edge5=e3, edge3=e5, stericity=ster, color='blue', thickness=1) #We add now non canonical edges       
        else:
            v.add_aux_BP(i, j, edge5=e5, edge3=e3, stericity=ster, color='blue', thickness=1) #We add now non canonical edges
    
        #We now color the part that correspond to the mapping using different styles
    style1 = varnaapi.param.BasesStyle(fill="#FF0000", label="#FFFFFF")
    style2 = varnaapi.param.BasesStyle(fill="#FFFFFF", outline="#000000")
    mapped = [nodes_target.index(j) for i,j in mapping]
    for (i,j) in mapping:
        if DEBUG:
            print(str(i), nodes_target.index(j))
        v.add_annotation(BaseAnnotation(str(i), nodes_target.index(j) + 1, color="#FF0000", size=10))
    v.add_bases_style(style1, mapped)
    v.add_bases_style(style2, [i for i in range(len(nodes_target)) if i not in mapped])
    return v
    
    

def print_mapping_on_target_graph(GP, GT, mapping = [], output_format = "pdf", name_file = "", show=1, L=0, E=0, G=0, maxGAPdistance=3, nb_samples=1000, respect_injectivity=1, D = 5, Distancer_preprocessed = {}, nb_procs = 1):
    """
    Input : - Two graphs, the pattern graph GP and the target graph GT. 
            - A mapping between GP and GT, if this mapping is not specified, it is calculated here.
            - output_format is the type of file that we want to save, specify None for no save.
            - name_file is the name of the output graph if specified.
            - show is set to 1 if we want to observe the output with directly.
            - L, E, G are fuzzyness options specified for the main function if mapping is not given.
    Output : Build the GT graph with a red colored part for each node in GP mapped in GT, can save it as a file or plot it depending on the options selected.
    """
    nodes_target = list(GT.nodes())
    nodes_target.sort()
    if mapping == []:
        mapping = main(GP, GT, L, E, G, maxGAPdistance=maxGAPdistance, nb_samples=nb_samples, respect_injectivity=respect_injectivity, D = D, Distancer_preprocessed = Distancer_preprocessed, nb_procs = nb_procs)
        mapping = [map for (map, length) in mapping if length == len(GP.nodes())]
        mapping = mapping[0] #We take only the first result
    v = draw_varna(GT, nodes_target, mapping)
    if output_format:
        if name_file == "":
            name = "GPintoGT" + "." + output_format
        else:
            name = name_file + "." + output_format
        v.savefig(name)
    if show:
        v.show()




