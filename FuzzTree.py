import infrared as ir
import infrared.rna as rna
import matplotlib.pyplot as plt

import re
import collections
import math
import sys
import os

#!export PYTHONPATH="${PYTHONPATH}:/home/uqamportable/miniconda3/lib/python3.9/site-packages/infrared/"
from infrared import def_constraint_class, def_function_class

import pickle
import networkx as nx
import varnaapi


def interaction_to_number(interact_char):
    """ 
    Input : A non-canonical interaction/label as a string, interact_char.
    Output : The corresponding index of this interaction in the Distance_matrix/ IDI_matrix 
            and -1 if it is a backbone.
    """
    if interact_char == 'B53':
        return -1
    if interact_char[2] == 'W':
        interact_char = interact_char[0] + interact_char[2] + interact_char[1]
    if interact_char[1] == 'S' and interact_char[2] == 'H':
        interact_char = interact_char[0] + 'HS'
    traduction = ['CHH', 'TWH', 'CWW', 'THS', 'CWS', 'CSS', 'CWH', 'CHS','TWS','TSS','TWW','THH']
    return traduction.index(interact_char)
    
def DistanceEdgesLabel(x, y, Distance_matrix, len_edges_pattern):
    """
    Input : - Two edges labels x and y respectively from GP and GT.
               - A Distance Matrix between non-canonical interactions/labels (for instance the IDI matricx for isostericity distance).
               - len_edges_patterns to scale the impact of errors with the size of the graph pattern.
    Output : A value between 0 and 1 that indicates how much these labels are close from each other.
    """
    (xx, yy) = interaction_to_number(x), interaction_to_number(y)
    if xx == yy: #including the case where xx == -1 (is a backbone)
        return 1
    if xx == -1 or yy == -1 or Distance_matrix[xx][yy] > 10: 
        return 0 #Too different or backbone not respected
    else:
        return 1 - Distance_matrix[xx][yy]/(20*len_edges_pattern)
    #TODO : carry on by distinguishing cases depending on distance on label
    #TODO : Les nombres placés devront être de la forme 0.5 /len(edges_pattern) pour garder la logique de pourcentages 


def _EdgeLabelRespect(x, y, label_edge_ij, len_edges_pattern, nodes_target, label_edge_target, edges_target, IDI_matrix, oriented):
    """
    Input : - Two nodes x and y from the GT respectively the mappings of i and j from the GP.
            - The label of edge between i and j, label_edge_ij.
            - len_edges_patterns to scale the impact of errors with the size of the graph pattern.
            - The list of the nodes in GT nodes_target to map the index of infrared functions with real names of nodes in the graph.
            - The list of all label of the edges in GT, label_edge_target.
            - The list of all edges in GT, edges_target.
            - A IDI Matrix giving distances between non-canonical interactions/labels (for instance the IDI matricx for isostericity distance).
            - oriented specifies if the graph is oriented (1 if it is the case)
    Output : A value between 0 and 1 that indicates if the edge is present or not but also how much these labels are close from each other.
    """
    if x == y: #If neighbors are mapped to the same node in GT, we can already reject 
        return 0 
        #TODO: Etendre ce critère d'injectivité à + ou - 1 le long du backbone ?
    if (not oriented):
        x, y =  min(x, y), max(x, y)
    if (nodes_target[x], nodes_target[y]) in edges_target: #We check that this edge is present or not in GT
        return DistanceEdgesLabel(label_edge_ij,label_edge_target[edges_target.index((nodes_target[x], nodes_target[y]))], IDI_matrix, len_edges_pattern)       
    return 0 #currently 0 here simply eliminates this pattern when the edge is not present
    #TODO : introduced more granularity about missing edges 



class EdgeLabelRespect(ir.infrared.WeightedFunction):
    """
    Constrain complementarity mapping of any pair of nodes (i,j) that are an edge in graph_pattern.
    ```
    Label(i, j,label_edge_ij, len_edges_pattern, nodes_target, label_edge_target, edges_target, IDI_matrix, oriented).
    ```
    The constraint is satisfied if label of the mapped edge (i,j) is an edge and if its label near or equal to the label of (i, j).
    """
def_function_class('EdgeLabelRespect', lambda i, j, label_edge_ij, len_edges_pattern, nodes_target, label_edge_target, edges_target, IDI_matrix, oriented: [i, j],
            lambda x, y, label_edge_ij, len_edges_pattern, nodes_target, label_edge_target, edges_target, IDI_matrix, oriented: _EdgeLabelRespect(x, y, label_edge_ij, len_edges_pattern, nodes_target, label_edge_target, edges_target, IDI_matrix, oriented),
            module=__name__)


def main(GP, GT, nb_samples=1000): 
    """
    Input : - Two graphs, the pattern graph GP and the target graph GT.
    	    - nb_samples, maximum numbers of samples allowed to look for a pattern in GT.
    Output : List of "pre"mapping between GP and GT (as injectivity is not respected yet) with the numbers of nodes of GP that are correctly mapped
            The mapping is correct if this number is equals to n_pattern.
    """

    #The matrix for the isostericity between non-canonical interactions
    IDI_matrix = [[8.9, 12.0, 14.7, 14.0, 13.7, 12.7, 15.1, 14.7, 16.2, 16.6, 16.2, 14.0], 
    [12.0, 2.6, 10.6, 9.7, 14.3, 15.6, 11.2, 15.2, 13.8, 15.4, 11.9, 11.4], 
    [14.7, 10.6, 4.1, 8.2, 9.2, 13.1, 14.5, 16.0, 12.4, 11.3, 11.1, 15.5], 
    [14.0, 9.7, 8.2, 2.1, 7.0, 12.7, 12.0, 12.1, 10.0, 11.9, 13.1, 15.8], 
    [13.7, 14.3, 9.2, 7.0, 3.5, 7.4, 14.9, 12.3, 10.9, 10.8, 14.6, 17.7], 
    [12.7, 15.6, 13.1, 12.7, 7.4, 1.3, 15.8, 12.9, 13.8, 12.0, 17.1, 19.0], 
    [15.1, 11.2, 14.5, 12.0, 14.9, 15.8, 3.2, 8.8, 8.4, 11.5, 10.6, 10.8], 
    [14.7, 15.2, 16.0, 12.1, 12.3, 12.9, 8.8, 2.4, 7.9, 11.2, 14.7, 14.9], 
    [16.2, 13.8, 12.4, 10.0, 10.9, 13.8, 8.4, 7.9, 3.4, 6.4, 9.6, 14.4], 
    [16.6, 15.4, 11.3, 11.9, 10.8, 12.0, 11.5, 11.2, 6.4, 2.2, 9.0, 14.4], 
    [16.2, 11.9, 11.1, 13.1, 14.6, 17.1, 10.6, 14.7, 9.6, 9.0, 3.8, 9.0],
    [14.0, 11.4, 15.5, 15.8, 17.7, 19.0, 10.8, 14.9, 14.4, 14.4, 9.0, 4.0]]
    
    n_pattern = len(GP.nodes)
    n_target = len(GT.nodes)

    edges_pattern = list(GP.edges())
    edges_target = list(GT.edges())

    nodes_pattern = list(GP.nodes())
    nodes_target = list(GT.nodes())

    label_edge_pattern = [t['label'] for (_, _, t) in GP.edges.data()]
    label_edge_target = [t['label'] for (_, _, t) in GT.edges.data()]
    
    #No label on vertices for now, can be introduced if needed.
    
    print("Graphs Infos\n n_pattern:", n_pattern, "\n n_target:", 
    n_target,"\n edges_pattern:", edges_pattern,"\n edges_target:", 
    edges_target,"\n nodes_pattern:", nodes_pattern,"\n nodes_target:", 
    nodes_target,"\n label_edge_pattern:", label_edge_pattern,
    "\n label_edge_target:", label_edge_target)

   
    model = ir.Model(n_pattern, n_target)


    #We define Label respect function on each edges in the pattern to check if the label in the target correspond to it.
    model.add_functions([EdgeLabelRespect(nodes_pattern.index(i), nodes_pattern.index(j), label_edge_pattern[k], len(edges_pattern), nodes_target, label_edge_target, edges_target, IDI_matrix, oriented=1) for k, (i,j) in enumerate(edges_pattern)], 'EdgeLabelRespect')
    #len(edges_pattern) here will allow us to scale the error depending on the size of the motif without touching the tolerancy indicated in the set_target part here.
    
    #We now take some samples of results given the built model.
    sampler = ir.Sampler(model)

    #Next we ensure that edges in the pattern are represented correctly for now (or enough in the future) in the target graph.
    sampler.set_target(1*len(edges_pattern), 0.999 , 'EdgeLabelRespect') 
    #With the tolerancy that is not scale on length of edges_pattern here, we ensure namely that if we put 0 as a return for EdgeLabelRespect then the pattern is automatically eliminated.
    #Other fuzzy "mistakes" must return values higher than 0 and "close" to 1 to stay in the long run as an admissible matching.
    
    samples = [sampler.targeted_sample() for _ in range(nb_samples)]
    resu = [([(nodes_pattern[k], nodes_target[x]) for k,x in enumerate(sample.values())], len(list(set(sample.values())))) for sample in samples]
    
    return resu






def print_mapping_on_target_graph(GP, GT, mapping = [], output_format = "pdf", name_file = "", show=1):
    """
    Input : - Two graphs, the pattern graph GP and the target graph GT 
            - A mapping between GP and GT, if this mapping is not specified, it is calculated here
            - output_format is the type of file that we want to save, specify None for no save.
            - name_file is the name of the output graph if specified.
            - show is set to 1 if we want to observe the output with matplotlib.
    Output : Build the GT graph with a red colored part for each node in GP mapped in GT, can save it as a file or plot it depending on the optons selected.
    """
    if mapping == []:
        mapping = main(GP, GT)
        mapping = [map for (map, length) in mapping if length == len(GP.nodes())]
        mapping = mapping[0] #We take only the first result
    GG = GT.copy()
    pos = nx.nx_agraph.graphviz_layout(GG, prog="sfdp") #The layout can be modified is the graph is too "packed" visually
    #pos = nx.spring_layout(GG, k=0.2)  #An alternative layout
    
    nx.draw(GG, pos, width=1, linewidths=1)
    
    # nodes drawing
    options = {"edgecolors": "tab:black", "node_size": 500, "alpha": 0.9}
    mapped = [j for (_, j) in mapping]
    nx.draw_networkx_nodes(GG, pos, nodelist=mapped, node_color="tab:red")
    nx.draw_networkx_nodes(GG, pos, nodelist=[i for i in list(GG.nodes()) if i not in mapped])

    # edges drawing
    nx.draw_networkx_edges(GG, pos)

    # labels drawing
    labels = {}
    for (i, j) in mapping:
        labels[j] = i
    nx.draw_networkx_labels(GG, pos, labels, font_color="whitesmoke")
    labels_edges = {}
    for (i, j, t) in GT.edges.data():
        labels_edges[(i,j)] = t['label']
    nx.draw_networkx_edge_labels(GG, pos, edge_labels = labels_edges, font_size=5)

    plt.axis("off")
    if output_format:
        if name_file == "":
            name = "GPintoGT" + "." + output_format
        else:
            name = name_file + "." + output_format
        plt.savefig(name) 
    if show:
        plt.show()



#TESTS PART FOR NOW

with open("rin147.pickle",'rb') as fP:
    GP = pickle.load(fP)

with open("6v9d.pickle",'rb') as fT:
    GT = pickle.load(fT)


#mapping = main(GP, GT)
#print(mapping)
#sur 10 samples on retrouve : ([(1, ('B', 7)), (2, ('B', 8)), (3, ('B', 12)), (4, ('B', 13)), (5, ('B', 17)), (6, ('B', 18)), (7, ('B', 23)), (8, ('B', 24))], 8)

print_mapping_on_target_graph(GP, GT, 
mapping = [(1, ('B', 7)), (2, ('B', 8)), (3, ('B', 12)), 
(4, ('B', 13)), (5, ('B', 17)), (6, ('B', 18)), (7, ('B', 23)), 
(8, ('B', 24))], name_file="rin147into6v9d", output_format="png")

#En chronométrant à la main le nombre de samples (sachant qu'ils sont print) ne semblent pas beaucoup impacter la complexité (2 min 20 sec pour 10 samples et 3min12 sec pour 1000  samples)
# En particulier sur 1000 samples fuzzy sur les labels on arrive à retrouver des samples corrects (8 exactement), dont certains qu'on ne trouver pas avant sur la strand E

print_mapping_on_target_graph(GP, GT, 
name_file="rin147into6v9d", output_format="png")
