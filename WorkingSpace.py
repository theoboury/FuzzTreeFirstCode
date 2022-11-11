import os, glob, pickle
import networkx as nx
from FuzzTree import main
from VarnaDrawing import print_mapping_on_target_graph
from TestFuzzTree import first_test_mapping, first_test_varna_with_mapping, first_test_varna_without_mapping, test_graph_where_pattern_is_detected
import time
from multiprocessing import Process, Queue
from RIN import import_rin


#name="1Y27-1B53misslabeledCWW1edgesmissing1falselabelCHSintoCWH.nxpickle"
#with open("1Y27.nxpickle",'rb') as fG:
    #GG = pickle.load(fG)
    #print([(i,j) for (i,j,t) in GG.edges.data() if i == ('X', 24) or j == ('X', 20)])
    #GG.remove_edge(('X', 53), ('X', 54))

    #GG.remove_edge(('X', 48), ('X', 52))
    #GG.remove_edge(('X', 52), ('X', 53))
    #GG.remove_edge(('X', 24), ('X', 20))
    #print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 50) and j == ('X', 49)])
    #GG.remove_edge(('X', 52), ('X', 53))
    #GG.add_edge(('X', 52), ('X', 53), label='CWW', near=False, long_range=False)
    #print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 52) and j == ('X', 53)])
    #print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 50) and j == ('X', 49)])
    #GG.remove_edge(('X', 50), ('X', 49))
    #GG.add_edge(('X', 50), ('X', 49), label='CWH', near=False, long_range=False)
    #print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 50) and j == ('X', 49)])
#with open(name, 'wb') as ff:
#    pickle.dump(GG, ff)

#name="1Y27-4falselabelsCHSintoCWHandTHWintoTHSandTHWintoTHSandcWWintoTHS.nxpickle"
#with open("1Y27.nxpickle",'rb') as fG:
#    GG = pickle.load(fG)

#    print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 50) and j == ('X', 49)])
#    GG.remove_edge(('X', 50), ('X', 49))
#    GG.add_edge(('X', 50), ('X', 49), label='CWH', near=False, long_range=False)
#    print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 50) and j == ('X', 49)])

#    print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 53) and j == ('X', 20)])
#    GG.remove_edge(('X', 53), ('X', 20))
#    GG.add_edge(('X', 53), ('X', 20), label='THS', near=False, long_range=False)
#    print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 53) and j == ('X', 20)])

#    print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 20) and j == ('X', 53)])
#    GG.remove_edge(('X', 20), ('X', 53))
#    GG.add_edge(('X', 20), ('X', 53), label='THS', near=False, long_range=False)
#    print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 20) and j == ('X', 53)])

#    print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 46) and j == ('X', 54)])
#    GG.remove_edge(('X', 46), ('X', 54))
#    GG.add_edge(('X', 46), ('X', 54), label='THS', near=False, long_range=False)
#    print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 46) and j == ('X', 54)])
#with open(name, 'wb') as ff:
#    pickle.dump(GG, ff)

#name="1Y27-1gapinside1edgesmissing4falselabelsCHSintoCWHandTHWintoTHSandTHWintoTHSandcWWintoTHS.nxpickle"
#name="1Y27-1gapinside.nxpickle"
#with open("1Y27.nxpickle",'rb') as fG:
#    GG = pickle.load(fG)
#print(GG.nodes[('X', 51)])
#from FuzzTree import distance

#GG.remove_edge(('X', 24), ('X', 25))
#GG.add_node(('X', 69), pdb_position = 'unkown', nt= 'U', real_nt = 'U', atoms = GG.nodes[('X', 25)]['atoms'])

#GG.add_edge(('X', 24),('X', 69), label='B53', long_range=False, near=False)
#GG.add_edge(('X', 69),('X', 25), label='B53', long_range=False, near=False)
#print(distance(('X', 24), ('X', 69),GG))
#print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 50) and j == ('X', 49)])
#GG.remove_edge(('X', 50), ('X', 49))
#GG.add_edge(('X', 50), ('X', 49), label='CWH', near=False, long_range=False)
#print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 50) and j == ('X', 49)])

#print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 53) and j == ('X', 20)])
#GG.remove_edge(('X', 53), ('X', 20))
#GG.add_edge(('X', 53), ('X', 20), label='THS', near=False, long_range=False)
#print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 53) and j == ('X', 20)])

#print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 20) and j == ('X', 53)])
#GG.remove_edge(('X', 20), ('X', 53))
#GG.add_edge(('X', 20), ('X', 53), label='THS', near=False, long_range=False)
#print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 20) and j == ('X', 53)])

#print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 46) and j == ('X', 54)])
#GG.remove_edge(('X', 46), ('X', 54))
#GG.add_edge(('X', 46), ('X', 54), label='THS', near=False, long_range=False)
#print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 46) and j == ('X', 54)])
#print(GG.edges.data())

#GG.remove_edge(('X', 27), ('X', 17))
    #GG.remove_edge(('X', 52), ('X', 53))
#GG.remove_edge(('X', 24), ('X', 20))
#print(GG.nodes[('X', 69)])
#with open(name, 'wb') as ff:
#    pickle.dump(GG, ff)


#name="1Y27-1B53misslabeledCWW1falselabelCHSintoCWH.nxpickle"
#with open("1Y27.nxpickle",'rb') as fG:
#    GG = pickle.load(fG)
#    print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 52) and j == ('X', 53)])
#    GG.remove_edge(('X', 52), ('X', 53))
#    GG.add_edge(('X', 52), ('X', 53), label='CWW', near=False, long_range=False)
#    print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 52) and j == ('X', 53)])
#    print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 50) and j == ('X', 49)])
#    GG.remove_edge(('X', 50), ('X', 49))
#    GG.add_edge(('X', 50), ('X', 49), label='CWH', near=False, long_range=False)
#    print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 50) and j == ('X', 49)])
#with open(name, 'wb') as ff:
#    pickle.dump(GG, ff)

#with open("1Y27.nxpickle",'rb') as fG:
#    GG = pickle.load(fG)
##with open("1Y27-4falselabelsCHSintoCWHandTHWintoTHSandTHWintoTHSandcWWintoTHS.nxpickle",'rb') as fT:
#with open("1Y27-2edgesmissing.nxpickle",'rb') as fT:    
#    GT = pickle.load(fT)
#Glist = [(i,j, t['label']) for (i,j,t) in GG.edges.data()]
#Tlist = [(i,j, t['label']) for (i,j,t) in GT.edges.data()]
#temp =[(elem1, elem2) for elem in Glist if elem != Glist[k]]
#print([i for i in Glist if i not in Tlist])
#from FuzzTree import distance
#for (i,j, t) in [i for i in Glist if i not in Tlist]:
#    print("dist", distance(i,j, GG))
#import_rin(163)

#first_test_mapping()
first_test_varna_with_mapping()
first_test_varna_without_mapping()


#test_graph_where_pattern_is_detected(GPpath = "rin163.pickle", GTlistfolder = "RNAstorage", BooleanListGPinGT = [0,0,0,0,0,0,0,0,0,0,0], nb_samples=100)


#([(1, ('X', 52)), (2, ('X', 46)), (3, ('X', 47)), (4, ('X', 48)), (5, ('X', 49)), (6, ('X', 50)), (7, ('X', 51)), (8, ('X', 20)), (9, ('X', 53)), (10, ('X', 22)), (11, ('X', 24)), (12, ('X', 25)), (13, ('X', 21)), (14, ('X', 54))], 14)
#found for rin163 into 1Y27 in 506.2300956249237 secondes

