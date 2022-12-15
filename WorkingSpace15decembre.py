import os, glob, pickle
import networkx as nx
from FuzzTree import main
from VarnaDrawing import print_mapping_on_target_graph
from TestFuzzTree import test_perfect_mapping_multiprocess_oneRNA_sliced, first_test_mapping, first_test_varna_with_mapping, first_test_varna_without_mapping, test_graph_where_pattern_is_detected, test_perfect_mapping, test_perfect_mapping_multiprocess,test_perfect_mapping_multiprocess_multiple_occurences, test_perfect_mapping_multiprocess_multiple_occurences_sliced
import time
from multiprocessing import Process, Queue
from RIN import import_rin
from Extractor import extractor, csv_parse
from TestFuzzTree import bar_graph, bar_graph2
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
    #with open("RNAstorage/1Y27-4falselabelsCHSintoCWHandTHWintoTHSandTHWintoTHSandcWWintoTHS.nxpickle",'rb') as fG:
    #    GG = pickle.load(fG)
    #print([(i, t) for (i,t) in GG.nodes.data()][0])
    #print([(i,j, t['label']) for (i,j,t) in GG.edges.data() if i == ('X', 50) and j == ('X', 49)])
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
    #name="1Y27-1hugegapinside.nxpickle"
    #with open("1Y27.nxpickle",'rb') as fG:
    #    GG = pickle.load(fG)
    #print(GG.nodes[('X', 51)])
    #from FuzzTree import distance

    #GG.remove_edge(('X', 24), ('X', 25))
    #GG.add_node(('X', 69), pdb_position = 'unkown', nt= 'U', real_nt = 'U', atoms = GG.nodes[('X', 24)]['atoms'])
    #GG.add_node(('X', 70), pdb_position = 'unkown', nt= 'U', real_nt = 'U', atoms = GG.nodes[('X', 25)]['atoms'])
    #GG.add_edge(('X', 24),('X', 69), label='B53', long_range=False, near=False)
    #GG.add_edge(('X', 69),('X', 70), label='B53', long_range=False, near=False)
    #GG.add_edge(('X', 70),('X', 25), label='B53', long_range=False, near=False)
    #print(distance(('X', 24), ('X', 69),GG))
    #print(distance(('X', 69), ('X', 70),GG))
    #print(distance(('X', 70), ('X', 25),GG))
    #print(distance(('X', 24), ('X', 70),GG))
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
    #first_test_varna_with_mapping()
    #first_test_varna_without_mapping()

    #test_graph_where_pattern_is_detected(GPpath = "rin163.pickle", GTlistfolder = "RNAstorage", E=80 , B=2, A=4, maxGAPdistance=3, nb_samples=100, remove_near=True,timeout = 200)# timeout=800)
    #test_graph_where_p attern_is_detected(GPpath = "rin163.pickle", GTlistfolder = "RNAstorage", BooleanListGPinGT = [0,0,0,0,0,0,0,0,0,0,0], nb_samples=100)


    #([(1, ('X', 52)), (2, ('X', 46)), (3, ('X', 47)), (4, ('X', 48)), (5, ('X', 49)), (6, ('X', 50)), (7, ('X', 51)), (8, ('X', 20)), (9, ('X', 53)), (10, ('X', 22)), (11, ('X', 24)), (12, ('X', 25)), (13, ('X', 21)), (14, ('X', 54))], 14)
    #found for rin163 into 1Y27 in 506.2300956249237 secondes



    #VERY FIRST TEST WITH EDGE MISSING

    #list_nodes = [('H', 11), ('H', 12), ('H', 16), ('H', 17), ('H', 18), ('H', 43), ('H', 44), ('H', 45), ('H', 47), ('H', 48), ('H', 49), ('H', 50)]
    #cutting_edges = [(('H', 18), ('H', 43))]
    #extractor("bigRNAstorage/3RW6.nxpickle", "IL_3RW6_002", list_nodes, cutting_edges)

    ##list_nodes = [('A', 242), ('A', 243), ('A', 245), ('A', 246), ('A', 247), ('A', 277), ('A', 278), ('A', 279), ('A', 281), ('A', 282), ('A', 283), ('A', 50)]
    ##cutting_edges = [(('A', 247), ('A', 277))]
    ##extractor("bigRNAstorage/4LFB.nxpickle", "IL_4LFB_011", list_nodes, cutting_edges)

    #list_nodes = [('X', 1272), ('X', 1273), ('X', 1274), ('X', 1275), ('X', 1276), ('X', 1246), ('X', 1247), ('X', 1248), ('X', 1250), ('X', 1251), ('X', 1252), ('X', 1253)]
    #cutting_edges = [(('X', 1276), ('X', 1246))]
    #extractor("bigRNAstorage/4WF9.nxpickle", "IL_4WF9_046", list_nodes, cutting_edges)

    #list_nodes = [('X', 1247), ('X', 1248), ('X', 1249), ('X', 1250), ('X', 1251), ('X', 1221), ('X', 1222), ('X', 1223), ('X', 1225), ('X', 1226), ('X', 1227), ('X', 1228)]
    #cutting_edges = [(('X', 1251), ('X', 1221))]
    #extractor("bigRNAstorage/7A0S.nxpickle", "IL_7A0S_043", list_nodes, cutting_edges)

    #list_nodes = [('C', 17), ('C', 18), ('C', 19), ('C', 20), ('C', 21), ('C', 31), ('C', 32), ('C', 33), ('C', 35), ('C', 36), ('C', 37), ('C', 38)]
    #cutting_edges = [(('C', 21), ('C', 31))]
    #extractor("bigRNAstorage/3V7E.nxpickle", "IL_3V7E_002", list_nodes, cutting_edges)

    #test_graph_where_pattern_is_detected(GPpath = "rin163.pickle", GTlistfolder = "RNAstorage", E=80 , B=2, A=4, maxGAPdistance=3, nb_samples=100, remove_near=True,timeout = 300)# timeout=800)
    #test_graph_where_pattern_is_detected(GPpath = "kinkturnpattern/IL_4WF9_046.pickle", GTlistfolder = "kinkturntarget", E=0 , B=5, A=0, maxGAPdistance=3, nb_samples=100, remove_near=True,timeout = 300, D = 5)# timeout=800)

    #Optimal parameter above E = 0, B = 4, A = 0, D = 5 car D = 3 n'est pas suffisant
def work(test = 10):
    if test == 1:
    #-----ETAPE 1------ 2 kink turns similaires a GAP près IL_4LCK_006 et IL_5XTM_007 respectivement dans les familles IL_51265.1 et IL_74051.1 

    #-----ETAPE 1.1------ Paramétrage optimal dans le cas le + proche
        list_nodes = [('F', 6), ('F', 7), ('F', 8), ('F', 10), ('F', 11), ('F', 12), ('F', 94), ('F', 95), ('F', 96), ('F', 97)]
        cutting_edges = [(('F', 12), ('F', 94))]
        extractor("bigRNAstorage/4LCK.nxpickle", "IL_4LCK_006", list_nodes, cutting_edges)
    #WARNING 20ZB is very strange for no reason same for 3SIV    ?

        list_nodes = [('D', 16), ('D', 17), ('D', 18), ('D', 20), ('D', 21), ('D', 22), ('D', 23), ('D', 30), ('D', 31), ('D', 32), ('D', 33)]
        cutting_edges = [(('D', 23), ('D', 30))]
        extractor("bigRNAstorage/5XTM.nxpickle", "IL_5XTM_007", list_nodes, cutting_edges)

        test_graph_where_pattern_is_detected(GPpath = "kinkturnpattern/IL_4LCK_006.pickle", GTlistfolder = "kinkturntarget", E=0, B=0, A=5, maxGAPdistance=7, nb_samples=100, remove_near=True,timeout = 300, D = 5)# timeout=800)


    #-----ETAPE 1.2------ Dessin avec Varna
        with open("kinkturnpattern/IL_4LCK_006.pickle",'rb') as fP:
            GP = pickle.load(fP)
        with open("kinkturntarget/IL_5XTM_007.nxpickle",'rb') as fT:
            GT = pickle.load(fT)
        mapping =   [(1, ('D', 16)),
        (2, ('D', 17)),
        (10, ('D', 33)),
        (3, ('D', 18)),
        (4, ('D', 20)),
        (5, ('D', 21)),
        (9, ('D', 32)),
        (6, ('D', 23)),
        (8, ('D', 31)),
        (7, ('D', 30))]
        print([(i,j, t['label']) for (i,j,t) in GT.edges.data() if t['label'] != 'B53'])
        print_mapping_on_target_graph(GP, GT, mapping = mapping, output_format = "png", name_file = "IL_4LCK_006intoIL_5XTM_007")
    #-----ETAPE 1.3------ Extension de recherche dans toute la famille IL_74051.1 #Pas d'extensions sur le cote de l'autre famille car grapes malformés ? Labellés Kink-turn from U4 ?

        list_nodes = [('F', 6), ('F', 7), ('F', 8), ('F', 10), ('F', 11), ('F', 12), ('F', 94), ('F', 95), ('F', 96), ('F', 97)]
        cutting_edges = [(('F', 12), ('F', 94))]
        extractor("bigRNAstorage/4LCK.nxpickle", "IL_4LCK_006", list_nodes, cutting_edges)

        list_nodes = [('D', 16), ('D', 17), ('D', 18), ('D', 20), ('D', 21), ('D', 22), ('D', 23), ('D', 30), ('D', 31), ('D', 32), ('D', 33)]
        cutting_edges = [(('D', 23), ('D', 30))]
        extractor("bigRNAstorage/5XTM.nxpickle", "IL_5XTM_007", list_nodes, cutting_edges)

        list_nodes = [('B', 16), ('B', 17), ('B', 18), ('B', 20), ('B', 21), ('B', 22), ('B', 23), ('B', 30), ('B', 31), ('B', 32), ('B', 33)]
        cutting_edges = [(('B', 23), ('B', 30))]
        extractor("bigRNAstorage/5XTM.nxpickle", "IL_5XTM_003", list_nodes, cutting_edges)

        list_nodes = [('D', 18), ('D', 19), ('D', 20), ('D', 22), ('D', 23), ('D', 24), ('D', 25), ('D', 32), ('D', 33), ('D', 34), ('D', 35)]
        cutting_edges = [(('D', 25), ('D', 32))]
        extractor("bigRNAstorage/5DCV.nxpickle", "IL_5DCV_008", list_nodes, cutting_edges)

        test_graph_where_pattern_is_detected(GPpath = "kinkturnpattern/IL_4LCK_006.pickle", GTlistfolder = "kinkturntarget", E=2, B=2, A=5, maxGAPdistance=7, nb_samples=100, remove_near=True,timeout = 300, D = 5)# timeout=800)
    #We have a tHS that becomes tSS in 5DCV, we need so fuzzy label 2 seems to be enough here, in addition one edge is missing in 5DCV, as it is in the two direction we need in parmeter 2 edge missing allowed, D = 5 for edge missing is required here too
    
    if test == 2:
    #-----ETAPE 2------ 2 kink turns similaires à 1 label près et 1 edge de différence IL_5J7L_035 et IL_4CS1_002 IL_90922.1 et IL_68780.2

    #-----ETAPE 2.1------ Paramétrage optimal dans le cas le + proche

        list_nodes = [('AA', 684), ('AA', 685), ('AA', 686), ('AA', 687), ('AA', 688), ('AA', 699), ('AA', 700), ('AA', 701), ('AA', 703), ('AA', 704),('AA', 705), ('AA', 706)]
        cutting_edges = [(('AA', 688), ('AA', 699))]
        extractor("bigRNAstorage/5J7L.nxpickle", "IL_5J7L_035", list_nodes, cutting_edges)

        list_nodes = [('H', 3), ('H', 4), ('H', 5), ('H', 6), ('H', 7), ('H', 26), ('H', 27), ('H', 28), ('H', 30), ('H', 31), ('H', 32), ('H', 33)]
        cutting_edges = [(('H', 7), ('H', 26))]
        extractor("bigRNAstorage/5FJ4.nxpickle", "IL_5FJ4_002", list_nodes, cutting_edges)

        test_graph_where_pattern_is_detected(GPpath = "kinkturnpattern/IL_5J7L_035.pickle", GTlistfolder = "kinkturntarget", E=15, B=2, A=0, maxGAPdistance=7, nb_samples=100, remove_near=True,timeout = 300, D = 3)# timeout=800)
    #WARNING SHOULD BE TUNED A BIT MORE E = 12 not OK D = 5 OK D = 3 OK too

    #-----ETAPE 2.2------ Dessin avec Varna

        with open("kinkturnpattern/IL_5J7L_035.pickle",'rb') as fP:
            GP = pickle.load(fP)
        with open("kinkturntarget/IL_5FJ4_002.nxpickle",'rb') as fT:
            GT = pickle.load(fT)
        mapping =   [(1, ('H', 3)),
    (2, ('H', 4)),
    (12, ('H', 33)),
    (3, ('H', 5)),
    (11, ('H', 32)),
    (4, ('H', 6)),
    (10, ('H', 31)),
    (5, ('H', 7)),
    (9, ('H', 30)),
    (6, ('H', 26)),
    (7, ('H', 27)),
    (8, ('H', 28))]

        print_mapping_on_target_graph(GP, GT, mapping = mapping, output_format = "png", name_file = "IL_5J7L_035intoIL_5FJ4_002")

    #-----ETAPE 2.3------ Extension de recherche dans toute la famille IL_74051.1 #Pas d'extensions sur le cote de l'autre famille car famille à 1 élément
        #TODO : end modifying requestiong nodes H must still be A at some places
        list_nodes = [('AA', 684), ('AA', 685), ('AA', 686), ('AA', 687), ('AA', 688), ('AA', 699), ('AA', 700), ('AA', 701), ('AA', 703), ('AA', 704),('AA', 705), ('AA', 706)]
        cutting_edges = [(('AA', 688), ('AA', 699))]
        extractor("bigRNAstorage/5J7L.nxpickle", "IL_5J7L_035", list_nodes, cutting_edges)

        list_nodes = [('H', 18), ('H', 19), ('H', 20), ('H', 21), ('H', 22), ('H', 68), ('H', 69), ('H', 70), ('H', 73), ('H', 74), ('H', 75), ('H', 76)]
        cutting_edges = [(('H', 22), ('H', 68))]
        extractor("bigRNAstorage/6DVK.nxpickle", "IL_6DVK_002", list_nodes, cutting_edges)

        list_nodes = [('A', 3), ('A', 4), ('A', 5), ('A', 6), ('A', 7), ('A', 16), ('A', 17), ('A', 18), ('A', 20), ('A', 21), ('A', 22), ('A', 23)]
        cutting_edges = [(('A', 7), ('A', 16))]
        extractor("bigRNAstorage/4BW0.nxpickle", "IL_4BW0_001", list_nodes, cutting_edges)

        list_nodes = [('H', 3), ('H', 4), ('H', 5), ('H', 6), ('H', 7), ('H', 26), ('H', 27), ('H', 28), ('H', 30), ('H', 31), ('H', 32), ('H', 33)]
        cutting_edges = [(('H', 7), ('H', 26))]
        extractor("bigRNAstorage/5FJ4.nxpickle", "IL_5FJ4_002", list_nodes, cutting_edges)

        list_nodes = [('A', 13), ('A', 14), ('A', 15), ('A', 16), ('A', 17), ('A', 3), ('A', 4), ('A', 5), ('A', 7), ('A', 8), ('A', 9), ('A', 10)]
        cutting_edges = [(('A', 17), ('A', 3))]
        extractor("bigRNAstorage/4CS1.nxpickle", "IL_4CS1_002", list_nodes, cutting_edges)

        test_graph_where_pattern_is_detected(GPpath = "kinkturnpattern/IL_5J7L_035.pickle", GTlistfolder = "kinkturntarget", E=15, B=4, A=0, maxGAPdistance=7, nb_samples=100, remove_near=True,timeout = 300, D = 18)# timeout=800)
        #D = 18 is necessary here as D=15 is not enough, but warning the source of the problem 4CS1 is a strange pattern with strange symmetry !!!!
    if test == 3:
        list_nodes = [('C', 28), ('C', 29), ('C', 30), ('C', 32), ('C', 33), ('C', 34), ('C', 42), ('C', 43), ('C', 44), ('C', 45)]
        cutting_edges = [(('C', 34), ('C', 42))]
        extractor("bigRNAstorage/2OZB.nxpickle", "IL_2OZB_002", list_nodes, cutting_edges)
    if test == 4:
        list_nodes = [('C', 28), ('C', 29), ('C', 30), ('C', 32), ('C', 33), ('C', 34), ('C', 42), ('C', 43), ('C', 44), ('C', 45)]
        with open("bigRNAstorage/2OZB.nxpickle",'rb') as fG:
            G = pickle.load(fG)
        print("Gpath", "bigRNAstorage/2OZB.nxpickle")
        print ("list_nodes", list_nodes, "\nedges with one in list_nodes :")
        for (i,t) in G.nodes.data():
            print("A NODE", i, t)
        for (i,j,t) in G.edges.data():
            if i in list_nodes or j in list_nodes:
                print("NOT AN EDGE", i, j, t)
            if i in list_nodes and j in list_nodes:
                print("GOT ONE EDGE",i, j, t)
    if test == 5:
        csv_parse("IL_29549.9", [(5,6)])
    if test == 6:
        csv_parse("IL_74051.1", [(7,8)])
    if test == 7:
        csv_parse("IL_51265.1", [(6,7)])
    if test == 8:
        #launch test = 7 , deplace pattern, laucn test = 6 to have target, launch test = 8
            test_graph_where_pattern_is_detected(GPpath = "kinkturnpattern/4IL_51265.1into4LCK.pickle", GTlistfolder = "kinkturntarget", E=2, B=2, A=5, maxGAPdistance=7, nb_samples=100, remove_near=True,timeout = 300, D = 5)# timeout=800)
    if test == 9:
        csv_parse("IL_68780.2", [(5,6)])
        csv_parse("IL_29549.9", [(5,6)])
        #test_graph_where_pattern_is_detected(GPpath = "ALLkinkturnpattern/0IL_68780.2into4BW0.pickle", GTlistfolder = "ALLkinkturntarget", E=20, B=6, A=10, maxGAPdistance=10, nb_samples=100, remove_near=True,timeout = 20, D = 20)# timeout=800)
        test_graph_where_pattern_is_detected(GPpath = "ALLkinkturnpattern/0IL_68780.2into4BW0.pickle", GTlistfolder = "ALLkinkturntarget", E=0, B=4, A=0, maxGAPdistance=7, nb_samples=100, remove_near=True,timeout = 20, D = 5)# timeout=800)
        #9IL_29549.9into7A0S not found if only 4 edge missing insteaf of 6
        #In addtiion, 12IL_29549.9into4V88.nxpickle; 27IL_29549.9into4LFB.nxpickle', '25IL_29549.9into4V9F and 3IL_29549.9into1T0K 28IL_29549.9into6CZR,5IL_29549.9into4AOB, 31IL_29549.9into4V88 and 1IL_29549.9into4GXY.nxpickle if only 2 edges missing
        #The right to do gap and to do fuzzy label is not important here.
        #it is because it is almost the same family but with one cWW that is double in one family and simple in the other family.
    if test == 10:
        csv_parse("IL_68780.2", [(5,6)])
        csv_parse("IL_68057.1", [(7,8)])
        #test_graph_where_pattern_is_detected(GPpath = "ALLkinkturnpattern/0IL_68780.2into4BW0.pickle", GTlistfolder = "ALLkinkturntarget", E=20, B=6, A=10, maxGAPdistance=10, nb_samples=100, remove_near=True,timeout = 20, D = 20)# timeout=800)
        test_graph_where_pattern_is_detected(GPpath = "ALLkinkturnpattern/0IL_68780.2into4BW0.pickle", GTlistfolder = "ALLkinkturntarget", E=13, B=2, A=0, maxGAPdistance=10, nb_samples=100, remove_near=True,timeout = 20, D = 5)# timeout=800)
        #When E = 0 3IL_68780.2into6DVK is not working anymore the rest is fine


        #with open("kinkturnpattern/IL_5J7L_035.pickle",'rb') as fP:
        #    GP = pickle.load(fP)
        with open("ALLkinkturntarget/2IL_68057.1into3Q3Z.nxpickle",'rb') as fT:
            GT = pickle.load(fT)
        mapping =   [(1, ('V', 54)), 
        (2, ('V', 55)), 
        (12, ('V', 23)), 
        (3, ('V', 56)), 
        (11, ('V', 22)), 
        (4, ('V', 57)), 
        (10, ('V', 21)), 
        (5, ('V', 58)), 
        (9, ('V', 20)), 
        (6, ('V', 16)), 
        (7, ('V', 17)), 
        (8, ('V', 18))]

        print_mapping_on_target_graph([], GT, mapping = mapping, output_format = "png", name_file = "IL_68780.2into4BW05intoIL_68057.1into3Q3Z")
    if test == 11:
        perfect_mapping = csv_parse("IL_29549.9", [(5,6)])
        csv_parse("IL_68780.2", [(5,6)])
        print("perfect mapping", perfect_mapping)
        perfect_mapping = [perfect_mapping[i] for i in range(len(perfect_mapping)) if perfect_mapping[i][0] in ['5Y7M', '7RQB', '4V88', '5J7L', '6CZR', '5J7L', '4V88']]
        print(len(perfect_mapping))
        #test_perfect_mapping(perfect_mapping, GPpath = "ALLkinkturnpattern/7IL_29549.9into6UFG.pickle", E=0 , B=0, A=0, maxGAPdistance = 0, nb_samples=10, remove_near=True, timeout=800, D = 5)
        test_perfect_mapping(perfect_mapping, GPpath = "ALLkinkturnpattern/0IL_68780.2into4BW0.pickle", E=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, remove_near=True, timeout=3600, D = 5)
            #TODO : do the same test with reduced gap allowed maxGap = 7 and A = 10 for instance and see if we mis some specific instances. 
    #LE MAPPING EN 1T0K implique une chaine D !!!!!
    #REMARQUES :

    #DANS IL_74051.1 : on a un 5-9  qui devient un 5 - 11 sorte de gap mais pas vraiment ? Il faudrait pouvoir replier deux nucléotides sr un même nucléotide ?

#nb_nodes_GT_before 102 nb_edges_GT_before 172
#nb_nodes_GT_after 51 nb_edges_GT_after 85
#mapping_first ([(1, ('D', 29)), (2, ('D', 30)), (12, ('D', 22)), (3, ('D', 31)), (11, ('D', 21)), (4, ('D', 32)), (10, ('D', 20)), (5, ('D', 34)), (9, ('D', 17)), (6, ('D', 14)), (7, ('D', 15)), (8, ('D', 16))], 12)
#filename, proportion, time ('5Y7M.nxpickle', 0.0, 12.564391613006592)
    #Il n'y a aucun kink turn où un label change au sein de la famille
    if test == 12:
        #IN REALITY WE CAN DIRECTLY OBTAINED RESU AS THE RESULTS FROM test_perfect_mapping 
        resu = [('6UFM.nxpickle', 0.88, 527.328111410141), ('4GXY.nxpickle', 0.95, 1999.8474061489105),('7RQB.nxpickle', -1, 0.08730363845825195), ('1T0K.nxpickle', 0.58, 13.746737957000732), ('4V9F.nxpickle', -1, 3600.0996148586273), ('4AOB.nxpickle', 0.8, 312.53955245018005), ('6UFH.nxpickle', 0.95, 2243.756167411804)
        , ('6UFG.nxpickle', 0.9, 2067.5632150173187), ('5Y7M.nxpickle', 0.0, 37.807730197906494), ('7A0S.nxpickle', -1, 3600.1246066093445), ('5G4U.nxpickle', 0.33, 24.88226819038391), ('4KQY.nxpickle', 0.76, 797.4127531051636), ('4V88.nxpickle', -1, 0.05709958076477051)
        , ('5G4T.nxpickle', 0.16, 4.063181400299072), ('5G4U.nxpickle', 0.32, 36.06558680534363),('5FJC.nxpickle', 0.08, 100.25284171104431), ('4C4W.nxpickle', 0.64, 31.49244499206543), ('4WF9.nxpickle', -1, 3600.0999608039856), ('3V7E.nxpickle', 0.89, 847.1687686443329), ('5TBW.nxpickle', -1, 3600.1048793792725) ,
        ('7RQB.nxpickle', -1, 0.09167170524597168), ('5J7L.nxpickle', -1, 0.08423399925231934), ('4V9F.nxpickle', -1, 3600.122519016266), ('7A0S.nxpickle', -1, 3600.1033158302307),  ('4V9F.nxpickle', -1, 3600.1037492752075), ('4WF9.nxpickle', -1, 3600.112297773361), ('4LFB.nxpickle', -1, 3600.101831674576), ('6CZR.nxpickle', -1, 0.06897544860839844)
        , ('3RW6.nxpickle', 0.35, 77.90302109718323), ('5J7L.nxpickle', -1, 0.06272149085998535), ('4V88.nxpickle', -1, 0.061377763748168945)]
        print(len(resu))
        #fatal_list = ['5Y7M', '7RQB', '4V88', '5J7L', '6CZR', '5J7L', '4V88']
        resu1 = resu[:10] #[(i,j,k) for (i, j, k) in resu[:14] if i[:4] not in fatal_list]
        resu2 = resu[10:20] #[(i,j,k) for (i, j, k) in resu[14:] if i[:4] not in fatal_list]
        resu3 = resu[20:]
        bar_graph(resu1, "0IL_68780.2into4BW0--familyIL29549.9\nE=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, timeout=3600, D = 5")
        bar_graph(resu2, "0IL_68780.2into4BW0--familyIL29549.9\nE=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, timeout=3600, D = 5")
        bar_graph(resu3, "0IL_68780.2into4BW0--familyIL29549.9\nE=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, timeout=3600, D = 5")
    if test == 13:
        csv_parse("IL_68780.2", [(5,6)])
        perfect_mapping = csv_parse("IL_68057.1", [(7,8)])
        print(len(perfect_mapping))
        #resu = test_perfect_mapping(perfect_mapping, GPpath = "ALLkinkturnpattern/0IL_68780.2into4BW0.pickle", E=13, B=2, A=20, maxGAPdistance=10, nb_samples=100, remove_near=True,timeout = 50, D = 5)
        #ABOVE WHAT IS MANDATORY IN RESU(not known about gaps)
        motifs_mapping=[(1,8), (2,9), (3, 10), (4, 11), (5, 12), (6, 1), (7, 2), (8, 3), (9, 4), (10, 5), (11, 6), (12, 7)]
        resu, mapping = test_perfect_mapping_multiprocess(perfect_mapping, GPpath = "ALLkinkturnpattern/0IL_68780.2into4BW0.pickle", E=13, B=0, A=10, maxGAPdistance=7, nb_samples=100, remove_near=True,timeout = 50, D = 5, motifs_mapping=motifs_mapping)
        #with 10000 samples does not converge
        print(resu)
        bar_graph(resu, "0IL_68780.2into4BW0--familyIL_68057.1\nE=13, B=2, A=20, maxGAPdistance=10, nb_samples=100, timeout = 3600, D = 5")
    if test == 14:
        perfect_mapping = csv_parse("IL_29549.9", [(5,6)])
        print("perfect mapping", perfect_mapping)
        #perfect_mapping = [perfect_mapping[i] for i in range(len(perfect_mapping)) if perfect_mapping[i][0] in ['5Y7M', '7RQB', '4V88', '5J7L', '6CZR', '5J7L', '4V88']]
        #perfect_mapping = [perfect_mapping[i] for i in range(len(perfect_mapping)) if i > 16]
        print(len(perfect_mapping))
        #test_perfect_mapping(perfect_mapping, GPpath = "ALLkinkturnpattern/7IL_29549.9into6UFG.pickle", E=0 , B=0, A=0, maxGAPdistance = 0, nb_samples=10, remove_near=True, timeout=800, D = 5)
        resu, mapping = test_perfect_mapping_multiprocess(perfect_mapping, GPpath = "ALLkinkturnpattern/22IL_29549.9into5J7L.pickle", E=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, remove_near=True, timeout=3600 * 4, D = 5)
        print("resu", resu)
        print("mapping", mapping)
    if test == 15:
        resu = [('5G4T,B,A', 0.13, 5.378903388977051), ('1T0K,C,D', 0.67, 19.369725704193115), ('4C4W,H', 0.63, 44.13654041290283), ('5G4U,J,I', 0.26, 58.97392225265503), ('5Y7M,D', 0.0, 63.21999931335449), ('5G4U,I,J', 0.31, 63.92445158958435), ('3U4M,B', 0.08, 112.59432363510132), ('5FJC,A', 0.08, 166.35041213035583), ('3RW6,H', 0.33, 162.470632314682), ('4AOB,A', 0.79, 488.4759957790375), ('6UFM,B', 0.85, 693.7601215839386), ('4KQY,A', 0.78, 1088.923911333084), ('3V7E,C', 0.94, 1262.6511182785034), ('6UFG,A', 0.93, 2970.834854364395), ('6UFH,A', 0.95, 3161.3698856830597), ('4GXY,A', 0.93, 3207.8569464683533), ('4V88,A6', -1, 3600.1291987895966), ('7RQB,1A', -1, 3600.125324487686), ('4V9F,0', -1, 3600.143296480179), ('7A0S,X', -1, 3600.148289680481), ('4WF9,X', -1, 3600.145581483841), ('5TBW,1', -1, 3600.1269764900208), ('7RQB,1A', -1, 3600.097424507141), ('5J7L,DA', -1, 3600.1047308444977), ('4V9F,0', -1, 3600.146371603012), ('4LFB,A', -1, 3600.155631542206), ('6CZR,1a', -1, 3600.2748243808746), ('7A0S,X', -1, 3600.1278727054596), ('5J7L,AA', -1, 3600.103734254837), ('4V88,A6', -1, 3600.169770717621), ('4V9F,0', -1, 3600.103699207306), ('4WF9,X', -1, 3600.157623052597)]
        resu1 = resu[:10] #[(i,j,k) for (i, j, k) in resu[:14] if i[:4] not in fatal_list]
        resu2 = resu[10:20] #[(i,j,k) for (i, j, k) in resu[14:] if i[:4] not in fatal_list]
        resu3 = resu[20:]
        bar_graph(resu1, "22IL_29549.9into5J7L--familyIL29549.9\nE=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, timeout=3600, D = 5")
        bar_graph(resu2, "22IL_29549.9into5J7L--familyIL29549.9\nE=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, timeout=3600, D = 5")
        bar_graph(resu3, "22IL_29549.9into5J7L--familyIL29549.9\nE=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, timeout=3600, D = 5")
    if test == 16:
        perfect_mapping = csv_parse("IL_29549.9", [(5,6)])
        print("perfect mapping", perfect_mapping)
        #perfect_mapping = [perfect_mapping[i] for i in range(len(perfect_mapping)) if perfect_mapping[i][0] in ['5Y7M', '7RQB', '4V88', '5J7L', '6CZR', '5J7L', '4V88']]
        #perfect_mapping = [perfect_mapping[i] for i in range(len(perfect_mapping)) if i > 16]
        print(len(perfect_mapping))
        #test_perfect_mapping(perfect_mapping, GPpath = "ALLkinkturnpattern/7IL_29549.9into6UFG.pickle", E=0 , B=0, A=0, maxGAPdistance = 0, nb_samples=10, remove_near=True, timeout=800, D = 5)
        resu, mapping = test_perfect_mapping_multiprocess_multiple_occurences(perfect_mapping, GPpath = "ALLkinkturnpattern/22IL_29549.9into5J7L.pickle", E=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, remove_near=True, timeout=3600 * 24 * 5, D = 5, strong_mapping=1)
        print("resu", resu)
        print("mapping", mapping)
    if test == 17:
        #TOY EXAMPLE HERE
        resu = [('5G4T,A,B', [0.09, 0.3], 10.983399629592896), ('6UFM,B', [0.4, 0.07, 0.2], 20.046091318130493), ('4GXY,A', [], 20.076128482818604), ('1T0K,C,D', [], 20000.205880165100098), ('6UFH,A', [], 20.07821297645569), ('6UFG,A', [], 20.058774709701538), ('4AOB,A', [], 20.106066942214966), ('5Y7M,D', [], 20.17974543571472), ('4KQY,A', [], 20.068326950073242), ('5G4U,I,J', [], 20.30578112602234), ('5FJC,A', [], 20.082794666290283), ('4C4W,H', [], 20.10532236099243), ('4V88,A6', [], 20.079351663589478), ('7RQB,1A', [], 20.062228202819824), ('3U4M,B', [], 20.09431767463684), ('3V7E,C', [], 20.080447673797607), ('4V9F,0', [], 20.087060928344727), ('7A0S,X', [], 20.057557106018066), ('4LFB,A', [], 20.072710275650024), ('3RW6,H', [], 20.107783794403076), ('4WF9,X', [], 20.037036895751953), ('6CZR,1a', [], 20.103374481201172), ('5J7L,AA', [], 20.046199798583984), ('5TBW,1', [], 20.068098783493042), ('5J7L,DA', [], 20.076056480407715)]
        resu1 = resu[:10] #[(i,j,k) for (i, j, k) in resu[:14] if i[:4] not in fatal_list]
        resu2 = resu[10:20] #[(i,j,k) for (i, j, k) in resu[14:] if i[:4] not in fatal_list]
        resu3 = resu[20:]
        bar_graph2(resu1, "22IL_29549.9into5J7L--familyIL29549.9\nE=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, timeout=20, D = 5")
        bar_graph2(resu2, "22IL_29549.9into5J7L--familyIL29549.9\nE=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, timeout=20, D = 5")
        bar_graph2(resu3, "22IL_29549.9into5J7L--familyIL29549.9\nE=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, timeout=20, D = 5")
    if test == 18:
        resu = [('5G4T,B,A', 0.06, 3.864753246307373), ('1T0K,C,D', 0.29, 25.263774394989014), ('4C4W,H', 0.32, 50.28081154823303), ('5Y7M,D', 0.0, 62.246861934661865), ('5G4U,I,J', 0.22, 68.29138660430908), ('5G4U,J,I', 0.18, 70.35977911949158), ('3U4M,B', 0.01, 102.1136040687561), ('5FJC,A', 0.06, 155.3191819190979), ('3RW6,H', 0.07, 179.20753979682922), ('4AOB,A', 0.52, 467.89054322242737), ('6UFM,B', 0.61, 779.5730571746826), ('4KQY,A', 0.42, 1260.3879652023315), ('3V7E,C', 0.66, 1510.5184276103973), ('6UFG,A', 0.61, 3424.0232150554657), ('6UFH,A', 0.58, 3621.7064003944397), ('4GXY,A', 0.93, 3694.1355867385864), ('4V88,A6', -1, 14400.089645385742), ('7RQB,1A', -1, 14400.283421993256), ('4V9F,0', -1, 14400.472038507462), ('7A0S,X', -1, 14400.492266893387), ('4WF9,X', -1, 14400.558336019516), ('5TBW,1', -1, 14400.457538366318), ('7RQB,1A', -1, 14400.18429017067), ('5J7L,DA', -1, 14400.297923326492), ('4V9F,0', -1, 14400.447693824768), ('4LFB,A', -1, 14400.129451036453), ('6CZR,1a', -1, 14400.145420312881), ('7A0S,X', -1, 14400.486805438995), ('5J7L,AA', -1, 14400.109685659409), ('4V88,A6', -1, 14400.134433746338), ('4V9F,0', -1, 14400.506151914597), ('4WF9,X', -1, 14400.446256637573)]
        resu1 = resu[:10] #[(i,j,k) for (i, j, k) in resu[:14] if i[:4] not in fatal_list]
        resu2 = resu[10:20] #[(i,j,k) for (i, j, k) in resu[14:] if i[:4] not in fatal_list]
        resu3 = resu[20:]
        bar_graph(resu1, "22IL_29549.9into5J7L--familyIL29549.9\nE=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, timeout=3600, D = 5")
        bar_graph(resu2, "22IL_29549.9into5J7L--familyIL29549.9\nE=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, timeout=3600, D = 5")
        bar_graph(resu3, "22IL_29549.9into5J7L--familyIL29549.9\nE=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, timeout=3600, D = 5")
    if test == 19:
        #TOY EXAMPLE HERE
        resu = [('5G4T,A,B', [0.04], 4.4998390674591064), ('1T0K,C,D', [0.39], 19.38454580307007), ('4C4W,H', [0.32], 40.12394309043884), ('5G4U,I,J', [0.2], 47.487282514572144), ('5Y7M,D', [0.0], 48.74524903297424), ('3U4M,B', [0.0], 82.06036257743835), ('5FJC,A', [0.01], 131.43501663208008), ('3RW6,H', [0.1], 135.13832187652588), ('4AOB,A', [0.44], 442.8611431121826), ('6UFM,B', [0.6], 722.5020854473114), ('4KQY,A', [0.5], 1224.7470755577087), ('3V7E,C', [0.6], 1480.1191973686218), ('6UFG,A', [0.57], 3028.474091053009), ('6UFH,A', [0.56], 3567.6298830509186), ('4GXY,A', [0.92], 3591.281092405319), ('4V88,A6', [], 43200.13017988205), ('7RQB,1A', [], 43200.231367349625), ('4V9F,0', [], 43200.12150168419), ('7A0S,X', [], 43200.14021587372), ('4WF9,X', [], 43200.1250333786), ('4LFB,A', [], 43200.274493932724), ('6CZR,1a', [], 43200.11847472191), ('5TBW,1', [], 43200.13448238373), ('5J7L,DA', [], 43200.08169221878), ('5J7L,AA', [], 43200.28419327736)]        
        resu1 = resu[:10] #[(i,j,k) for (i, j, k) in resu[:14] if i[:4] not in fatal_list]
        resu2 = resu[10:20] #[(i,j,k) for (i, j, k) in resu[14:] if i[:4] not in fatal_list]
        resu3 = resu[20:]
        bar_graph2(resu1, "22IL_29549.9into5J7L--familyIL29549.9\nE=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, timeout=43200, D = 5, fuzzy = parameterover2")
        bar_graph2(resu2, "22IL_29549.9into5J7L--familyIL29549.9\nE=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, timeout=43200, D = 5, fuzzy = parameterover2")
        bar_graph2(resu3, "22IL_29549.9into5J7L--familyIL29549.9\nE=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, timeout=43200, D = 5, fuzzy = parameterover2")
    if test == 20:
        perfect_mapping = csv_parse("IL_29549.9", [(5,6)])
        print("perfect mapping", perfect_mapping)
        #perfect_mapping = [perfect_mapping[i] for i in range(len(perfect_mapping)) if perfect_mapping[i][0] in ['5Y7M', '7RQB', '4V88', '5J7L', '6CZR', '5J7L', '4V88']]
        #perfect_mapping = [perfect_mapping[i] for i in range(len(perfect_mapping)) if i > 16]
        print(len(perfect_mapping))
        #test_perfect_mapping(perfect_mapping, GPpath = "ALLkinkturnpattern/7IL_29549.9into6UFG.pickle", E=0 , B=0, A=0, maxGAPdistance = 0, nb_samples=10, remove_near=True, timeout=800, D = 5)
        resu, mapping = test_perfect_mapping_multiprocess_multiple_occurences_sliced(perfect_mapping, GPpath = "ALLkinkturnpattern/22IL_29549.9into5J7L.pickle", pattern_name= "22IL_29549.9into5J7L", E=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, remove_near=True, timeout=3600 * 24 * 5, D = 5, strong_mapping=1)
        print("resu", resu)
        print("mapping", mapping)

    if test == 21:
        perfect_mapping = csv_parse("IL_29549.9", [(5,6)])
        print("perfect mapping", perfect_mapping)
        perfect_mapping = [perfect_mapping[i] for i in range(len(perfect_mapping)) if perfect_mapping[i][0] in ['5Y7M']]
        #perfect_mapping = [perfect_mapping[i] for i in range(len(perfect_mapping)) if i > 16]
        print(len(perfect_mapping))
        #test_perfect_mapping(perfect_mapping, GPpath = "ALLkinkturnpattern/7IL_29549.9into6UFG.pickle", E=0 , B=0, A=0, maxGAPdistance = 0, nb_samples=10, remove_near=True, timeout=800, D = 5)
        resu = test_perfect_mapping_multiprocess_oneRNA_sliced(perfect_mapping[0], GPpath = "ALLkinkturnpattern/22IL_29549.9into5J7L.pickle", pattern_name= "22IL_29549.9into5J7L", E=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, remove_near=True, timeout=3600 * 24 * 5, D = 5, strong_mapping=1)
        print("resu", resu)
    if test == 22:
        perfect_mapping = csv_parse("IL_29549.9", [(5,6)])
        print("perfect mapping", perfect_mapping)
        perfect_mapping = [perfect_mapping[i] for i in range(len(perfect_mapping)) if perfect_mapping[i][0] in ['4V88']]
        #perfect_mapping = [perfect_mapping[i] for i in range(len(perfect_mapping)) if i > 16]
        print(len(perfect_mapping))
        #test_perfect_mapping(perfect_mapping, GPpath = "ALLkinkturnpattern/7IL_29549.9into6UFG.pickle", E=0 , B=0, A=0, maxGAPdistance = 0, nb_samples=10, remove_near=True, timeout=800, D = 5)
        resu = test_perfect_mapping_multiprocess_oneRNA_sliced(perfect_mapping[0], GPpath = "ALLkinkturnpattern/22IL_29549.9into5J7L.pickle", pattern_name= "22IL_29549.9into5J7L", E=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, remove_near=True, timeout=3600 * 24 * 5, D = 5, strong_mapping=1)
        print("resu", resu)
#list of studied RNA files ['/home/uqamportable/Documents/FuzzTreeFirstCode/bigRNAstorage/3NVI.nxpickle', '/home/uqamportable/Documents/FuzzTreeFirstCode/bigRNAstorage/3NMU.nxpickle', '/home/uqamportable/Documents/FuzzTreeFirstCode/bigRNAstorage/3Q3Z.nxpickle']
#size of near removal 0
#chains ['F']
#nb_nodes_GT_before 48 nb_edges_GT_before 78
#nb_nodes_GT_after 24 nb_edges_GT_after 39
#mapping_first ([(1, ('F', 3)), (2, ('F', 4)), (12, ('F', 21)), (3, ('F', 5)), (11, ('F', 20)), (4, ('F', 6)), (10, ('F', 19)), (5, ('F', 7)), (9, ('F', 18)), (6, ('F', 14)), (7, ('F', 15)), (8, ('F', 16))], 12)
#filename, proportion, time ('3NVI.nxpickle', 0.0, 1.468977928161621)
#size of near removal 0
#chains ['E']
#nb_nodes_GT_before 94 nb_edges_GT_before 134
#nb_nodes_GT_after 34 nb_edges_GT_after 47
#"mapping_first ([(1, ('E', 13)), (2, ('E', 14)), (12, ('E', 31)), (3, ('E', 15)), (11, ('E', 30)), (4, ('E', 16)), (10, ('E', 29)), (5, ('E', 17)), (9, ('E', 28)), (6, ('E', 24)), (7, ('E', 25)), (8, ('E', 26))], 12)
#filename, proportion, time ('3NMU.nxpickle', 0.0, 2.7047994136810303)
#size of near removal 0
#chains ['V']
#nb_nodes_GT_before 150 nb_edges_GT_before 280
#nb_nodes_GT_after 75 nb_edges_GT_after 139
#mapping_first ([(1, ('V', 47)), (2, ('V', 49)), (12, ('V', 18)), (3, ('V', 50)), (11, ('V', 16)), (4, ('V', 51)), (10, ('V', 15)), (5, ('V', 52)), (9, ('V', 14)), (6, ('V', 10)), (7, ('V', 11)), (8, ('V', 12))], 12)
#filename, proportion, time ('3Q3Z.nxpickle', 0.0, 16.299317836761475)
#[('3NVI.nxpickle', 0.0, 1.468977928161621), ('3NMU.nxpickle', 0.0, 2.7047994136810303), ('3Q3Z.nxpickle', 0.0, 16.299317836761475)]
#work(test = 19)
work(test = 21)

