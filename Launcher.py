# FuzzTree
# Copyright (C) 2023 THEO BOURY 

import csv
import pickle
import argparse
import matplotlib.pyplot as plt
from random import random
import seaborn as sns
import pandas as pd
from FuzzTree import main
from VarnaDrawing import print_mapping_on_target_graph
from TestFuzzTree import test_mapping, test_varna
from TestFuzzTree import test_GP_into_multiples_GT
from RIN import import_rin
from Extractor import csv_parse
from TestFuzzTree import bar_graph_time_by_filename
from ExtractAndSearchGeometry import full_metrics
from Cartography import fromorigincartograph
from RNAalignusage import full_treatement, mean_RMSD

parser = argparse.ArgumentParser()
parser.add_argument('--task', type=str, required=True)
parser.add_argument('--number', type=str, required=False)
parser.add_argument('--procs', type=str, required=False)
parser.add_argument('--timeout', type=str, required=False)
parser.add_argument('--near', type=str, required=False)
parser.add_argument('--pattern', type=str, required=False)
parser.add_argument('--target', type=str, required=False)
parser.add_argument('--thresholdbigGT', type=str, required=False)
parser.add_argument('--L', type=str, required=False)
parser.add_argument('--E', type=str, required=False)
parser.add_argument('--G', type=str, required=False)
parser.add_argument('--Dedge', type=str, required=False)
parser.add_argument('--Dgap', type=str, required=False)
parser.add_argument('--samples', type=str, required=False)
parser.add_argument('--patterncut', type=str, required=False)
args = parser.parse_args()
if args.task == "create_patterns_and_targets":
    print("Launching storage of Pattern and Target graphs!")
    csv_parse("kink_turn", -1, RNAstorage = "bigRNAstorage/", csvlocation = "RNAcsv/", pattern_place="ALLkinkturnpattern/", target_place ="ALLkinkturntarget/", withgaps = 0)
    print("Pattern and Target graphs partially extracted and stored!")
    csv_parse("kink_turn", -1, RNAstorage = "bigRNAstorage/", csvlocation = "RNAcsv/", pattern_place="ALLkinkturnpatternwithgaps/", target_place ="ALLkinkturntargetwithgaps/", withgaps = 1)
    print("Pattern and Target graphs completely extracted and stored!")
if args.task == "launch_sanity_test":
    print("Launching a simple test of the FuzzTree method and of Varna visualization!")
    test_mapping("ALLkinkturnpattern/53kink_turninto3NVI.pickle", "bigRNAstorage/3NVI.nxpickle", L=10, E=4, G=20, maxGAPdistance=10, nb_samples=1000, D = 5, nb_procs = 1)
    print("FuzzTree method is setted up!")
    test_varna("SanityCheck","ALLkinkturnpattern/53kink_turninto3NVI.pickle", "bigRNAstorage/3NVI.nxpickle", show=1, output_format="png", L = 10, E = 4, G = 20, maxGAPdistance=10, nb_samples=1000, D = 5, nb_procs = 1)
    print("FuzzTree method + Varna drawing are setted up! You can see SanityCheck.png for the mapping drawn with Varna.")
if args.task == "launch_sanity_test_novarna":
    print("Launching a simple test of the FuzzTree method and of Varna visualization!")
    test_mapping("ALLkinkturnpattern/53kink_turninto3NVI.pickle", "bigRNAstorage/3NVI.nxpickle", L=10, E=4, G=20, maxGAPdistance=10, nb_samples=1000, D = 5, nb_procs = 1)
    print("FuzzTree method is setted up!")
if args.task == "launch_fuzztree":
    if args.pattern:
        pattern = "ALLkinkturnpattern/" + args.pattern + ".pickle"
    else:
        pattern = "ALLkinkturnpattern/53kink_turninto3NVI.pickle"

    rm_near = True
    GTlistfolder = "bigRNAstorage"
    if args.near:
        if int(args.near) == 0:
            rm_near = False
            GTlistfolder = "bigRNAstoragenear"
    if args.target:
        target = GTlistfolder + "/" + args.target + ".nxpickle"
    else:
        target = GTlistfolder + "/3NVI.nxpickle"
    L = 20
    E = 4
    G = 20
    Dedge = 5
    Dgap = 10
    nb_samples=1000
    nb_procs = 64
    if args.L:
        L = int(args.L)
    if args.E:
        E = int(args.E)
    if args.G:
        G = int(args.G)
    if args.samples:
        nb_samples = int(args.samples)
    if args.procs:
        nb_procs = int(args.procs)
    if args.Dedge:
        Dedge = int(args.Dedge)
    if args.Dgap:
        Dgap = int(args.Dgap)
    mappings = test_mapping(pattern, target, L=L, E=E, G=G, maxGAPdistance=Dgap, nb_samples=nb_samples, D = Dedge, nb_procs = nb_procs)
    print("Found mappings", mappings)
if args.task == "launch_varna_mapping":
    if args.pattern:
        pattern = "ALLkinkturnpattern/" + args.pattern + ".pickle"
    else:
        pattern = "ALLkinkturnpattern/53kink_turninto3NVI.pickle"

    rm_near = True
    GTlistfolder = "bigRNAstorage"
    if args.near:
        if int(args.near) == 0:
            rm_near = False
            GTlistfolder = "bigRNAstoragenear"
    if args.target:
        target = GTlistfolder + "/" + args.target + ".nxpickle"
    else:
        target = GTlistfolder + "/3NVI.nxpickle"
    L = 20
    E = 4
    G = 20
    Dedge = 5
    Dgap = 10
    nb_samples=1000
    nb_procs = 64
    if args.L:
        L = int(args.L)
    if args.E:
        E = int(args.E)
    if args.G:
        G = int(args.G)
    if args.samples:
        nb_samples = int(args.samples)
    if args.procs:
        nb_procs = int(args.procs)
    if args.Dedge:
        Dedge = int(args.Dedge)
    if args.Dgap:
        Dgap = int(args.Dgap)
    test_varna("VarnaMapping", pattern, target, show=1, output_format="png",L=L, E=E, G=G, maxGAPdistance=Dgap, nb_samples=nb_samples, D = Dedge, nb_procs = nb_procs)
    print("See VarnaMapping.png for the mapping drawn with Varna.")
if args.task == "launch_custom_test":
    if args.pattern:
        pattern = "ALLkinkturnpattern/" + args.pattern + ".pickle"
    else:
        pattern = "ALLkinkturnpattern/20kink_turninto5TBW.pickle"
    rm_near = True
    GTlistfolder = "bigRNAstorage"
    if args.near:
        if int(args.near) == 0:
            rm_near = False
            GTlistfolder = "bigRNAstoragenear"
    threshold_bigGT = 500
    if args.thresholdbigGT:
        threshold_bigGT = int(args.thresholdbigGT)
    L = 20
    E = 4
    G = 20
    Dedge = 5
    Dgap = 10
    nb_samples=1000
    nb_procs = 64
    if args.L:
        L = int(args.L)
    if args.E:
        E = int(args.E)
    if args.G:
        G = int(args.G)
    if args.samples:
        nb_samples = int(args.samples)
    if args.procs:
        nb_procs = int(args.procs)
    if args.Dedge:
        Dedge = int(args.Dedge)
    if args.Dgap:
        Dgap = int(args.Dgap)
    if args.number:
        pid = int(args.number)
    else:
        pid = 0
        print("Not working on a distributed system, computation done only on small RNAs")
    filelist= ["smallRNA", "4LFB", "4V9F", "4V88", "4WF9", "5J7L", "5TBW", "6CZR", "7A0S", "7RQB"] #'4CS1' skipped due to symmetry
    file = filelist[pid]
    if pid == 0:
        timeout = 36000 * 3
    else:
        timeout = 2000 #3600
    if args.timeout:
        timeout = int(args.timeout)
    perfect_mapping = csv_parse(file, -1, csvlocation="RNAcsv/byRNA/")
    resu = test_GP_into_multiples_GT(pattern, GTlistfolder = GTlistfolder, threshold_bigGT = threshold_bigGT, strong_mapping = 0.8, respect_injectivity=1, L=L , E=E, G=G, maxGAPdistance = Dgap, nb_samples=nb_samples, remove_near=rm_near, timeout= timeout, D = Dedge, nb_procs = nb_procs, perfect_mapping=perfect_mapping)
    print("\nresu", resu)
if args.task == "compute_metrics_example":
    #Compute metrics, we recommand 32 procs to do that 
    near = True
    GTlistfolder = "bigRNAstorage"
    if args.near:
        if int(args.near) == 0:
            near = False
            GTlistfolder = "bigRNAstoragenear"
    if near:
        #We import the results that we have for the near computation after the FuzzTree method task "launch_custom_test"
        #Computation was launched on pattern 20kink_turninto5TBW on a distributed system for gain of time.
        from results_storage import resuwithnear
        list_resu = resuwithnear()
    else:
        #We import the results that we have after the FuzzTree method task "launch_custom_test"
        #Computation was launched on pattern 20kink_turninto5TBW on a distributed system for gain of time.
        from results_storage import resuwithoutnear
        list_resu = resuwithoutnear()
    nb_procs = 32
    if args.procs:
        nb_procs = int(args.procs)
    if nb_procs < 32:
        print("WARNING: We recommand to launch this test with at least 32 procs. \n")
    pattern_cut = 5
    if args.patterncut:
        pattern_cut = int(args.patterncut)
    dicto = {}
    for (name, blub1, blub2, mappings) in list_resu:
        if name in dicto.keys():
            dicto[name] +=mappings
        else:
            dicto[name] = mappings
    full_metrics(dicto, GTlistfolder = GTlistfolder, nb_procs = nb_procs, cutting_edge = [pattern_cut])
if args.task == "compute_cartography_kink_turns":
    #Cartography computation for the Kink-Turn family, we recommand 32 procs to do that.
    pattern = "20kink_turninto5TBW"
    if args.pattern:
        pattern = args.pattern
    nb_procs = 32
    if args.procs:
        nb_procs = int(args.procs)
    if nb_procs < 32:
        print("WARNING : We recommand to launch this test with at least 32 procs. \n")
    fromorigincartograph(pattern, Lmax = 50, Emax = 12, Gmax = 50, Dedgemax = 20, Dgapmax = 20, nb_procs = nb_procs)
if args.task == "draw_cartography_example":
    li = [(48, (-1, -1, -1)), (49,(-1, -1, -1)), (61,(-1, -1, -1)), (15, (0, 0, 0)), (53, (13, 0, 0)), (16, (0, 4, 0)), (5, (0, 4, 0)), (17, (0, 4, 0)), (8, (0, 2, 0)), (13, (0, 2, 0)), (0, (0, 2, 0)), (14, (0, 2, 0)), (44, (0, 2, 0)), (54, (13, 2, 0)), (66, (-1, -1, -1)), (18, (0, 2, 0)), (23, (0, 0, 0)), (9, (0, 8, 0)), (69, (-1, -1, -1)), (4, (0, 2, 0)), (19, (0, 2, 0)), (68, (0, 6, 0)), (3, (0, 4, 0)), (12, (0, 6, 0)), (71, (16, 2, 0)), (26, (0, 4, 0)), (22, (0, 2, 0)), (7, (0, 4, 0)), (51, (0, 2, 0)), (62, (0, 8, 0)), (21, (0, 2, 0)), (24, (0, 4, 0)), (6, (0, 4, 0)), (55, (13, 0, 0)), (65, (24, 0, 8)), (2, (0, 2, 0)), (11, (0, 4, 0)), (20, (0, 0, 0)), (43, (0, 2, 0)), (1, (0, 6, 0)), (10, (0, 4, 0)), (70, (14, 6, 2)), (34, (0, 0, 10)), (35, (0, 0, 10)), (33, (0, 0, 8)), (32, (0, 0, 8)), (28, (0, 6, 2)), (30, (0, 6, 2)), (25, (0, 6, 2)), (37, (0, 0, 8)), (27, (0, 6, 2)), (36, (0, 0, 8)), (31, (0, 6, 2)), (41, (0, 6, 5)), (64, (16, 4, 6)), (46, (0, 4, 2)), (40, (0, 4, 5)), (29, (0, 4, 2)), (39, (0, 4, 5)), (63, (16, 2, 9)), (42, (0, 4, 5)), (57, (0, 0, 5)), (38, (0, 4, 5)), (56, (0, 0, 5)), (47, (0, 2, 6)), (58, (0, 2, 5)), (50, (28, 2, 9)), (67, (25, 4, 11)), (52, (0, 2, 11)), (60, (31, 2, 20)), (59, (31, 2, 19))]
    li.sort()
    li = [j for (i, j) in li]
    for indi1,(i1, j1, k1) in enumerate(li):
        for indi2, (i2, j2, k2) in enumerate(li):
            if indi2 > indi1:
                if i1 == i2 and j1 == j2 and k1 == k2 and i1 != -1:
                    li[indi2] = (i2 + (0.5 - random())*5, j2 + (0.5 - random())*2, k2 + (1 - random())*2) #avoid number to be exactly at the same position
                    (i3, j3, k3)= li[indi2]  
                    li[indi2] = (max(i3, 0), max(j3, 0), max(k3, 0))
    color = [('IL_29549.9', 32), ('Others families', 11), ('IL_68780.2', 3), ('Others families', 25)]
    full_color = []
    for (fam, num) in color:
        full_color+= [fam]*num
    full_li = []
    nb = 0
    for k in range(len(li)):
        (a1, a2, a3) = li[k]
        if a1 != -1:
            full_li.append((a1, a2, a3, full_color[k]))
        else:
            nb += 1
    with open('carte1.csv', 'w') as f:
        writer = csv.writer(f)
        headers = ["Label Difference", "Edge Difference", "Gap difference", "Family"]
        writer.writerow(headers)
        for (a1, a2, a3, a4) in full_li:
            data = [a1, a2, a3, a4]
            writer.writerow(data)
    carte = pd.read_csv("carte1.csv")   
    print("nb_kink_turns_not_cartografied", nb) 
    plot = sns.pairplot(data=carte, hue="Family")
    fig = plot.fig
    plt.savefig("Kink_Turns_cartography"+ '.png', format='png')
    plt.savefig("Kink_Turns_cartography"+ '.pdf', format='pdf')
if args.task == "compute_new_motifs_RMSD_example":
    print("WARNING: This computation can take time!")
    full_treatement("closest")
    full_treatement("full neighborhood")
    full_treatement("farthest")
if args.task == "compute_new_motifs_mean_RMSD_example":
    print("WARNING: This computation can take time!")
    liRMSD = mean_RMSD()
    from results_storage import return_connex_list
    liconnex = return_connex_list()
    new_li = []
    for k,li in enumerate(liconnex):
        new_li.append((liRMSD[k], li))
    print(new_li)
if args.task == "plot_metrics_example":
    if args.near:
        #We import the results that we have after the postprocessing task "compute_metrics_example"
        from results_storage import postprocessresuwithnear
        resu = postprocessresuwithnear()
    else:
        #We import the results that we have for the near computation after the postprocessing task "compute_metrics_example"
        from results_storage import postprocessresuwithoutnear 
        resu = postprocessresuwithoutnear()
    resu = resu[1:]
    spec_list = []
    name_list = []
    sens_list = []
    num= 0
    for elem in resu:
        (filename, chains, (precision, specificity, sensitivity, F, nb_found_motifs, found_motifs)) = elem
        if sensitivity >= 0.75:
            num+=1
        spec_list.append(specificity)
        sens_list.append(sensitivity)
        file = filename
        for cha in chains:
            file = file + "_" + cha
        name_list.append(file)
    fig, ax1 = plt.subplots()
    ax1.set(xlim=(-0.1, 1.1), ylim=(-0.1, 1.1))
    ax1.set_xlabel("Sensitivity", color='black')
    ax1.set_ylabel("Specificity", color='black')
    ax1.plot(sens_list, spec_list, color = 'blue', marker = 'o', markerfacecolor='None',linestyle = 'None')
    offset = 0
    for k in range(len(name_list)):
        if sens_list[k] != 1:
            addon = 0
            if sens_list[k] < 0.1 and spec_list[k] < 0.1:
                addon = offset
            elif name_list[k] == '4V88_A6': #Necessary to avoid text to be above each others
                addon = -20
            label = name_list[k]
            plt.annotate(label, # this is the text
            (sens_list[k],spec_list[k]), # these are the coordinates to position the label
            textcoords="offset points", # how to position the text
            xytext=(0, 7 + addon), # distance from text to points (x,y)
            ha='center') 
        if sens_list[k] < 0.1 and spec_list[k] < 0.1:
            offset = offset + 11
    ax1.set_facecolor(color='white')
    fig.set_facecolor(color='white')
    title = "Sensitivity and specificity of found mappings for the Kink Turn family"
    if args.near:
        title = "Sensitivity and specificity of found mappings for the Kink Turn family with near"
    plt.title(title)
    plt.savefig(title + '.png', format='png')
    plt.savefig(title + '.pdf', format='pdf')
if args.task == "time_graphs_example":
    if args.near:
        entry_raw = [('1E7K', 1.4278004169464111, [[0.0, 0.175, 0.699]]), ('6HCT', 1.5956840515136719, [[0.0, 0.203, 0.762]]), ('5G4T', 1.6446595191955566, [[0.115, 0.115, 0.62]]), ('3NVI', 1.7362689971923828, [[0.518, 0.518, 1.0]]), ('3SIU', 2.3870444297790527, [[0.0, 0.184, 0.788]]), ('4BW0', 2.5716769695281982, [[0.181, 0.181, 0.663]]), ('1T0K', 3.274966239929199, [[0.128, 0.128, 0.557]]), ('3NMU', 3.298027753829956, [[0.257, 0.257, 1.0]]), ('5FJ4', 4.195765018463135, [[0.1, 0.1, 0.587]]), ('4C4W', 4.322925806045532, [[0.139, 0.139, 0.625]]), ('6HCT', 5.505043268203735, [[0.0, 0.157, 0.356], [0.0, 0.125, 0.42]]), ('5G4U', 5.8360772132873535, [[0.081, 0.081, 0.349], [0.096, 0.096, 0.295]]), ('2OZB', 6.400434494018555, [[0.0, 0.0, 0.0]]), ('5XTM', 8.633862018585205, [[0.0, 0.0, 1.0]]), ('5XTM', 8.646195411682129, [[0.0, 0.0, 0.0]]), ('5DCV', 9.876028060913086, [[0.0, 0.0, 1.0]]), ('1U63', 11.58642053604126, [[0.0, 0.0, 0.54]]), ('3SIV', 12.602856397628784, [[0.0, 0.0, 0.0], [0.0, 0.155, 0.671]]), ('5Y7M', 13.279045343399048, [[0.0, 0.0, 0.0]]), ('2HW8', 24.179459810256958, [[0.0, 0.0, 0.0]]), ('5D8H', 24.81542706489563, [[0.0, 0.125, 0.619]]), ('3U4M', 31.55271315574646, [[0.113, 0.113, 0.607]]), ('3Q3Z', 31.67681050300598, [[0.199, 0.199, 0.641]]), ('2VPL', 46.498624324798584, [[0.0, 0.0, 0.0]]), ('5FJC', 46.54033088684082, [[0.118, 0.118, 0.609]]), ('6UFM', 52.27826380729675, [[0.084, 0.084, 0.506]]), ('4LCK', 54.75152587890625, [[0.0, 0.222, 0.786]]), ('6DVK', 58.71999716758728, [[0.0, 0.0, 0.243]]), ('4AOB', 61.546157360076904, [[0.0, 0.0, 0.365]]), ('4KQY', 84.89300441741943, [[0.157, 0.157, 0.65]]), ('3RW6', 95.82019829750061, [[0.077, 0.077, 0.453]]), ('3V7E', 128.75165915489197, [[0.0, 0.0, 0.529]]), ('6UFG', 271.3213586807251, [[0.094, 0.094, 0.515]]), ('6UFH', 360.8544452190399, [[0.251, 0.251, 0.767]]), ('4GXY', 366.75551652908325, [[0.228, 0.228, 0.746]]), ('2R8S', 533.8074650764465, [[0.0, 0.0, 0.0]]), ('1U6B', 1053.009474515915, [[0.0, 0.0, 0.0]]), ('6SY6', 1612.3541376590729, [[0, 0, 0]]), ('6SY4', 2101.450353384018, [[0, 0, 0]])]
        filestorage = "Timegraphbrut/timewithnear"
        name = "Time graph FuzzTree method with near "
    else:
        entry_raw = [('1E7K', 0.9891200065612793, [[0.0, 0.157, 0.709]]), ('5G4T', 1.1437277793884277, [[0.118, 0.118, 0.593]]), ('6HCT', 1.1279723644256592, [[0.0, 0.0, 0.492]]), ('3NVI', 1.203505277633667, [[0.48, 0.48, 1.0]]), ('3SIU', 1.604506254196167, [[0.0, 0.173, 0.751]]), ('4BW0', 1.7208876609802246, [[0.113, 0.113, 0.617]]), ('3NMU', 2.1028857231140137, [[0.3, 0.3, 1.0]]), ('1T0K', 2.616589069366455, [[0.0, 0.0, 0.422]]), ('5FJ4', 2.818898916244507, [[0.091, 0.091, 0.554]]), ('4C4W', 2.8795292377471924, [[0.0, 0.0, 0.255]]), ('2OZB', 4.376844882965088, [[0.0, 0.0, 0.0]]), ('5G4U', 4.467447519302368, [[0.0, 0.0, 0.0], [0.0, 0.0, 0.358]]), ('6HCT', 4.600372314453125, [[0.0, 0.0, 0.526], [0.0, 0.0, 0.0]]), ('5XTM', 5.748464822769165, [[0.0, 0.0, 1.0]]), ('3SIV', 8.21737289428711, [[0.0, 0.0, 0.0], [0.0, 0.165, 0.637]]), ('1U63', 9.44894528388977, [[0.0, 0.0, 0.835]]), ('5XTM', 13.136950254440308, [[0.0, 0.0, 1.0]]), ('3U4M', 20.403923988342285, [[0.06, 0.06, 0.487]]), ('5D8H', 20.6782124042511, [[0.0, 0.0, 0.279]]), ('3Q3Z', 21.08718967437744, [[0.223, 0.223, 0.682]]), ('5FJC', 30.549811601638794, [[0.103, 0.103, 0.586]]), ('4LCK', 35.4049026966095, [[0.0, 0.145, 0.728]]), ('6DVK', 38.29049301147461, [[0.0, 0.0, 0.228]]), ('6UFM', 43.587788581848145, [[0.0, 0.0, 0.556]]), ('3RW6', 62.94833707809448, [[0.084, 0.084, 0.444]]), ('4KQY', 71.06677460670471, [[0.0, 0.0, 0.416]]), ('3V7E', 84.64838099479675, [[0.0, 0.0, 0.577]]), ('6UFG', 176.6623773574829, [[0.0, 0.0, 0.444]]), ('4AOB', 182.92915296554565, [[0.311, 0.311, 0.899]]), ('6UFH', 187.40633749961853, [[0.0, 0.0, 0.308]]), ('2R8S', 344.2688467502594, [[0.0, 0.0, 0.0]]), ('1U6B', 398.7175359725952, [[0.0, 0.0, 0.0]]), ('4GXY', 447.5979459285736, [[0.0, 0.0, 0.408]]), ('6SY6', 1051.7494266033173, [[0, 0, 0]]), ('6SY4', 1430.0863449573517, [[0, 0, 0]]), ('2HW8', 3270.106005191803, [[0, 0, 0]]), ('5DCV', 3953.3911316394806, [[0, 0, 0]]), ('5Y7M', 4198.523247003555, [[0, 0, 0]]), ('2VPL', 6039.295374155045, [[0, 0, 0]])]
        filestorage = "Timegraphbrut/timewithoutnear"
        name = "Time graph FuzzTree method "
    entry = []
    for (name, time, prop) in entry_raw:
        entry.append((name, time, "smallRNA"))
    li_big_RNA = ["blub", "4LFB", "4V9F", "4V88", "4WF9", "5J7L", "5TBW", "6CZR", "7A0S", "7RQB"]
    for i in range(1, 10):
        booleen = 1
        resu = []
        val = 0
        with open(filestorage + str(i) + ".txt", "r") as f:
            for line in f.readlines():
                #print("line", line)
                if len(line) >= 3:
                    if booleen:
                        resu.append(float(line[:-2]))
                    booleen = (booleen + 1) % 2
        for num in resu:
            val+= num
        entry.append((li_big_RNA[i], val, "bigRNA"))
    small1 = [(i,j) for (i, j, k) in entry if k == "smallRNA"]
    cut = (len(small1) + 1)/2
    small2 = [elem for k, elem in enumerate(small1) if k >= cut]
    small1 = [elem for k, elem in enumerate(small1) if k < cut]
    bar_graph_time_by_filename(small1, name + "part 1 out of 3", bar_length = 0.3)
    bar_graph_time_by_filename(small2, name + "part 2 out of 3", bar_length = 0.3)
    bar_graph_time_by_filename([(i,j) for (i, j, k) in entry if k == "bigRNA"], name + "with slicing in sphere part 3 out of 3", bar_length = 0.3)
if args.task == "connexity_graphs_creation_example":
    #Return connexity list with near
    #We import the results that we have for the near computation after the postprocessing task "compute_metrics_example"
    from results_storage import postprocessresuwithnear
    perfect_mapping = csv_parse("kink_turn", -1)
    resu = postprocessresuwithnear()[1:]
    new_resu = []
    for (RNA, chains, (blub1, blub2, blub3, blub4, blub5, mappings)) in resu:
        new_resu_loc = []
        loc_perfect_mapping = [perfect_mapping[i] for i in range(len(perfect_mapping)) if perfect_mapping[i][0] in [RNA]]
        local_nodes_pdb = []
        for (blub, li) in loc_perfect_mapping:
            with open("bigRNAstorage/" + RNA + ".nxpickle",'rb') as f:
                G = pickle.load(f)
            for (a,(i,j)) in li:
                b = [(ii, jj) for ((ii,jj), tt) in G.nodes.data() if ii == i and tt['pdb_position'] == j][0]
                local_nodes_pdb.append(b)
        for mapp in mappings:
            booleen = 1
            for (a,b) in mapp:
                if b in local_nodes_pdb:
                    booleen = 0
                    break
            if booleen:
                new_resu_loc.append(mapp)
        new_resu.append((RNA, new_resu_loc.copy()))
    def purge(li, elem_li):
        new_li = li.copy()
        for elem in elem_li:
            ind = new_li.index(elem)
            new_li = new_li[:ind] + new_li[(ind + 1):]
        return new_li
    connex_by_RNA = []
    for (RNA, mapping) in new_resu:
        with open("bigRNAstorage/" + RNA + ".nxpickle",'rb') as f:
            G = pickle.load(f)  
        mapping_unfold = []
        for mapp in mapping:
            for (num,e) in mapp:
                if e not in mapping_unfold:
                    mapping_unfold.append(e)  
        connex_graphs = []
        while mapping_unfold:
            elem = mapping_unfold[0]
            mapping_unfold = purge(mapping_unfold, [elem])
            new_connex = [elem]
            old_connex = []
            while new_connex != old_connex:
                old_connex = new_connex.copy()
                for elem1 in new_connex:
                    pred = [i for i in G.predecessors(elem1)]
                    predsucc = pred + [i for i in G.successors(elem1) if i not in pred]
                    adj = [e for e in predsucc if e in mapping_unfold and e not in new_connex]
                    mapping_unfold = purge(mapping_unfold, adj)
                    new_connex += adj
            connex_graphs.append(new_connex.copy())
        connex_by_RNA.append((RNA, connex_graphs.copy()))
    print("\nconnex_by_RNA", connex_by_RNA)
if args.task == "connexity_graphs_creation_with_pdb_example":
    #same but with pdb
    #We import the results that we have for the near computation after the postprocessing task "compute_metrics_example"        
    from results_storage import postprocessresuwithnear
    perfect_mapping = csv_parse("kink_turn", -1)
    resu = postprocessresuwithnear()[1:]
    new_resu = []
    for (RNA, chains, (blub1, blub2, blub3, blub4, blub5, mappings)) in resu:
        new_resu_loc = []
        loc_perfect_mapping = [perfect_mapping[i] for i in range(len(perfect_mapping)) if perfect_mapping[i][0] in [RNA]]
        local_nodes_pdb = []
        with open("bigRNAstorage/" + RNA + ".nxpickle",'rb') as f:
            G = pickle.load(f)
        for (blub, li) in loc_perfect_mapping:
            for (a,(i,j)) in li:
                b = [(ii, jj) for ((ii,jj), tt) in G.nodes.data() if ii == i and tt['pdb_position'] == j][0]
                local_nodes_pdb.append(b)
        for mapp in mappings:
            booleen = 1
            for (a,b) in mapp:
                if b in local_nodes_pdb:
                    booleen = 0
                    break
            if booleen:
                new_resu_loc.append(mapp)
        new_resu.append((RNA, new_resu_loc.copy()))
    def purge(li, elem_li):
        new_li = li.copy()
        for elem in elem_li:
            ind = new_li.index(elem)
            new_li = new_li[:ind] + new_li[(ind + 1):]
        return new_li
    connex_by_RNA = []
    for (RNA, mapping) in new_resu:
        with open("bigRNAstorage/" + RNA + ".nxpickle",'rb') as f:
            G = pickle.load(f)  
        mapping_unfold = []
        for mapp in mapping:
            for (num,e) in mapp:
                if e not in mapping_unfold:
                    mapping_unfold.append(e)  
        connex_graphs = []
        while mapping_unfold:
            elem = mapping_unfold[0]
            mapping_unfold = purge(mapping_unfold, [elem])
            new_connex = [elem]
            old_connex = []
            while new_connex != old_connex:
                old_connex = new_connex.copy()
                for elem1 in new_connex:
                    pred = [i for i in G.predecessors(elem1)]
                    predsucc = pred + [i for i in G.successors(elem1) if i not in pred]
                    adj = [e for e in predsucc if e in mapping_unfold and e not in new_connex]
                    mapping_unfold = purge(mapping_unfold, adj)
                    new_connex += adj
            pdb_connex = []
            for (i_m, j_m) in new_connex:
                b_m = [(ii, tt['pdb_position']) for ((ii,jj), tt) in G.nodes.data() if ii == i_m and jj == j_m][0]
                pdb_connex.append(b_m)
            connex_graphs.append(pdb_connex.copy())
        connex_by_RNA.append((RNA, connex_graphs.copy()))
    print("\nconnex_by_RNA_with_pdb", connex_by_RNA)