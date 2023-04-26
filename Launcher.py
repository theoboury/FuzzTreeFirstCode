# FuzzTree
# Copyright (C) 2023 THEO BOURY 
    
import argparse
from FuzzTree import main
from VarnaDrawing import print_mapping_on_target_graph
from TestFuzzTree import test_mapping, test_varna
from TestFuzzTree import test_GP_into_multiples_GT
from RIN import import_rin
from Extractor import csv_parse
from TestFuzzTree import bar_graph_3proportions_1time_by_filename, bar_graph_time_by_filename
from ExtractAndSearchGeometry import full_metrics
from Cartography import fromorigincartograph, plot_cartography

parser = argparse.ArgumentParser()
parser.add_argument('--task', type=str, required=True)
parser.add_argument('--number', type=str, required=False)
parser.add_argument('--procs', type=str, required=False)
parser.add_argument('--timeout', type=str, required=False)
parser.add_argument('--near', type=str, required=False)
parser.add_argument('--pattern', type=str, required=False)
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
    csv_parse("kink_turn", -1, RNAstorage = "bigRNAstorage/", csvlocation = "RNAcsv/", pattern_place="ALLkinkturnpattern/", target_place ="ALLkinkturntarget/", withgaps = 0)
    csv_parse("kink_turn", -1, RNAstorage = "bigRNAstorage/", csvlocation = "RNAcsv/", pattern_place="ALLkinkturnpatternwithgaps/", target_place ="ALLkinkturntargetwithgaps/", withgaps = 1)
if args.task == "launch_sanity_test":
    test_mapping("ALLkinkturnpattern/53kink_turninto3NVI.pickle", "bigRNAstorage/3NVI.nxpickle", E=10, B=4, A=20, maxGAPdistance=10, nb_samples=1000, D = 5, nb_procs = 1)
    print("FuzzTree method is setted up!")
    test_varna("SanityCheck","ALLkinkturnpattern/53kink_turninto3NVI.pickle", "bigRNAstorage/3NVI.nxpickle", show=1, output_format="png", E = 10, B = 4, A = 20, maxGAPdistance=10, nb_samples=1000, D = 5, nb_procs = 1)
    print("FuzzTree method + Varna drawing are setted up! You can see SanityCheck.png for the mapping drawn with Varna.")
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
        print("Not working in a distributed system, computation done only on small RNAs")
    filelist= ["smallRNA", "4LFB", "4V9F", "4V88", "4WF9", "5J7L", "5TBW", "6CZR", "7A0S", "7RQB"] #'4CS1' skipped due to symmetry
    file = filelist[pid]
    if pid == 0:
        timeout = 36000 * 3
    else:
        timeout = 2000 #3600
    if args.timeout:
        timeout = int(args.timeout)
    perfect_mapping = csv_parse(file, -1, csvlocation="RNAcsv/byRNA/")
    resu = test_GP_into_multiples_GT(pattern, GTlistfolder = GTlistfolder, threshold_bigGT = threshold_bigGT, strong_mapping = 0.8, respect_injectivity=1, E=L , B=E, A=G, maxGAPdistance = Dgap, nb_samples=nb_samples, remove_near=rm_near, timeout= timeout, D = Dedge, nb_procs = nb_procs, perfect_mapping=perfect_mapping)
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
        #from example2 import example3 
        #list_resu = example3()
        from example2 import example4
        list_resu = example4()
    else:
        from example2 import example2
        list_resu = example2()
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
if args.task == "compute_cartography_example":
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



