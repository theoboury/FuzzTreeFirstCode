from Extractor import csv_parse
from FuzzTree import main
import os, glob, pickle
from TestFuzzTree import open_graph
from multiprocessing import Pool
csv_parse("kink_turn", -1, RNAstorage = "bigRNAstorage/", csvlocation = "RNAcsv/", pattern_place="ALLkinkturnpatternwithgaps/", target_place ="ALLkinkturntargetwithgaps/", withgaps = 1)



def wrapper(GPsrc_GTsrc_Lmax_Emax_Gmax_Dedgemax_Dgapmax):
    (GPsrc, GTsrc, Lmax, Emax, Gmax, Dedgemax, Dgapmax) = GPsrc_GTsrc_Lmax_Emax_Gmax_Dedgemax_Dgapmax
    GP, GT = GPsrc, GTsrc
    if len(GPsrc.nodes()) > len(GTsrc.nodes()):
        GP, GT = GTsrc, GPsrc
    print("lenGP", len(GPsrc.nodes()), "lenGT" , len(GTsrc.nodes()))
    mappings = main(GP, GT, Lmax, Emax, Gmax, maxGAPdistance=Dgapmax, D = Dedgemax)
    if len(mappings) == 0:
        return (-1,-1,-1)
    else:
        left = 0
        right = Emax
        while left < right:
            print("E", left, right)
            mappings = main(GP, GT, Lmax, left + int((right - left)/2), Gmax, maxGAPdistance=Dgapmax, D = Dedgemax)
            if len(mappings) == 0:
                left = left + int((right - left)/2) + 1
            else: 
                right = left + int((right - left)/2)
        print("E", left, right)
        Efinal = right
        left = 0
        right = Lmax
        while left < right:
            print("L", left, right)
            mappings = main(GP, GT, left + int((right - left)/2), Efinal, Gmax, maxGAPdistance=Dgapmax, D = Dedgemax)
            if len(mappings) == 0:
                left = left + int((right - left)/2) + 1
            else: 
                right = left + int((right - left)/2)
        print("L", left, right)
        Lfinal = right
        left = 0
        right = Gmax
        while left < right:
            print("G", left, right)
            mappings = main(GP, GT, Lfinal, Efinal, left + int((right - left)/2), maxGAPdistance=Dgapmax, D = Dedgemax)
            if len(mappings) == 0:
                left = left + int((right - left)/2) + 1
            else: 
                right = left + int((right - left)/2)
        print("G", left, right)
        Gfinal = right
        print("param", (Lfinal, Efinal, Gfinal))
        return (Lfinal, Efinal, Gfinal)

def fromorigincartograph(graphname, Lmax, Emax, Gmax, Dedgemax, Dgapmax, nb_procs=32):
    GTlist = []
    resu = []
    GPsrc = open_graph("ALLkinkturnpatternwithgaps/" + graphname + ".pickle")
    path = os.path.abspath(os.getcwd()) + "/" + "ALLkinkturntargetwithgaps"
    path_list = glob.glob(os.path.join(path, '*.nxpickle'))
    for filename in path_list:
        GTlist.append(open_graph(filename))
    entry = []
    for GTsrc in GTlist:
        entry.append((GPsrc, GTsrc, Lmax, Emax, Gmax, Dedgemax, Dgapmax))
    with Pool(nb_procs) as pool:
        resu = list(pool.imap_unordered(wrapper, entry))
    print("resu", resu)
    return resu
fromorigincartograph("20kink_turninto5TBW", 50, 12, 50, 20, 20, 32)

