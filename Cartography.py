# FuzzTree
# Copyright (C) 2023 THEO BOURY 

from FuzzTree import main
import os, glob, pickle
from TestFuzzTree import open_graph
from multiprocessing import Pool 


def wrapper(GPsrc_GTsrc_Lmax_Emax_Gmax_Dedgemax_Dgapmax_filename):
    """
    A wrapper to compute cartography parameters for each couple of graphs.
    """
    (GPsrc, GTsrc, Lmax, Emax, Gmax, Dedgemax, Dgapmax, filename) = GPsrc_GTsrc_Lmax_Emax_Gmax_Dedgemax_Dgapmax_filename
    GP, GT = GPsrc, GTsrc
    if len(GPsrc.nodes()) > len(GTsrc.nodes()):
        GP, GT = GTsrc, GPsrc
    mappings = main(GP, GT, Lmax, Emax, Gmax, maxGAPdistance=Dgapmax, D = Dedgemax)
    if len(mappings) == 0:
        print("Parameters not retrieved for ", filename)
        return (-1,-1,-1)
    else:
        left = 0
        right = Emax
        while left < right:
            mappings = main(GP, GT, Lmax, left + int((right - left)/2), Gmax, maxGAPdistance=Dgapmax, D = Dedgemax)
            if len(mappings) == 0:
                left = left + int((right - left)/2) + 1
            else: 
                right = left + int((right - left)/2)
        Efinal = right
        left = 0
        right = Lmax
        while left < right:
            mappings = main(GP, GT, left + int((right - left)/2), Efinal, Gmax, maxGAPdistance=Dgapmax, D = Dedgemax)
            if len(mappings) == 0:
                left = left + int((right - left)/2) + 1
            else: 
                right = left + int((right - left)/2)
        Lfinal = right
        left = 0
        right = Gmax
        while left < right:
            mappings = main(GP, GT, Lfinal, Efinal, left + int((right - left)/2), maxGAPdistance=Dgapmax, D = Dedgemax)
            if len(mappings) == 0:
                left = left + int((right - left)/2) + 1
            else: 
                right = left + int((right - left)/2)
        Gfinal = right
        print("Parameters retrieved for ", filename, (Lfinal, Efinal, Gfinal))
        return (filename, (Lfinal, Efinal, Gfinal))

def fromorigincartograph(graphname, Lmax, Emax, Gmax, Dedgemax, Dgapmax, nb_procs=32):
    """
    Input : - graphname, the name of the graph that serves as origin of the cartography.
            - Lmax, the maxium threshold in term of sum of isostericity allowed.
            - Emax, the maximum threshold in term of number of missing edges allowed.
            - Gmax, the maximum threshold in term of sum in Angstrom of gaps allowed.
            - Dedgemax the maximal distance after which we are not looking foredge anymore.
            - Dgapmax the maximal distance after which we are not looking for gap anymore.
    	    - nb_samples, maximum numbers of samples allowed to look for a pattern in GT.
            - respect_injectivity is a boolean, if set to true we filter patterns that are not injective by postprocessing
    Output : List of (filename, (Lfinal, Efinal, Gfinal)) with filename a Kink-Turn and (Lfinal, Efinal, Gfinal) the optimal parameters that are required according to the cartography.
    
    """
    GTlist = []
    resu = []
    GPsrc = open_graph("ALLkinkturnpatternwithgaps/" + graphname + ".pickle")
    path = os.path.abspath(os.getcwd()) + "/" + "ALLkinkturntargetwithgaps"
    path_list = glob.glob(os.path.join(path, '*.nxpickle'))
    for filename in path_list:
        GTlist.append((filename, open_graph(filename)))
    entry = []
    for filename, GTsrc in GTlist:
        entry.append((GPsrc, GTsrc, Lmax, Emax, Gmax, Dedgemax, Dgapmax, filename))
    with Pool(nb_procs) as pool:
        resu = list(pool.imap_unordered(wrapper, entry))
    print("Parameters for each file", resu)
    return resu


