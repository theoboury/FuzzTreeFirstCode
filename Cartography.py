from Extractor import csv_parse
from FuzzTree import main
import os, glob, pickle
from TestFuzzTree import open_graph
from multiprocessing import Pool
import matplotlib.pyplot as plt  


#csv_parse("kink_turn", -1, RNAstorage = "bigRNAstorage/", csvlocation = "RNAcsv/", pattern_place="ALLkinkturnpatternwithgaps/", target_place ="ALLkinkturntargetwithgaps/", withgaps = 1)



def wrapper(GPsrc_GTsrc_Lmax_Emax_Gmax_Dedgemax_Dgapmax_filename):
    (GPsrc, GTsrc, Lmax, Emax, Gmax, Dedgemax, Dgapmax, filename) = GPsrc_GTsrc_Lmax_Emax_Gmax_Dedgemax_Dgapmax_filename
    GP, GT = GPsrc, GTsrc
    if len(GPsrc.nodes()) > len(GTsrc.nodes()):
        GP, GT = GTsrc, GPsrc
    print("lenGP", len(GPsrc.nodes()), "lenGT" , len(GTsrc.nodes()))
    mappings = main(GP, GT, Lmax, Emax, Gmax, maxGAPdistance=Dgapmax, D = Dedgemax)
    if len(mappings) == 0:
        print("BLUB")
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
        return (filename, (Lfinal, Efinal, Gfinal))

def fromorigincartograph(graphname, Lmax, Emax, Gmax, Dedgemax, Dgapmax, nb_procs=32):
    GTlist = []
    resu = []
    GPsrc = open_graph("ALLkinkturnpatternwithgaps/" + graphname + ".pickle")
    path = os.path.abspath(os.getcwd()) + "/" + "ALLkinkturntargetwithgaps"
    path_list = glob.glob(os.path.join(path, '*.nxpickle'))
    for filename in path_list[:3]:
        GTlist.append((filename, open_graph(filename)))
    entry = []
    for filename, GTsrc in GTlist:
        entry.append((GPsrc, GTsrc, Lmax, Emax, Gmax, Dedgemax, Dgapmax, filename))
    with Pool(nb_procs) as pool:
        resu = list(pool.imap_unordered(wrapper, entry))
    print("resu", resu)
    return resu
#fromorigincartograph("20kink_turninto5TBW", 50, 12, 50, 20, 20, 3)

from random import random

def plot(resu, title, precolor, xlabel = "", ylabel = "", zlabel = "", param_plotted = [True, True, True]):
    """
    Input: - resu, the list of triplets (L, E, G) obtained during the cartography.
           - title, a title for the plot.
           - param_plotted choice of param to plot.
    Output: returns a dot plot of the different family from an origin.
    """
    for indi1,(i1, j1, k1) in enumerate(resu):
        for indi2, (i2, j2, k2) in enumerate(resu):
            if indi2 > indi1:
                if i1 == i2 and j1 == j2 and k1 == k2:
                    resu[indi2] = (i2 + (0.5 - random())*5, j2 + (0.5 - random())*2, k2 + (1 - random())*2) #avoid number to be exactly at the same position
                    (i3, j3, k3)= resu[indi2]  
                    resu[indi2] = (max(i3, 0), max(j3, 0), max(k3, 0))
    ind = [i for i in range(len(param_plotted)) if param_plotted[i] == True]
    if len(ind) == 1:
        print("1D not supported for now")
        x = [triplet[ind[0]] for triplet in resu]
        return
    elif len(ind) == 2:
        x = [triplet[ind[0]] for triplet in resu]
        y = [triplet[ind[1]] for triplet in resu]
    else:
        print("3D not supported for now")
        x = [triplet[ind[0]] for triplet in resu]
        y = [triplet[ind[1]] for triplet in resu]
        z = [triplet[ind[2]] for triplet in resu]
        return

    plt.rc('xtick', labelsize=7)
    color = []
    family = []
    marker_list = []
    for (fam, col, nb, form) in precolor:
        color += [col]*nb
        family += [fam]*nb
        marker_list += [form]*nb
    print("len", len(resu), len(color))
    fig, ax1 = plt.subplots()
    ax1.set_xlabel(xlabel, color='black')
    ax1.set_ylabel(ylabel, color='black')
    ax1.set(xlim=(-0.5, max(x) + 0.5 + 20), ylim=(-0.5, max(y) + 0.5))
    old_family = ''
    for i in range(len(color)):
        if old_family != family[i]:
            lab = family[i]
            old_family = lab
        else:
            lab = None
        ax1.plot(x[i], y[i], color = color[i],
            markerfacecolor='None',linestyle = 'None', marker=marker_list[i], label= lab)

    ax1.set_facecolor(color='white')
    fig.set_facecolor(color='white')
    ax1.legend()
    plt.title(title)
    plt.show()

resu = [(-1, -1, -1), (-1, -1, -1), (-1, -1, -1), (0, 0, 0), (13, 0, 0), (0, 4, 0), (0, 2, 0), (0, 2, 0), 
(0, 4, 0), (0, 4, 0), (0, 6, 0), (0, 2, 0), (0, 2, 0), (0, 2, 0), (0, 2, 0), (0, 0, 0), (0, 2, 0), (-1, -1, -1), 
(-1, -1, -1), (0, 8, 0), (13, 2, 0), (0, 6, 0), (0, 2, 0), (0, 4, 0), (16, 2, 0), (0, 2, 0), (0, 4, 0), (0, 4, 0), 
(0, 2, 0), (0, 8, 0), (24, 0, 8), (0, 4, 0), (0, 2, 0), (0, 4, 0), (13, 0, 0), (0, 2, 0), (0, 4, 0), (0, 0, 0), (0, 2, 0), (0, 6, 0), (0, 4, 0), (14, 6, 2), (0, 0, 10), (0, 0, 8), (0, 0, 10), (0, 6, 2), (0, 0, 8), (0, 6, 2), (0, 0, 8), (0, 0, 8), (0, 6, 2), (0, 6, 2), (0, 6, 5), (0, 6, 2), (16, 4, 6), (0, 4, 2), (0, 4, 5), (16, 2, 9), (0, 4, 5), (0, 4, 2), (0, 4, 5), (0, 0, 5), (0, 4, 5), (0, 0, 5), (0, 2, 6), (0, 2, 5), (27, 2, 9), (25, 4, 11), (0, 2, 11), (31, 2, 20), (31, 2, 19)]
color = [('IL_29549.9', 'orange', 32), ('IL_51265.1', 'blue', 6), ('IL_90538.3', 'green', 5), ('IL_68780.2', 'pink', 3),
('IL_35904.1', 'brown', 3), ("IL_64900.1", 'cyan', 3), ('IL_68057.1', 'olive', 3), ('IL_74051.1', 'purple', 3), ("IL_53581.1", 'sienna', 2),
('IL_16456.1', 'darkblue', 2), ("IL_37722.1", 'darkkhaki', 2), ("IL_93094.1", 'grey', 1), ("IL_22436.1", 'violet', 1),
('IL_34780.1', "bisque", 1), ("IL_65876.1", "lightgreen", 1), ("IL_76900.1", "dodgerblue", 1), ('IL_95993.1', 'indigo', 1), ("IL_90922.1", "lime", 1)]

color = [('IL_29549.9', 'orange', 32, "o"), ('IL_29549.9', 'blue', 11, "x"), ('IL_68780.2', 'green', 3, "s"),
('Others families', 'blue', 25, "x")]
#'IL_68780.2' number of instance reduced by one due to symmetry

#plot([(0,1,2), (1,2,3), (3,5,6)], "Essai", [("IL32", 'red', 1), ("IL33", 'blue', 2)], param_plotted = [True, False, True])

plot(resu, "Label and gap distances cartography from IL_5TBW_059", color, xlabel="Label distance (Isostericity)", ylabel = "Gap distance (Angstrom)", param_plotted=[True, False, True])
plot(resu, "Label and edge distances cartography from IL_5TBW_059", color, xlabel="Label distance (Isostericity)", ylabel = "Edge distance (Number)", param_plotted=[True, True, False])
plot(resu, "Edge and gap distances cartography from IL_5TBW_059", color, xlabel = "Edge distance (Number)", ylabel = "Gap distance (Angstrom)", param_plotted=[False, True, True])