from FuzzTree import main
import os, glob, pickle
from TestFuzzTree import open_graph
from multiprocessing import Pool
import matplotlib.pyplot as plt  
from mpl_toolkits.mplot3d import Axes3D
from random import random


def wrapper(GPsrc_GTsrc_Lmax_Emax_Gmax_Dedgemax_Dgapmax_filename):
    (GPsrc, GTsrc, Lmax, Emax, Gmax, Dedgemax, Dgapmax, filename) = GPsrc_GTsrc_Lmax_Emax_Gmax_Dedgemax_Dgapmax_filename
    GP, GT = GPsrc, GTsrc
    if len(GPsrc.nodes()) > len(GTsrc.nodes()):
        GP, GT = GTsrc, GPsrc
    print(GP.nodes())
    print(GP.edges())
    mappings = main(GP, GT, Lmax, Emax, Gmax, maxGAPdistance=Dgapmax, D = Dedgemax)
    if len(mappings) == 0:
        print("Not found", filename)
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
        print("param", filename, (Lfinal, Efinal, Gfinal))
        return (filename, (Lfinal, Efinal, Gfinal))

def fromorigincartograph(graphname, Lmax, Emax, Gmax, Dedgemax, Dgapmax, nb_procs=32):
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
    print("resu", resu)
    return resu



def plot_cartography(resu, title, precolor, xlabel = "", ylabel = "", zlabel = "", param_plotted = [True, True, True]):
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
        threeD = 0
    else:
        x = [triplet[ind[0]] for triplet in resu]
        y = [triplet[ind[1]] for triplet in resu]
        z = [triplet[ind[2]] for triplet in resu]
        threeD = 1

    plt.rc('xtick', labelsize=7)
    color = []
    family = []
    marker_list = []
    for (fam, col, nb, form) in precolor:
        color += [col]*nb
        family += [fam]*nb
        marker_list += [form]*nb
    print("len", len(resu), len(color))
    if threeD:
        fig = plt.figure()
        ax1 = fig.add_subplot(111, projection='3d')
        #fig, ax1 = plt.subplots(projection='3d')  
        ax1.set_zlabel(zlabel, color='black')
        ax1.set(xlim=(-0.5, max(x) + 0.5), ylim=(-0.5, max(y) + 0.5), zlim=(-0.5, max(z) + 0.5))
    else:
        fig, ax1 = plt.subplots()
        ax1.set(xlim=(-0.5, max(x) + 0.5 + 20), ylim=(-0.5, max(y) + 0.5))
    ax1.set_xlabel(xlabel, color='black')
    ax1.set_ylabel(ylabel, color='black')
    old_family = ''
    for i in range(len(color)):
        if old_family != family[i]:
            lab = family[i]
            old_family = lab
        else:
            lab = None
        if threeD:
            ax1.plot(x[i], y[i], z[i], color = color[i], markerfacecolor='None',linestyle = 'None', marker=marker_list[i], label= lab)
        else:
            ax1.plot(x[i], y[i], color = color[i], markerfacecolor='None',linestyle = 'None', marker=marker_list[i], label= lab)

    ax1.set_facecolor(color='white')
    fig.set_facecolor(color='white')
    if threeD:
        fig.legend(loc = 'lower right')
    else:
        ax1.legend()
    plt.title(title)
    plt.show()

