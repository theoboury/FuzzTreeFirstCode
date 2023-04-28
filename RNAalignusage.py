# FuzzTree
# Copyright (C) 2023 THEO BOURY 

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import glob
from TestFuzzTree import open_graph
from Extractor import csv_parse
import os
import csv

def visualise_distribution(eRMSD_kinkturn, eRMSD_external, eRMSD_external_vs_kinkturn, computation):
    """
    Input : - eRMSD_kinkturn, the list of RMSDs between the Kink-Turns
            - eRMSD_external, the list of RMSDs between the new patterns found with FuzzTree method that were not identified as Kink-Turns
            - eRMSD_external_vs_kinkturn, the list of RMSDs between the Kink-Turns and the new patterns found with FuzzTree method that were not identified as Kink-Turns
            - computation can be 
                * "closest", when only minimal RMSD between a new motif vs all the Kink-Turns are retrieved.
                * "farthest", when only maximal RMSD between a new motif vs all the Kink-Turns are retrieved.
                * "full neighborhood", when all RMSDs between a new motif vs all the Kink-Turns are retrieved.
    Output : Return a KDE plots of the RMSD for each of the RMSD lists above.
    """
    with open('eRMSD_distribution.csv', 'w') as f:
        writer = csv.writer(f)
        headers = ["RMSD", "Family"]
        writer.writerow(headers)
        for a in eRMSD_kinkturn:
            data = [a, "RMSD between Kink-Turns"]
            writer.writerow(data)
        for a in eRMSD_external:
            data = [a, "RMSD between new patterns"]
            writer.writerow(data)
        for a in eRMSD_external_vs_kinkturn:
            data = [a, "RMSD between new patterns and their " + computation + " Kink-Turns"]
            writer.writerow(data)
    distrib = pd.read_csv("eRMSD_distribution.csv")   
    plot = sns.displot(data=distrib, x="RMSD", hue="Family", fill=True, common_norm = False, kind="kde") #
    fig = plot.fig
    plt.savefig("RMSDs_Kink_Turns_vs_new_found_motifs_" + computation + '.png', format='png')
    plt.savefig("RMSDs_Kink_Turns_vs_new_found_motifs_" + computation + '.pdf', format='pdf')
def compute_eRMSD(computation):
    """
    Input : - computation can be 
                * "closest", when only minimal RMSD between a new motif vs all the Kink-Turns are retrieved.
                * "farthest", when only maximal RMSD between a new motif vs all the Kink-Turns are retrieved.
                * "full neighborhood", when all RMSDs between a new motif vs all the Kink-Turns are retrieved.
    Output : Compute RMSD for chosen computation method for all Kink-Turns vs the new motifs found with FuzzTree.
    """
    nb_kink_turn = 71
    nb_new_pattern = 198
    li_internal = []
    li_external = []
    li_internal_external = []
    for i in range(nb_kink_turn):
        for j in range(nb_kink_turn):
            if j > i:
                mini = 100.0
                query = "pdb_kink_files/" + str(i) + ".pdb" 
                pdb = "pdb_kink_files/" + str(j) + ".pdb" 
                resu_to_parse = os.system("./RNAalign -A " + query + " -B " + pdb + " > temp_kink_pdb.txt")
                with open("temp_kink_pdb.txt", "r") as f:
                    resu_to_parse = f.readlines()
                results = []
                for line in resu_to_parse:
                    if line[0:7] == "Aligned":
                        line = line.split(", ")
                        for l in line:
                            if l[0:5] == "RMSD=":
                                results.append(float(l[5:]))
                for k in range(len(results)):
                    mini = min(mini, results[k])
                li_internal.append(mini)
    max_kink_turn = max(li_internal)
    for i in range(nb_new_pattern):
        for j in range(nb_new_pattern):
            if j > i:
                mini = 100.0
                query = "pdb_files/" + str(i) + ".pdb" 
                pdb = "pdb_files/" + str(j) + ".pdb" 
                resu_to_parse = os.system("./RNAalign -A " + query + " -B " + pdb + " > temp_kink_pdb.txt")
                with open("temp_kink_pdb.txt", "r") as f:
                    resu_to_parse = f.readlines()
                results = []
                for line in resu_to_parse:
                    if line[0:7] == "Aligned":
                        line = line.split(", ")
                        for l in line:
                            if l[0:5] == "RMSD=":
                                results.append(float(l[5:]))
                for k in range(len(results)):
                    mini = min(mini, results[k])
                li_external.append(mini)
    for j in range(nb_new_pattern):
        if computation == "closest":
            mini = 100.0
        if computation == "farthest":
            mini = 0.0
        for i in range(nb_kink_turn):
            if computation == "full neighborhood":
                mini = 100.0 
            query = "pdb_kink_files/" + str(i) + ".pdb" 
            pdb = "pdb_files/" + str(j) + ".pdb"
            resu_to_parse = os.system("./RNAalign -A " + query + " -B " + pdb + " > temp_kink_pdb.txt")
            with open("temp_kink_pdb.txt", "r") as f:
                resu_to_parse = f.readlines()
            results = []
            for line in resu_to_parse:
                if line[0:7] == "Aligned":
                    line = line.split(", ")
                    for l in line:
                        if l[0:5] == "RMSD=":
                            results.append(float(l[5:]))
            for k in range(len(results)):
                if computation == "farthest":
                    mini = max(mini, results[k])
                else:
                    mini = min(mini, results[k])
            if computation == "full neighborhood":
                li_internal_external.append(mini) 
        if computation != "full neighborhood":
            li_internal_external.append(mini) 
            if computation == "farthest":
                if mini <= max_kink_turn:
                    print("New motif: " + str(j) + ", RMSD: " + str(mini) + ", below Kink-Turn RMSD threshold: " + str(max_kink_turn))
    return (li_internal, li_external, li_internal_external)

def full_treatement(computation):
    """
    Input : - computation can be 
                * "closest", when only minimal RMSD between a new motif vs all the Kink-Turns are retrieved.
                * "farthest", when only maximal RMSD between a new motif vs all the Kink-Turns are retrieved.
                * "full neighborhood", when all RMSDs between a new motif vs all the Kink-Turns are retrieved.
    Output : Compute RMSD for chosen computation method for all Kink-Turns vs the new motifs found with FuzzTree and output the corresponding KDE plot.
    """
    (eRMSD_kinkturn, eRMSD_external, eRMSD_external_vs_kinkturn) = compute_eRMSD(computation)
    visualise_distribution(eRMSD_kinkturn, eRMSD_external, eRMSD_external_vs_kinkturn, computation)
    return

def mean_RMSD():
    """
    Output :Works similarly to compute_eRMSD, but retrieve only for a given 
    new motif the mean of the RMSD between this motif and all the known Kink-Turns.
    """
    nb_kink_turn = 71
    nb_new_pattern = 198
    li_internal_external = []
    for j in range(nb_new_pattern):
        mean = 0.0
        for i in range(nb_kink_turn):
            mini = 100.0
            query = "pdb_kink_files/" + str(i) + ".pdb" 
            pdb = "pdb_files/" + str(j) + ".pdb"
            resu_to_parse = os.system("./RNAalign -A " + query + " -B " + pdb + " > temp_kink_pdb.txt")
            with open("temp_kink_pdb.txt", "r") as f:
                resu_to_parse = f.readlines()
            results = []
            for line in resu_to_parse:
                if line[0:7] == "Aligned":
                    line = line.split(", ")
                    for l in line:
                        if l[0:5] == "RMSD=":
                            results.append(float(l[5:]))
            for k in range(len(results)):
                mini = min(mini, results[k])
            mean+=mini
        li_internal_external.append(mean/nb_kink_turn)
    return li_internal_external

