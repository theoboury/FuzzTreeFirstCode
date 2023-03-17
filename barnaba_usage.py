# FuzzTree
# Copyright (C) 2023 THEO BOURY 

import seaborn as sns
import pandas as pd
import barnaba as bb
import glob


def visualise_distribution(eRMSD_kinkturn, eRMSD_external_vs_kinkturn):
    with open('eRMSD_distribution.csv', 'w') as f:
        writer = csv.writer(f)
        headers = ["eRMSD", "Family"]
        writer.writerow(headers)
        for a in eRMSD_kinkturn:
            data = [a, "eRMSD between Kink-Turns"]
            writer.writerow(data)
        for a in eRMSD_external_vs_kinkturn:
            data = [a, "eRMSD between new patterns and their closer Kink-Turns"]
            writer.writerow(data)
    distrib = pd.read_csv("eRMSD_distribution.csv")   
    sns.displot(data=distrib, x="eRMSD", hue="Family", kind="kde", fill=True)
    
# find all SARCIN motifs in H.Marismortui large ribosomal subunit (PDB 1S72)
#    query = "pdb_files/2.pdb" #pdb_files/0.pdb" 
#    pdb = "pdb_files/1.pdb"#"pdb_files/1.pdb" 
#  results = bb.ds_motif(query,pdb,l1=10,l2=11,bulges=0,threshold=0.7,out='sarcin_motif')
def compute_eRMSD(threshold, bulges):
    li_internal = []
    li_external = []
    nb_kink_turn = 72
    li_l1_l2 = [] #TODO: by hand number of nucleotides on each kink turn strand
    nb_new_pattern = 197
    for i in range(nb_kink_turn):
        for j in range(nb_kink_turn):
            if j > i:
                mini = 100.0
                query = "pdb_kink_files/" + i + ".pdb" #TODO: watch out for real name of these files and locations
                pdb = "pdb_kink_files/" + j + ".pdb"
                (l1, l2) = li_l1_l2[i]
                results = bb.ds_motif(query,pdb,l1=li1,l2=li2,bulges=bulges,threshold=threshold,out='pdb_resu/kink_internal_motif')
                pdbs = glob.glob("pdb_resu/kink_internal*.pdb")
                for k in range(len(results)):
                    mini = min(mini, results[k][1])
                li_internal.append(mini)
    for j in range(nb_new_pattern):
        mini = 100.0 #min_of_mini = 100.0
        for i in range(nb_kink_turn):
            query = "pdb_kink_files/" + i + ".pdb" #TODO: watch out for real name of these files and locations
            pdb = "pdb_files/" + j + ".pdb"
            #mini = 100.0
            (l1, l2) = li_l1_l2[i]
            results = bb.ds_motif(query,pdb,l1=li1,l2=li2,bulges=bulges,threshold=threshold,out='pdb_resu/kink_external_motif')
            pdbs = glob.glob("pdb_resu/kink_external*.pdb")
            for k in range(len(results)):
                mini = min(mini, results[k][1])
            #min_of_mini = min(min_of_mini, mini)
        li_external.append(mini) #(min_of_mini)
    return (li_internal, li_external)

def full_barnaba_treatement(threshold=0.7, bulges=0):
    (eRMSD_kinkturn, eRMSD_external_vs_kinkturn) = compute_eRMSD(threshold, bulges)
    visualise_distribution(eRMSD_kinkturn, eRMSD_external_vs_kinkturn)
    return

full_barnaba_treatement(threshold=0.7, bulges=0)