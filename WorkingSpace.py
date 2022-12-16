
from FuzzTree import main
from VarnaDrawing import print_mapping_on_target_graph
from TestFuzzTree import test_mapping, test_varna
from TestFuzzTree import test_GP_into_multiples_GT, test_perfect_mapping_multiprocess_oneRNA_sliced
from RIN import import_rin
from Extractor import  csv_parse
from TestFuzzTree import bar_graph_3proportions_1time_by_filename
 
def work(test = 1):
    if test == 1:
        csv_parse("IL_29549.9", [(5,6)])
        resu = test_GP_into_multiples_GT("ALLkinkturnpattern/22IL_29549.9into5J7L.pickle", GTlistfolder = "ALLkinkturntarget", threshold_bigGT = 500, strong_mapping = 1, respect_injectivity=1, E=80 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, remove_near=True, timeout=-1, D = 7, nb_procs = 4)
        bar_graph_3proportions_1time_by_filename(resu, "turboessai")
    if test == 2:
        csv_parse("IL_29549.9", [(5,6)])
        resu = test_GP_into_multiples_GT("ALLkinkturnpattern/22IL_29549.9into5J7L.pickle", GTlistfolder = "ALLkinkturntarget", threshold_bigGT = 5, strong_mapping = 1, respect_injectivity=1, E=80 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, remove_near=True, timeout=-1, D = 7, nb_procs = 4)
        bar_graph_3proportions_1time_by_filename(resu, "turboessai2")
    if test == 3:
        perfect_mapping = csv_parse("IL_29549.9", [(5,6)])
        #perfect_mapping = [perfect_mapping[i] for i in range(len(perfect_mapping)) if perfect_mapping[i][0] in ['4V88']]
        resu = test_GP_into_multiples_GT("ALLkinkturnpattern/22IL_29549.9into5J7L.pickle", GTlistfolder = "bigRNAstorage", threshold_bigGT = 500, strong_mapping = 1, respect_injectivity=1, E=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, remove_near=True, timeout= 3600, D = 5, nb_procs = 32, perfect_mapping=perfect_mapping)
        print("resu", resu)
        resu1 = resu[:10] 
        resu2 = resu[10:20] 
        resu3 = resu[20:]
        bar_graph_3proportions_1time_by_filename(resu1, "22IL_29549.9into5J7L--familyIL29549.9\nE=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, timeout=3600, D = 5")
        bar_graph_3proportions_1time_by_filename(resu2, "22IL_29549.9into5J7L--familyIL29549.9\nE=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, timeout=3600, D = 5")
        bar_graph_3proportions_1time_by_filename(resu3, "22IL_29549.9into5J7L--familyIL29549.9\nE=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, timeout=3600, D = 5")
    if test == 4:
        perfect_mapping = csv_parse("IL_29549.9", [(5,6)])
        perfect_mapping = [perfect_mapping[i] for i in range(len(perfect_mapping)) if perfect_mapping[i][0] in ['4V88']]
        resu = test_GP_into_multiples_GT("ALLkinkturnpattern/22IL_29549.9into5J7L.pickle", GTlistfolder = "bigRNAstorage", threshold_bigGT = 500, strong_mapping = 1, respect_injectivity=1, E=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, remove_near=True, timeout= -1, D = 5, nb_procs = 32, perfect_mapping=perfect_mapping)
        print("resu", resu)

work(test = 3)



