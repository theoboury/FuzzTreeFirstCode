
from FuzzTree import main
from VarnaDrawing import print_mapping_on_target_graph
from TestFuzzTree import test_mapping, test_varna
from TestFuzzTree import test_GP_into_multiples_GT, test_perfect_mapping_multiprocess_oneRNA_sliced
from RIN import import_rin
from Extractor import  csv_parse
from TestFuzzTree import bar_graph_3proportions_1time_by_filename
 
def work(test = 1):
    if test == 1:
        test_GP_into_multiples_GT("rin163.pickle", GTlistfolder = "RNAstorage", threshold_bigGT = 500, strong_mapping = 1, respect_injectivity=1, E=80 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, remove_near=True, timeout=-1, D = 7, nb_procs = 4)
    if test == 2:
        perfect_mapping = csv_parse("IL_29549.9", [(5,6)])
        print("perfect mapping", perfect_mapping)
        perfect_mapping = [perfect_mapping[i] for i in range(len(perfect_mapping)) if perfect_mapping[i][0] in ['4V88']]
        #perfect_mapping = [perfect_mapping[i] for i in range(len(perfect_mapping)) if i > 16]
        print(len(perfect_mapping))
        #test_perfect_mapping(perfect_mapping, GPpath = "ALLkinkturnpattern/7IL_29549.9into6UFG.pickle", E=0 , B=0, A=0, maxGAPdistance = 0, nb_samples=10, remove_near=True, timeout=800, D = 5)
        resu = test_perfect_mapping_multiprocess_oneRNA_sliced(perfect_mapping[0], GPpath = "ALLkinkturnpattern/22IL_29549.9into5J7L.pickle", pattern_name= "22IL_29549.9into5J7L", E=0 , B=4, A=20, maxGAPdistance = 10, nb_samples=100, remove_near=True, timeout=3600 * 24 * 5, D = 5, strong_mapping=1)
        print("resu", resu)

work(test = 1)



