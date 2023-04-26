TODO : introduce a limitation on the set target quickly and explain that it nededs a filter (already done the filter)
TODO : have zip files or link to bigRNAstorage(near) datas.
TODO : intergrate the ipynb in launcher (as some scripts as useless and others required multiprocessing) and remove all the remaining DEBUG messages
TODO : add somewhere Vladimir's script to create the pdb files from connex list.
TODO : update README
TODO: make clean ExtractandsearchGeometry.py and copy it in FuzzTree folder
TODO: make clean RMalignusage.py.py, computation of function in it should be put in launcher.py and copy it in FuzzTree folder
TODO : make sure that example2 is renamed, cleaned and whatch out about te proportion to dont make seems hard to understand for nothing, cleanly explained what are all of these retrieve values and from what.
TODO ; iimport_connex_list_of_list_of_RNA_chain_pdbposition.py should be merged with example2
TODO: test the code in a new folder with none of the pattern/target folder built
# FuzzTree
# Copyright (C) 2023 THEO BOURY 

# FuzzTree Project

This project present FyzzTree, a new tool to solve the Subgraph Isomorphism Problem in its "Fuzzy" version, that allowed to take into account and search for neighborhoods of the requested pattern.

## Prerequisite and launch

### Dependencies 
The project uses the following non-standard libraries:

* infrared (see https://hal.archives-ouvertes.fr/hal-03711828/document by Hua-Ting Yao, 2022 for more details.
* matplotlib (for vizualizations).
* pickle (for graph imports). (natively installed for python above 3.7 ?)
* networkx (for graph vizualizations).
* varnaapi (for graph vizualizations). 


### Files and repositories
The project is entirely available on the current archive. 

For Target RNA graphs, the used graphs were downloaded at https://uqam-my.sharepoint.com/:u:/g/personal/reinharz_vladimir_uqam_ca/Eb9rc26hI-FPgBuRH8EHR3gByiQmcgCgMu6hAdlK0-EVyA?e=R9FocY. (1.5 GB) for Our recommandation is to extract them directly into the folder bigRNAstorage

The different files are:
* FuzzTree.py: launches the tool to search for Fuzzy Subgraph Isomorphism
* Extractor.py to extract the RNA pattern from csv to put them in pickle format files
* VarnaDrawing.py: contains necessary wrappers to launch Varna on graphs and their mapping.
* TestFuzzTree.py: contains the test framework 
* TODO : to remove for final version -- WorkingSpace.py: to write and import functions to import rin, build graphs and launch some tests 

The different folders are:
* bigRNAstorage mandatory to build it): To store the Target RNA graphs


### Get started

First, this project will require Anaconda (or at least Miniconda), make sure that conda is installed or you can have a look at https://conda.io/en/latest/miniconda.html for the installation

Make sure next that conda path is known with environement variable:
```bash
export PATHTOCONDA=Your_path_to_conda/conda
```

First init conda locally and link it to avoid to replace your python by conda :

```bash
conda config --set auto_activate_base False
ln -s your/path/to/conda/conda conda
```

Make sure to create a virtual environnement for conda called FuzzTreeEnv :

```bash
make env
```

To install the dependencies package you can next type :

```bash
make dependencies
```

TODO: not sure we need it anymore
```bash
conda activate fish
```

You can now launch the test by using: #TODO: for now the real test is to laucn the WorkingSpace.py file directly, different tests can be done by switching the parameter "test" in entry of function work at the end of the file.

```bash
make tests
```
If no mistake appears here, you are good to go.

#### Examples of executions

TODO : propose examples

## Option for Fuzzyness

Parameters can be modified to customize the searched Fuzziness, they are taken in input of the main function from Fuzztree.py
Here are the different options :

- E, specifies the threshold on the sum of distances in isostericity between labels in the searched pattern and corresponding labels in the found pattern. 
- B, specifies the threshold on the number of interactions that we allow to be missing (number of edges missing).
- A, specifies the threshold on the sum of geometrical distances covered by gaps.
- maxGapAllowed, specifies the geometric distance at which we allow to look for potential gap, also distance for which we look for missing edges.


## Contributors

* Théo Boury
* Vladimir Reinharz
* Yann Ponty






