# FuzzTree Project

This project aims to present a new tool to solve the Subgraph Isomorphism Problem in its Fuzzy version.

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

For Target RNA graphs, the used graphs were downloaded at https://uqam-my.sharepoint.com/:u:/g/personal/reinharz_vladimir_uqam_ca/Eb9rc26hI-FPgBuRH8EHR3gByiQmcgCgMu6hAdlK0-EVyA?e=R9FocY. (1.5 GB) Your recommandation is to extract them directly into the folder RNAstorage

The different files are:
* FuzzTree.py: launches the tool to search for Fuzzy Subgraph Isomorphism
* FuzzynessParameters.py: contains the parameters than can be modified to customize the searched Fuzziness
* TestFuzzTree.py: contains the test framework 

The different folders are:
* RNAstorage: To store the Target RNA graphs


### Get started

First, this project will require Anaconda (or at least Miniconda), make sure that conda is installed or you can have a look at https://conda.io/en/latest/miniconda.html for the installation

Make sure next that conda is activate

```bash
conda activate  fish
```

To install the dependencies pacckage you can next type :

```bash
make dependencies
```

You can now launch the test by using:

```bash
make tests
```
If no mistake appears here, you are good to go.

#### Examples of executions


## Option for Fuzzyness

These options can be modified directy in the file FuzzynessParameters.py

TODO : introduce the different options


## Contributors

* Théo Boury
* Vladimir Reinharz
* Yann Ponty






