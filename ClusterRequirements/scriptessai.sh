#!/bin/bash
#SBATCH --time=00:15:00
#SBATCH --account=def-tboury
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=2009G

module load cmake
module load gcc/11.3.0
module load gengetopt
module load python/3.10.2
module load python-build-bundle/2022a
source ENV/bin/activate
python3 WorkingSpace.py
deactivate
