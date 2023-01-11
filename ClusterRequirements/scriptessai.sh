#!/bin/bash
#SBATCH --time=00:05:00
#SBATCH --account=def-vreinhar
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=25G

module load cmake
module load gcc/11.3.0
module load gengetopt
module load python/3.10.2
module load python-build-bundle/2022a
source ../../ENV/bin/activate
export PYTHONPATH="${PYTHONPATH}:/home/tboury/Infrared-master/build/lib.linux-x86_64-3.10/infrared"
python3 -u WorkingSpace.py
deactivate
