#!/bin/bash
#SBATCH --array=0,1,2,3,4,5,6,7,8,9,10
#SBATCH --time=48:00:00
#SBATCH --account=def-vreinhar
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=2009G

module load cmake
module load gcc/11.3.0
module load gengetopt
module load python/3.10.2
module load python-build-bundle/2022a
source ../../ENV/bin/activate
export PYTHONPATH="${PYTHONPATH}:/home/tboury/Infrared-master/src/infrared/"
python3 -u WorkingSpace.py --number $SLURM_ARRAY_TASK_ID
deactivate
