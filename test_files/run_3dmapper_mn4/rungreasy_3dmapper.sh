#!/bin/bash
#SBATCH -o /gpfs/scratch/bsc08/bsc08927/3Dmapper_analysis/debug/output_%A.txt
#SBATCH -e /gpfs/scratch/bsc08/bsc08927/3Dmapper_analysis/debug/errors_%A.txt
#SBATCH -J 3dmapper
#SBATCH --qos=debug
#SBATCH --ntasks=20
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=6

# Como tengo que correr 3dmapper con 6 data studies diferentes, cada uno está en carpetas diferentes para no mezclar resultados
study="$1"
#######################################################################
# List of tasks 
#######################################################################
FILE=/gpfs/scratch/bsc08/bsc08927/3Dmapper_analysis/$study/input/tasks.txt

echo $FILE
#######################################################################
# Greasy log file
#######################################################################
export GREASY_LOGFILE=/gpfs/scratch/bsc08/bsc08927/3Dmapper_analysis/debug/log_3dmapper_$study.txt

#######################################################################
# Run greasy
#######################################################################
# load necessary modules
module load greasy
module load python/3.6.1

# activate environment
. /gpfs/home/bsc08/bsc08927/PDBmapper-master/venv/bin/activate

# indicate pythonpath (esto es porque el paquete dask no lo encontraba o algo así y lo instalé manualmente, no lo recuerdo bien. Puede que tú no tengas el mismo problema).
PYTHONPATH=/gpfs/home/bsc08/bsc08927/PDBmapper-master/venv/lib:/gpfs/home/bsc08/bsc08927/PDBmapper-master/venv/lib/dask_jobqueue-0.7.2+6.ge7fcb31-py3.6.egg

#run greasy
greasy $FILE

