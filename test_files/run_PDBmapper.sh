#!/bin/bash
#SBATCH -o /gpfs/home/bsc51/bsc51927/pdbmapper/debug/output_%A.txt
#SBATCH -e /gpfs/home/bsc51/bsc51927/pdbmapper/debug/errors_%A.txt
#SBATCH -J pdbmapper
#SBATCH --qos=debug
#SBATCH --ntasks=672
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
# activate virual env?
module load python/3.6.1

. /gpfs/home/bsc51/bsc51927/PDBmapper-master/venv/bin/activate

PYTHONPATH=/gpfs/home/bsc51/bsc51927/PDBmapper-master/venv/lib:/gpfs/home/bsc51/bsc51927/PDBmapper-master/venv/lib/dask_jobqueue-0.7.2+6.ge7fcb31-py3.6.egg

pdbmapper -pid "$1" -psdb "$2" -vdb "$3" -csv --force -l -o "$4"
