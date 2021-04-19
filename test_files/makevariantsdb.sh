#!/bin/bash
#SBATCH -o /gpfs/home/bsc51/bsc51927/pdbmapper/debug/output_%A.txt
#SBATCH -e /gpfs/home/bsc51/bsc51927/pdbmapper/debug/errors_%A.txt
#SBATCH -J mkvardb
#SBATCH --qos=debug
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=48

module load python/3.6.1

. /gpfs/home/bsc51/bsc51927/PDBmapper-master/venv/bin/activate

makevariantsdb -vf "$1" -o "$2" --force
