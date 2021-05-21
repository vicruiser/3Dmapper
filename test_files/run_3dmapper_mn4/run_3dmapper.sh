#!/bin/bash
#SBATCH -o /gpfs/scratch/bsc08/bsc08927/3Dmapper_analysis/debug/output_%A.txt
#SBATCH -e /gpfs/scratch/bsc08/bsc08927/3Dmapper_analysis/debug/errors_%A.txt
#SBATCH -J 3dmapper
#SBATCH --qos=debug
#SBATCH --ntasks=20
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=6

# load necessary modules
module load python/3.6.1

# activate environment
. /gpfs/home/bsc08/bsc08927/PDBmapper-master/venv/bin/activate

# indicate pythonpath (esto es porque el paquete dask no lo encontraba o algo así y lo instalé manualmente, no lo recuerdo bien. Puede que tú no tengas el mismo problema).
PYTHONPATH=/gpfs/home/bsc08/bsc08927/PDBmapper-master/venv/lib:/gpfs/home/bsc08/bsc08927/PDBmapper-master/venv/lib/dask_jobqueue-0.7.2+6.ge7fcb31-py3.6.egg

# run 3dmapper
wd=/gpfs/scratch/bsc51/bsc51927/3Dmapper_analysis
wd2=$wd/sars2
pdbmapper -pid E_protein_SARS-CoV-2 -psdb $wd2/DBs/psdb -vdb $wd2/DBs/varDB -o $wd2/3dmapper_results -dic $wd/needed_input_data/dict_geneprot_sars2.txt -csv -a -l
