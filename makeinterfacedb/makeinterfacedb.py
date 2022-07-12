# -*- coding: utf-8 -*-
import subprocess
import os
import glob
from joblib import Parallel, delayed, parallel_backend
import multiprocessing as mp7
import sys

#sys.path.append('/gpfs/scratch/bsc08/bsc08927/3dmapper_v2/3Dmapper/makeinterfacedb/')
#from memory_profiler import profile
from .parse_argv import parse_commandline
from .input_isfile import isfile

def parallel(parallel, njobs): 
    if parallel is True:
        if njobs != 1:
            num_cores = int(njobs)
        else:
            num_cores = mp.cpu_count()-1
    else:
        num_cores = njobs
    return(num_cores)

# by @Six in stackoverflow
def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    # from whichcraft import which
    from shutil import which

    return which(name) is not None

#fp=open('memory_profiler.log','w+')
#@profile(stream=fp)
def pipeline(f): 
    # check if R and BLAST are available
    if is_tool('R') is False: 
        sys.exit("Error: R not found or not executable.") 
    if is_tool('blastp') is False:
        sys.exit("Error: BLAST not found or not executable.") 
    # parse command line options
    args = parse_commandline()
    pdbf = os.path.basename(f)
    pdbid = pdbf.strip('.gz')

    chainseqs_outdir = os.path.join(args.out, 'pdb_chainseqs')
    blast_outdir = os.path.join(args.out, 'blast_results')
    interfaces_outdir = os.path.join(args.out, 'predicted_interfaces')
    mapped_outdir = os.path.join(args.out, 'structuralDB')
    # Extract PDB fasta sequences
    #chainseqs_outdir = os.path.join(args.out, 'pdb_chainseqs')
    #if not os.path.exists(chainseqs_outdir):
    #    os.makedirs(chainseqs_outdir)    

    root = os.path.join(os.path.dirname (os.path.abspath (__file__)), 'rscripts')
    rscript1 = os.path.join(root, "extract_pdb_chain_seqs.R")
    rscript2 = os.path.join(root, "filter_blast_output.R")
    rscript3 = os.path.join(root, "predict_PDB_interfaces_main.R")
    rscript4 = os.path.join(root, "map_interfaces.R")

    process1 = subprocess.Popen(["Rscript %s %s %s" % (rscript1, f, chainseqs_outdir) ], shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    #process1 = subprocess.Popen(["R --slave --quiet --no-restore --file=%s --args %s %s" % (rscript1, f, chainseqs_outdir) ], shell=True)
    process1.wait()
    out, err = process1.communicate()
    process1.terminate()
    # Run BLAST
    
    filtered_blast = []
    for file in glob.glob(os.path.join(chainseqs_outdir  ,(pdbid + "*_chain*.fasta"))): 
        process2 = subprocess.Popen(["blastp -query %s -db %s -outfmt '6 qseqid qlen sseqid slen qstart qend sstart send evalue length pident nident qseq sseq gaps' \
        -out %s.blast" % (file, args.blastdb, os.path.join(blast_outdir , os.path.basename(file)))] , shell =True, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
        process2.wait()
        out2, err2 = process2.communicate() 
        process2.terminate()
        # filter blast results
        if err2 is None and out2 == b'':
            process3 = subprocess.Popen(["Rscript %s %s %s %s %s %s %s.filtered" %(rscript2, os.path.join(blast_outdir , os.path.basename(file)+ '.blast'),
            args.pident, args.coverage, args.evalue, blast_outdir, os.path.basename(file))+ '.blast'], shell =True, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
            process3.wait()
            out3, err3 = process3.communicate() 
            if os.path.isfile(os.path.join(blast_outdir , os.path.basename(file)+ '.filtered.blast')):
                filtered_blast.append(True)
            #else: 
            #    filtered_blast.append(False)
            process3.terminate()
        else: 
            print(out2.decode('utf8'))
    # Predict interfaces # To be included? ideally hidrophobicity, ss type and full structure
    if any(filtered_blast) is True : 
        
        process4 = subprocess.Popen(["Rscript %s %s %s %s %s %s %s %s" % (rscript3, root, f, interfaces_outdir, args.dist, args.type, args.int, args.biolip)] , shell =True, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
        process4.wait()
        out4, err4 = process4.communicate() 
        process4.terminate()
        # Map PDB and protein information
        
        process5 = subprocess.Popen(["Rscript %s %s %s %s %s %s" % (rscript4, root, pdbid, blast_outdir, interfaces_outdir, mapped_outdir)], shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
        process5.wait()
        out5, err5 = process5.communicate() 
        process5.terminate()

def main():
    # parse command line options
    args = parse_commandline()
    # Extract PDB fasta sequences

    chainseqs_outdir = os.path.join(args.out, 'pdb_chainseqs')
    if not os.path.exists(chainseqs_outdir):
        os.makedirs(chainseqs_outdir, exist_ok=True)

    blast_outdir = os.path.join(args.out, 'blast_results')
    if not os.path.exists(blast_outdir):
        os.makedirs(blast_outdir, exist_ok=True)
    
    interfaces_outdir = os.path.join(args.out, 'predicted_interfaces')
    if not os.path.exists(interfaces_outdir):
        os.makedirs(interfaces_outdir, exist_ok=True)
    
    mapped_outdir = os.path.join(args.out, 'structuralDB')
    if not os.path.exists(mapped_outdir):
        os.makedirs(mapped_outdir, exist_ok=True)

    if args.parallel is True: 
        num_cores = parallel(args.parallel, args.njobs)
        if isfile(args.pdb) == "is_file":
            Parallel(n_jobs=num_cores)(delayed(pipeline)(f) for f in args.pdb)
        elif  isfile(args.pdb) == "list_files": 
            for f in args.pdb: 
                with open(f) as list_int_files:
                    int_f = list_int_files.read().splitlines()
                    Parallel(n_jobs=num_cores)(delayed(pipeline)(int_file) for int_file in int_f)
    else:
        if isfile(args.pdb) == 'list_files':
            for f in args.pdb:
                with open(f) as list_int_files:
                            int_f = list_int_files.read().splitlines()
                            for int_infile in int_f:
                                pipeline(int_infile)
        elif isfile(args.pdb) == 'is_file':
            for f in args.pdb:
                pipeline(f)
        elif isfile(args.pdb) == 'file_not_recognized':
            print(['Error: Not such file: %s',  (f)])
            exit(-1)

    
    

    

