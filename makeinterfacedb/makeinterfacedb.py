# -*- coding: utf-8 -*-
import subprocess
import os
import glob
from .parse_argv import parse_commandline



def main():

    # parse command line options
    args = parse_commandline()

    pdbf = os.path.basename(args.pdb)
    pdbid = pdbf.strip('.gz')

    # Extract PDB fasta sequences
    chainseqs_outdir = os.path.join(args.out, 'PDBchainseqs')
    if not os.path.exists(chainseqs_outdir):
        os.makedirs(chainseqs_outdir)    

    root = os.path.join(os.path.dirname (os.path.abspath (__file__)), 'rscripts')
    rscript1 = os.path.join(root, "extract_pdb_chain_seqs.R")
    rscript2 = os.path.join(root, "filter_blast_output.R")
    rscript3 = os.path.join(root, "predict_PDB_interfaces_main.R")
    rscript4 = os.path.join(root, "map_interfaces.R")

    #subprocess.call ("Rscript " + rscript1 +" " + args.pdb + " " + chainseqs_outdir, shell=True)

    process1 = subprocess.Popen(["Rscript %s %s %s" % (rscript1, args.pdb, chainseqs_outdir) ], shell=True)
    process1.wait()
    # Run BLAST
    blast_outdir = os.path.join(args.out, 'BLASTresults')
    if not os.path.exists(blast_outdir):
        os.makedirs(blast_outdir)
    for file in glob.glob(os.path.join(chainseqs_outdir  ,"*_chain*.fasta")): 

        # subprocess.call("blastp -query " + file + \
        # " -db " + args.blastdb + \
        # " -outfmt '6 qseqid qlen sseqid slen qstart qend sstart send evalue length pident nident qseq sseq gaps' \
        # -out " + os.path.join(blast_outdir , os.path.basename(file)+ '.blast') , shell =True)

        process2 = subprocess.Popen(["blastp -query %s -db %s -outfmt '6 qseqid qlen sseqid slen qstart qend sstart send evalue length pident nident qseq sseq gaps' \
        -out %s.blast" % (file, args.blastdb, os.path.join(blast_outdir , os.path.basename(file)))] , shell =True)
        process2.wait()
    
        # filter blast results
        process3 = subprocess.Popen(["Rscript %s %s %s %s %s %s %s.filtered" %(rscript2, os.path.join(blast_outdir , os.path.basename(file)+ '.blast'),
        args.pident, args.evalue, args.coverage, blast_outdir, os.path.basename(file))], shell =True)
        process3.wait()
  
    # Predict interfaces
    interfaces_outdir = os.path.join(args.out, 'interfacesDB')
    if not os.path.exists(interfaces_outdir):
        os.makedirs(interfaces_outdir)
    process4 = subprocess.Popen(["Rscript %s %s %s %s %s %s" % (rscript3, root, args.pdb, interfaces_outdir, args.dist, args.type)] , shell =True)
    process4.wait()

    # Map PDB and protein information
    process5 = subprocess.Popen(["Rscript %s %s %s %s %s %s" % (rscript4, root, pdbid, blast_outdir, interfaces_outdir, args.out)], shell=True)
    process5.wait()

    

