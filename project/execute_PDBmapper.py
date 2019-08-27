#!/usr/bin/python3

import parse_argv as pa

# main function:
def main():
    # get command line options
    options = pa.parse_commandline()
    # read sequences from fasta file, and catch error reading file
    try:
        sequences = readSequences(open(options.fasta))
    except IOError:
        print "ERROR: cannot open or read fasta input file:", fastafile
        exit(-1)
    #print 'mean length of sequences: %s'%(str(retrieve_info(sequences)))
    families = ['Rhodanese']
    desired = store_desired_sequences(sequences,families,'test')
    write_sequences_fasta(options.fasta,desired,'test')


if __name__ == "__main__":
    main()

# last line