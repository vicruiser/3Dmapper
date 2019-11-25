
def f():
    # PDBmapper accepts single or multiple protein ids
    # as input as well as prot ids stored in a file
    for pid in args.protid:
        open("jeje")
        print("BUENAS CHAVALES")
        # check if input is a file
        try:

            with open(pid) as f:
                lines = f.read().splitlines()

                input = "file"
                # for every prot id

                # get geneID

        except:
            input = "not_file"

        if input == "file":
           # for pids in lines:
            ensemblIDs = translate_ensembl(pids)
            geneID = ensemblIDs['geneID']
            # run PDBmapper
            try:
                PDBmapper(pids, geneID, int_db_dir,
                          vcf_db_dir, args.out, args.pident)

                # error handling
            except IOError:
                log = open(out_dir + '/log_ensembl.File', 'a')
                log.write('Warning: ' + pids +
                          ' has no ENGS.\n')
                continue
        # input is not a file but one or more protein ids
        # given in command line
        elif input == "not_file":
                # for prot id get the gene id
            try:
                ensemblIDs = translate_ensembl(pid)
                geneID = ensemblIDs['geneID']
                # run PDBmapper
                try:
                    PDBmapper(pid, geneID, int_db_dir,
                              vcf_db_dir, args.out, args.pident)
            # error handling
                except IOError:
                    next
            except IOError:
                next
        else:
            print("wrong input!!")


# original function

def f():
            # PDBmapper accepts single or multiple protein ids
            # as input as well as prot ids stored in a file
    for pid in args.protid:
        open("jeje")
        print("BUENAS CHAVALES")
        # check if input is a file
        try:

            with open(pid) as f:
                lines = f.read().splitlines()
                # for every prot id
                for pids in lines:
                   # get geneID
                    try:
                        ensemblIDs = translate_ensembl(pids)
                        geneID = ensemblIDs['geneID']
                        # run PDBmapper
                        try:
                            PDBmapper(pids, geneID, int_db_dir,
                                      vcf_db_dir, args.out, args.pident)
                            del()
                            gc.collect()
                     # error handling
                        except IOError:
                            continue
                    except IOError:
                        log = open(out_dir + '/log_ensembl.File', 'a')
                        log.write('Warning: ' + pids +
                                  ' has no ENGS.\n')
                        continue
        # input is not a file but one or more protein ids
        # given in command line
        except:
            print("SE HA PASAO DE ROSCA")
            # for prot id get the gene id
            try:
                ensemblIDs = translate_ensembl(pid)
                geneID = ensemblIDs['geneID']
                # run PDBmapper
                try:
                    PDBmapper(pid, geneID, int_db_dir,
                              vcf_db_dir, args.out, args.pident)
            # error handling
                except IOError:
                    next
            except IOError:
                next
