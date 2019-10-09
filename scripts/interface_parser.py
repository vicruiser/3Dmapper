# Extract the info corresponding to the prot ID (Interface parse)
def interfaceParse(interfacesDB, protID):
    '''Parse input interfaces database to put it in the right format.
    Parameters
    ----------
    protID : str
        Ensemble protein id 
    interfacesDB_filepath : str
        DESCRIPTION MISSING!!
    Returns
    -------
    subset_interfaces_db
        DESCRIPTION MISSING!!
    '''
    # read file
    fp = glob.glob(interfacesDB + '/' + protID + '*.csv')[0] 
    prot_interface = pd.read_csv(fp, sep=' ', header=0)

    # store subspace
    subset_prot_interface = prot_interface[['pdb.id',
                                            'ensembl.prot.id',
                                            'temp.chain',
                                            'int.chain',
                                            'interaction',
                                            'resid_sseq',
                                            'mapped.real.pos',
                                            'pdb.pos']]
   
    # put it into right format
    subset_prot_interface.columns = \
        subset_prot_interface.columns.str.replace('\\.', '_')

    for col in ('resid_sseq', 'mapped_real_pos', 'pdb_pos'):
        subset_prot_interface.loc[:, col] = \
            subset_prot_interface[col].str.replace('-', ',')
        subset_prot_interface.loc[:, col] = \
            subset_prot_interface[col].str.split(',')

    subset_prot_interface = mt.explode(subset_prot_interface,
                                       ['resid_sseq',
                                        'mapped_real_pos',
                                        'pdb_pos'])
    subset_prot_interface.rename(columns={'mapped_real_pos':
                                          'Protein_position'}, inplace=True)
    # create region id for setID file
    subset_prot_interface['region_id'] = subset_prot_interface['pdb_id'] + \
        '_' + subset_prot_interface['ensembl_prot_id'] + '_' + \
        subset_prot_interface['temp_chain'] + '_' + \
        subset_prot_interface['int_chain'] + '_' + \
        subset_prot_interface['interaction']

    return subset_prot_interface
