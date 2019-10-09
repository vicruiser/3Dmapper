def vcfParser(VCF_dir, geneID, sep,  *args):
    '''Parse vcf file and put it in the right format if necessary.

    Parameters
    ----------
    VEP_dir : str
        Where th VEP file is
    geneID : str
        Ensemble gene ID corresponding to the translated protein ID
    
    Returns
    -------
    VCF_subset
        Data frame containing subset information vep
    '''
    # read the vep file
    VCF_file = glob.glob(VCF_dir + '/' + geneID + '*.csv')[0]
    print(VCF_file)
    VCF_subset = pd.read_csv(VCF_file, sep=sep, header=0)
    # when the thingy is varmap
    if args[0] == 'varmap':
        # drop columns
        VCF_subset = VCF_subset[['CHROMOSOME',
                                 'COORDS',
                                 'USER_BASE',
                                 'USER_VARIANT',
                                 'ENSEMBL_BASE',
                                 'VEP_CODING_BASE',
                                 'GENE',
                                 'GENE_ACC',
                                 'TRANSCRIPT',
                                 'CODON_CHANGE',
                                 'VEP_AA',
                                 'UNIPROT_AA',
                                 'AA_CHANGE',
                                 'CHANGE_TYPE',
                                 'RES_NAME',
                                 'RES_NUM']]
        VCF_subset.drop_duplicates()
        VCF_subset = VCF_subset.rename(columns={'RES_NUM': 'Protein_position'})
        VCF_subset['#Uploaded_variation'] = \
            VCF_subset['CHROMOSOME'].map(str) +\
            '_' + VCF_subset['COORDS'].map(str) + \
            '_' + VCF_subset['USER_BASE'].map(str) + \
            '_' + VCF_subset['USER_VARIANT'].map(str)
    # return the loaded subset
    return VCF_subset