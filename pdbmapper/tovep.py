mafFields = ['Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build',
             'Chromosome', 'Start_position', 'End_position', 'Strand',
             'Variant_Classification', 'Variant_Type', 'Reference_Allele',
             'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS', 'dbSNP_Val_Status',
             'Tumor_Sample_Barcode', 'Normal_Norm_SAmple_Barcode',
             'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele2', 'Tumor_Validation_Allele1',
             'Tumor_Validation_Allele2', 'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2',
             'Verification_Status', 'Validation_Status', 'Mutation_Status', 'Sequencing_Phase', 'Sequence_Source',
             'Validation_Method', 'Score', 'BAM_File', 'Sequencer', 'Tumor_Sample_UUID', 'Matched_Norm_Sample_UUID']

maf = "/home/vruizser/PhD/2018-2019/Immunity_interfaces_analysis/raw_data/mc3.v0.2.8.PUBLIC.maf"


with open(maf) as f:
    # get col names
    cols = f.readline().split(',')
    # parse file
    for line in f:
        pass
