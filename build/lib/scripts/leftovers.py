import pandas as pd
from timeit import default_timer as timer

class Variant:
    
    def __init__(self, geneID, crossref_filepath, VEPs_dir):
        self.gene = geneID
        self.cr = crossref_filepath
        self.VEPdir = VEPs_dir
        self.VEP_filename = None
        aa_variants = None
    
    def retrieve(self, geneID, crossref_filepath, VEPs_dir): 
        crossref_df = pd.read_csv(self.cr, sep = "\t")
        self.VEP_filename = crossref_df[crossref_df['ids'].str.contains(self.gene)]["VEPfile"].values[0] 
        # extract the info for that region id from the VEP file  
        VEP_filepath = self.VEPdir + '/' + self.VEP_filename
        # read the selected VEP file
        VEP_file = pd.read_csv(VEP_filepath, sep = "\t", skiprows = 42)
        # subset data frame regarding the selected gene id. 
        # only subset variants affecting the amino acid sequence 
        self.aa_variants = VEP_file[VEP_file['Gene'].str.contains(self.gene) & (VEP_file['Amino_acids'] != "-")] 
        
# Put in Script "execute_PDBmapper.py"
    # geneID = sys.argv[1]
    # crossref_file = sys.argv[2]
    # VEPs_dir = sys.argv[3]

# for each of the regions ids of interest (from the interfaces DB) see where is the match in the VEP file. 
geneID = 'ENSG00000171502'
crossref_filepath = "/home/vruizser/PhD/2018-2019/git/PDBmapper/project/geneids_VEPfiles_crossref.txt"
VEPs_dir = "/home/vruizser/PhD/2018-2019/PanCancer_data/vep_output"

# execute class code
start = timer()
var = Variant.retrieve(geneID, crossref_filepath, VEPs_dir)
end = timer()
print(end - start) 
#var.retrieve()
#print(var.VEP_filename)



# Put in Script "execute_PDBmapper.py"
    # geneID = sys.argv[1]
    # crossref_file = sys.argv[2]
    # VEPs_dir = sys.argv[3]

# for each of the regions ids of interest (from the interfaces DB) see where is the match in the VEP file. 
geneID = 'ENSG00000171502'
crossref_filepath = "/home/vruizser/PhD/2018-2019/git/PDBmapper/project/geneids_VEPfiles_crossref.txt"
VEPs_dir = "/home/vruizser/PhD/2018-2019/PanCancer_data/vep_output"

def GetVariantsInfo (geneID, crossref_filepath, VEPs_dir):
    # read the crossref file
    crossref_df = pd.read_csv(crossref_filepath, sep = "\t")

    # get the vep filename associated to the gene ID
    VEP_filename = crossref_df[crossref_df['ids'].str.contains(geneID)]["VEPfile"].values[0] 
    print(VEP_filename)
    # extract the info for that region id from the VEP file  
    VEP_filepath = VEPs_dir + '/' + VEP_filename

    # read the selected VEP file
    VEP_file = pd.read_csv(VEP_filepath, sep = "\t", skiprows = 42)

    # subset data frame regarding the selected gene id. 
    # only subset variants affecting the amino acid sequence 
    aa_variants = VEP_file[VEP_file['Gene'].str.contains(geneID) & (VEP_file['Amino_acids'] != "-")] 
    return(aa_variants)
    
start = timer()
jeje = GetVariantsInfo(geneID, crossref_filepath, VEPs_dir)
end = timer()
print(end - start) 

   # get the translated id. If gene then prot and viceversa. 
    def translator(self):
    
    def VEPinfo(self, crossref_filepath, VEPs_dir): 
        
        crossref_df = pd.read_csv(crossref_filepath, sep = "\t")
        # get the vep filename associated to the gene ID
        VEP_filename = crossref_df[crossref_df['ids'].str.contains(geneID)]["VEPfile"].values[0] 