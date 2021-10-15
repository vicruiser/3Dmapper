[![build status](http://img.shields.io/travis/vicruiser/3Dmapper/master.svg?style=flat)](https://travis-ci.com/username/vicruiser/3Dmapper) [![Join the chat at https://gitter.im/pdbmapper/community](https://badges.gitter.im/pdbmapper/community.svg)](https://gitter.im/pdbmapper/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

# Overview

3Dmapper is a command line tool based on R and Python programming languages that maps annotated genomic variants or positions to protein structures.

<p align="center">

<img src="./docs/pdbmapper_methods.png" width="800" heigh="400"/>

</p>

# Dependencies

## Installed by user

-   BLAST standalone software version >= 2.6. Follow these [instructions](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) to download and use the command line tool.
-   Python > 3.6
-   R version > 3.5

## Automatically installed

-   BCFtools 1.13

-   HTSlib 1.13

-   R packages:

    -   seqinr 4.2-4
    -   reshape2 1.4.4
    -   data.table 1.13.6
    -   dplyr 1.0.4
    -   plyr 1.8.6
    -   bio3d 2.4-2
    -   stringr 1.4.0
    -   tidyr 1.1.2\
    -   veriNA3d 1.0.3

# Install

In a Linux OS open a terminal and enter the following:

``` bash
git clone https://github.com/vicruiser/3Dmapper.git
cd 3Dmapper
pip install . 
```

# Tutorial

## 1. Generation of a local protein structures database with `makeinterfacedb`

### Overview

For each of the considered PDB files, `makeinterfacedb` automatically will: 
1) Extract the PDB chain sequences. 
2) BLAST PDB chain sequences (query) against the target proteome of interest (subject). 
3) Retrieve structural data of hits passing the selected homology filtering.

### Input files

-   A set of *PDB or CIF* files of interest (either real structures or models).
-   A target proteome BLAST database. More details on how to do this can be found in the example below.

### Example

In this example, we are going to map variants or positions to human protein structures.

1)  Retrieve the human proteome in FASTA format. This can be done using a public repository such as [UniProt](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz) or [Ensembl](http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz).

2)  Build a BLAST protein database for the BLAST search: 

``` markdown
makeblastdb -in target_proteome.fasta -dbtype prot -out target_proteome_db
```
This command will generate three files with name "target_proteome_db" and extensions `.phr`, `.pin` and `.psq`.

3)  Download a set of PDB files of interest. If you would like to retrieve structural data for proteins that currently do not have structure in PDB, you can download all the files in PDB and rely on sequence homology to find structural homologs.

4)  Build the structural database. You can input the PDB data in three different ways:

  - As a list of files directly specified in the command line:
``` markdown
makeinterfacedb -pdb pdb_dir/* --blast_db target_proteome_db
```
-   As a plain text file (.txt) with one PDB file path per line: 
``` markdown
makeinterfacedb -pdb list_pdbs.txt --blast_db target_proteome_db
``` 
-   As individual tasks (one per PDB file). This setting is useful when running jobs in parallel using job arrays or greasy. The results will be appended to the same output file.

``` markdown
makeinterfacedb -pdb file1.pdb --blast_db target_proteome_db
makeinterfacedb -pdb file2.pdb --blast_db target_proteome_db
makeinterfacedb -pdb file3.pdb --blast_db target_proteome_db
...
makeinterfacedb -pdb fileN.pdb --blast_db target_proteome_db
```
Note: Crystallographic artifacts can be removed with the option `-b`, based on the [BioLiP artifact ligand list](https://zhanggroup.org/BioLiP/ligand_list).

### Output interfaceDB.txt

The output is a 22 column tab-delimited file. In all cases, **"PDB chain"** refers to the extracted PDB chain or query sequence from each PDB file and **"Protein"** refers to the hit sequence found with the BLAST search against the target proteome. A more detailed description of the meaning of each column ID is specified in the table below.


| Column name                              | Notes                                                                                                                                           |
|:-----------------------------------------|:--------------------------------------------------------------------------------------------------------------------|
| **Protein_accession**  | Target protein ID Length of the target protein sequence   |
| **Protein_length**     | Length of protein sequence  |
| **Protein_position**                     | Positions relative to the target protein sequence  |
| **Protein_aa**                           | Amino acids corresponding to the target protein positions                                                                                           |
| **PDB_code**                             | PDB ID                                                                                                                                              |
| **PDB_chain**                            | ID of the template PDB protein chain                                                                                                                |
| **PDB_chain_length**                     | Length of PDB chain sequence                                                                                                                    |
| **PDB_3D_position**                      | Position in the PDB chain **structure**                                                                                                           |
| **PDB_seq_position**                     | Position in the PDB chain **sequence**                                                                                                              |
| **PDB_aa**                               | Amino acids corresponding to both the PDB sequence and 3D positions                                                                                 |
| **Evalue**                               | E-value of the alignment between the query or PDB chain sequence and the target protein                                                             |
| **Pident**                               | Identity percent between the query (PDB chain) and the target sequence (protein).                                                                   |
| **Protein_coverage**                     | Coverage (%) of the target protein by the PDB chain sequence                                                                                        |
| **Length_alignment**                     | Total length of the alignment between the query or PDB chain sequence and the target protein                                                        |
| **Interaction_type**                     | Type of interface interaction: "protein","nucleic" or "ligand". NA means no interaction which represents the positions of the rest of the structure |
| **PDB_interacting_chain**                | Interacting PDB chain ID with the template PDB chain. NA means no interaction which represents the positions of the rest of the structure           |
| **PDB_interacting_3D_position**          | Position in the interacting PDB chain **structure**                                                                                                 |
| **PDB_interacting_aa**                   | Amino acids corresponding to the interacting PDB structure positions                                                                                |
| **Interface_min_distance**               | Minimum existing distance between the pair of selected positions participating in the interface                                                     |
| **PDB_B\_factor**                        | Minimum B factor (or pLDDT in the case of AF2 models) observed in each PDB 3D position                                                              |
| **PDB_interacting_B\_factor**            | Minimum B factor (or pLDDT in the case of AF2 models) observed in each PDB interacting 3D position                                                  |
| **Protein_alignment_start**              | Alignment start position in target protein sequence                                                                                                 |
| **Protein_alignment_end**                | Alignment end position in target protein sequence                                                                                                   |
| **PDB_alignment_start**                  | Alignment start position in PDB chain protein sequence                                                                                              |
| **PDB_alignment_end**                    | Alignment end position in PDB chain protein sequence                                                                                                |
| **Structure_feature_id**                 | As `PDB_code`\_`Protein_accession`\_`PDB_chain`\_`PDB_interacting_chain`\_`Interaction_type`                                                        |

## 2. Split structural data DB

To reduce the computational workload during the mapping process, the structural data set generated in the previous step is divided by protein IDs into individual files by executing the following command:

``` markdown
makepsdb -psdb interfaces/interfacesDB.txt
```

The following directories and files are generated:

    |__DBs
         |_____makepsdb.log
         |_____makepsdb.report
         |_____psdb
                |_____prot_ID1.txt
                |_____prot_ID2.txt        
                |_____...        
                |_____prot_IDn.txt        

Files makepsdb.log and makepsdb.report report the progress of the executed command and then folder psdb contains allthe splitted files.

## 3. Split variants / annotated positions files

Similar to the previous step, we will perform a splitting of the variants or annotated positions files.

### Input files

#### Variants file

The input annotated genomic variants file must be either in [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format), [VEP](https://www.ensembl.org/info/docs/tools/vep/vep_formats.html#defaultout) or [MAF](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) default format. Additionally, a VEP-like format is admissible. This is in question same as VEP but not all the files are needed:

-   Uploaded_variation
-   Gene
-   Feature
-   Consequence
-   Protein_position
-   Amino_acids

#### Positions file

We can create a positions file using the the same format as for a variant file

### Example

``` markdown
makevariantsdb -vf variants.vep 
```

## 4. Map variants

### Example CSV output

``` markdown
mapper -pid ENSP00000356150 -psdb DBs/psdb -vdb DBs/varDB/ -ids dict_geneprot_GRCh38.txt  -f -csv 
```

### Example hdf output

``` markdown
mapper -pid ENSP00000356150 -psdb DBs/psdb -vdb DBs/varDB/ -ids dict_geneprot_GRCh38.txt  -f -hdf
```

# 3D visualization with ChimeraX

Results can be visualized running

``` markdown
makechimera xxxx
```

# Paralellization

3Dmapper can be parallelized parallelization. While `makepsdb` and `makevariantsdb` run with GNU parallel \[ref\], `mapper` uses the python module joblib.

An alternative to parallelaize 3Dmapper is to is to give as input the protein ids in individual tasks to perform a job array in a cluster computer?. The first task should write the initial files. The rest set the option `-force n` to prevent repeating innecesary steps.
