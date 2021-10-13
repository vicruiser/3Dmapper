[![build status](
  http://img.shields.io/travis/vicruiser/3Dmapper/master.svg?style=flat)](
 https://travis-ci.com/username/vicruiser/3Dmapper)
 [![Join the chat at https://gitter.im/pdbmapper/community](https://badges.gitter.im/pdbmapper/community.svg)](https://gitter.im/pdbmapper/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

# Overview

3Dmapper is a Python and R tool to map annotated genomic variants or positions to protein structures.
 <p align="center">
<img src= "./docs/pdbmapper_methods.png" width = "800" heigh = "400">
</p>

# Install 

```bash
git clone https://github.com/vicruiser/3Dmapper.git
cd 3Dmapper
pip install . 
```

# Generation of local interfaces database

## Requirements 
 - A set of *PDB or CIF* files of interest (either real structures or models).  
 - BLAST standalone software. Follow these [instructions](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) to download and use the command line tool.
 - A target proteome built running the BLAST command `makeblastdb` with the set of protein sequences of interest.
 - R version > 3.0.0
 
## Overview

For each of the PDB files we are going to consider, `makeinterfacedb` automatically will do the following:  
  1) Extract its PDB chain sequences. 
  2) BLAST PDB chain sequences against the proteome of interest. 
  3) Predict interfaces of hits passing the selected homology filtering. 
  4) Merge PDB and proteome data. 

## Example

If we were interested in mapping variants or positions to the protein structures of human, we have to download the human proteome from [UniProt](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz) or [Ensembl](http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz) for example. 

Next, we build the target proteome database for the BLAST search executing the following command:
```markdown
makeblastdb -in target_proteome.fasta -dbtype protein -out target_proteome
```

Them 

```markdown
makeinterfacedb -pdb file.pdb --blast-db proteins_db  -b
```

## Results
The output interfaces database is a 22 column tab-delimited file. In all cases, **"PDB chain"** refers to the extracted PDB chain or query sequence from each PDB file and **"Protein"** refers to the hit sequence found with the Blast search against the target proteome. A more detailed description of the meaning of each column ID is specified in the table below.

| Column name               | Notes                                                                                                                                          |
| :------------------------ | :--------------------------------------------------------------------------------------------------------------------------------------------- |
| **Protein_accession**     | Target protein ID |
| **Protein_length**        | Length of the target protein sequence |  
| **Protein_position**      | Positions relative to the target protein sequence |
| **Protein_aa**            | Amino acids corresponding to the target protein positions |                                  
| **PDB_code**              | PDB ID        |
| **PDB_chain**             | ID of the template PDB protein chain |
| **PDB_chain_length**      | Length of the PDB chain sequence |
| **PDB_3D_position**       | Position in the PDB chain **structure** |
| **PDB_seq_position**      | Position in the PDB chain **sequence** | 
| **PDB_aa**      | Amino acids corresponding to both the PDB sequence and 3D positions |  
| **Evalue**      | E-value of the alignment between the query or PDB chain sequence and the target protein |  
| **Pident**      | Identity percent between the query (PDB chain) and the target sequence (protein). | 
| **Protein_coverage**      | Coverage (%) of the target protein by the PDB chain sequence |     
| **Length_alignment**      | Total length of the alignment between the query or PDB chain sequence and the target protein |    
| **Interaction_type**           | Type of interface interaction: “protein”,”nucleic” or  “ligand”. NA means no interaction which represents the positions of the rest of the structure |
| **PDB_interacting_chain** | Interacting PDB chain ID with the template PDB chain. NA means no interaction which represents the positions of the rest of the structure |
| **PDB_interacting_3D_position** | Position in the interacting PDB  chain **structure**  |
| **PDB_interacting_aa** | Amino acids corresponding to the interacting PDB structure positions                              |
| **Interface_min_distance** | Minimum existing distance between the pair of selected positions participating in the interface|
| **PDB_B_factor**     | Minimum B factor (or pLDDT in the case of AF2 models) observed in each PDB 3D position                                                                           |
| **PDB_interacting_B_factor**       | Minimum B factor (or pLDDT in the case of AF2 models) observed in each PDB interacting 3D position                                                                               |
| **Protein_alignment_start**          | Alignment start position in target protein sequence                   |
| **Protein_alignment_end**                | Alignment end position in target protein sequence                   |
| **PDB_alignment_start**           | Alignment start position in PDB chain protein sequence|
| **PDB_alignment_end**                | Alignment end position in PDB chain protein sequence |
| **Structure_feature_id**               | As `PDB_code`\_`Protein_accession`\_`PDB_chain`\_`PDB_interacting_chain`\_`Interaction_type` |


# Split variants / position files
```markdown
makevariantsdb -vf variants.vep 
```

# Split interface DB
```markdown
makepsdb -psdb interfaces/interfacesDB.txt -s
```

The default input format is a simple whitespace-separated format (columns may be separated by space or tab characters), containing six required columns. Any extra columns will be included in the results. Empty values are denoted by '-'.

- PROTEIN_ACCESION: Ensembl or Uniprot accesion identifier. 
- TRANSCRIPT_ACCESION: Ensembl identifier. 
- PROTEIN_POSITION: residue position on the protein. 
- PDB_CODE
- PDB_CHAIN
- PROTEIN_FEATURE_ID : structural identifier made up by the user, e.g.: PROTEIN_POSITION + PDB_CODE + PDB_CHAIN

Optional but recommendable for downstream analysis:
- AMINO_ACIDS: original residue. 
- PDB_AMINO_ACIDS: original residue in the PDB file. 
- PIDENT: sequence identity percentage between the protein sequence and the PDB chain. When available, this can be moludated as a filter parameter with the option `--pident`. 

Example

`ENSP0000023123  ENST0000023123  123    1wqs    A   ep300_1wqs_A`

# Map variants
```markdown
mapper -pid ENSP00000356150 -psdb DBs/psdb -vdb DBs/varDB/ -ids dict_geneprot_GRCh38.txt  -f 
```

The input annotated genomic variants file must be either in [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format), [VEP](https://www.ensembl.org/info/docs/tools/vep/vep_formats.html#defaultout) or [MAF] (https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) default format. Additionally, a VEP-like format is admisible. This is in question same as VEP but not all the files are needed: 

- Uploaded_variation
- Gene
- Feature
- Consequence
- Protein_position
- Amino_acids (optional but recommedable)


# 3D visualization with ChimeraX

Results can be visualized running

```markdown
makechimera xxxx
```

# Dependencies
- samtools

# Paralellization

3Dmapper can be parallelized parallelization. While `makepsdb` and `makevariantsdb` run with GNU parallel [ref], `mapper` uses the python module joblib. 

An alternative to parallelaize 3Dmapper is to is to give as input the protein ids in individual tasks to perform a job array in a cluster computer?. The first task should write the initial files. The rest set the option `-force n` to prevent repeating innecesary steps. 
