<img src= "./docs/LogoPDBcopy2.png" width = "600" heigh = "300">

# Overview

3Dmapper is a tool to map annotated genomic positions to protein structures.

# Install and documentation

## Just variants mapping

## Generation of local interfaces database

# Running PDBmapper

```markdown
makepsdb -psdb p53_ep300_intdb.dat --out out/path/dir 
makevariantsdb -vcf PDBmapper/pdbmapper/data/p53_ep300_ExACvariants.vep -out ./test_pdbmapper/
pdbmapper -vardb test_pdbmapper/DBs/varDB/ -intdb test_pdbmapper/DBs/intDB/ -protid ENSP00000263253 ENSP00000482258 -out test_pdbmapper/ 
```
### Test

<!-- More details in ./test -->

## Reference

| <!--    | Setup                   | Command | Notes |
| :------ | :---------------------- | :------ |
| install | `pip install pdbmapper` | -->     |
<!-- 
| <!--   | Creating a CLI         | Command                                   | Notes |
| :----- | :--------------------- | :---------------------------------------- |
| import | `import fire`          |
| Call   | `fire.Fire()`          | Turns the current module into a Fire CLI. |
| Call   | `fire.Fire(component)` | Turns `component` into a Fire CLI. -->    |
<!-- 
| <!-- Using a CLI                                | Command                                 | Notes                                                    |
| :---------------------------------------------- | :-------------------------------------- | :------------------------------------------------------- |
| [Help](docs/using-cli.md#help-flag)             | `command --help` or `command -- --help` |
| [REPL](docs/using-cli.md#interactive-flag)      | `command -- --protid`                   | Protein id ensembl.                                      |
| [Separator](docs/using-cli.md#separator-flag)   | `command -- --separator=X`              | Sets the separator to `X`. The default separator is `-`. |
| [Completion](docs/using-cli.md#completion-flag) | `command -- --completion [shell]`       | Generates a completion script for the CLI.               |
| [Trace](docs/using-cli.md#trace-flag)           | `command -- --trace`                    | Gets a Fire trace for the command.                       |
| [Verbose](docs/using-cli.md#verbose-flag)       | `command -- --verbose`                  | -->                                                      | --> --> |

### Paralellization

PDBmapper has an option to speed up the running time by means of parallelization. While `makepsdb` and `makevariantsdb` run with GNU parallel [ref], `pdbmapper` has an algorithm in python. 

An alternative to parallelaize pdbmapper is to is to give as input the protein ids in individual tasks to perform a job array in a cluster computer?. The first task should write the initial files. The rest set the option `-force n` to prevent repeating innecesary steps. 

**Explain better or implemen a parallel option** 

# PDBmapper databases

## Protein structural features input format

The default format is a simple whitespace-separated format (columns may be separated by space or tab characters), containing six required columns. Any extra columns will be included in the results. Empty values are denoted by '-'.

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

Put an example here:

`ENSP0000023123  ENST0000023123  123    1wqs    A   ep300_1wqs_A`



## Interfaces database 

For our own research, a pre-computed database was used. We did generate this by following the protocol described in [ref to paper]. The format is a  10 column tab-delimited file: 

| Column name               | Notes                                                                                                                                          |
| :------------------------ | :--------------------------------------------------------------------------------------------------------------------------------------------- |
| **PDB_code**              | Biological Assembly PDB ID                                                                                                                     |
| **Protein_accession**     | Ensembl protein ID                                                                                                                             |
| **PDB_chain**             | Template chain (only protein)                                                                                                                  |
| **PDB_interacting_chain** | Interacting chain (protein, ligand or DNA)                                                                                                     |
| **Protein_length**        | Length of the Ensembl protein sequence (without including the gaps of the alignment)                                                           |
| **Protein_start_pos**     | Start of alignment in the Ensembl protein sequence (template chain)                                                                            |
| **Protein_end_pos**       | End of alignment in the Ensembl protein sequence (template chain)                                                                              |
| **Length_align**          | Length of the alignment between Ensembl protein sequence (subject) and PDB chain sequence (query)                                              |
| **Pident**                | Sequence Identity percent. This parameter serves as threshold. Only results with pident equal or higher to 50% are included                    |
| **Interaction**           | Type of interfacial interaction, i.e., “protein”,”nucleic” or “ligand”                                                                         |
| **PDB_aa**                | residues of aligned query (PDB chain) sequence (i.e.: only aligned and includes gaps.)                                                         |
| **Protein_aa**            | residues of aligned subject (Ensembl) sequence (i.e.: only aligned and includes gaps.)                                                         |
| **PDB_pos**               | index position of each residue of the aligned query (PDB chain) sequence starting from 1                                                       |
| **Protein_pos**           | index position of each residue of the aligned subject (Ensembl) sequence starting from 1                                                       |
| **PDB_align_pos**         | real index position of each residue of the aligned query (PDB chain) sequence                                                                  |
| **Protein_align_pos**     | Position of the interfacial residues on the Ensembl protein sequence. This is the column of your interest! (Protein position in the MC3 file). |
| **PDB_pos**               | Corresponding position of the interfacial residues on the PDB chain sequence                                                                   |

## Variant annotated files

The input annotated genomic variants file must be either in [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format), [VEP](https://www.ensembl.org/info/docs/tools/vep/vep_formats.html#defaultout) or [MAF] (https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) default format. Additionally, a VEP-like format is admisible. This is in question same as VEP but not all the files are needed: 

- Uploaded_variation
- Gene
- Feature
- Consequence
- Protein_position
- Amino_acids (optional but recommedable)


## Custom

Para poner tus propias databases han de cumplir con los siguientes requisitos de formato. 


# 3D visualization with ChimeraX

soon available

# FAQs

is not intended as a interface generator. 


You can use the [editor on GitHub](https://github.com/vicruiser/PDBmapper/edit/master/README.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/vicruiser/PDBmapper/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://help.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.

