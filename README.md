<img src= "./pics/logo/image4.png" width = "600" heigh = "300">

# Introduction

## Motivation
PDB mapper is a tool to map kamsldkmalskdmaksdmalkdm because aoskdakdosakdpoaksd. 

## Overview
The PDBmaper program supports three different search methods:

    Users can ssh into the instance and run stand-alone BLAST+ at the command-line.
    The instance has a simple webpage for search submission.
    The instance supports the NCBI-BLAST Common URL API interface.



# Usage

## Running PDBmapper

```markdown
python PDBmapper -protid ensmblprotid 
```

## Reference

| Setup   | Command             | Notes
| :------ | :------------------ | :---------
| install | `pip install fire`  |

| Creating a CLI | Command                | Notes
| :--------------| :--------------------- | :---------
| import         | `import fire`          |
| Call           | `fire.Fire()`          | Turns the current module into a Fire CLI.
| Call           | `fire.Fire(component)` | Turns `component` into a Fire CLI.

Using a CLI                                     | Command                                 | Notes
:---------------------------------------------- | :-------------------------------------- | :----
[Help](docs/using-cli.md#help-flag)             | `command --help` or `command -- --help` |
[REPL](docs/using-cli.md#interactive-flag)      | `command -- --protid`                   | Protein id ensembl.
[Separator](docs/using-cli.md#separator-flag)   | `command -- --separator=X`              | Sets the separator to `X`. The default separator is `-`.
[Completion](docs/using-cli.md#completion-flag) | `command -- --completion [shell]`       | Generates a completion script for the CLI.
[Trace](docs/using-cli.md#trace-flag)           | `command -- --trace`                    | Gets a Fire trace for the command.
[Verbose](docs/using-cli.md#verbose-flag)       | `command -- --verbose`                  |

# PDBmapper databases

## Interfaces database

## Variant annotated files

### ClinVar

## Custom

Para poner tus propias databases han de cumplir con los siguientes requisitos de formato. 




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

