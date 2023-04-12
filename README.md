[[![Build Status](https://app.travis-ci.com/vicruiser/3Dmapper.svg?branch=master)](https://app.travis-ci.com/vicruiser/3Dmapper) [![Join the chat at https://gitter.im/3d_mapper/community](https://badges.gitter.im/pdbmapper/community.svg)](https://gitter.im/3d_mapper/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

## Overview

3Dmapper is a command line tool based on R and Python programming languages that maps annotated genomic variants or positions to protein structures.

## Dependencies

-   BLAST standalone software version >= 2.6. Follow these [instructions](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) to download and use the command line tool.
-   Python > 3.6
-   R version > 3.5
-   [GNU parallel](https://www.gnu.org/software/parallel/)

## Quick install

Execute the following code in Terminal:

``` bash
git clone https://github.com/vicruiser/3Dmapper.git
cd 3Dmapper
pip install . 
sh r_dependencies.sh
```

## Tutorial
To learn how to use 3Dmapper, please refer to our [Wiki](https://github.com/vicruiser/3Dmapper/wiki). 
