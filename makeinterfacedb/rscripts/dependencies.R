#! /usr/bin/Rscript

# https://stackoverflow.com/questions/50378639/how-to-install-package-with-pip-that-has-access-to-an-r-script
warn.conflicts = FALSE
library(veriNA3d)
requiredPackages = c('tidyr', 'stringr', 'bio3d',  'plyr','dplyr', 'data.table','reshape2','seqinr') #'parallel',
for (p in requiredPackages) {
  if (!require(p, character.only = TRUE))
    install.packages(p)
  suppressMessages(library(p, character.only = TRUE))
}


#   - seqinr 4.2-4
# - reshape2 1.4.4
# - data.table 1.13.6 
# - dplyr 1.0.4 
# - plyr 1.8.6
# - bio3d 2.4-2
# - stringr 1.4.0
# - tidyr 1.1.2      
# - veriNA3d 1.0.3 


if(!require("veriNA3d")){
  system('wget mmb.irbbarcelona.org/gitlab/dgallego/veriNA3d-dev/repository/archive.zip?ref=master -O ~/PDBmapper-master/veriNA3d_0.99.0.zip')
  system('unzip ~/PDBmapper-master/veriNA3d_0.99.0.zip')
  system('mv ~/PDBmapper-master/veriNA3d-*master* ~/PDBmapper-master/veriNA3d_0.99.0')
  system('R CMD build ~/PDBmapper-master/veriNA3d_0.99.0 --no-build-vignettes')
  system('R CMD INSTALL ~/PDBmapper-master/veriNA3d*.tar.gz')
  
}
invisible(suppressPackageStartupMessages(library("veriNA3d", character.only = TRUE)))



