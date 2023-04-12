if ! [ -x "$(command -v git)" ]; then
  echo 'Error: git is not installed.' >&2
  exit 1
else


mkdir -p r_dependencies
cd r_dependencies

echo Installing seqinr...
wget https://cran.r-project.org/src/contrib/Archive/seqinr/seqinr_4.2-4.tar.gz
R CMD INSTALL seqinr_4.2-4.tar.gz
echo Installing seqinr...Done!

echo Installing reshape2...
wget https://cran.r-project.org/src/contrib/Archive/reshape2/reshape2_1.4.3.tar.gz
R CMD INSTALL reshape2_1.4.3.tar.gz
echo Installing reshape2...Done!

echo Installing data.table...
wget https://cran.r-project.org/src/contrib/Archive/data.table/data.table_1.13.6.tar.gz
R CMD INSTALL data.table_1.13.6.tar.gz
echo Installing data.table...Done!

echo Installing dplyr...
wget https://cran.r-project.org/src/contrib/Archive/dplyr/dplyr_1.0.4.tar.gz
R CMD INSTALL dplyr_1.0.4.tar.gz
echo Installing dplyr...Done!

echo Installing plyr...
wget https://cran.r-project.org/src/contrib/Archive/plyr/plyr_1.8.5.tar.gz
R CMD INSTALL plyr_1.8.5.tar.gz
echo Installing plyr...Done!

echo Installing bio3d...
wget https://cran.r-project.org/src/contrib/Archive/bio3d/bio3d_2.4-1.tar.gz
R CMD INSTALL bio3d_2.4-1.tar.gz
echo Installing bio3d...Done!

echo Installing stringr...
wget https://cran.r-project.org/src/contrib/Archive/stringr/stringr_1.3.1.tar.gz
R CMD INSTALL stringr_1.3.1.tar.gz
echo Installing stringr...Done!

echo Installing tidyr...
wget https://cran.r-project.org/src/contrib/Archive/tidyr/tidyr_1.1.2.tar.gz
R CMD INSTALL tidyr_1.1.2.tar.gz
echo Installing tidyr...Done!

echo Installing flock...
wget https://cran.r-project.org/src/contrib/Archive/flock/flock_0.5.tar.gz
R CMD INSTALL flock_0.5.tar.gz
echo Installing flock...Done!

echo Installing GlobalOptions...
wget https://cran.rstudio.com/src/contrib/GlobalOptions_0.1.2.tar.gz
R CMD INSTALL GlobalOptions_0.1.2.tar.gz
echo Installing GlobalOptions...Done!

echo Installing shape...
wget https://cran.r-project.org/src/contrib/Archive/shape/shape_1.4.5.tar.gz
R CMD INSTALL shape_1.4.5.tar.gz
echo Installing shape...Done!

echo Installing circlize...
wget https://cran.r-project.org/src/contrib/Archive/circlize/circlize_0.4.14.tar.gz
R CMD INSTALL circlize_0.4.14.tar.gz
echo Installing circlize...Done!

echo Installing misc3d...
wget https://cran.rstudio.com/src/contrib/misc3d_0.9-1.tar.gz
R CMD INSTALL misc3d_0.9-1.tar.gz
echo Installing misc3d...Done!

echo Installing plot3D...
wget https://cran.rstudio.com/src/contrib/plot3D_1.4.tar.gz
R CMD INSTALL plot3D_1.4.tar.gz
echo Installing plot3D...Done!

echo Installing RANN...
wget https://cran.rstudio.com/src/contrib/RANN_2.6.1.tar.gz
R CMD INSTALL RANN_2.6.1.tar.gz
echo Installing RANN...Done!

echo Installing veriNA3d...
wget https://mmb.irbbarcelona.org/gitlab/dgallego/veriNA3d-dev/-/archive/master/veriNA3d-dev-master.zip -O ./veriNA3d_0.99.0.zip
unzip veriNA3d_0.99.0.zip
mv veriNA3d-*master* veriNA3d_0.99.0
R CMD build veriNA3d_0.99.0 --no-build-vignettes
R CMD INSTALL veriNA3d*.tar.gz
echo Installing veriNA3d...Done!

cd ..
fi
