## ----setup, include = FALSE---------------------------------------------------

knitr::opts_chunk$set(
    collapse   = TRUE,
    comment    = "#>",
    message    = FALSE,
    warning    = FALSE,
    fig.width  = 6,
    fig.height = 2,
    fig.wide   = TRUE
)


## ----?cif_accessors-----------------------------------------------------------
library(veriNA3d)
?cif_accessors

## ----cifParser(pdbID)---------------------------------------------------------
## To parse a local mmCIF file:
# cif <- cifParser("your-file.cif")
## To download from PDB directly:
cif <- cifParser("1bau")
cif
## To see the coordinates:
coords <- cifAtom_site(cif)
head(coords)

## ----cifAsPDB(cif)------------------------------------------------------------
pdb <- cifAsPDB(cif)
pdb

## ----?queryFunctions----------------------------------------------------------
?queryFunctions

## ----"queryTechnique(pdbID, verbose=TRUE)"------------------------------------
## Run a query for the first time, which will access the API
tech <- queryTechnique("4KQX", verbose=TRUE)

## Run the same query for the second time, which will get it from memory
tech <- queryTechnique("4KQX", verbose=TRUE)

## See result
tech

## ----"Example1: queryAPI(ID, API, string1, string2)"--------------------------
atpsummary <- queryAPI(ID="ATP", API="ebi", 
                        string1="pdb/compound/summary/", string2="")
str(atpsummary$ATP)

## ----"Example2: queryAPI(ID, API, string1, string2)"--------------------------
PISAsummary <- queryAPI(ID="3gcb", API="ebi", 
                        string1="pisa/noofinterfaces/", string2="0")
str(PISAsummary$"3gcb")

## ----"eRMSD(structure1, structure2)"------------------------------------------
## Parse cif file
cif <- cifParser("2d18")
## Select a couple of models
model1 <- selectModel(cif=cif, model=1)
model3 <- selectModel(cif=cif, model=3)

## Calculate the eRMSD
eRMSD(cif1=model1, cif2=model3)

## The RMSD can also be calculated easily
RMSD(cif1=model1, cif2=model3)


## ----"trimsphere(structure, chain)"-------------------------------------------
## Parse human ribosome - takes around 12 seconds in R-3.5
cif <- cifParser("6ek0")

## Query entities and check them
ent <- queryEntities("6ek0")
head(ent[, c("entity_id", "molecule_name", "in_chains")])

## Generate a smaller pdb with the 60S ribosomal protein L8
chain <- "LA"
protL8 <- trimSphere(cif, chain=chain, cutoff=0)
protL8

## The same command with the argument file would save it directly:
trimSphere(cif, chain=chain, cutoff=0, file="output.pdb")

## ----"trimsphere(structure, selection)"---------------------------------------
## Load bio3d library
library(bio3d)

## Get pdb object from CIF
pdb <- cifAsPDB(cif)

## Get list of ligands in the human ribosome 6EK0
queryLigands("6ek0", onlyligands=T)

## Get the atomic index for a desired ligand
HMTligand_inds <- pdb$atom$eleno[which(pdb$atom$resid == "HMT")]

## Use bio3d function to select the ligand using its atom indices
sel <- atom.select(pdb, eleno=HMTligand_inds)

## Get substructure and sorroundings at 10 Angstroms
HTMligand <- trimSphere(pdb, sel=sel, cutoff=5)

## And generate file to visualize it
trimSphere(pdb, sel=sel, file="output2.pdb", cutoff=5)

## ----"trimsphere(structure, selection2)"--------------------------------------
## Parse another pdb for this example
pdb <- cifAsPDB("1nyb")

## Find region of interaction between RNA and protein
data <- findBindingSite(pdb, select="RNA", byres=TRUE)

## Get atom indices from interacting region molecules
eleno <- append(data$eleno_A, data$eleno_B)

## Select using bio3d
sel <- atom.select(pdb, eleno=eleno)

## Get substructure
trimSphere(pdb, sel=sel, file="interacting_site.pdb", verbose=FALSE)

## ----measureElenoDist(structure, refeleno, eleno, cutoff)---------------------
## Parse pdb
pdb <- cifAsPDB("1nyb")

## Select P atoms by element number (eleno)
ind <- which(pdb$atom$elety == "P")
eleno <- pdb$atom$eleno[ind]

## Count number of phosphates
total <- length(eleno)

## Execute function to measure the distances
P_distances <- measureElenoDist(pdb=pdb, refeleno=eleno, eleno=eleno, 
                                n=total, cutoff=100)

## ----getRNAList(release, threshold)-------------------------------------------
## Get non-redundant list from Leontis website
rnalist <- getRNAList(release="3.47", threshold="2.0A")
head(rnalist)

## ----'getAltRepres(rnalist, type="protRNA")'----------------------------------
## Set progressbar=TRUE to see the progress
protrna <- getAltRepres(rnalist=rnalist, type="protRNA", progressbar=FALSE)
head(protrna)

## ----represAsDataFrame(nrlist)------------------------------------------------
nrlist <- represAsDataFrame(protrna)
head(nrlist)

## ----applyToPDB(listpdb=nrlist, FUN=hasHetAtm, hetAtms="MG")------------------
## Set progressbar=TRUE to see the progress
nrlist_mg <- applyToPDB(FUN=hasHetAtm, listpdb=nrlist$pdb, 
                        hetAtms="MG", progressbar=FALSE)
nrlist <- cbind(nrlist, Mg=nrlist_mg[, 2])
head(nrlist)

## To see only the structures containing Mg use
head(nrlist[nrlist$Mg == TRUE, ])

## ----pipeNucData(pdblist, chainlist, cores)-----------------------------------
## After the download is finished, this dataset is analysed in less than 2 min
## (single core, intel i5 2.3Ghz). Set progressbar=TRUE to see the progress.
ntinfo <- pipeNucData(pdbID=nrlist$pdb, 
                        model=nrlist$model,
                        chain=nrlist$chain, 
                        progressbar=FALSE, cores=2)
str(ntinfo)

## ----pipeNucData(pdblist, chainlist, cores, path, extension), eval=FALSE------
#  ntinfo <- pipeNucData(pdbID=nrlist$pdb,
#                          model=nrlist$model,
#                          chain=nrlist$chain,
#                          progressbar=FALSE, cores=2,
#                          path="/your/path/to/the/dataset/", extension=".cif.gz")

## ----pipeProtNucData(pdblist, chainlist, cores)-------------------------------
## After the download is finished, this dataset is analysed in less than 1 min
## (single core, intel i5 2.3Ghz). Set progressbar=TRUE to see the progress.
aantinfo <- pipeProtNucData(pdbID=nrlist$pdb, 
                            model=nrlist$model,
                            chain=nrlist$chain, 
                            progressbar=FALSE, cores=2)
str(aantinfo)

## ----dssr(pdbID)--------------------------------------------------------------
## Execute dssr, the wrapper will donwload the (mmCIF) file if necessary
rna <- dssr("1bau")

## The contents of the `rna` object can be seen with
names(rna)

## Then, the contents of the list model inside `rna`
names(rna$model)

## And, inside the sublist parameters
names(rna$models$parameters)

## The number 1 herein selects the model of the structure, and is needed even
## in XRAY structures.
class(rna$models$parameters$pairs[[1]])

## Finally, to get the data about base pairs, use:
rna$models$parameters$pairs[[1]][, c("index", "nt1", "nt2", "LW")]


