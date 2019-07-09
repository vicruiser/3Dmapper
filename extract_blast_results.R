#!/~/Rscript
#indicate library path
.libPaths(c(.libPaths(), "/home/bsc51/bsc51927/R/x86_64-pc-linux-gnu-library"))
#Load necessary packages
library("RFLPtools")
library(stringr)

#######################################################################################
#   extractBLASTResults: transforms jackhmmer's tabular output into csv format.            #
#######################################################################################
#Input: 
#blastOutput                   .txt file containing the names of the BLAST output files (tabular format)
#                              that have been generated per query sequence
# outputFile=output            name of the output file. Defaul = "output"
# outputPath                   path indicating where to write the function output files.
#Output:
# 2 csv files that gather all the tblout and domtblout results of all the queries

extractBLASTResults<-function(InputPath, blastOutput, outputName, outputPath){
  #Create empty data frame. 
  blast_df<-data.frame() 
  #for loop to generate a data frame containing all queries and their respectives homologs.
  for(i in 1:length( blastOutput)){
    print(i)
    fileName <- as.character(file.path(InputPath,blastOutput[i]))
    print(fileName)
    #read_tblout() function from rhmmer package
    readBlast <- try(read.table(file = fileName, sep="\t"))
    head(readBlast)
    colnames(readBlast) <- c("qseqid", "qlen", "sseqid", "slen", "qstart",
                             "qend", "sstart", "send", "evalue", "length", "pident", "nident")
  }
  #write data frame in csv format
  outputName <- sub(".txt", "", basename(fileName))
  write.csv(readBlast, paste(outputPath,"/", outputName,".csv", sep=""))
}





#########################
# Pass input arguments  #
#########################
InputPath <- commandArgs(TRUE)[1]
blastOutput <-commandArgs(TRUE)[2] 
outputName <- commandArgs(TRUE)[3]
outputPath <- commandArgs(TRUE)[4]


#execute function
extractBLASTResults(InputPath,  blastOutput, outputName, outputPath)