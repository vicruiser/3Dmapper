#load necessary packages
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
# load file with appris isoforms ids
appris_isoforms = fread("/home/vruizser/PhD/2018-2019/Immunity_interfaces_analysis/raw_data/APPRIS_isoforms_ensembl_ids.txt")

# change its colnames
colnames(appris_isoforms) = c("APPRIS_annotation", "ENSP", "Hugo_Symbol")

##
unique_interfaces_overlap = fread( "/home/vruizser/PhD/2018-2019/Immunity_interfaces_analysis/raw_data/mpi_v94_principal_isoforms_overlap_labels.txt")
###############################################################################
# select only principal isoform from APPRIS
filtered_appris_isoforms = subset(appris_isoforms, APPRIS_annotation == "principal1" | APPRIS_annotation == "principal2" )

unique_interfaces_overlap

####################################################
# Coverage ENSP/ Genes - percentage of identity
####################################################
####################################################
# Plots. Coverage by alignment
####################################################
# We have a protein sequence (ENSP)
# that it is aligned against a PDB sequence. To calculate the
# coverage of the interface, we calculate it, relative
# to the length of the alignment 
# E.g.: 
#   xxxxxxxxxxxxxxxIIIxxxxIxxxxIIIxxxxxx ENSP
#   xx---------xxxxxxxxx-xx---xxxxxxxxxx PDB
# 
# length of ENSP = 35
# length of PDB = 23
# length of ali = 23
# interface length = 7
# coverage of interface = 7/23 ???

ENSP_vs_pident = unique(unique_interfaces_overlap[, c("ENSP","pident")])
ENSP_vs_pident = subset(ENSP_vs_pident, !is.na(pident))
ENSP_vs_pident = ENSP_vs_pident[order(ENSP_vs_pident$pident, decreasing = TRUE),]  

# same for gene ids
gene_vs_pident = unique(unique_interfaces_overlap[, c("Hugo_Symbol","pident")])
gene_vs_pident = subset(gene_vs_pident, !is.na(pident))
gene_vs_pident = gene_vs_pident[order(gene_vs_pident$pident, decreasing = TRUE),] 

# same for number of interfacial residues
 int_vs_pident = unique(unique_interfaces_overlap[, c("interfaces_ids", "mapped.real.pos", "length.ali")])
# # int_vs_pident = int_vs_pident  %>% 
# #   mutate(mapped.real.pos = strsplit(as.character(mapped.real.pos), "\\|")) %>% 
# #   unnest(mapped.real.pos)
# # int_vs_pident= unique(int_vs_pident)
# int_vs_pident = subset(int_vs_pident, !is.na(pident))
# int_vs_pident = int_vs_pident[order(int_vs_pident$pident, decreasing = TRUE),] 

# target sequence coverage
coverage_vs_pident = unique(unique_interfaces_overlap[, c("ENSP","pident", "temp.length", "length.ali")])
coverage_vs_pident = subset(coverage_vs_pident, !is.na(pident))
coverage_vs_pident = coverage_vs_pident[order(coverage_vs_pident$pident, decreasing = TRUE),] 

# split into list of data frames
ENSP_vs_pident_list = split(ENSP_vs_pident, f = ENSP_vs_pident$pident)
gene_vs_pident_list = split(gene_vs_pident, f = gene_vs_pident$pident)
coverage_vs_pident_list = split(coverage_vs_pident, f = coverage_vs_pident$pident)
#int_vs_pident_list = split(int_vs_pident, f = int_vs_pident$pident)

cum_sum_ENSP =  list()
cum_sum_gene = list()
cum_sum_coverage = list()
#cum_sum_int = list()

for(i in 1:length(ENSP_vs_pident_list)){
  print(i)
  cum_sum_ENSP[[i]] = length(unique(rbindlist(ENSP_vs_pident_list [1:i])$ENSP))
  cum_sum_gene[[i]] = length(unique(rbindlist(gene_vs_pident_list [1:i])$Hugo_Symbol))
  
  coverage_vs_pident_df = rbindlist(coverage_vs_pident_list [1:i])
  cum_sum_coverage[[i]] = (sum(coverage_vs_pident_df$length.ali )/ sum(coverage_vs_pident_df$temp.length))*100
  
  # int_vs_pident_df = rbindlist(int_vs_pident_list [1:i])
  # ddply(int_vs_pident_df, .(ENSP,id),summarize,sum=sum(amount),number=length(id))
  # cum_sum_int[[i]] = (sum(int_vs_pident_df$length_interface)/ sum(coverage_vs_pident_df$length.ali))*100
  
}


cum_sum_ENSP = unlist(cum_sum_ENSP)
cum_sum_gene = unlist(cum_sum_gene)
cum_sum_coverage = unlist(cum_sum_coverage)
#cum_sum_int = unlist(cum_sum_int)

cumsum_plot = rbind(data.frame(var = cum_sum_ENSP,
                               pident = sort(as.numeric(names(ENSP_vs_pident_list)), decreasing = TRUE),
                               id ="ENSP" ),
                    data.frame(var = cum_sum_gene,
                               pident = sort(as.numeric(names(ENSP_vs_pident_list)),
                                            decreasing = TRUE),
                               id = "geneID"))

cumsum_coverage_plot = rbind(data.frame(coverage_percent = cum_sum_coverage,
                               pident = sort(as.numeric(names(ENSP_vs_pident_list)), decreasing = TRUE)
                               ))

library("wesanderson")
pal <- wes_palette("Zissou1", 1, type = "continuous")


ggplot(cumsum_plot , aes(fill=id, y=var, x=pident, color=id)) + 
  geom_line()

ggplot(cumsum_coverage_plot , aes( y=coverage_percent, x=pident)) + 
  geom_line()

plot(sort(as.numeric(names(ENSP_vs_pident_list)), decreasing = TRUE), cum_sum_coverage)


####################################################
# Size distribution of
####################################################
  plot_df2 = unique(unique_interfaces_overlap[, c("ENSP","temp.chain","pdb.id","pident")])
plot_df2$pident = as.numeric(plot_df2$pident)
p2 = ggplot(plot_df2, aes(x=pident)) +
  geom_histogram(binwidth=10, color="black", fill="white")
p2

plot_df %>%
  group_by(ENSP, pident) %>%
  summarise(counts = n())

