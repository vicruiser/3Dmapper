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
ENSP_vs_pident = unique(unique_interfaces_overlap[, c("ENSP","pident")])
ENSP_vs_pident = subset(ENSP_vs_pident, !is.na(pident))
ENSP_vs_pident = ENSP_vs_pident[order(ENSP_vs_pident$pident, decreasing = TRUE),]  

# same for gene ids
gene_vs_pident = unique(unique_interfaces_overlap[, c("Hugo_Symbol","pident")])
gene_vs_pident = subset(gene_vs_pident, !is.na(pident))
gene_vs_pident = gene_vs_pident[order(gene_vs_pident$pident, decreasing = TRUE),] 

# same for number of interfacial residues
int_vs_pident = unique(unique_interfaces_overlap[, c("ENSP","pident", "interfaces_ids", "length_interface", "length.ali")])
int_vs_pident = subset(int_vs_pident, !is.na(pident))
int_vs_pident = int_vs_pident[order(int_vs_pident$pident, decreasing = TRUE),] 

# target sequence coverage
coverage_vs_pident = unique(unique_interfaces_overlap[, c("ENSP","pident", "interfaces_ids", "temp.length", "length.ali")])
coverage_vs_pident = subset(coverage_vs_pident, !is.na(pident))
coverage_vs_pident = coverage_vs_pident[order(coverage_vs_pident$pident, decreasing = TRUE),] 

# split into list of data frames
ENSP_vs_pident_list = split(ENSP_vs_pident, f = ENSP_vs_pident$pident)
gene_vs_pident_list = split(gene_vs_pident, f = gene_vs_pident$pident)
coverage_vs_pident_list = split(coverage_vs_pident, f = coverage_vs_pident$pident)
int_vs_pident_list = split(int_vs_pident, f = int_vs_pident$pident)

cum_sum_ENSP =  list()
cum_sum_gene = list()
cum_sum_coverage = list()
cum_sum_int = list()

for(i in 1:length(ENSP_vs_pident_list)){
  print(i)
  cum_sum_ENSP[[i]] = length(unique(rbindlist(ENSP_vs_pident_list [1:i])$ENSP))
  cum_sum_gene[[i]] = length(unique(rbindlist(gene_vs_pident_list [1:i])$gene))
  
  coverage_vs_pident_df = unique(rbindlist(coverage_vs_pident_list [1:i]))
  cum_sum_coverage[[i]] = (sum(coverage_vs_pident_df$length.ali )/ sum(coverage_vs_pident_df$temp.length))*100
  
  int_vs_pident_df = unique(rbindlist(int_vs_pident_list [1:i]))
  cum_sum_int[[i]] = (sum(int_vs_pident_df$length_interface)/ sum(coverage_vs_pident_df$length.ali))*100
  
}


cum_sum_ENSP = unlist(cum_sum_ENSP)
cum_sum_gene = unlist(cum_sum_gene)

cumsum_plot = rbind(data.frame(var = cum_sum_ENSP,  pident= sort(as.numeric(names(ENSP_vs_pident_list)), decreasing = TRUE) , id ="ENSP" ),
                    data.frame(var = cum_sum_gene,  pident= sort(as.numeric(names(ENSP_vs_pident_list)), decreasing = TRUE), id="geneID"))

library("wesanderson")
pal <- wes_palette("Zissou1", 1, type = "continuous")

p = ggplot(df, aes(x=pident, y=ENSP)) + geom_bar(stat = "identity", fill=rep (pal, 6))+
  theme_minimal() #+   labels()

ggplot(cumsum_plot , aes(fill=id, y=var, x=pident)) + 
  geom_line()+
  geom_point

####################################################
# Plots. target-sequence coverage - Nali/Ntot *100 
####################################################
plot_df = unique(unique_interfaces_overlap[, c("Hugo_Symbol","pident")])
plot_df= subset(plot_df, !is.na(pident))
plot_df$val = 1
# setorder(plot_df, by = "pident")
# plot_df.b = plot_df[!duplicated(plot_df$ENSP),]
plot_df$pident = as.numeric(plot_df$pident)

pdf= plot_df %>%
  # group_by(MemID) %>% 
  summarise(c50= length(unique(Hugo_Symbol[pident >=50])),
            c60= length(unique(Hugo_Symbol[pident >=60])),
            c70= length(unique(Hugo_Symbol[pident >=70])),
            c80= length(unique(Hugo_Symbol[pident >=80])),
            c90= length(unique(Hugo_Symbol[pident >=90])),
            c100= length(unique(Hugo_Symbol[pident >=100])))

df = data.frame(pident = names(pdf),
                Hugo_Symbol = t(pdf)[,1])
df$pident = factor(df$pident, levels = c("c50","c60","c70","c80","c90","c100"))

library("wesanderson")
pal <- wes_palette("Zissou1", 1, type = "continuous")

p = ggplot(df, aes(x=pident, y=ENSP)) + geom_bar(stat = "identity", fill=rep (pal, 6))+
  theme_minimal() #+   labels()
p


####################################################
# Plots. Coverage by interfaces
####################################################
  plot_df2 = unique(unique_interfaces_overlap[, c("ENSP","temp.chain","pdb.id","pident")])
plot_df2$pident = as.numeric(plot_df2$pident)
p2 = ggplot(plot_df2, aes(x=pident)) +
  geom_histogram(binwidth=10, color="black", fill="white")
p2

plot_df %>%
  group_by(ENSP, pident) %>%
  summarise(counts = n())

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


df_plot3 = unique(unique_interfaces_overlap[, c("ENSP",
                                                "temp.chain",
                                                "pdb.id",
                                                "pident",
                                                "temp.length",
                                                "length.ali")])
df_plot3 = subset(df_plot3, !is.na(temp.length))

pdf3= df_plot3 %>%
  # group_by(MemID) %>% 
  summarise(c50= (sum(length.ali[pident >=50]) / sum(temp.length[pident >=50]))*100,
            c60= (sum(length.ali[pident >=60]) / sum(temp.length[pident >=60]))*100,
            c70= (sum(length.ali[pident >=70]) / sum(temp.length[pident >=70]))*100,
            c80= (sum(length.ali[pident >=80]) / sum(temp.length[pident >=80]))*100,
            c90= (sum(length.ali[pident >=90]) / sum(temp.length[pident >=90]))*100,
            c100= (sum(length.ali[pident >=100]) / sum(temp.length[pident >=100]))*100)

df3 = data.frame(pident = names(pdf3),
                Coverage = t(pdf3)[,1])
df3$pident = factor(df3$pident, levels = c("c50","c60","c70","c80","c90","c100"))

pal <- wes_palette("Zissou1", 1, type = "continuous")

p = ggplot(df3, aes(x=pident, y=Coverage)) + geom_bar(stat = "identity", fill=rep (pal, 6))+
  theme_minimal() #+   labels()
p

####################################################
# Plots. Coverage by alignment
####################################################
df_plot4 = unique(unique_interfaces_overlap[, c("ENSP",
                                                "temp.chain",
                                                "int.chain",
                                                "pdb.id",
                                                "pident",
                                                "interaction",
                                                "temp.length",
                                                "length.ali",
                                                "length_interface")])
pdf3= df_plot3 %>%
  # group_by(MemID) %>% 
  summarise(c50= (sum(length.ali[pident >=50]) / sum(temp.length[pident >=50]))*100,
            c60= (sum(length.ali[pident >=60]) / sum(temp.length[pident >=60]))*100,
            c70= (sum(length.ali[pident >=70]) / sum(temp.length[pident >=70]))*100,
            c80= (sum(length.ali[pident >=80]) / sum(temp.length[pident >=80]))*100,
            c90= (sum(length.ali[pident >=90]) / sum(temp.length[pident >=90]))*100,
            c100= (sum(length.ali[pident >=100]) / sum(temp.length[pident >=100]))*100)
