#######################################################################################
#  diagrams of results comparisons
#######################################################################################
# Load necessary packages
library(data.table)
library(eulerr)
library(plyr)
library(dplyr)
library(stringr)

# load summary of new and reference mapped interfaces datasets 
new_mapped_interfaces <-
  fread("/home/vruizser/PhD/2018-2019/Mapped_interfaces/mpi_v89_5_vs_v94.csv")

# ref_mapped_interfaces<- fread("/home/vruizser/PhD/2018-2019/Eduard/interfaces_mapped_to_v89_5_reformatted.txt")
# inter_domains <- which(str_detect(ref_mapped_interfaces$chain.1, substr(ref_mapped_interfaces$pdb.id.x, 1, 4)))
# ref_mapped_interfaces <- ref_mapped_interfaces[-inter_domains,]
# write.table(ref_mapped_interfaces[,-1],
#             "/home/vruizser/PhD/2018-2019/Eduard/interfaces_mapped_to_v89_5_reformatted_without_interdomains.txt",
#             quote = FALSE,
#             row.names = T,
#             col.names = T)
ref_mapped_interfaces <- fread("/home/vruizser/PhD/2018-2019/Eduard/interfaces_mapped_to_v89_5_reformatted_without_interdomains.txt", drop=1)

# number of total biounit ids and biounits in new mapped interfaces
new_mapped_interfaces_pdbid_biounits <- unique(new_mapped_interfaces$id)
(length(new_mapped_interfaces_pdbid_biounits)) # 97163 biounits

new_pdb_biounits = new_mapped_interfaces_pdbid_biounits[which(str_detect(new_mapped_interfaces_pdbid_biounits, "\\.pdb"))] # 95212
new_cif_biounits = new_mapped_interfaces_pdbid_biounits[which(str_detect(new_mapped_interfaces_pdbid_biounits, "assembly"))] # 999
(length(new_pdb_biounits)) + length(new_cif_biounits) #should be 97163

# number of total pdb ids and biounits in new mapped interfaces
new_mapped_interfaces_pdbid <- str_remove(new_mapped_interfaces_pdbid_biounits, "-assembly*")
new_mapped_interfaces_pdbid <- unique(str_remove(new_mapped_interfaces_pdbid, "\\.pdb[0-9]*"))
length(new_mapped_interfaces_pdbid) # 65707 pdb files

# number of total pdb ids and biounits in referemce mapped interfaces
ref_mapped_interfaces_pdbid_biounits<- unique(ref_mapped_interfaces$pdb.id.x)
(length(ref_mapped_interfaces_pdbid_biounits)) #71440 (without interdomain)

no<- ref_mapped_interfaces_pdbid_biounits [which(! ref_mapped_interfaces_pdbid_biounits %in% new_mapped_interfaces_pdbid_biounits)]
#ref_mapped_interfaces <- ref_mapped_interfaces_pdbid_biounits[-which(! ref_mapped_interfaces_pdbid_biounits %in% new_mapped_interfaces_pdbid_biounits),]
# number of total pdb ids and biounits in ref mapped interfaces
ref_mapped_interfaces_pdbid <- unique(str_remove(ref_mapped_interfaces_pdbid_biounits, "\\.pdb[0-9]*"))
length(ref_mapped_interfaces_pdbid) # 47954 biounits  47604


# summarize all data for  diagrams
# here only the new pdb ids are taken into account. Those pdb that present results in ref but not in new
# are not considered but (I HAVE TO CALCULATE THAT AND ADD IT TO THE TABLE CAUSE IT SHOULD BE EASY)
# yes in the pdb ids  diagram
unique_new_mapped_interfaces <- unique(new_mapped_interfaces)
comparison_ref_new_table <- ddply(unique_new_mapped_interfaces, .(data.id), summarize,
	chain = sum(chain),
	chain.1 = sum(chain.1),
	ensembl.prot.id = sum(ensembl.prot.id),
	sstart = sum(sstart),
	send = sum(send),
	mapped.real.pos = sum(mapped.real.pos)
)


#  Diagrams
fit_chain <- euler(setNames(comparison_ref_new_table[1:3,2], comparison_ref_new_table[1:3,1]))
fit_chain.1 <- euler(setNames(comparison_ref_new_table[1:3,3], comparison_ref_new_table[1:3,1]))
fit_ensembl.prot.id <- euler(setNames(comparison_ref_new_table[1:3,4], comparison_ref_new_table[1:3,1]))
fit_sstart <- euler(setNames(comparison_ref_new_table[1:3,5], comparison_ref_new_table[1:3,1]))
fit_send <- euler(setNames(comparison_ref_new_table[1:3,6], comparison_ref_new_table[1:3,1]))
fit_real.pos <- euler(setNames(comparison_ref_new_table[1:3,7], comparison_ref_new_table[1:3,1]))

#Biounit
(intersect_biounit <- length(which(ref_mapped_interfaces_pdbid_biounits %in% new_mapped_interfaces_pdbid_biounits )))
(diff_ref_biounit <- length(which(! ref_mapped_interfaces_pdbid_biounits %in% new_mapped_interfaces_pdbid_biounits )))
(diff_new_biounit <- length(which(!new_mapped_interfaces_pdbid_biounits %in% ref_mapped_interfaces_pdbid_biounits )))

comparison_biounit <- c(diff_new_biounit, diff_ref_biounit, intersect_biounit)
fit_pdbid_biounit <- euler(setNames(comparison_biounit,comparison_ref_new_table[1:3,1]))

# PDB id
(intersect_pdbid <- length(which(ref_mapped_interfaces_pdbid %in% new_mapped_interfaces_pdbid )))
(diff_ref_pdbid <- length(which(! ref_mapped_interfaces_pdbid %in% new_mapped_interfaces_pdbid )))
(diff_new_pdbid <- length(which(! new_mapped_interfaces_pdbid %in% ref_mapped_interfaces_pdbid )))

comparison_pdbids<- c(diff_new_pdbid,diff_ref_pdbid, intersect_pdbid)
fit_pdbid <- euler(setNames(comparison_pdbids,comparison_ref_new_table[1:3,1]))

###############################################################################
# PLOTS
###############################################################################
library(VennDiagram)

grid.newpage()
venn.plot_chain <-
		draw.pairwise.venn(
				area1           = (comparison_ref_new_table[1,2] + comparison_ref_new_table[3,2]),
				area2           = comparison_ref_new_table[2,2] + comparison_ref_new_table[3,2],
				cross.area      = comparison_ref_new_table[3,2],
				fill = c("red", "darkblue"),
				alpha = c(0.6, 0.6),
				cex = 1.4,
				cat.fontface = 2,
				main = "Template chain",
				cat.just =list(c(1,1) , c(0,0))

		)

grid.newpage()
venn.plot_chain.1 <-
		draw.pairwise.venn(
				area1           = (comparison_ref_new_table[1,3] + comparison_ref_new_table[3,3]),
				area2           = comparison_ref_new_table[2,3] + comparison_ref_new_table[3,3],
				cross.area      = comparison_ref_new_table[3,3],
				fill = c("red", "darkblue"),
				alpha = c(0.6, 0.6),
				cex = 1.4,
				cat.fontface = 2,
				main = "Template chain",
				cat.just =list(c(1,1) , c(0,0))
		)

grid.newpage()
venn.plot_ensembl.id <-
		draw.pairwise.venn(
				area1           = (comparison_ref_new_table[1,4] + comparison_ref_new_table[3,4]),
				area2           = comparison_ref_new_table[2,4] + comparison_ref_new_table[3,4],
				cross.area      = comparison_ref_new_table[3,4],
				fill = c("red", "darkblue"),
				alpha = c(0.6, 0.6),
				cex = 1.4,
				cat.fontface = 2,
				main = "Template chain",
				cat.just =list(c(1,1) , c(0,0))
		)

grid.newpage()
venn.plot_sstart <-
		draw.pairwise.venn(
				area1           = (comparison_ref_new_table[1,5] + comparison_ref_new_table[3,5]),
				area2           = comparison_ref_new_table[2,5] + comparison_ref_new_table[3,5],
				cross.area      = comparison_ref_new_table[3,5],
				fill = c("red", "darkblue"),
				alpha = c(0.6, 0.6),
				cex = 1.4,
				cat.fontface = 2,
				main = "Template chain",
				cat.just =list(c(1,1) , c(0,0))
		)

grid.newpage()
venn.plot_send <-
		draw.pairwise.venn(
				area1           = (comparison_ref_new_table[1,6] + comparison_ref_new_table[3,6]),
				area2           = comparison_ref_new_table[2,6] + comparison_ref_new_table[3,6],
				cross.area      = comparison_ref_new_table[3,6],
				fill = c("red", "darkblue"),
				alpha = c(0.6, 0.6),
				cex = 1.4,
				cat.fontface = 2,
				main = "Template chain",
				cat.just =list(c(1,1) , c(0,0))
		)

grid.newpage()
venn.plot_mapped_pos <-
		draw.pairwise.venn(
				area1           = (comparison_ref_new_table[1,7] + comparison_ref_new_table[3,7]),
				area2           = comparison_ref_new_table[2,7] + comparison_ref_new_table[3,7],
				cross.area      = comparison_ref_new_table[3,7],
				fill = c("red", "darkblue"),
				alpha = c(0.6, 0.6),
				cex = 1.4,
				cat.fontface = 2,
				main = "Template chain",
				cat.just =list(c(1,1) , c(0,0))
		)

grid.newpage()
venn.plot_pdbid <-
		draw.pairwise.venn(
				area1           = diff_new_pdbid +  intersect_pdbid,
				area2           = diff_ref_pdbid + intersect_pdbid,
				cross.area      = intersect_pdbid,
				fill = c("red", "darkblue"),
				alpha = c(0.6, 0.6),
				cex = 1.4,
				cat.fontface = 2,
				main = "Template chain",
				cat.just =list(c(1,1) , c(0,0))
		)

grid.newpage()
venn.plot_biounit <-
		draw.pairwise.venn(
				area1           = diff_new_biounit +  intersect_biounit,
				area2           = diff_ref_biounit + intersect_biounit,
				cross.area      = intersect_biounit,
				fill = c("red", "darkblue"),
				alpha = c(0.6, 0.6),
				cex = 1.4,
				cat.fontface = 2,
				main = "Template chain",
				cat.just =list(c(1,1) , c(0,0))
		)

new_mapped_interfacesnew_mapped_interfaces
require(gridExtra)
library(cowplot)

cols <- c("#FF6666", "#6666B9")
lg <- legendGrob(labels=c("New (2019)","Ref (2017)"), pch=rep(19,length(c("New (2019)","Ref (2017)"))),
																	gp=gpar(col=cols,  fill="gray"),
																	nrow=1,
																	byrow=FALSE)

png("/home/vruizser/Dropbox/Interfaces_mapping/task4_Map_prot_interfaces/results/comparison_mapped_interfaces_results.png",
				width = 1000,
				height = 1800)

cowplot::plot_grid(
		gTree(children=venn.plot_pdbid),
		NULL,
		gTree(children=venn.plot_biounit),
		gTree(children=venn.plot_mapped_pos),
		NULL,
		gTree(children=venn.plot_ensembl.id),
		gTree(children= venn.plot_chain),
		NULL,
		gTree(children=venn.plot_chain.1),
		gTree(children=venn.plot_sstart),
		NULL,
		gTree(children=venn.plot_send),
		NULL,
		lg,
		labels=c("PDB IDs","", "Biounit IDs",
											"Mapped residues", "", "Ensembl ID",
											"Template chain","", "Interacting chain",
											"Alignment start position","", "Alignment end position"),
		rel_widths = c(1, 0.3, 1),
		rel_heights = c(1, 1, 1, 1, .1),
		ncol = 3

)
dev.off()
###############################################################################
# Plots with UpSetR                                                           #
###############################################################################

# library(UpSetR)
# mutations <- read.csv( system.file("extdata", "mutations.csv", package = "UpSetR"), header=T, sep = ",")
# upset(mutations, sets = c("PTEN", "TP53", "EGFR", "PIK3R1", "RB1"), sets.bar.color = "#56B4E9",
#       order.by = "freq", empty.intersections = "on")
# 
# row.names(comparison_ref_new_table)<- comparison_ref_new_table[,1]
# j<- t(comparison_ref_new_table[,-1])
#  upset(comparison_ref_new_table, sets = c("chain", "chain.1"), sets.bar.color = "#56B4E9",
#        order.by = "freq", empty.intersections = "on")

###############################################################################
# overlapping percentage per structure                                        #
###############################################################################
nmpi = fread("/home/vruizser/PhD/2018-2019/Mapped_interfaces/mpi_v89_5_vs_v94.csv", stringsAsFactors = T)
#nmpi_bystruct <- split(nmpi, f=nmpi$id)
#overlapping_percentage <- lapply(nmpi_bystruct, function(x) sum(x$mapped.real.pos))

total_nmpi =  unique(nmpi) %>%
     group_by(id)%>%
     summarize(total = sum(as.numeric(mapped.real.pos)))%>%
     filter(id != "3lo3.pdb1") %>%
     filter(id != "6bgn.pdb1")

overlap <- subset(unique(nmpi), data.id == "Reference(2017)&New(2019)")


nmpi_overlap = merge(total_nmpi, overlap, by ="id")
nmpi_overlap$percentage<- (nmpi_overlap$mapped.real.pos/nmpi_overlap$total) *100

#head(nmpi_overlap)
#nrow(subset(nmpi_overlap, percentage > 50))


library(ggplot2)
ggplot(nmpi_overlap, aes(x=percentage))+
  geom_histogram(color="darkblue", fill="lightblue")+
  labs(title="Ref vs. New overlapping percentage",x="Overlap percentage", y = "Number of biounits")


# filtering out no overlapping
nmpi_overlap_filtered<- subset(nmpi_overlap, percentage > 0)

library(ggplot2)
ggplot(nmpi_overlap_filtered, aes(x=percentage))+
  geom_histogram(color="darkblue", fill="lightblue")+
  labs(title="Ref vs. New overlapping percentage (f1)",x="Overlap percentage", y = "Number of biounits")

# filter so Comparison there is at least one residue. 
total_nmpi =  unique(nmpi) %>%
  group_by(id)%>%
  summarize(total = sum(as.numeric(mapped.real.pos)))%>%
  filter(id != "3lo3.pdb1") %>%
  filter(id != "6bgn.pdb1")

ggplot(nmpi_overlap, aes(x=percentage))+
  geom_histogram(color="darkblue", fill="lightblue")+
  labs(title="Ref vs. New overlapping percentage",x="Overlap percentage", y = "Number of biounits")

head(subset(nmpi_overlap, percentage < 10))

########################################################3
# different
#########################################################
 all_new_mpi <- fread("/home/vruizser/PhD/2018-2019/Mapped_interfaces/mpi_v94_all.csv",
                      verbose = F,
                      fill=TRUE)
all_new_mpi <- subset(all_new_mpi, interaction != "interaction")

#subset new results to overlapping
overlaping_ensemblID_new <-
  all_new_mpi[which( all_new_mpi$ensembl.prot.id %in% ref_mapped_interfaces$ensembl.prot.id ),]

overlaping_ensemblID_and_biounit_new <-
  unique(overlaping_ensemblID_new[which(overlaping_ensemblID_new$pdb.id %in% ref_mapped_interfaces$pdb.id.x),])

# subset old results to overlapping
overlaping_ensemblID_ref <- 
  ref_mapped_interfaces[which(ref_mapped_interfaces$ensembl.prot.id %in% overlaping_ensemblID_and_biounit_new$ensembl.prot.id), ]

overlaping_ensemblID_and_biounit_ref <- unique(overlaping_ensemblID_ref[which(overlaping_ensemblID_ref$pdb.id.x %in% overlaping_ensemblID_and_biounit_new$pdb.id), ])

which(! overlaping_ensemblID_and_biounit_new$ensembl.prot.id %in% overlaping_ensemblID_and_biounit_ref$ensembl.prot.id)

library(tidyr)
subset_overlaping_ensemblID <- separate_rows(overlaping_ensemblID, mapped.real.pos, sep="-")
