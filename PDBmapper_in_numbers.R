##
unique_interfaces_overlap = fread( "/home/vruizser/PhD/2018-2019/Immunity_interfaces_analysis/raw_data/mpi_v94_principal_isoforms_overlap_labels.txt")
###############################################################################
# select only principal isoform from APPRIS
filtered_appris_isoforms = subset(appris_isoforms, APPRIS_annotation == "principal1" | APPRIS_annotation == "principal2" )

# load input file of all estimated coefficients for 3D regions
input_file = "/home/vruizser/PhD/2018-2019/Immunity_interfaces_analysis/git/ICR/results/lm_3D_ICR/estimated_coeffs_lm_Mut3DLoc_ICR.txt"
estimated_coeffs_lm3D <- fread(input_file, fill = T)

# set the corresponding column names
colnames(estimated_coeffs_lm3D) = c("Covariate","Intercept", "Std.Error", "t.value", "p.value", "Interface.ID")

# subset estimated values of the covariate of interest (mutation 3D location)
estimated_coeffs_Mut3DLoc <- subset(estimated_coeffs_lm3D, str_detect(Covariate, "Mutation_3D_LocationInterface"))

# set pvalue column into correct format (numeric)
estimated_coeffs_Mut3DLoc$p.value <- as.numeric(estimated_coeffs_Mut3DLoc$p.value)

# extract ENSP id from feature id for later analysis
estimated_coeffs_Mut3DLoc$ENSP = str_split_fixed(estimated_coeffs_Mut3DLoc$Interface.ID, "_", 3)[,2]

# combine APPRIS information and the one from Ensembl
estimated_coeffs_Mut3DLoc_isoforms = right_join(estimated_coeffs_Mut3DLoc, filtered_appris_isoforms, by = "ENSP")
estimated_coeffs_Mut3DLoc_isoforms = subset(estimated_coeffs_Mut3DLoc, !is.na(Covariate))

fwrite(estimated_coeffs_Mut3DLoc_isoforms,
       "/home/vruizser/PhD/2018-2019/Immunity_interfaces_analysis/git/ICR/results/lm_3D_ICR/estimated_coeffs_lm_Mut3DLoc_ICR_isoforms.txt")


# 

overlap_file =  "/home/vruizser/PhD/2018-2019/Immunity_interfaces_analysis/data_ready_for_analysis/MC3_overlap_interfaces.txt"
overlaps =  fread(overlap_file)

###########################################33
unique_interfaces_overlap
library(ggplot2)

plot_df = unique(unique_interfaces_overlap[, c("ENSP","pident")])
plot_df= subset(plot_df, !is.na(pident))
plot_df$val = 1
# setorder(plot_df, by = "pident")
# plot_df.b = plot_df[!duplicated(plot_df$ENSP),]
plot_df$pident = as.numeric(plot_df$pident)

pdf= plot_df %>%
  # group_by(MemID) %>% 
  summarise(c50= length(unique(ENSP[pident >=50])),
            c60= length(unique(ENSP[pident >=60])),
            c70= length(unique(ENSP[pident >=70])),
            c80= length(unique(ENSP[pident >=80])),
            c90= length(unique(ENSP[pident >=90])),
            c100= length(unique(ENSP[pident >=100])))
plot(t(pdf))
#plot_df.b$pident = as.numeric(plot_df.b$pident)
p = ggplot(plot_df, aes(x=pident, y=cumsum(val))) + geom_line() + geom_point()
p

appris = 
  
  plot_df2 = unique(unique_interfaces_overlap[, c("ENSP","temp.chain","pdb.id","pident")])
plot_df2$pident = as.numeric(plot_df2$pident)
p2 = ggplot(plot_df2, aes(x=pident)) +
  geom_histogram(binwidth=1, color="black", fill="white")
p2

plot_df %>%
  group_by(ENSP, pident) %>%
  summarise(counts = n())


####################################################
# Plots
####################################################
