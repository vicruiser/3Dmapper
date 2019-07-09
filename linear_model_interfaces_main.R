###############################################################################
# Load necessary packages
###############################################################################
#indicate library path
.libPaths(c(.libPaths(),"/gpfs/home/bsc51/bsc51927/R/x86_64-pc-linux-gnu-library")) 

library(data.table)
library(dplyr)
library(nlme)

###############################################################################
# Input arguments                                                             #      
###############################################################################
interface_id          = as.character(commandArgs(TRUE)[1])
input_dir             = as.character(commandArgs(TRUE)[2])
interfaces_file       = as.character(commandArgs(TRUE)[3])
LF_CT_file            = as.character(commandArgs(TRUE)[4])
TB_file               = as.character(commandArgs(TRUE)[5])
output_dir            = as.character(commandArgs(TRUE)[6])

# interface_id = "5o1h.pdb2_ENSP00000269305_B_GOL"
# input_dir = "/home/vruizser/PhD/2018-2019/Immunity_interfaces_analysis/data_ready_for_analysis"
# interfaces_file = "enriched_interfaces_MC3_data.txt"
# LF_CT_file = "LeuFrac_CancerType_per_patient.txt"
# TB_file = "tumor_burden_per_patient.txt"
# output_dir = "/home/vruizser/PhD/2018-2019/Immunity_interfaces_analysis/results"

###############################################################################
# Load data                                                                   #
###############################################################################
interface_data <- fread(
  cmd = paste("awk -F \",\" 'NR==1; $9 == \"",
              interface_id,
              "\"' ",
              file.path(input_dir,
              interfaces_file),
              sep = ""
              )
  )

LF_CT <- fread(file.path(input_dir, LF_CT_file))
TB <-  fread(file.path(input_dir, TB_file))

###############################################################################
# Change colnames                                                             #
###############################################################################
names(LF_CT)[names(LF_CT)== "Study"] <- "Cancer_Type"
names(LF_CT)[names(LF_CT)== "leukocyte_fraction"] <- "Leukocyte_Fraction"
names(TB)[names(TB)== "total_number_mutations"] <- "Tumor_Mutation_Burden"
names(interface_data)[names(interface_data)== "interfaces_id"] <- "Mutated_Interface"
interface_data = unique(interface_data[,c("ParticipantBarcode", "Mutated_Interface")])

###############################################################################
# Define formula Linear model                                                 #
###############################################################################
lm_formula = formula(Leukocyte_Fraction ~ 1 + Cancer_Type + Tumor_Mutation_Burden + Mutated_Interface)

###############################################################################
# Prepare data (model data matrix)                                            #
###############################################################################
lm_data_matrix <- left_join(LF_CT, interface_data, by ="ParticipantBarcode")
lm_data_matrix <- left_join(lm_data_matrix, TB , by ="ParticipantBarcode")

lm_data_matrix$Mutated_Interface[is.na(lm_data_matrix$Mutated_Interface)] <- "mutationless"
lm_data_matrix$Tumor_Mutation_Burden[is.na(lm_data_matrix$Tumor_Mutation_Burden)] <- 0

lm_data_matrix$Mutated_Interface <- as.factor(lm_data_matrix$Mutated_Interface)
lm_data_matrix$Mutated_Interface <- relevel(lm_data_matrix$Mutated_Interface, ref = "mutationless")

###############################################################################
# Fit LM                                                                      #
###############################################################################
lm <- gls(lm_formula, lm_data_matrix)

###############################################################################
# Retrieve results                                                            #
###############################################################################
estimated_coefficients <- data.table(coefficients(summary(lm)))
estimated_coefficients$Mutated_Interface <- interface_id
rownames(estimated_coefficients) <- rownames(coefficients(summary(lm)))

model_fit_quality_param = data.frame(AIC = AIC(lm),
                                     BIC = BIC(lm),
                                     logLik = lm$logLik,
                                     interface_id = interface_id)
###############################################################################
# Write them into a file                                                      #
###############################################################################
fwrite(estimated_coefficients,
       file = file.path(output_dir, "estimated_coefficients_immunity_interfaces_analysis.txt"), 
       append = TRUE, 
       row.names = TRUE,
       col.names = !(file.exists(file.path(output_dir, "estimated_coefficients_immunity_interfaces_analysis.txt"))))

fwrite(model_fit_quality_param, 
       file = file.path(output_dir, "model_fit_quality_param.txt"), 
       append = TRUE, 
       col.names = !(file.exists(file.path(output_dir, "model_fit_quality_param.txt"))))
