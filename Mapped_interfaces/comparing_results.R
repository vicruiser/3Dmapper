library(bio3d)
pdb_file_3lea_pdb1 <- read.pdb("/home/vruizser/PhD/2018-2019/Mapped_interfaces/3lea.pdb1")

prot_chains_3lea <- atom.select(pdb_file_3lea_pdb1, string="protein")
lig_3lea <- atom.select(pdb_file_3lea_pdb1, string="ligand")

prot_residues_3lea <- pdb_file_3lea_pdb1$atom[prot_chains_3lea$atom,]
ligs_3lea <- pdb_file_3lea_pdb1$atom[lig_3lea$atom,]

head(prot_residues_3lea)
ligs_3lea_ZN <- subset(ligs_3lea, resid =="ZN")
ligs_3lea_Z93 <- subset(ligs_3lea, resid =="Z93")

interaction_chainA_ZN_3lea <- dist.xyz(prot_residues_3lea[,c("x","y","z")],
                                       ligs_3lea_ZN[,c("x","y","z")])

interaction_chainA_Z93_3lea <- dist.xyz(prot_residues_3lea[,c("x","y","z")],
                                        ligs_3lea_Z93[,c("x","y","z")])

resid313<- subset(prot_residues_3lea, resno == 313)
resid351<- subset(prot_residues_3lea, resno == 351)
interaction_chainA_Z93_3lea[as.numeric(rownames(resid313)),]
interaction_chainA_Z93_3lea[as.numeric(rownames(resid351)),]


pos_ZN<- which(interaction_chainA_ZN_3lea <= 5, arr.ind = T)
prot_residues_3lea[pos_ZN[,"row"],]



pos_Z93<- which(interaction_chainA_Z93_3lea <= 5, arr.ind = T)
prot_residues_3lea[pos_Z93[,"row"],]
sort(unique(prot_residues_3lea[pos_Z93[,"row"],]$resno))
