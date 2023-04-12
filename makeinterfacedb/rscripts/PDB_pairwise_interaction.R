#! /usr/bin/Rscript
################################################################################
options(echo = FALSE, verbose = F,warn = -1) 
#' Title:  INTERACTION PAIRS: select the option for pairs
#'  between interaction to be find
#'
#' @param pdb_file
#' @param comparison
#' @param atom_select
#'
#' @return
#' @export
#'
#' @examples
#'
PDB_pairwise_interaction <-
  function(pdb_file,
           comparison,
           atom_select) {
    # Select information of all the atoms of the pdb file
    all_chains <- pdb_file$atom
    
    # Select position of protein chains in PDB file
    pdb_protChains <-
    combine.select(
      combine.select(
        atom.select(pdb_file, ATOMS_INTERACTION, verbose = F),
        atom.select(pdb_file, "protein", verbose = F),
        verbose =F
      ), atom.select(pdb_file, resid= "UNK"),
      operator = "OR", verbose = F)
    # Select protein chains
    protein_chains <- all_chains[pdb_protChains$atom , ]
    if (atom_select == "calpha") {
      protein_chains = subset(protein_chains , elety == "CA")
    }
    if (atom_select == "cbeta") {
      df = rbind(
        subset(protein_chains, resid == "GLY" & elety == "CA"),
        subset(protein_chains , elety == "CB")
      )
      protein_chains = df[order(df$eleno), ]
    }

    # Add residue index for the mapping process:
    # Every chain sequence will start by 1 and if any residue coordinates are missing
    # and as a result the amino acid is not in the protein sequence (only in SEQRES)
    # then the counting index will skip this amino acid as well. I.e., as if it were
    # a continous protein sequence.
    # The same process is followed for nucleic and ligand chains
    protein_chains <-
      split(protein_chains, f = protein_chains$chain)
    protein_chains <-
      lapply(protein_chains, function(x) {
        respos <-
          unique(x[, c("resno", "insert")])  # take into account insertions of residues
        df <-
          data.frame(respos, real.pos = 1:nrow(respos))
        left_join(x, df, by = c("resno", "insert"))
      })
    protein_chains <- ldply(protein_chains, data.frame)
    
    # Count number of protein chains
    n_prot_chains <- length(unique(protein_chains$chain))

    # Select the type of interaction between pairs: prot-prot, prot-dna, prot-ligand
    ### PROTEIN - PROTEIN
    if (comparison == "protein") {
      # Stop if number of chains is not enoguh
      if (n_prot_chains == 0) {
        stop("This PDB file does not contain any protein structure.")
        
      } else if (n_prot_chains == 1) {
        stop(
          "This PDB file does contain only 1 protein chain.
          There are no possible inter-atomic interactions."
        )
        
      } else {
        # Print out number of protein chains
        #print(paste("This PDB file contains", n_prot_chains, "protein chains"))
        # Return two sets containing all protein chains
        return(list(SetChains1 = protein_chains, SetChains2 = protein_chains))
      }
    }
    
    ### PROTEIN -DNA
    else if (comparison == "nucleic") {
      # Select position in PDB file of nucleic chains
      pdb_nucleicChains <- atom.select(pdb_file, string = "nucleic", verbose = F)
      # If there are nucleic chains in the PDB file
      if (length(pdb_nucleicChains$atom) > 0) {
        # Select nucleic chains
        nucleic_chains <- all_chains[pdb_nucleicChains$atom , ]
        # Count and print out number of protein and DNA chains
        n_nucl_chains <- length(unique(nucleic_chains$chain))
        #print(
        #  paste(
        #    "This PDB file contains",
        #    n_prot_chains,
        #    "protein chains and",
        #    n_nucl_chains,
        #    "nucleic chains."
        #  )
        #)
        # Return two sets of chains, one proteins and nucleics the other
        return(list(SetChains1 = protein_chains, SetChains2 = nucleic_chains))
        
      } # STOP if no nucleic chains are found
      else {
        stop("No DNA molecule is present in this PDB file.")
      }
    }
    
    ### PROTEIN - LIGAND
    else if (comparison == "ligand") {
      # Select position in PDB file of ligand chains
      pdb_ligandChains <- atom.select(pdb_file, string = "ligand", verbose = F)
      # If there are ligand chains in the PDB file
      if (length(pdb_ligandChains$atom) > 0) {
        # Select ligand chains
        ligand_chains <- all_chains[pdb_ligandChains$atom , ]
        # Count and print out number of protein and DNA chains
        n_ligand_chains <- length(unique(ligand_chains$chain))
        # print(
        #   paste(
        #     "This PDB file contains",
        #     n_prot_chains,
        #     "protein chains and",
        #     n_ligand_chains,
        #     "ligand chains."
        #   )
        # )
        # Return two sets of chains, one containing proteins and ligands the other
        return(list(SetChains1 = protein_chains, SetChains2 = ligand_chains))
      }# STOP if no nucleic chains are found
      else {
        stop("No ligands are present in this PDB file.")
      }
    }# STOP if wrong type of comparison
    else{
      stop(
        "Wrong type of comparison. Please select one comparison:
        'protein', 'nucleic' or 'ligand'."
      )
    }
}