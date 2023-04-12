#' An S4 class to represent a structure parsed from its mmCIF file.
#'
#' All mmCIF files in the PDB at date 2018-Feb-19th contain (just) 14 common
#' attributes, which are represented in the CIF objects herein with the same
#' names as found in mmCIF documentation
#' [http://mmcif.wwpdb.org/](http://mmcif.wwpdb.org/).
#'
#' @slot entry The ID code
#' @slot audit_conform 'mmCIF' dictionary version
#' @slot database_2 Cross-reference ID codes to other databases
#' @slot pdbx_database_status Deposition data
#' @slot audit_author Authors
#' @slot entity Entities (molecules & ions) in the structure
#' @slot chem_comp Residues (ATOM & HETATM) in the structure
#' @slot exptl Experimental technique
#' @slot struct Author description of the structure
#' @slot struct_keywords Author description key words
#' @slot struct_asym Chain-entity equivalences
#' @slot atom_sites Details about the crystallographic cell
#' @slot atom_type Details about the atoms in structure
#' @slot atom_site The atomic coordinates
#'
#' @seealso To create CIF objects use [cifParser()]
#'
#' @author Diego Gallego
#'
CIF <- setClass("CIF",
                slots = list(
                                entry                = "character",
                                audit_conform        = "character",
                                database_2           = "data.frame",
                                pdbx_database_status = "character",
                                audit_author         = "data.frame",
                                entity               = "data.frame",
                                chem_comp            = "data.frame",
                                exptl                = "data.frame", 
                                struct               = "character",
                                struct_keywords      = "character",
                                struct_asym          = "data.frame",
                                atom_sites           = "character",
                                atom_type            = "data.frame",
                                atom_site            = "data.frame")
)

setMethod(
    "show",
    signature="CIF",
    definition=function(object) {
        cat("\n-- mmCIF with ID: ", object@entry,
            " ------------------------------------------------------\n\n",
            sep="")
        description <- paste("Author description: ", object@struct[2], sep="")
        description <- gsub('(.{1,80})(\\s|$)', '\\1\n', description)
        cat(description, "\n")
        cat("mmCIF version:       ", object@audit_conform[2], "\n\n", sep="")
        cat("To extract coordinates and other data use accessor functions",
            "\n(type ?cif_accessors for details)", "\n\n", sep="")
    }
)
