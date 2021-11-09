#!/usr/bin/python

# standard library imports
import argparse
import os
import pandas as pd
import re
import sys
  

def add_details(main_details: dict, args) -> dict:
    """
    Adds all remaining details for ChimeraX options provided through the command
    line arguments.

    Args:
        main_details: details pertaining to the mapped files.
        args: arguments passed to this tool through the terminal.

    Returns:
        dict: details necessary to fill the script_template (keys match the 
            template's placeholders).
    """           
    details = {
        'pdb_code': main_details['pdb'],
        'assembly': main_details['asmbl'], 
        'lighting': args.lighting,
        'bg': args.bg,
        'silhouettes': 'false' if args.sil else 'true',
        'mol_style': args.mol_style, 
        'interfaces': ''.join(main_details['interfaces']), 
        'int_style': args.itf_style,     
        'variants_itf': ''.join(main_details['itf_variants']),
        'variants_str': ''.join(main_details['str_variants'])
        }    
    
    return details


def check_args(args: argparse.Namespace):
    """
    Checks the arguments passed to the tool through the command line for 
    potential issues. Raises appropriate error if an issue is found.

    Args:
        args (argparse.Namespace): arguments obtained through the argparse 
            package standard process.
    """
    # PDB CODE
    # ========
    
    # if PDBs are passed as a list within a file, read it and overwrite args.pdb
    if args.pdb_list:
        try:
            with open(args.pdb[0], 'r') as f:
                args.pdb = f.read().strip().split()
        except IOError:
            sys.exit("Error: path to file containing the PDBs is incorrect or file does not exist.")
    
    # check PDB code matches standard format with or without bioassembly
    for pdb in args.pdb:
        if not re.match(r'^[a-zA-Z0-9_]{4}\.pdb[0-9]+$', pdb):
            if not re.match(r'^[a-zA-Z0-9_]{4}$', pdb):
                sys.exit(f'Error: the PDB code "{pdb}" does not follow standard format (If you are trying to pass a file containing PDB codes to the tool, add the parameter --pdb_list).')
       
        
    # MAPPED FILES
    # ============
    
    # check that at least one of the two mapped files is provided
    if all(a is None for a in (args.interface_file, args.structure_file)):
        sys.exit('Error: file containing either mapped-to-interface or mapped-to-structure variants must be provided (or both).')
    
    # Mapped-to-interface file
    # ==============================
    if args.interface_file:
        
        # check that mapped-to-interface file exists
        if not os.path.isfile(args.interface_file):
            sys.exit('Error: path to mapped-to-interface file provided is incorrect or file does not exist.')
            
        # check that mapped-to-interface file contains necessary columns
        necess_cols = ["PDB_code", "PDB_chain", "PDB_3D_position", 
                    "Interaction_type", "PDB_interacting_chain", 
                    "Chimera_3D_position"] 
        
        with open(args.interface_file, 'r') as f:
            mapped = pd.read_csv(f)
            
        if not pd.Series(necess_cols).isin(mapped.columns).all():
            sys.exit('Error: mapped-to-interface file does not contain necessary columns to write the script.')
    
    # Mapped-to-structure file
    # ========================
    if args.structure_file:
        
        # check that mapped-to-structure file exists
        if not os.path.isfile(args.structure_file):
            sys.exit('Error: path to mapped-to-structure file provided is incorrect or file does not exist.')
            
        # check that mapped-to-structure file contains necessary columns
        necess_cols = ["PDB_code", "PDB_chain", "PDB_3D_position"] 
        
        with open(args.structure_file, 'r') as f:
            mapped = pd.read_csv(f)
        
        if not pd.Series(necess_cols).isin(mapped.columns).all():
            sys.exit('Error: mapped-to-structure file does not contain necessary columns to write the script.')
    
    
    # OUTPUT FOLDER
    # =============

    if args.output:
        
        # check that output folder exists
        if not os.path.isdir(args.output):
            sys.exit('Error: path to output folder provided is incorrect or folder does not exist.')
            
        # check that user has permission to write (and execute) in output folder
        if not os.access(args.output, os.W_OK | os.X_OK):
            sys.exit('Error: you do not have permission to write or execute files in the output folder.')
    
    
    # BASE FILE NAME
    # ==============
    
    if args.name:
        if len(args.name) > 200:
            sys.exit('Error: length of base name for output file exceeds 200 characters.')
        
    # FILTERS
    # =======
    
    if args.filter_it:
        
        # check that interaction type provided is one of the available options
        if not args.filter_it in ('protein', 'ligand', 'nucleic'):
            sys.exit('Error: interaction type provided is not one of the available options: protein or ligand.')
    
    
    # CHIMERA OPTIONS
    # ===============
    
    # check that lighting option provided is one of the available options
    if not args.lighting in ('full', 'soft', 'simple'):
        sys.exit('Error: lighting option provided is not one of the available options: full, soft or simple.')

    # check that background option provided is one of the available options
    if not args.bg in ('white', 'black'):
        sys.exit('Error: background option provided is not one of the available options: white or black.')

    # check that molecule style option provided is one of the available options
    if not args.mol_style in ('ball', 'sphere', 'stick'):
        sys.exit('Error: "molecule style" global option provided is not one of the available options: ball, sphere or stick.')
    
    # check that lighting option provided is one of the available options
    if not args.itf_style in ('ball', 'sphere', 'stick'):
        sys.exit('Error: "molecule style of protein interface" option provided is not one of the available options: ball, sphere or stick.')
    

def check_path_available(path: str):
    """
    Checks if file in the specified path already exists. If it does, it 
    displays an error message and exits the program.

    Args:
        path (str): path to file to be checked.
    """    
    if os.path.isfile(path):
        sys.exit("Error: ChimeraX script already exists. To overwrite the existing script add the option '--force'")


def get_assemblies(pdb: str, data: pd.DataFrame) -> tuple:
    """
    Checks if bioassembly is given with PDB code and if not, finds all mapped 
    assemblies and returns them.

    Args:
        pdb (str): PDB code.
        data (pd.DataFrame): mapped file.

    Returns:
        tuple: PDB code and its available bioassemblies.
    """
    if '.pdb' in pdb:
        return pdb.split('.pdb')
    
    bioassemblies = sorted([code[8:] for code in data.PDB_code.unique().tolist()])
    bioassemblies = [b.split('.')[0] for b in bioassemblies]
    return (pdb[:4], bioassemblies)
  

def get_interfaces(data: pd.DataFrame) -> list:
    """
    Returns all unique interfaces from mapped-to-interface data. 

    Args:
        data (pd.DataFrame): variants mapped to interfaces.

    Returns:
        list: unique interfaces.
    """
    interfaces = []
    chains = data.PDB_chain.unique().tolist()
    
    for chain in chains:
        coords = data[data.PDB_chain == chain].Chimera_3D_position.unique().tolist()
        merged_coords = sorted(list(set('-'.join([coord for coord in coords]).split('-'))))
        interfaces.append(f'/{chain}:{",".join(merged_coords)}')
    
    return interfaces
    

def get_interface_variants(data: pd.DataFrame) -> list:
    """
    Returns all unique variants from the mapped-to-interface data.

    Args:
        data (pd.DataFrame): variants mapped to interfaces.

    Returns:
        list: unique variant locations.
    """
    variants = []
    chains = data.PDB_chain.unique().tolist()
        
    for chain in chains:        
        positions = sorted(data[data.PDB_chain == chain].PDB_3D_position.unique().tolist())
        variants.append(f'/{chain}:{",".join(str(pos) for pos in positions)}') 
        
    return variants


def get_structure_variants(data: pd.DataFrame) -> list:
    """
    Returns all unique variants from the mapped-to-structure data.

    Args:
        data (pd.DataFrame): variants mapped to protein structure.

    Returns:
        list: unique variant locations.
    """
    variants = []
    chains = data.PDB_chain.unique().tolist()
    
    for chain in chains:
        positions = sorted(data[data.PDB_chain == chain].PDB_3D_position.unique().tolist())
        variants.append(f'/{chain}:{",".join(str(pos) for pos in positions)}')
    
    return variants

        
def filter_data(pdb: str, data: pd.DataFrame, args) -> pd.DataFrame:
    """
    Filters dataframe based on filters provided through the arguments passed
    to the tool.

    Args:
        pdb (str): PDB code.
        data (pd.DataFrame): mapped data.
        args: arguments passed to this tool through the terminal.

    Returns:
        pd.DataFrame: filtered mapped data.
    """  
    subset = data  
        
    if args.filter_it:
        subset = data[data['Interaction_type'] == args.filter_it]
        if subset.empty:
            print(f'PDB "{pdb}" does not have any interactions of type {args.filter_it}.')
    
    return subset


def read_interface_data(map_file: str, cols: list = None) -> pd.DataFrame:
    """
    Reads a CSV file with the variants mapped to interfaces from 3Dmapper results.

    Args:
        map_file (str): path to file with variants mapped to interfaces.
        cols (list, optional): list of columns to be selected from the file. 
            Defaults to None.

    Returns:
        pd.DataFrame: dataframe containing only the columns with the necessary 
            details to produce the ChimeraX scripts, unless a list of columns 
            is specified.
    """ 
    if not cols:
        cols = ["PDB_code", "PDB_chain", "PDB_3D_position", "Interaction_type",
                "PDB_interacting_chain", "Chimera_3D_position"]  
    
    with open(map_file, 'r') as f:
        data = pd.read_csv(f)
        return data[cols]
    
    
def read_structure_data(map_file: str, cols: list = None) -> pd.DataFrame:
    """
    Reads a CSV file with the variants mapped to structure from 3Dmapper results.

    Args:
        map_file (str): path to file with variants mapped to structure.
        cols (list, optional): list of columns to be selected from the file. 
            Defaults to None.

    Returns:
        pd.DataFrame: dataframe containing only the columns with the necessary 
            details to produce the ChimeraX scripts, unless a list of columns 
            is specified.
    """ 
    if not cols:
        cols = ["PDB_code", "PDB_chain", "PDB_3D_position"] 
    
    with open(map_file, 'r') as f:
        data = pd.read_csv(f)
        return data[cols] 
    
    
def write_script(file: str, template: str, details: dict):  
    """
    Outputs contents of template filled with details into given file.

    Args:
        file (str): path to file that will contain output (if does not exist, 
            will be created)
        template (str): CXC template with placeholders for string formatting.
        details (dict): dictionary containing values for the placeholders of 
            the template.
    """ 
    with open(file, "w") as f:
        f.write(template.format(**details))
    
            
    
    
    
