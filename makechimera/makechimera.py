#!/usr/bin/python

# standard library imports
import datetime
import os
import pandas as pd
import sys
import time

# third-party packages
from halo import Halo

# local application imports
from makechimera import helpers
from makechimera import parsers
from makechimera import templates
from makechimera import logger


def main():
    """
    Gets the show on the road.
    """ 
    # display CLI aesthetics
    print(templates.description)
    print(templates.epilog)
    
    # initialize spinner decorator
    spinner = Halo(text='Loading', spinner='dots12', color="cyan")
    spinner.start()
    
    # parse arguments and check for errors  
    args = parsers.parse_args()
    helpers.check_args(args)
    
    # set up log file
    outdir = args.output if args.output else os.getcwd()
    log = logger.get_logger('main', outdir)
    
    # set up report
    report = open(os.path.join(outdir, 'makechimera.report'), 'w')
    report.write(templates.description)
    report.write(templates.epilog)
    report.write('''
    Command line input:
    -------------------
    \n''')
    progname = os.path.basename(sys.argv[0])
    report.write(progname + ' ' + " ".join(sys.argv[1:]) + '\n' + '\n' + '\n')
    
    # start timer
    time_format = '[' + time.ctime(time.time()) + ']'
    start = time.time()
    
    # load mapped-to-interface data
    if args.interface_file:
        interface_data = helpers.read_interface_data(args.interface_file)
    
    # load mapped-to-structure data
    if args.structure_file:
        structure_data = helpers.read_structure_data(args.structure_file)
    
    for pdb in args.pdb: 
        
        # set up default values 
        pdb_interface = pdb_structure = pd.DataFrame()
        assemblies_interface = assemblies_structure = []
        
        if args.interface_file:
            # subset interface data to match PDB code
            pdb_interface = helpers.filter_data(pdb, interface_data[interface_data.PDB_code.str.contains(pdb)], args)

            # find all available bioassemblies if one is not provided
            pdb, assemblies_interface = helpers.get_assemblies(pdb, pdb_interface)
        
        if args.structure_file:
            # subset structure data to match PDB code
            pdb_structure = structure_data[structure_data.PDB_code.str.contains(pdb)]
            
            # find all available bioassemblies if one is not provided
            pdb, assemblies_structure = helpers.get_assemblies(pdb, pdb_structure)
            
        # if PDB code is not found on either file, move on to next PDB    
        if pdb_interface.empty and pdb_structure.empty:
            print(f'PDB code "{pdb}" was not found within the data with selected parameters.')
            continue
        
        # create union between assemblies from both files (if provided)
        bioassemblies = list(set().union(assemblies_structure, assemblies_interface))  
             
        
        for bioassembly in bioassemblies:
            
            # set up default values
            asmbl_data_interface = asmbl_data_structure = pd.DataFrame()
            
            # subset interface data to match bioassembly
            if args.interface_file:
                asmbl_data_interface = pdb_interface[pdb_interface.PDB_code.str.contains(bioassembly)]
            
            # subset structure data to match bioassembly
            if args.structure_file:
                asmbl_data_structure = pdb_structure[pdb_structure.PDB_code.str.contains(bioassembly)]
            
            # if subsets generated are empty (bioassembly not in the data), move on to next assembly
            if asmbl_data_interface.empty and asmbl_data_structure.empty:
                continue
            
            # set up main details (model number 1000 has been chosen because the selection in ChimeraX 
            # will always come up empty)
            main_details = {
                'pdb': pdb,
                'asmbl': bioassembly,
                'interfaces': '#1000',
                'itf_variants': '#1000',
                'str_variants': '#1000',
            }
            
            # add interface details
            if args.interface_file:
                # add interfaces
                main_details['interfaces'] = helpers.get_interfaces(asmbl_data_interface)
                
                # add variant locations for interface data
                main_details['itf_variants'] = helpers.get_interface_variants(asmbl_data_interface)
            
            # add structure details
            if args.structure_file:
                # add variant locations for structure data
                main_details['str_variants'] = helpers.get_structure_variants(asmbl_data_structure) 
            
            # add remaining details (ChimeraX options)
            details = helpers.add_details(main_details, args)
            
            # define parts of script name
            filter_it = f'_it-{args.filter_it}' if args.filter_it else ''
            base_name = args.name if args.name else ''
            
            # define script name and path
            script_name = base_name + f'{pdb}.pdb{bioassembly}' + filter_it + '.cxc'
            script_path = os.path.join(args.output, script_name) if args.output else os.path.join(os.getcwd(), script_name)

            # check if file exists and overwrite flag not given
            if not args.overwrite:
                helpers.check_path_available(script_path)
                
            # write ChimeraX script based on template
            helpers.write_script(script_path, templates.script_template, details)
            log.info(f'{script_name} script was written successfully.')
            report.write(f'{time_format} {script_name} script was written successfully.\n')
    
    # end timer
    end = time.time()
    
    # write success statement to report and close report
    report.write(f'{time_format} ChimeraX scripts have been successfully writen to {outdir}. ' + 
                 f'Total time: {str(datetime.timedelta(seconds=round(end-start)))}\n')
    
    report.close()
    
    # stop spinner
    spinner.stop()
    
    # print success statement to command line
    spinner.stop_and_persist(symbol='\U0001F4CD',
        text=' makechimera process finished. Total time: ' + 
        str(datetime.timedelta(seconds=round(end-start))))

    
if __name__ == '__main__':
    main()
