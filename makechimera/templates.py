#!/usr/bin/python

"""
Contains:
- ChimeraX script templates with placeholders to facilitate string formatting.
- CLI aesthetics.
"""

# CHIMERAX TEMPLATES
# ==================

script_template = """\

# Chimera X script to visualize the variants mapped to the interfaces of {pdb_code}.pdb{assembly} in 3D.


# NOTE: to manually customize the details of the script, refer to the sections labeled with the OPTIONS keyword and modify
#     the line of code below accordingly. the placeholder '#1000' is used when data pertaining to either mapped-to-interface 
#     or mapped-to-structure variants is missing because it will always point out to an inexistent PDB model, which allowed 
#     us to simplify the script template.


# PDB code selection and assembly
# =====================================================

# select protein by PDB code
open {pdb_code}

# select biological assembly
sym #1 assembly {assembly}

# include only chosen biological assembly
close #1,3-end



# General Settings
# =====================================================

# set lighting (OPTIONS: full, soft, simple)
lighting {lighting}

# set background (OPTIONS: white, black)
set bg {bg}

# set up silhoutte outline (OPTIONS: true, false)
set silhouettes {silhouettes}

# set molecule style (OPTIONS: ball, sphere, stick)
style {mol_style}

# display ribbons if not already displayed (OPTIONS: show, hide)
show ribbons

# hide atoms if not already hidden (OPTIONS: show, hide)
hide atoms



# Protein Complex Customization
# =====================================================

# color protein chains
color bychain



# Interface Customization
# =====================================================

# display atoms from the affected interface and show them in molecule style (OPTIONS: sphere, ball, stick)
select {interfaces}; show sel atoms; style sel {int_style}; ~select



# Variant Customization (interface)
# =====================================================

# color variant residues mapped to interfaces (OPTIONS: https://www.rbvi.ucsf.edu/chimerax/docs/user/commands/colornames.html)
color {variants_itf} red

# label varying residues mapped to interfaces
label {variants_itf} height 1.5 size 68



# Variant Customization (structure)
# =====================================================

# color variant residues mapped to structure (OPTIONS: https://www.rbvi.ucsf.edu/chimerax/docs/user/commands/colornames.html)
color {variants_str} red

# label varying residue mapped to structure
label {variants_str} height 1.5 size 68

"""

test_template = """\
pdb_code: {pdb_code}
assembly: {assembly}
main_chain: {main_chain}
second_chain: {second_chain}
main_resid: {main_resid}
second_resid: {second_resid}
var_pos: {var_pos}
"""



# CLI AESTHETICS
# ==============

description = '''\n
    ----------------------------------------- Welcome to ----------------------------------------------
    $$$$$$$\  $$$$$$$\  $$$$$$$\  
    $$  __$$\ $$  __$$\ $$  __$$\  
    $$ |  $$ |$$ |  $$ |$$ |  $$ |$$$$$$\$$$$\   $$$$$$\   $$$$$$\   $$$$$$\   $$$$$$\   $$$$$$\  
    $$$$$$$  |$$ |  $$ |$$$$$$$\ |$$  _$$  _$$\  \____$$\ $$  __$$\ $$  __$$\ $$  __$$\ $$  __$$\  
    $$  ____/ $$ |  $$ |$$  __$$\ $$ / $$ / $$ | $$$$$$$ |$$ /  $$ |$$ /  $$ |$$$$$$$$ |$$ |  \__|
    $$ |      $$ |  $$ |$$ |  $$ |$$ | $$ | $$ |$$  __$$ |$$ |  $$ |$$ |  $$ |$$   ____|$$ |
    $$ |      $$$$$$$  |$$$$$$$  |$$ | $$ | $$ |\$$$$$$$ |$$$$$$$  |$$$$$$$  |\$$$$$$$\ $$ |
    \__|      \_______/ \_______/ \__| \__| \__| \_______|$$  ____/ $$  ____/  \_______|\__|
                                                            $$ |      $$ |
                                                            $$ |      $$ |
                                                            \__|      \__|
          ---------------  Map genomic variants to protein data in 3D. -----------------
    \n'''
    
epilog = \
        '''
          ------------------------------------------------------------------------------------
         |  Copyright (c) 2019 Victoria Ruiz --                                               |
         |  victoria.ruizserrra@bsc.es -- https://www.bsc.es/ruiz-serra-victoria-isabel       |
          ------------------------------------------------------------------------------------
        '''
        


    






    

