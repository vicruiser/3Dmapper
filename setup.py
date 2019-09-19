#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
        name         = 'PDBmapper',
        version      = '1.0',
        author       = 'Victoria Ruiz-Serra',
        author_email = 'victoria.ruizserra@bsc.es',
        # ext_modules  = [pytadbit_module, pytadbit_module_old,
        #                 eqv_rmsd_module, centroid_module,
        #                 consistency_module, aligner3d_module,
        #                 squared_distance_matrix_module],
        # package_dir  = {'pytadbit': PATH + '/_pytadbit'},
        packages     = find_packages(),
        #['pytadbit', 'pytadbit.parsers', 'pytadbit.tools',
        #                 'pytadbit.boundary_aligner', 'pytadbit.utils',
        #                 'pytadbit.tad_clustering', 'pytadbit.modelling',
        #                 'pytadbit.mapping'],
        # py_modules   = ["pytadbit"],
        platforms = "OS Independent",
        license = "GPLv3",
        description  = 'Map rare variants to protein interfaces data in 3D.',
        long_description = (open("README.rst").read() +
                            open("docs/source/install.rst").read()),
        #classifiers  = TAGS,
        #provides     = ["pytadbit"],
        #keywords     = ["testing"],
        url          = 'https://github.com/vicruiser/PDBmapper',
        download_url = 'https://github.com/vicruiser/PDBmapper/tarball/master',
        scripts      = ['scripts/PDBmapper.py', 'scripts/mapping_tools.py',
                        'scripts/VEPcrossref.py', 'scripts/parse_argv.py'],
        #data_files   = [(path.expanduser('~'),
        #                 ['extras/.bash_completion'])]
)