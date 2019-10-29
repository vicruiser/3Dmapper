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
        package_dir  = {'pdbmapper': 'pdbmapper'}, #PATH + '/_pytadbit'},
        packages     = find_packages(),
        #['pytadbit', 'pytadbit.parsers', 'pytadbit.tools',
        #                 'pytadbit.boundary_aligner', 'pytadbit.utils',
        #                 'pytadbit.tad_clustering', 'pytadbit.modelling',
        #                 'pytadbit.mapping'],
        # py_modules   = ["pytadbit"],
        platforms = "OS Independent",
        license = "GPLv3",
        description  = 'Map annotated genomic variants to protein interfaces data in 3D.',
        long_description = (open("README.rst").read() +
                            open("docs/source/install.rst").read()),
        #classifiers  = TAGS,
        #provides     = ["pytadbit"],
        #keywords     = ["testing"],
        url          = 'https://github.com/vicruiser/PDBmapper',
        download_url = 'https://github.com/vicruiser/PDBmapper/tarball/master',
        scripts      = ['scripts/PDBmapper.py',
                        'scripts/mapping_tools.py',
                        'scripts/VEPcrossref.py',
                        'scripts/parse_argv.py'],
        #data_files   = [(path.expanduser('~'),
        #                 ['extras/.bash_completion'])]
)


import io
from setuptools import setup, find_packages  # pylint: disable=no-name-in-module,import-error


def dependencies(file):
    with open(file) as f:
        return f.read().splitlines()


with io.open("README.md", encoding='utf-8') as infile:
    long_description = infile.read()

setup(
    name='halo',
    packages=find_packages(exclude=('tests', 'examples')),
    version='0.0.28',
    license='MIT',
    description='Beautiful terminal spinners in Python',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='Manraj Singh',
    author_email='manrajsinghgrover@gmail.com',
    url='https://github.com/manrajgrover/halo',
    keywords=[
        "console",
        "loading",
        "indicator",
        "progress",
        "cli",
        "spinner",
        "spinners",
        "terminal",
        "term",
        "busy",
        "wait",
        "idle"
    ],
    install_requires=dependencies('requirements.txt'),
    tests_require=dependencies('requirements-dev.txt'),
    include_package_data=True,
    extras_require={
        'ipython': [
            'IPython==5.7.0',
            'ipywidgets==7.1.0',
        ]
    }
)