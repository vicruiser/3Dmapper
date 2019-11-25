#!/usr/bin/env python3

from setuptools import setup, find_packages  # pylint: disable=no-name-in-module,import-error
import io

setup(
    name='pdbmapper',
    version='0.1.0',
    author='Victoria Ruiz-Serra',
    author_email='victoria.ruizserra@bsc.es',
    url="https://github.com/vicruiser/PDBmapper/master",
    platforms="OS Independent",
    license="GPLv3",
    description='Map annotated genomic variants to protein interfaces data in 3D.',
    long_description=(open("README.rst").read() +
                      open("docs/source/install.rst").read()),
    # For a list of valid classifiers, see https://pypi.org/classifiers/
    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        #'Intended Audience :: Developers',
        #'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish
        #'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        # These classifiers are *not* checked by 'pip install'. See instead
        # 'python_requires' below.
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    keywords=["pdbmapper", "gene variation", "protein 3D"],

    # When your source code is in a subdirectory under the project root, e.g.
    # `src/`, it is necessary to specify the `package_dir` argument.
    package_dir={'': 'scripts'},  # Optional

    # Specify which Python versions you support. In contrast to the
    # 'Programming Language' classifiers above, 'pip install' will check this
    # and refuse to install the project if the version does not match. If you
    # do not support Python 2, you can simplify this to '>=3.5' or similar, see
    # https://packaging.python.org/guides/distributing-packages-using-setuptools/#python-requires
    python_requires='>=3, <4',

    # This field lists other packages that your project depends on to run.
    # Any package you put here will be installed by pip when your project is
    # installed, so they must be valid existing projects.
    #
    # For an analysis of "install_requires" vs pip's requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['halo', 'vcfpy'],  # Optional

    #package_data={'PDBmapper': ['data/splitted_interfaces_db/*']},


    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # `pip` to create the appropriate form of executable for the target
    # platform.
    #
    # For example, the following would provide a command called `sample` which
    # executes the function `main` from this package when invoked:
    entry_points={
        "console_scripts": ['pdbmapper=scripts.execute_pdbmapper:main']
    },

    # List additional URLs that are relevant to your project as a dict.


    download_url='https://github.com/vicruiser/PDBmapper/tarball/master',
    packages=find_packages(where='scripts'),  # Required,

    project_urls={  # Optional
        'Bug Reports': 'https://github.com/vicruiser/PDBmapper/issues',
        # 'Funding': 'https://donate.pypi.org',
        # 'Say Thanks!': 'http://saythanks.io/to/example',
        'Source': 'https://github.com/vicruiser/PDBmapper/'
    },
    # data_files   = [(path.expanduser('~'),
    #                 ['extras/.bash_completion'])]
)
