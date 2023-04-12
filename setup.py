"""A setuptools based setup module.
See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""
# Always prefer setuptools over distutils
# from distutils.core import setup, Extension
from setuptools.command.install import install as DistutilsInstall
from setuptools import find_packages
from setuptools import setup  # pylint: disable=no-name-in-module,import-error
#from setuptools.command.build_ext import build_ext
from urllib.request import urlopen
import subprocess
import io
from os import path
import os
import sys

here = path.abspath(path.dirname(__file__))
# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Get the requirements of the packages
with open('dependencies.txt') as f:
    requirements = f.read().splitlines()


class git_clone_external(DistutilsInstall):

    def internet_on(self):
        try:
            urlopen('http://python.org/', timeout=1)
            return True
        except urlopen.URLError as err:
            return False

    # def run_command(self, command, cwd):
    #     # print("about to execute command `{0}`. sources = {1}"
    #     # .format(command, sources))
    #     proc = subprocess.Popen(command,cwd = cwd,  stderr=subprocess.STDOUT, shell=True)
    #     output, stderr = proc.communicate(input)
    #     status = proc.wait()
    #     if status:
    #         raise Exception("command = {0} failed with output = {1} status {2:d}\n"
    #                         .format(command, output, status))

    
    def run(self):

        #connection = self.internet_on()
        #if connection is True:
            d = os.path.dirname(os.path.realpath(__file__))
            bcftools_dir = os.path.dirname(
                os.path.realpath(__file__)) + "/bcftools"
            htslib_dir = os.path.dirname(os.path.realpath(__file__))+" /htslib"
            verina_dir = os.path.dirname(os.path.realpath(__file__))+ "/verina"
            
            if not os.path.exists(htslib_dir):
                cmd1 = ['git', 'clone',
                            'https://github.com/samtools/htslib.git']

                #self.run_command( cmd1, d)
                subprocess.call(cmd1, cwd=os.path.dirname(
                    os.path.realpath(__file__)))

            if not os.path.exists(path.join(bcftools_dir, 'bcftools')):
                cmd2 = ['git', 'clone',
                            'https://github.com/samtools/bcftools.git']

                #self.run_command(self, cmd2, d)
                subprocess.call(
                    cmd2, cwd=os.path.dirname(os.path.realpath(__file__)))
                
            # if not os.path.exists(verina_dir):
            #     cmd3 = ['git','clone',
            #                 'https://mmb.irbbarcelona.org/gitlab/dgallego/veriNA3d-dev.git']

            #     #self.run_command( cmd1, d)
            #     subprocess.call(cmd3, cwd=os.path.dirname(
            #         os.path.realpath(__file__)))

            #self.run_command(self, ["make", "clean"], bcftools_dir)
            #self.run_command(self, ["make"], bcftools_dir)
            subprocess.call(
                ["make", "clean"], cwd=bcftools_dir)

            subprocess.call(['make'], cwd=bcftools_dir)
            #self.run_command(self, ["cp", "bcftools", "%s/bin/" % os.environ.get('VIRTUAL_ENV', '/usr/local/')], bcftools_dir)
            #self.run_command(self, ["cp", "plugins/split-vep.so", "%s/bin/" % os.environ.get('VIRTUAL_ENV', '/usr/local/')], bcftools_dir)
            subprocess.call(
                ["cp", "bcftools", "%s/bin/" % os.environ.get('VIRTUAL_ENV', '/usr/local/')], cwd=bcftools_dir)

            subprocess.call(
                ["cp", "plugins/split-vep.so", "%s/bin/" % os.environ.get('VIRTUAL_ENV', '/usr/local/')], cwd=bcftools_dir)

            # subprocess.call(
            #     ["R","CMD", "build", "%s" % os.environ.get('VIRTUAL_ENV', '/usr/local/') , "--no-build-vignettes"], cwd=verina_dir)

            # subprocess.call(
            #     ["R","CMD", "build", "%s" % os.environ.get('VIRTUAL_ENV', '/usr/local/') , "--no-build-vignettes"], cwd=verina_dir)
            
            #subprocess.call(
            #    ["tar", "-zxf" ,"stride.tar.gz"], cwd=os.path.dirname(
            #        os.path.realpath(__file__)))
            
            #subprocess.call(
            #    ["make"], cwd=stride_dir)
            
            
            DistutilsInstall.run(self)


# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.
setup(

    # This is the name of your project. The first time you publish this
    # package, this name will be registered for you. It will determine how
    # users can install this project, e.g.:
    #
    # $ pip install sampleproject
    #
    # And where it will live on PyPI: https://pypi.org/project/sampleproject/
    #
    # There are some restrictions on what makes a valid project name
    # specification here:
    # https://packaging.python.org/specifications/core-metadata/#name
    name='3dmapper',  # Required
    
    # Versions should comply with PEP 440:
    # https://www.python.org/dev/peps/pep-0440/
    #
    # For a discussion on single-sourcing the version across setup.py and the
    # project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version='1.0.0',  # Required

    # This should be your name or the name of the organization which owns the
    # project.
    author='Victoria Ruiz-Serra',  # Optional

    # This should be a valid email address corresponding to the author listed
    # above.
    author_email='victoria.ruiz.serra@gmail.es',  # Optional

    # This should be a valid link to your project's main homepage.
    #
    # This field corresponds to the "Home-Page" metadata field:
    url="https://github.com/vicruiser/3Dmapper",  # Optional

    # You have previously uploaded your project to your github repository. Now,
    #  we create a new release version of your project on github. This release
    #  will then be downloaded by anyone that runs the “pip install YourPackage” command.
    # download_url='https://github.com/vicruiser/PDBmapper/archive/v_01.tar.gz',
    # This is a one-line description or tagline of what your project does. This
    # corresponds to the "Summary" metadata field:
    description='Map annotated genomic positions to protein interfaces data in 3D.',  # Optional

    # This field corresponds to the "Description" metadata field:
    # https://packaging.python.org/specifications/core-metadata/#description-optional
    long_description=long_description,  # Optional

    # Denotes that our long_description is in Markdown; valid values are
    # text/plain, text/x-rst, and text/markdown
    #
    # Optional if long_description is written in reStructuredText (rst) but
    # required for plain-text or Markdown; if unspecified, "applications should
    # attempt to render [the long_description] as text/x-rst; charset=UTF-8 and
    # fall back to text/plain if it is not valid rst" (see link below)
    #
    # This field corresponds to the "Description-Content-Type" metadata field:
    # https://packaging.python.org/specifications/core-metadata/#description-content-type-optional
    # Optional (see note above)
    long_description_content_type='text/markdown',

    # For a list of valid classifiers, see https://pypi.org/classifiers/
    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        # Pick your license as you wish
        'License :: GPLv3',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        # These classifiers are *not* checked by 'pip install'. See instead
        # 'python_requires' below.
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],

    # This field adds keywords for your project which will appear on the
    # project page. What does your project relate to?
    #
    # Note that this is a string of words separated by whitespace, not a list.
    keywords=["mapping","CLI", "genetic positions", "protein structure", "PDB"],  # Optional

    # When your source code is in a subdirectory under the project root, e.g.
    # `src/`, it is necessary to specify the `package_dir` argument.
    # package_dir={'': 'pdbmapper'},  # Optional

    # You can just specify package directories manually here if your project is
    # simple. Or you can use find_packages().
    #
    # Alternatively, if you just want to distribute a single Python file, use
    # the `py_modules` argument instead as follows, which will expect a file
    # called `my_module.py` to exist:
    #
    #   py_modules=["my_module"],
    #

    packages=['mapper', #'makepsdb',
              'makevariantsdb', 'makeinterfacedb', 'makechimera' ],  # Required!!!!!

    # Specify which Python versions you support. In contrast to the
    # 'Programming Language' classifiers above, 'pip install' will check this
    # and refuse to install the project if the version does not match. If you
    # do not support Python 2, you can simplify this to '>=3.5' or similar, see
    # https://packaging.python.org/guides/distributing-packages-using-setuptools/#python-requires
    python_requires='>=3.6, <4',  # Optional

    # This field lists other packages that your project depends on to run.
    # Any package you put here will be installed by pip when your project is
    # installed, so they must be valid existing projects.
    #
    # For an analysis of "install_requires" vs pip's requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=requirements,
    setup_requires=['wheel'],
    # External requirements
    cmdclass={'install': git_clone_external},
    # ext_modules=[module1],
    # Data requirements
    package_data={
        # And include any *.dat files found in the 'data' subdirectory
        # of the 'mypkg' package, also:
        '': ['data/*','rscripts/*']
    },
    include_package_data = True, 

    # If your project depends on packages that don’t exist on PyPI, you may
    # still be able to depend on them, as long as they are available for
    #  download as:
    #   an egg, in the standard distutils sdist format,
    #   a single .py file, or
    #   a VCS repository(Subversion, Mercurial, or Git).

    # The dependency_links option takes the form of a list of URL strings. For
    # example, this will cause a search of the specified page for eggs or
    # source distributions, if the package’s dependencies aren’t already
    # installed
    # dependency_links=[                              # Optional
    #    'https://github.com/samtools/bcftools@develop#egg=bcftools-1.8.0'
    # ],
    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # `pip` to create the appropriate form of executable for the target
    # platform.
    #
    # For example, the following would provide a command called `sample` which
    # executes the function `main` from this package when invoked:
    entry_points={
        "console_scripts": ['mapper=mapper.__main__:main',
                            #'makepsdb=makepsdb.__main__:main',
                            'makevariantsdb=makevariantsdb.__main__:main',
                            'makevisualization=makechimera.__main__:main',
                           'makestructuraldb=makeinterfacedb.__main__:main']

    },

    # List additional URLs that are relevant to your project as a dict.
    #
    # This field corresponds to the "Project-URL" metadata fields:
    # https://packaging.python.org/specifications/core-metadata/#project-url-multiple-use
    #
    # Examples listed include a pattern for specifying where the package tracks
    # issues, where the source is hosted, where to say thanks to the package
    # maintainers, and where to support the project financially. The key is
    # what's used to render the link text on PyPI.
    project_urls={  # Optional
        'Bug Reports': 'https://github.com/vicruiser/3Dmapper/issues',
        # 'Funding': 'https://donate.pypi.org',
        # 'Say Thanks!': 'http://saythanks.io/to/example',
        'Source': 'https://github.com/vicruiser/3Dmapper/',
    },


)
