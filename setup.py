#!/usr/bin/env python

import sys
from setuptools import setup
from setuptools.extension import Extension
from setuptools.command.test import test as TestCommand

class PyTest(TestCommand):
    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []
    def run(self):
        import pytest
        errno = pytest.main()
        sys.exit(errno)

psfgenfiles=[
   "./src/charmm_file.c",
   "./src/charmm_parse_topo_defs.c",
   "./src/extract_alias.c",
   "./src/hash.c",
   "./src/hasharray.c",
   "./src/memarena.c",
   "./src/pdb_file.c",
   "./src/pdb_file_extract.c",
   "./src/psf_file.c",
   "./src/psf_file_extract.c",
   "./src/python_psfgen.c",
   "./src/stringhash.c",
   "./src/topo_defs.c",
   "./src/topo_mol.c",
   "./src/topo_mol_output.c"]

psfext = Extension('_psfgen',
                   define_macros=[('PSFGENTCLDLL_EXPORTS', '1')],
                   include_dirs=["/usr/include", "./src"],
                   libraries=["tcl8.5"],
                   library_dirs=["/usr/lib"],
                   sources=psfgenfiles
                  )

setup(name="psfgen",
      version="1.7.0",
      description="Protein structure file generator",
      author="Robin Betz, Justin Gullingsrud and Jim Phillips",
      author_email="robin@robinbetz.com",
      url="https://psfgen.robinbetz.com",
      packages=['psfgen'],
      ext_modules=[psfext],
      tests_require=["pytest", "vmd-python"],
      cmdclass = {'test': PyTest}
     )

