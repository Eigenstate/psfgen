from setuptools import setup
from setuptools.extension import Extension


psfgenfiles = [
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
    "./src/topo_mol_output.c",
]

psfext = Extension(
    "_psfgen",
    define_macros=[("PSFGENTCLDLL_EXPORTS", "1")],
    include_dirs=["./src"],
    libraries=[],
    library_dirs=[],
    sources=psfgenfiles,
)

setup(
    packages=["psfgen"],
    ext_modules=[psfext],
    include_package_data=True,
)
