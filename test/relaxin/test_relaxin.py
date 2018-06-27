#/usr/bin/env python
"""
Tests a simple protein with many disulfides: disulfides between two identical
resids on different chains, disulfides on the same chain, and unrelated
disulfides on different chains.
"""
import os
from vmd import atomsel, molecule

dir = os.path.dirname(__file__)
topdir = os.path.abspath(os.path.join(dir, ".."))

#==============================================================================
def check_correctness(molid):
    """ Verifies molecule is sane """

    molecule.set_top(molid)

    # Check the protein is there with the correct capping groups
    assert len(atomsel("protein or resname ACE NMA NME")) == 828
    assert len(set(atomsel("protein").get("fragment"))) == 2
    assert len(set(atomsel("resname ACE NMA NME").get("residue"))) == 4

    # Check for 6 cysteines, 2 with same resid
    assert len(set(atomsel("resname CYS CYX").get("residue"))) == 6

    # Check connectivity between cysteines is correct
    for res in set(atomsel("resname CYS CYX").get("residue")):
        assert len(atomsel("residue %d" % res)) == 10
        assert len(atomsel("residue %d and name SG" % res)) == 1
        idxs = atomsel("residue %d and name SG" % res).bonds[0]
        assert set(atomsel("index %s"
                           % " ".join(str(i) for i in idxs)).get("name")) \
            == set(["CB", "SG"])

#==============================================================================

def test_single_chain(tmpdir):

    from psfgen import PsfGen
    p = str(tmpdir.mkdir("single_chain"))
    os.chdir(dir)

    gen = PsfGen(output="/dev/null")
    gen.read_topology(os.path.join(topdir,"top_all36_caps.rtf"))
    gen.read_topology(os.path.join(topdir,"top_all36_prot.rtf"))
    gen.read_topology(os.path.join(topdir,"top_water_ions.rtf"))

    # Read protein
    gen.add_segment(segid="P0", pdbfile="psf_protein_P0.pdb")
    gen.read_coords(segid="P0", filename="psf_protein_P0.pdb")

    gen.add_segment(segid="P1", pdbfile="psf_protein_P1.pdb")
    gen.read_coords(segid="P1", filename="psf_protein_P1.pdb")

    # Read waters, with 10k atoms per file to avoid PDB limitations
    gen.add_segment(segid="W0", pdbfile="psf_wat_0.pdb")
    gen.read_coords(segid="W0", filename="psf_wat_0.pdb")

    gen.add_segment(segid="W1", pdbfile="psf_wat_1.pdb")
    gen.read_coords(segid="W1", filename="psf_wat_1.pdb")

    # Read ions
    gen.add_segment(segid="I", pdbfile="psf_ions.pdb")
    gen.read_coords(segid="I", filename="psf_ions.pdb")

    # Add disulfides
    gen.patch(patchname="DISU", targets=[("P0","10"), ("P0","15")])
    gen.patch(patchname="DISU", targets=[("P0","24"), ("P1","23")])
    gen.patch(patchname="DISU", targets=[("P0","11"), ("P1","11")])

    # Regenerate
    gen.regenerate_angles()
    gen.regenerate_dihedrals()

    # Write
    os.chdir(p)
    gen.write_psf(filename="output.psf")
    gen.write_pdb(filename="output.pdb")

    # Load as a molecule with vmd-python and check it's correct
    m = molecule.load("psf", "output.psf", "pdb", "output.pdb")
    check_correctness(m)
    molecule.delete(m)

#==============================================================================
