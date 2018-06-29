.. _examples:

Documentation examples
======================

The documentation for TCL-psfgen features a few examples. Here, I have
implemented the same examples using the Python interface.


Preparing separate PDB files
----------------------------
This example uses multiple PDB files, containing water, ions, and two protein
fragments, and composites them into a single PSF file. Several disulfide bonds
are added with a patch. All files for this example can be found in
`test/relaxin` or `test/` in the case of the parameter files.

.. code-block:: python

    from psfgen import Psfgen
    gen = Psfgen()
    gen.read_topology("top_all36_caps.rtf")
    gen.read_topology("top_all36_prot.rtf")
    gen.read_topology("top_water_ions.rtf")

    # Read protein
    gen.add_segment(segid="P0", pdbfile="psf_protein_P0.pdb")
    gen.read_coords(segid="P0", filename="psf_protein_P0.pdb")

    gen.add_segment(segid="P1", pdbfile="psf_protein_P1.pdb")
    gen.read_coords(segid="P1", filename="psf_protein_P1.pdb")

    # Read waters, with 10k atoms per file to avoid PDB limitations
    # (this limitation may be fixed later)
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
    gen.write_psf(filename="relaxin.psf")
    gen.write_pdb(filename="relaxin.pdb")


Coming soon
-----------
More examples coming soon
