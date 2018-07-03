.. _examples:

Examples
========

Here are a few examples of using the Python psfgen interface, that should
cover most typical uses.

These examples form the basis of the psfgen-Python test suite, so all
input topology files and protein structures for these examples can be
found `here <https://github.com/Eigenstate/psfgen/tree/master/test>`_.


Combining separate PDB files
----------------------------
This example uses multiple PDB files, containing water, ions, and two protein
fragments, and composites them into a single PSF file. Several disulfide bonds
are added with a patch.

These proteins were prepared with "capping groups" of ACE and NMA on the
N- and C-termini that are represented as their own residues rather than the
Charmm style of capping patches. Topology information for these residues is
found in `top_all36_caps.rtf`, but this is just a convention used in my lab.

Although not required to parse the test case files, the common task of aliasing
different histidine protonation states to the same topology residue definition
is done.

.. code-block:: python

    from psfgen import PsfGen
    gen = PsfGen(output="/dev/null")  # Suppress output since there's too much
    gen.read_topology("top_all36_caps.rtf")
    gen.read_topology("top_all36_prot.rtf")
    gen.read_topology("top_water_ions.rtf")

    # Set up an alias for histidine protonation states
    gen.alias_residue(top_resname="HIS", pdb_resname="HIE")
    gen.alias_residue(top_resname="HIS", pdb_resname="HID")

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


Mutating a residue
------------------
Using the same proteins as the previous example, we'll load only the first chain
to keep the example simple. Then, we'll mutate Leu:2 to an Alanine and have
psfgen guess coordinates for new atoms.


.. code-block:: python

    from psfgen import PsfGen
    gen = PsfGen()
    gen.read_topology("top_all36_caps.rtf")
    gen.read_topology("top_all36_prot.rtf")

    # Read protein and do the mutation
    gen.add_segment(segid="P", pdbfile="psf_protein_P0.pdb",
                    mutate=[("2","ALA")])

    # Add disulfide bond
    gen.patch(patchname="DISU", targets=[("P","10"), ("P","15")])

    # Guess coordinates for ALA mutation
    gen.guess_coords()

    # Regenerate angles and dihedrals since we added a patch
    gen.regenerate_angles()
    gen.regenerate_dihedrals()

    # Save output
    gen.write_psf(filename="L2A.psf")
    gen.write_pdb(filename="L2A.pdb")


We keep the default output from the PsfGen object, so will see it in stdout (the
terminal window) and will see a few warnings about failing to set coordinates
for LEU:2 when calling `read_coords`. This is because we're trying to read in
coordinates for a residue that doesn't exist, as we have an alanine at that
location, and the warning can therefore be safely ignored.


Adding capping groups
---------------------

Let's say you have a protein without explicitly represented capping groups,
called `protein_nocaps.pdb` (again, this file can be found in the test
directory).

This example shows how the `add_segment` function can be used to
set the patches applied to the N- and C-termini. We'll apply the neutral NTER
patch to the N-terminus, and have a positively charged GLUP patch on the
C-terminus. We'll use the `residues` argument to add an ALA and a GLU to the
end of the protein chain (after resid 24) as well, as there isn't a glutamate
at the end of the input protein.

.. code-block:: python

    from psfgen import PsfGen
    gen = PsfGen()
    gen.read_topology("top_all36_prot.rtf")

    # Read protein, set the patches for the ends, and add ALA-GLU to the end
    gen.add_segment(segid="P", pdbfile="protein_nocaps.pdb",
                    first="NTER", last="GLUP",
                    residues=[("25","ALA"), ("26","GLU")])

    # Read in the coordinates we have and guess the remainder
    gen.read_coords(segid="P", filename="protein_nocaps.pdb")
    gen.guess_coords()

    # Save the result
    gen.write_psf(filename="patch_ends.psf")
    gen.write_pdb(filename="patch_ends.pdb")

You'll see in the output a lot of warnings about poorly guessed coordinates, but
this is a contrived example so this is okay.


Working with velocities
-----------------------

Psfgen can actually read and write NAMD binary files, including velocities.
Let's pretend we have a pre-equilibrated lipid membrane, and a simulation of
protein in solution. We'll take the protein and lipid and combine them while
preserving velocities from simulation, using vmd-python to figure out which
water molecules to delete to insert the membrane.

This is an extremely contrived example designed to show off how to use
:meth:`psfgen.PsfGen.read_psf` and :meth:`psfgen.PsfGen.write_namdbin`, and
isn't an actually recommended way to set up simulations. We'll pretend that
the coordinates of everything are aligned so the files can be combined, too,
and that the lipid will be exactly from z=5 to z=15.

.. code-block:: python

    from psfgen import PsfGen
    from vmd-python import atomsel, molecule

    gen = PsfGen()

    # Psfgen will attempt to load topologies listed in the PSF file, but to be
    # safe and to avoid problems with absolute vs. relative paths we explicitly
    # load them here.
    gen.read_topology("top_all36_prot.rtf")
    gen.read_topology("top_all3_lipid.rtf")

    # Read the protein segment, including coordinates and velocities.
    gen.read_psf(filename="protein_equil.psf",
                 namdbinfile="protein_equil.bin",
                 velnamdbinfile="protein_equil_vel.bin")

    # Get the segment name that was read in, will use for deleting later
    pseg = gen.get_segids()[0]

    # Now the lipid segment
    gen.read_psf(filename="popc_equil.psf",
                 namdbinfile="popc_equil.bin",
                 velnamdbinfile="popc_equil_vel.bin")

    # Load files in vmd-python to figure out which waters to delete
    pid = molecule.load("psf", "protein_equil.psf",
                        "namdbin", "protein_equil.bin")
    to_delete = set(atomsel("water and z > 5 and z < 15").get("resid"))
    delete(pid)

    # Delete overlapping waters in the water segment
    for resid in to_delete:
        gen.delete_atoms(segid=pseg, resid=resid)

    # Write the output files, including velocities
    gen.write_psf(filename="combined.psf")
    gen.write_namdbin(filename="combined.bin",
                      velocity_filename="combined_vel.bin")


