.. _all_commands:
.. currentmodule:: all_commands
.. highlight:: bash

TCL <=> Python equivalents
==========================

Unlike TCL, Python is much more naturally object-oriented. This means that all
psfgen functionality is performed on a PsfGen object that holds the state of
the system. For example:

.. code-block:: python

    from psfgen import PsfGen
    gen = PsfGen()

Then all commands are done on the PsfGen object. The syntax in the following
table refers to this generic PsfGen object I've chosen to name ``gen``.

Topology and naming functions
-----------------------------
.. list-table::
   :header-rows: 1
   :widths: 5 40 40

   * - Task
     - TCL command
     - Python class equivalent
   * - Load topology definitions
     - ``topology <file name>``
     - ``gen.read_topology(filename)`` :meth:`psfgen.PsfGen.read_topology`
   * - Provide alternate names for residues in topology file
     - ``topology alias <desired residue name> <topology residue name>``
     - ``gen.alias_residue(top_resname, pdb_resname)`` :meth:`psfgen.PsfGen.alias_residue`
   * - Provide alternate names for residues in pdb file
     - ``pdbalias residue <PDB residue name> <desired residue name>``
     - ``gen.alias_residue(top_resname, pdb_resname)`` :meth:`psfgen.PsfGen.alias_residue`
   * - Provide translations for atom names found in PDB files to those
       in topology files
     - ``pdbalias atom <residue name> <PDB atomname> <topology atomname>``
     - ``gen.alias_atom(resname, pdb_atomname, top_atomname)`` :meth:`psfgen.PsfGen.alias_atom`

Query functions
---------------
.. list-table::
   :header-rows: 1

   * - Task
     - TCL command
     - Python class equivalent
   * - List loaded topology files
     - ``topology list``
     - 
   * - List current segids
     - ``segment segids``
     - ``gen.get_segids()`` :meth:`psfgen.PsfGen.get_segids`
   * - List current resids in a segment
     - ``segment resids``
     - ``gen.get_resids(segid)`` :meth:`psfgen.PsfGen.get_resids`
   * - Get residue name given a resid and a segment
     - ``segment residue <segment ID> <resid>``
     - ``gen.get_resname(segid, resid)`` :meth:`psfgen.PsfGen.get_resname`
   * - List all currently applied patches, including default patches
     - ``patch listall``
     - 
   * - List all explicitly applied patches
     - ``patch list``
     - 
   * - Get the name of the patch applied to the beginning of a given segment
     - ``segment first <segment ID>``
     - ``gen.get_first(segid)`` :meth:`psfgen.PsfGen.get_first`
   * - Get the name of the patch applied to the end of a given segment
     - ``segment last <segment ID>``
     - ``gen.get_last(segid)`` :meth:`psfgen.PsfGen.get_last`
   * - Get a list of atoms given a resid and segment
     - ``segment atoms <segment ID> <resid>``
     - 
   * - Get x,y,z coordinates for a given atom
     - ``segment coordinates <segment ID> <resid>``
     - 

System building functions
-------------------------
.. list-table::
   :header-rows: 1

   * - Task
     - TCL command
     - Python class equivalent
   * - Guess coordinates of atoms that aren't explicitly set
     - ``guesscoord``
     -  
   * - Set coordinates of a given atom
     - ``coord <segment ID> <resid> <atomname> {x y z}``
     - ``set_coord(segid, resid, atomname, position)`` :meth:`psfgen.PsfGen.set_coord`
   * - Delete all atoms in a segment
     - ``delatom <segment ID>``
     - 
   * - Delete all atoms in a residue
     - ``delatom <segment ID> <resid>``
     - 
   * - Delete a single atom
     - ``delatom <segment ID> <resid> <atom name>``
     - 
   * - Create multiple images of a set of atoms for locally enhanced sampling
     - ``multiply <factor> <segid[:resid[:atomname]]> ...``
     - Not implemented
   * - Remove insertion codes and modify resids minimially for uniqueness
     - ``regenerate resids``
     - ``gen.regenerate_resids()`` :meth:`psfgen.PsfGen.regenerate_resids`
   * - Remove angles and regenerate them, after patching
     - ``regenerate angles``
     - ``gen.regenerate_angles()`` :meth:`psfgen.PsfGen.regenerate_angles`
   * - Remove dihedrals and regenerate them, after patching.
     - ``regenerate dihedrals``
     - ``gen.regenerate_dihedrals()`` :meth:`psfgen.PsfGen.regenerate_dihedrals`
   * - Apply a patch to one or more residues, as determined by the patch
     - ``patch <patchname> <segid:resid> [...]``
     - ``gen.patch(patch_name, targets)`` :meth:`psfgen.PsfGen.patch`

Adding a new segment
--------------------

All of the following TCL code to create a new segment is equivalent to the
following Python code:

.. code-block:: tcl

    segment <segment ID> {                    # Create a new segment with given ID
        first <patchname>                     # Set patch applied to first residue
        last  <patchname>                     # Set patch applied to last residue
        pdb <pdbfilename>                     # Extract sequence info from given PDB file
        auto [angles|dihedrals|none]          # Automatically generate these

        residue <resid> <resname> [<chain>]   # Add a residue to the end of segment
        mutate <resid> <resname>              # Mutate a residue in current segment to a new one
    }


All of this block is equivalent to the following function call in Python,
with the `segid` argument being mandatory, and all others optional.

.. code-block:: python

    gen.add_segment(segid=<segment ID>,
                    first=None,
                    last=None,
                    auto=None,
                    residue=None,
                    mutate=None,
                   )


I/O functions
-------------
.. list-table:: I/O functions
   :header-rows: 1

   * - Task
     - TCL command
     - Python class equivalent
   * - Write out structure to a PSF file
     - ``writepsf [charmm] [x-plor] [cmap|nocmap] <filename>``
     - ``gen.write_psf(filename, type)`` :meth:`psfgen.PsfGen.write_psf` 
   * - Write out structure to a PDB file
     - ``writepdb <filename>``
     - ``gen.write_pdb(filename)`` :meth:`psfgen.PsfGen.write_pdb`
   * - Write out a NAMD binary input file
     - ``writenamdbin <filename>``
     - 
   * - Read in a PSF file and add it to the current structure. Optionally
       read coordinates from a pdb file or a namdbin file.
     - ``readpsf <filename> [pdb] <pdbfilename>``
     - 
   * - Read coordinates in from a PDB file, matching segment, residue, and atom
       names
     - ``coordpdb <filename> [segid]``
     - ``gen.read_coords(filename, segid)`` :meth:`psfgen.PsfGen.read_coords`

Context functions
-----------------
The TCL version of psfgen has the concept of a "psf context", that is, a given
molecule or system that is being built. You can use the TCL methods in the
table below to generate new contexts or switch between them to build multiple
molecules at once.

However since Python is object oriented, all you need to do to have multiple
systems in memory is create multiple PsfGen objects. Each one is self-contained.

.. list-table::
   :header-rows: 1

   * - Task
     - TCL command
     - Python class equivalent
   * - Create a new context, with its own structure, topology definitions,
       and aliases, and set it to active
     - ``psfcontext new``
     - ``gen = new PsfGen()``
   * - Create a new context, but do not switch to it.
     - ``psfcontext create``
     - ``gen = new PsfGen()``
   * - Delete a psfcontext
     - ``psfcontext delete <context>``
     - ``del gen``
   * - Switch to a diferent psfcontext
     - ``psfcontext <id>``
     - Use multiple PsfGen objects
   * - Make context case sensitive (by default it is not)
     - ``psfcontext mixedcase``
     - 
   * - Make context case insensitive (default setting)
     - ``psfcontext allcaps``
     - 
   * - Clear the structure, topology definitions, and aliases
     - ``psfcontext reset``
     - Delete and make a new PsfGen object
   * - Evaluate commands in a given context, returning to the current one when
       done
     - ``psfcontext eval <context> { <commands> }``
     - Use multiple PsfGen objects
   * - Get total number of contexts created and destroyed so far
     - ``psfcontext stats``
     - Use multiple PsfGen objects
