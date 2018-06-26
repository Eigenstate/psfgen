.. _all_commands:
.. currentmodule:: all_commands
.. highlight:: bash

TCL <=> Python equivalents
==========================

Unlike TCL, Python is much more naturally object-oriented. This means that all
psfgen functionality is performed on a Psfgen object that holds the state of
the system. For example:

.. code-block:: python

    from psfgen import Psfgen
    gen = Psfgen()

Then all commands are done on the Psfgen object. The syntax in the following
table refers to this generic Psfgen object I've chosen to name ``gen``.

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
     - ``gen.read_topology(filename)``
   * - Provide alternate names for residues in topology file
     - ``topology alias <desired residue name> <topology residue name>``
     - ``gen.alias_residue(top_resname, pdb_resname)`` TODO CHECK
   * - Provide alternate names for residues in pdb file
     - ``pdbalias residue <PDB residue name> <desired residue name>``
     - ``gen.alias_residue(top_resname, pdb_resname)`` TODO CHECK
   * - Provide translations for atom names found in PDB files to those
       in topology files
     - ``pdbalias atom <residue name> <PDB atomname> <topology atomname>``
     - ``gen.alias_atom(resname, pdb_atomname, top_atomname)``

Query functions
---------------
.. list-table:: Query functions
   :header-rows: 1

   * - Task
     - TCL command
     - Python class equivalent
   * - List loaded topology files
     - ``topology list``
     - 
   * - List current segids
     - ``segment segids``
     - 
   * - List current resids
     - ``segment resids``
     - 
   * - List all currently applied patches, including default patches
     - ``patch listall``
     - 
   * - List all explicitly applied patches
     - ``patch list``
     - 
   * - Get residue name given a resid and a segment
     - ``segment residue <segment ID> <resid>``
     - 
   * - Get a list of atoms given a resid and segment
     - ``segment atoms <segment ID> <resid>``
     - 
   * - Get x,y,z coordinates for a given atom
     - ``segment coordinates <segment ID> <resid>``
     - 
   * - Get the name of the patch applied to the beginning of a given segment
     - ``segment first <segment ID>``
     - 
   * - Get the name of the patch applied to the end of a given segment
     - ``segment last <segment ID>``
     - 

System building functions
-------------------------
.. list-table:: Segment functions
   :header-rows: 1

   * - Task
     - TCL command
     - Python class equivalent
   * - Guess coordinates of atoms that aren't explicitly set
     - ``guesscoord``
     -  
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
     - 
   * - Remove angles and regenerate them, after patching
     - ``regenerate angles``
     - 
   * - Remove dihedrals and regenerate them, after patching.
     - ``regenerate dihedrals``
     - 
   * - Apply a patch to one or more residues, as determined by the patch
     - ``patch <patchname> <segid:resid> [...]``
     - 

Within-segment functions
------------------------
These methods appear within the ``segment {}`` declaration in TCL.

.. list-table:: Segment functions
   :header-rows: 1

   * - Task
     - TCL command
     - Python class equivalent
   * - Add a new segment to the molecule
     - ``segment <segment ID> {}``
     - 
   * - Mutate a residue in the current segment to a new one
     - ``mutate <resid> <resname>`` in ``segment {}`` block
     - 
   * - Extract sequence information from a PDB file when building segment.
     - ``pdb <pdbfilename>`` in ``segment {}`` block
     - 
   * - Add a singe residue to the end of current segment
     - ``residue <resid> <resname> [<chain>]``
     - 
   * - Override default patch for first residue in segment
     - ``first <patchname>``
     - 
   * - Override default patch for last residue in segment
     - ``last <patchname>``
     -
   * - Ovrride default topology settings for automatic generation of angles
       and dihedrals in segment
     - ``auto [angles] [dihedrals] [none]``
     - 


I/O functions
-------------
.. list-table:: I/O functions
   :header-rows: 1

   * - Task
     - TCL command
     - Python class equivalent
   * - Write out structure to a PSF file
     - ``writepsf [charmm] [x-plor] [cmap|nocmap] <filename>``
     - 
   * - Write out structure to a PDB file
     - ``writepdb <filename>``
     - 
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
     - 

Context functions
-----------------
.. list-table:: Context functions
   :header-rows: 1

   * - Task
     - TCL command
     - Python class equivalent
   * - Create a new context, with its own structure, topology definitions,
       and aliases, and set it to active
     - ``psfcontext new``
     - ``gen = new Psfgen()``
   * - Create a new context, but do not switch to it.
     - ``psfcontext create``
     - ``gen = new Psfgen()``
   * - Delete a psfcontext
     - ``psfcontext delete <context>``
     - ``del gen``
   * - Switch to a diferent psfcontext
     - ``psfcontext <id>``
     - No equivalent - use multiple Psfgen objects
   * - Make context case sensitive (by default it is not)
     - ``psfcontext mixedcase``
     - 
   * - Make context case insensitive (default setting)
     - ``psfcontext allcaps``
     - 
   * - Clear the structure, topology definitions, and aliases
     - ``psfcontext reset``
     - No equivalent - delete and make a new Psfgen object
   * - Evaluate commands in a given context, returning to the current one when
       done
     - ``psfcontext eval <context> { <commands> }``
     - No equivalent - use multiple Psfgen objects
   * - Get total number of contexts created and destroyed so far
     - ``psfcontext stats``
     - No equivalent - use multiple Psfgen objects
