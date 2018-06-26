.. _examples:
.. currentmodule:: examples
.. highlight:: python

Documentation examples
======================

The documentation for TCL-psfgen features a few examples. Here, I have
implemented the same examples using the Python interface.


Preparing separate PDB files
----------------------------
This example uses three PDB files, containing water and two protein fragments,
and composites them into a single PSF file.

.. code-block:: python

    from psfgen import Psfgen
    gen = Psfgen()


Deleting unwanted atoms
-----------------------
This example deletes some water molecues that are too far from the substrate.
The vmd-python interface is used here to select which waters.

.. code-block:: python

    from psfgen import Psfgen
    gen = Psfgen()


BPTI example
------------
Just like in the psfgen documentation, you need the charmm topology files and
the PDB file 6PTI.pdb from the Protein Data Bank.

.. code-block:: python

    from psfgen import Psfgen
    gen = Psfgen()


Building solvent around a protein
---------------------------------
Here we will use psfgen and vmd-python to prepare a system from the command line.
I would recommend using Dabble instead here.
