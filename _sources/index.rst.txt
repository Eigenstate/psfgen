.. psfgen documentation master file, created by
   sphinx-quickstart on Mon Jun 25 13:46:34 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

psfgen - python
===============
`Psfgen <http://www.ks.uiuc.edu/Research/vmd/plugins/psfgen/>`_ is a structure
building tool with a TCL interface commonly used as a VMD plugin to set up
simulations using NAMD or other molecular dynamics software. This package
provides Python bindings for psfgen.

Why a Python interface?
-----------------------
The usual way of interfacing with psfgen is through a series of TCL commands.
The only reason many researchers learn TCL is for using psfgen and VMD, but now
that VMD can be built with an `internal Python interpreter
<http://www.ks.uiuc.edu/Research/vmd/current/ug/node161.html>`_ or used as a
`standalone Python module <https://github.com/Eigenstate/vmd-python>`_, it makes
sense to also allow psfgen to listen to commands in Python.

As the internals of psfgen are all written in C, it wasn't too hard to write
a Python interface using Python's C API. I also think this kind of thing is fun
for some reason. I also have some selfish motivations, as having a Python
interface greatly simplifies how my system building software `Dabble
<https://dabble.robinbetz.com>`_ interfaces with psfgen.

Installation
------------
You need a working C compiler. Either Python 2 or 3 will work, but I'd recommend
switching to Python 3 as soon as possible. Building is easy:

.. code-block:: bash

    python setup.py install

Usage
-----
All interactions with psfgen are done by creating a PsfGen object that
represents an individual molecular system. A very truncated example:

.. code-block:: python

    from psfgen import PsfGen
    gen = PsfGen()
    # system building commands...
    gen.write_psf(filename="system.psf")

You may want to jump to some :doc:`examples </examples>`, or if you're familiar with the
TCL interface, view :doc:`equivalent function calls in Python </all_commands>`,
or read the whole :doc:`API spec </api>`.

.. toctree::
    :maxdepth: 2
    :hidden:

    all_commands
    examples
    api

