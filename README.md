# psfgen

This is psfgen 1.7.0 with some patches and a SWIG python interface.

## Why???

The idea for this is based on my finding of a swig interface file in the VMD
CVS plugins repository from 2002. 

I like the idea of a Python interface to psfgen, as many people learn Tcl only
for obscure MD preparation tools and they shouldn't have to.

This will also enable Dabble to avoid writing Tcl files and executing them,
which is generally more elegant.

## Installation

To build and install, you'll need a C compiler:

    python setup.py install

## Usage

Well, I'm still figuring that out.

## License

This software is classified as a "derivative work" of the original psfgen
code. It uses the main psfgen code and a swig interface file that are the work
of the original developers (Justin Gullingsrud and Jim Phillips). I have made
it installable, and am fixing warnings in the C code. I will probably rewrite
the Python interface at some point (this is early).

psfgen is licensed by the University of Illinois' Non-Exclusive, Non-Commercial
use license.

"This software includes code developed by the Theoretical and Computational
Biophysics Group in the Beckman Institute for Advanced Science and Technology
at the University of Illinois at Urbana-Champaign."

