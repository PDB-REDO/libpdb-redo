libpdb-redo
===========

This is the README file for libpdb-redo. This library contains code
shared by the various tools we develop at the NKI for the
[PDB-REDO](https://pdb-redo.eu/) project.

Requirements
------------

The code for this library was written in C++17. You therefore need a
recent compiler to build it. For the development gcc with version 9.4.0
and up and clang 9.0 have been used. On MS Windows you'll need at least
the 2019 version of MSVC.

Other requirements are:

- The clipper library, either the latest from CCP4 or version 2020-11-09
- [newuoa-cpp](https://github.com/elsid/newuoa-cpp), required to
  calculate atom radii
- [libcifpp](https://github.com/PDB-REDO/libcifpp.git), a library containing
  code to read and manipulate macro molecular models in mmCIF and PDB format.
- [gsl](https://www.gnu.org/software/gsl/), the GNU Scientific Library.
  Usually you can install this using a package manager on your OS. In
  Debian/Ubuntu the required package is libgsl-dev.
- [cmake](https://cmake.org)

