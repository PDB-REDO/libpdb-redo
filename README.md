libpdb-redo
===========

This is the README file for libpdb-redo. This library contains code
shared by the various tools we develop at the NKI for the
[PDB-REDO](https://pdb-redo.eu/) project.

Requirements
------------

The code for this library was written in C++17. You therefore need a
recent compiler to build it. For the development gcc 9.3 and clang 9.0
have been used.

Other requirements are:

- The clipper library, either the latest from CCP4 or version 2020-11-09
- GNU make version 4.1 or higher
- GSL, the GNU Scientific Library
- [newuoa-cpp](https://github.com/elsid/newuoa-cpp), required to
  calculate atom radii
- [libzeep](https://github.com/mhekkel/libzeep), a library that
  contains a full validating XML parser as well as a complete HTTP,
  SOAP and REST server implementation
- autoconf and the autoconf-archive

Building
--------

Make sure you install the libraries and tools in the list above first,
then simply run:

```
autoreconf -if
./configure
make
sudo make install
```
