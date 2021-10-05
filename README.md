libpdb-redo
===========

This is the README file for libpdb-redo. This library contains code
shared by the various tools we develop at the NKI for the
[PDB-REDO](https://pdb-redo.eu/) project.

Requirements
------------

The code for this library was written in C++17. You therefore need a
recent compiler to build it. For the development gcc 9.3 and clang 9.0
have been used. On MS Windows you'll need at least the 2019 version
of MSVC.

Other requirements are:

- The clipper library, either the latest from CCP4 or version 2020-11-09
- [newuoa-cpp](https://github.com/elsid/newuoa-cpp), required to
  calculate atom radii
- [libzeep](https://github.com/mhekkel/libzeep), a library that
  contains a full validating XML parser as well as a complete HTTP,
  SOAP and REST server implementation
- [cmake](https://cmake.org)


Building
--------

This library uses [cmake](https://cmake.org). The usual way of building
and installing is to create a `build` directory and run cmake there.

On linux e.g. you would issue the following commands:

```
	git clone https://github.com/PDB-REDO/libpdb-redo.git
	cd libpdb-redo
	mkdir build
	cd build
	cmake ..
	cmake --build . --config Release
	ctest -C Release
	cmake --install .
```

This checks out the source code from github, creates a new directory
where cmake stores its files. Run a configure, build the code and run
tests. And then it installs the library and auxiliary files.

The default is to install everything in `$HOME/.local` on Linux and
`%LOCALAPPDATA%` on Windows (the AppData/Local folder in your home directory).
You can change this by specifying the prefix with the
[CMAKE_INSTALL_PREFIX](https://cmake.org/cmake/help/v3.21/variable/CMAKE_INSTALL_PREFIX.html)
variable.

