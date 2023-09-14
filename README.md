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

Building
--------

_NOTE_: The following recipe is a bit out of date. My apologies.

This recipe assumes you're building on a Unix like operating system and will build only static libraries which will be installed in the `.local` folder in your home directory, which means you do not need super user powers to build and install. However, this also means you should have the `$HOME/.local/bin` path in your _PATH_ environment variable.

To set the path correctly, add the following line to your .bashrc file.

```bash
export PATH="$HOME/.local/bin:$PATH"
```

This tutorial is building on Ubuntu 18.04 LTS which comes with a very old build environment. We therefore `apt-get install g++-8` to make sure we have a usable compiler. If you have a different compiler you should replace the `g++-8` string with the appropriate compiler name.

cmake
-----

Building the libraries and tools is done using [CMake](https://cmake.org/). Please make sure you have at least version 3.16. If not, you will have to build and install cmake yourself and use it instead of the system provided version.

Building and installing cmake is done as follows:

```bash
wget https://github.com/Kitware/CMake/releases/download/v3.21.3/cmake-3.21.3.tar.gz
tar xf cmake-3.21.3.tar.gz
cd cmake-3.21.3
./configure --prefix=$HOME/.local
make
make install
```

After installing cmake, make sure it is found and works correctly. In the terminal type `cmake --version` and the result should be:

```bash
$ cmake --version
cmake version 3.21.3

CMake suite maintained and supported by Kitware (kitware.com/cmake).
```

boost
-----

So assuming you have a good compiler, we start by making sure you have the correct boost. This is perhaps the most tricky part and most often results in link errors or crashing applications. If your system comes with a packages boost version 1.70 or higher I suggest you use that and skip this step. If however the OS is bundled with an older version you will have to build boost yourself. Start by downloading the boost source code, extract, build and then install it. However, before building boost make sure you have the development packages for libz and bzip2 installed! For Ubuntu use `apt-get install zlib1g-dev libbz2-dev`.

(In this example I'm leaving out the python module of boost since it takes a long time to build and is not needed by any of my code.)

```bash
wget https://boostorg.jfrog.io/artifactory/main/release/1.77.0/source/boost_1_77_0.tar.bz2
tar xf boost_1_77_0.tar.bz2
cd boost_1_77
./bootstrap.sh
./b2 stage link=static --prefix=$HOME/.local --without-python
./b2 install link=static --prefix=$HOME/.local --without-python
```

That should install a recent boost in your .local folder. I've used boost version 1.77 in this example, but anything since 1.70 should work.

Please check the output of the `b2` command above to make sure the zlib and bzip2 components are included.

mrc
---

Although not strictly required, and not supported at all on MacOS, the tool [`mrc`](https://github.com/mhekkel/mrc.git) is highly recommended. It can be used to bundle data files into an executable making the executable more easily deployable. If you decide not to use this tool, you must make sure you have the data files for libcifpp available at run time.

```bash
git clone https://github.com/mhekkel/mrc.git
cd mrc
mkdir build
cd build
cmake .. -DCMAKE_CXX_COMPILER=g++-8 -DCMAKE_INSTALL_PREFIX=$HOME/.local
cmake --build .
ctest .
cmake --install .
```

Now run the tests:

```bash
ctest .
```

There should be no errors and all tests should pass.

```bash
cmake --install .
```

libcifpp
--------

We will build libcifpp with cmake as well. Every dependency should be up to date by now, so we can simply fetch, extract, configure and build the lib:

```bash
git clone https://github.com/PDB-REDO/libcifpp.git
cd libcifpp
mkdir build
cd build
cmake .. -DCMAKE_CXX_COMPILER=g++-8 -DBUILD_TESTING=ON
cmake --build .
```

Again, we've opted to build tests. Make sure they work without errors:

```bash
ctest .
```

If all is OK, install:

```bash
cmake --install .
```

An optional step here is to install the update script for CCD and pdbx dictionary files. To do this, you can add the -DINSTALL_UPDATE_SCRIPT=ON flag to the configure step:

```bash
cmake .. -DCMAKE_CXX_COMPILER=g++-8 -DBUILD_TESTING=ON -DINSTALL_UPDATE_SCRIPT=ON
```

If you then run `sudo cmake --install .` the cron script will be installed. However, you will need to have sudo rights to do so.

At this stage you can already build tools like [`DSSP`](https://github.com/PDB-REDO/dssp.git) and [`tortoize`](https://github.com/PDB-REDO/tortoize.git).

clipper
-------

The pdb-redo library depends on a recent version of clipper. If your system contains version 2.1.20201109 of clipper, you can simply install the system supplied development package for clipper and skip this step. On Debian 10 you would e.g. type `apt-get install libclipper-dev` to do this.

But we're using Ubuntu 18.04 here, so we need to build clipper ourselves.

Note that clipper depends on _ccp4_ and _fftw_. These libraries are hopefully supported by your operating system. On Ubuntu 18.04 I had to install

- fftw-dev
- sfftw-dev
- libccp4-dev

Unfortunately, this is the only step where I'm using `sudo`. If you cannot use sudo, you'll have to build the packages yourself installing them in `$(HOME)/.local`.

```bash
sudo apt-get install fftw-dev sfftw-dev libccp4-dev
```

Then the steps required were:

```bash
wget ftp://ftp.ccp4.ac.uk/opensource/clipper-2.1.20201109.tar.gz
tar xf clipper-2.1.20201109.tar.gz
cd clipper-2.1/
./configure --prefix=$HOME/.local --enable-ccp4
make
make install
```

newuoa
------

This library contains a C++ implementation of the NEWUOA algorithm. This is required to calculate atom radii.

The default package will build a shared library, but we prefer a static one:

```bash
git clone https://github.com/elsid/newuoa-cpp.git
cd newuoa-cpp
mkdir build
cd build
sed -e's/ SHARED / STATIC /g' -i.bak ../CMakeLists.txt
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/.local
cmake --build .
cmake --install .
```

libpdb-redo
-----------

Finally, we've arrived at building [`libpdb-redo`](https://github.com/PDB-REDO/libpdb-redo.git). The steps required now are perhaps familiar:

```bash
git clone https://github.com/PDB-REDO/libpdb-redo.git
cd libpdb-redo
mkdir build
cd build
cmake .. -DCMAKE_CXX_COMPILER=g++-8 -DBUILD_TESTING=ON
cmake --build .
ctest .
cmake --install .
```

Again, errors should not happen here.
