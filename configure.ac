AC_PREREQ([2.69])

AC_INIT([libpdb-redo], 1.0, [m.hekkelman@nki.nl])

dnl Switch to a decent C++ compiler, and check if it works.
AC_LANG(C++)
AX_CXX_COMPILE_STDCXX_17([noext])

AX_CHECK_COMPILE_FLAG([-fstandalone-debug],
	[
		CXXFLAGS="$CXXFLAGS -fstandalone-debug"
	] , , [-Werror])

AC_CONFIG_SRCDIR([src/AtomShape.cpp])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_HEADERS([src/config.hpp])
AC_CONFIG_MACRO_DIR([config/m4])

AC_PREFIX_DEFAULT(/usr/local)

AC_DEFUN([read_test], [AC_LANG_SOURCE(
	esyscmd(config/tools/m4esc.sh config-tests/$1))])

AC_PROG_INSTALL

LT_INIT([disable-shared pic-only])
AC_SUBST(LIBTOOL_DEPS)

dnl versioning, first for libtool
LIBPDB_REDO_CURRENT=1
LIBPDB_REDO_REVISION=0
LIBPDB_REDO_AGE=0

LIBPDB_REDO_LT_CURRENT="${LIBPDB_REDO_CURRENT}"
LIBPDB_REDO_LT_VERSION="${LIBPDB_REDO_CURRENT}:${LIBPDB_REDO_REVISION}:${LIBPDB_REDO_AGE}"

AC_SUBST(LIBPDB_REDO_LT_CURRENT)
AC_SUBST(LIBPDB_REDO_LT_VERSION)

dnl and now for the semantic version
LIBPDB_REDO_SEMANTIC_VERSION=1.0.0
AC_SUBST(LIBPDB_REDO_SEMANTIC_VERSION)

AC_ARG_VAR([DEBUG], [Build a debug version of the application])

AC_ARG_VAR([MRC], [Specify a location for the mrc executable])

PKG_PROG_PKG_CONFIG

dnl using resources?
USE_RSRC=0

if test "x$MRC" = "x"; then
	AC_PATH_PROG([MRC], [mrc])
fi

if test "x$MRC" = "x"; then
	AC_MSG_WARN([The mrc application was not found, not using resources.])
else
	AC_ARG_ENABLE(
		resources,
		[AS_HELP_STRING([--disable-resources], [Use mrc to store data in resources])])

	AS_IF([test "x$enable_resources" != "xno" ], [
		USE_RSRC=1
	])
fi

AC_SUBST([USE_RSRC], [$USE_RSRC])

AC_DEFINE_UNQUOTED([USE_RSRC], [$USE_RSRC], [Use mrc to store resources])

dnl revision numbering is something used internally at the NKI
AC_ARG_ENABLE(
	revision,
	[AS_HELP_STRING([--enable-revision], [Create a build number as revision])])

AS_IF([test "x$enable_revision" = "xyes" ], [
	UPDATE_REVISION=1
])

AC_SUBST([UPDATE_REVISION], [$UPDATE_REVISION])

AC_CHECK_LIB([m], [sin])

AX_PTHREAD

AX_BOOST_BASE([1.71], [], [
	AC_MSG_ERROR([Sorry, your boost is not found or not up-to-date.])
])
AX_BOOST_IOSTREAMS
AX_BOOST_DATE_TIME
AX_BOOST_REGEX


AC_CHECK_HEADER([filesystem], [], [AC_MSG_ERROR([The file <filesystem> is missing, perhaps you should install a more recent libstdc++ implementation.])])

dnl check if we need stdc++fs as library
AC_TRY_LINK(
	[#include <filesystem>],
	[(void)std::filesystem::current_path();],
	[],
	[
		LIBS="$LIBS -lstdc++fs"

		AC_TRY_LINK(
			[#include <filesystem>],
			[(void)std::filesystem::current_path();],
			[],
			[
				AC_MSG_ERROR([Could not link filesystem])
			]
		)
	]
)

AC_CHECK_HEADER([endian.h], [AC_DEFINE([HAVE_ENDIAN_H], [1], [Define to 1 if you have the <endian.h> header file.])], [AC_MSG_WARN([The file <endian.h> is missing, the code most likely will not build.])])

CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"

AC_ARG_WITH([cif++],
	AS_HELP_STRING([--with-cif++=@<:@location@:>@],
		[Use the cif++ library as specified.]),
		[
			AS_IF([test -d ${withval}/include], [], [
				AC_MSG_ERROR(['${withval}'' is not a valid directory for --with-cif++])
			])
			CIFPP_CFLAGS="-I ${withval}/include"
			CIFPP_LIBS="-L${withval}/.libs -lcifpp"
			CIFPP_RSRC="${withval}/rsrc"

			AC_SUBST([CIFPP_CFLAGS], [$CIFPP_CFLAGS])
			AC_SUBST([CIFPP_LIBS], [$CIFPP_LIBS])
		])

AC_SUBST([CIFPP_RSRC], [$CIFPP_RSRC])

AS_IF([test "x$CIFPP_LIBS" = "x"], [
	if test -x "$PKG_CONFIG"
	then
		AX_PKG_CHECK_MODULES([CIFPP], [libcifpp], [], [], [AC_MSG_ERROR([the required package libcifpp is not installed])])
	else
		AC_CHECK_HEADER(
	        [cif++/Cif++.hpp],
	        [
				dnl CIFPP_CFLAGS="-I ${withval}/include"
			],
			[AC_MSG_ERROR([
Can't find the libcif++ header, Config.hpp.  Make sure that it
is installed, and either use the --with-cif++ option or install
pkg-config.])])

			AX_CHECK_LIBRARY([CIFPP], [cif++/Cif++.hpp], [cifpp],
				 [
					 LIBS="-lcifpp $LIBS"
				 ],
                 [AC_MSG_ERROR([libcif++ not found])])
	fi
])

dnl check for GSL
if test -x "$PKG_CONFIG"
then
	AX_PKG_CHECK_MODULES([GSL], [gsl], [], [], [AC_MSG_ERROR([the required package gsl-dev is not installed])])
else
	AC_CHECK_HEADER(
		[gsl/gsl_version.h],
		[],
		[AC_MSG_ERROR([Can't find the gsl header, gsl/gsl_version.h.  Make sure that libgsl is installed.])])

		AX_CHECK_LIBRARY([GSL], [gsl/gsl_version.h], [gsl],
				[ GSL_LIBS="-lgsl -lgslcblas -lm" ],
				[AC_MSG_ERROR([libgsl not found])])
fi

dnl Check for libzeep

AC_ARG_WITH([zeep],
	AS_HELP_STRING([--with-zeep=@<:@location@:>@],
		[Use the libzeep library as specified.]),
		[
			AS_IF([test -d ${withval}/include], [], [
				AC_MSG_ERROR(['${withval}'' is not a valid directory for --with-zeep])
			])

			ZEEP_CFLAGS="-I ${withval}/include"
			ZEEP_LIBS="-L${withval}/.libs -lzeep"

			AC_SUBST([ZEEP_CFLAGS], [$ZEEP_CFLAGS])
			AC_SUBST([ZEEP_LIBS], [$ZEEP_LIBS])
		])

AS_IF([test "x$ZEEP_LIBS" = "x"], [
	if test -x "$PKG_CONFIG"
	then
		AX_PKG_CHECK_MODULES([ZEEP], [libzeep], [], [], [AC_MSG_ERROR([the required package libzeep is not installed])])
	else
		AC_CHECK_HEADER(
			[zeep/json/element.hpp],
			[],
			[AC_MSG_ERROR([
Can't find the libzeep header, zeep/json/element.hpp.  Make sure that libzeep
is installed, and either use the --with-zeep option or install
pkg-config.])])

			AX_CHECK_LIBRARY([ZEEP], [zeep/json/element.hpp], [zeep],
					[ LIBS="-lzeep $LIBS" ],
					[AC_MSG_ERROR([libzeep not found])])
	fi
])

AC_ARG_VAR([CCP4], [The location where CCP4 is installed])

dnl check if we need stdc++fs as library
AC_TRY_LINK(
	[#include <filesystem>],
	[(void)std::filesystem::current_path();],
	[],
	[
		LIBS="$LIBS -lstdc++fs"

		AC_TRY_LINK(
			[#include <filesystem>],
			[(void)std::filesystem::current_path();],
			[],
			[
				AC_MSG_ERROR([Could not link filesystem])
			]
		)
	]
)

AS_IF([test x"$CCP4" != x""],
	[
	CPPFLAGS="$CPPFLAGS -I ${CCP4}/include"
	CXXFLAGS="$CXXFLAGS -I ${CCP4}/include"
	LDFLAGS="$LDFLAGS -L${CCP4}/lib -Wl,-rpath=${CCP4}/lib"

	AC_CHECK_HEADER(
		[clipper/clipper.h],
		[],
		[AC_MSG_ERROR([
Can't find the libclipper header, clipper/clipper.h.  Make sure that libclipper
is installed, and either use the --with-clipper option or install
pkg-config.])])

	AX_CHECK_LIBRARY([CLIPPER], [clipper/clipper.h], [clipper-core],
		[
			CLIPPER_LIBS="-L${CCP4}/lib -lclipper-ccp4 -lclipper-cif -lclipper-minimol -lclipper-mmdb -lclipper-cns -lclipper-phs -lclipper-contrib -lclipper-core -lccp4c -lmmdb2 -lrfftw -lfftw"
		],
		[AC_MSG_ERROR([libclipper not found])])
	])

dnl Check for libclipper

AC_ARG_WITH([clipper],
	AS_HELP_STRING([--with-clipper=@<:@location@:>@],
		[Use the libclipper library as specified.]),
		[
			AS_IF([test -d ${withval}], [], [
				AC_MSG_ERROR(['${withval}'' is not a valid directory for --with-clipper])
			])

			CLIPPER_CFLAGS="-I ${withval}"
			CLIPPER_LIBS="-L${withval}/lib -lclipper-ccp4 -lclipper-cif -lclipper-minimol -lclipper-mmdb -lclipper-cns -lclipper-phs -lclipper-contrib -lclipper-core -lccp4c -lmmdb2"

			AC_SUBST([CLIPPER_CFLAGS], [$CLIPPER_CFLAGS])
			AC_SUBST([CLIPPER_LIBS], [$CLIPPER_LIBS])
		])

AS_IF([test "x$CLIPPER_LIBS" = "x"], [
	if test -x "$PKG_CONFIG"
	then
		AX_PKG_CHECK_MODULES([CLIPPER], [clipper], [], [], [AC_MSG_ERROR([the required package libclipper is not installed])])
	else
		AC_CHECK_HEADER(
			[clipper/clipper.h],
			[],
			[AC_MSG_ERROR([
Can't find the libclipper header, clipper/clipper.h.  Make sure that libclipper
is installed, and either use the --with-clipper option or install
pkg-config.])])

		AX_CHECK_LIBRARY([CLIPPER], [clipper.h], [clipper-core],
				[ LIBS="-lclipper-fortran -lclipper-ccp4 -lclipper-cif -lclipper-minimol -lclipper-mmdb -lclipper-cns -lclipper-phs -lclipper-contrib -lclipper-core -lccp4c -lmmdb2 $LIBS" ],
				[AC_MSG_ERROR([libclipper not found])])
	fi
])

AC_MSG_CHECKING([clipper version])
AC_COMPILE_IFELSE(
	[read_test(clipper-test.cpp)],
	[AC_MSG_RESULT([ok])],
	[AC_MSG_ERROR([The version of clipper is not up to date])])


AX_CHECK_LIBRARY([DLIB], [dlib/global_optimization.h], [dlib],
		[
			HAVE_DLIB=1
			DLIB_LIBS="-ldlib $LIBS"
		],
		[AC_MSG_WARN([libdlib-dev not found])])

AC_SUBST([DLIB_CFLAGS], [])
AC_SUBST([DLIB_LIBS], [$DLIB_LIBS])

if test -x "$PKG_CONFIG"
then
	AX_PKG_CHECK_MODULES([NEWUOA], [newuoa], [], [
		HAVE_NEWUOA=1
		AC_DEFINE([HAVE_NEWUOA], [1])
	], [AC_MSG_WARN([the required package libnewuoa-dev is not installed])])
else
	AC_CHECK_HEADER(
		[newuoa.h], [], [AC_MSG_ERROR([
Can't find the libcif++ header, Config.hpp.  Make sure that it
is installed, and either use the --with-cif++ option or install
pkg-config.])])

	AX_CHECK_LIBRARY([NEWUOA], [newuoa.h], [newuoa],
		[
			LIBS="-lnewuoa $LIBS"
		],
		[AC_MSG_WARN([libnewuoa not found])])
fi

AS_IF([ test "x${HAVE_DLIB}" = "x" && test "x${HAVE_NEWUOA}" = "x" ], [
	AC_MSG_ERROR([Either libnewuoa-dev or libdlib-dev should be installed])
])

dnl These are still needed outside the Debian environment
AX_CHECK_LIBRARY([LIBZ], [zlib.h], [z],
	[ LIBS="$LIBS -lz" ],
	[AC_MSG_ERROR([libz not found - compressed files not supported])])
AX_CHECK_LIBRARY([LIBBZ2], [bzlib.h], [bz2],
	[ LIBS="$LIBS -lbz2"],
	[AC_MSG_ERROR([libbz2 not found - compressed files not supported])])

AC_SUBST([LIBS], [$LIBS])

dnl Process Makefile.in to create Makefile
AC_OUTPUT([GNUmakefile
           libpdb-redo.pc])
