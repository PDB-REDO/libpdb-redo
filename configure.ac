AC_PREREQ([2.69])

AC_INIT([libpdb-redo], 1.0, [m.hekkelman@nki.nl])

dnl Switch to a decent C++ compiler, and check if it works.
AC_LANG(C++)
AX_CXX_COMPILE_STDCXX_17([noext])

AC_CONFIG_SRCDIR([src/AtomShape.cpp])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_HEADERS([src/config.hpp])
# AC_CONFIG_MACRO_DIR([config/m4])

AC_PREFIX_DEFAULT(/usr/local)

AC_DEFUN([read_test], [AC_LANG_SOURCE(
	esyscmd(config/tools/m4esc.sh config-tests/$1))])

AC_PROG_INSTALL

LT_INIT
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

AX_CHECK_COMPILE_FLAG([-fstandalone-debug], , , [-Werror])

AC_PATH_PROG([PKG_CONFIG], [pkg-config])

AC_CHECK_LIB([m], [sin])
AX_PKG_CHECK_MODULES([GSL], [gsl], [], [], [AC_MSG_ERROR([the required package gsl-dev is not installed])])

AX_PTHREAD

AC_CHECK_HEADER([filesystem], [], [AC_MSG_ERROR([The file <filesystem> is missing, perhaps you should install a more recent libstdc++ implementation.])])

AC_ARG_WITH([cif++],
	AS_HELP_STRING([--with-cif++=@<:@location@:>@],
		[Use the cif++ library as specified.]),
		[
			AS_IF([test -d ${withval}/include], [], [
				AC_MSG_ERROR(['${withval}'' is not a valid directory for --with-cif++])
			])
			CIFPP_CFLAGS="-I ${withval}/include"
			CIFPP_LIBS="-L${withval}/.libs -lcifpp"

			AC_SUBST([CIFPP_CFLAGS], [$CIFPP_CFLAGS])
			AC_SUBST([CIFPP_LIBS], [$CIFPP_LIBS])
		])

AS_IF([test "x$CIFPP_LIBS" = "x"], [
	if test -x "$PKG_CONFIG"
	then
		AX_PKG_CHECK_MODULES([CIFPP], [libcifpp], [], [], [AC_MSG_ERROR([the required package libcif++ is not installed])])
	else
		AC_CHECK_HEADER(
	        [cif++/Config.hpp],
	        [
				dnl CIFPP_CFLAGS="-I ${withval}/include"
			],
			[AC_MSG_ERROR([
Can't find the libcif++ header, Config.hpp.  Make sure that it
is installed, and either use the --with-cif++ option or install
pkg-config.])])

			AX_CHECK_LIBRARY([CIFPP], [cif++/Config.hpp], [cifpp],
				 [
					 LIBS="-lcifpp $LIBS"
				 ],
                 [AC_MSG_ERROR([libcif++ not found])])
	fi
])

AC_DEFINE([HAVE_CIFPP], [1], [Use CIFPP])

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

AC_DEFINE([HAVE_ZEEP], [1], [Use libzeep])

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

AC_DEFINE([HAVE_CLIPPER], [1], [Use libclipper])

AX_CHECK_LIBRARY([DLIB], [dlib/global_optimization.h], [dlib],
		[
			DLIB_LIBS="-ldlib $LIBS"
		],
		[AC_MSG_ERROR([libdlib-dev not found])])

AC_SUBST([DLIB_CFLAGS], [])
AC_SUBST([DLIB_LIBS], [$DLIB_LIBS])

dnl Process Makefile.in to create Makefile
AC_OUTPUT([GNUmakefile])
