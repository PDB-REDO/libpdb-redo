AC_PREREQ([2.69])

m4_define([libpdb_redo_version_major],[2])
m4_define([libpdb_redo_version_minor],[0])
m4_define([libpdb_redo_version_micro],[0])
m4_define([libpdb_redo_version_extra],[])
m4_define([libpdb_redo_version],[libpdb_redo_version_major().libpdb_redo_version_minor().libpdb_redo_version_micro()libpdb_redo_version_extra])

AC_INIT([libpdb-redo], [libpdb_redo_version], [m.hekkelman@nki.nl])

dnl Switch to a decent C++ compiler, and check if it works.
AC_LANG(C++)
AX_CXX_COMPILE_STDCXX_17([noext])

AX_CHECK_COMPILE_FLAG([-fstandalone-debug],
	[
		CXXFLAGS="$CXXFLAGS -fstandalone-debug"
	] , , [-Werror])

AC_CONFIG_SRCDIR([src/AtomShape.cpp])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([config/m4])
AC_CONFIG_HEADERS([src/config.hpp])

AC_PREFIX_DEFAULT(/usr/local)

AC_DEFUN([read_test], [AC_LANG_SOURCE(
	esyscmd(config/tools/m4esc.sh config/test/$1))])

AC_PROG_INSTALL
PKG_PROG_PKG_CONFIG

AX_PTHREAD

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
LIBPDB_REDO_SEMANTIC_VERSION=$libpdb_redo_version
AC_SUBST(LIBPDB_REDO_SEMANTIC_VERSION)

AC_ARG_VAR([DEBUG], [Build a debug version of the application])

AX_MRC

dnl revision numbering is something used internally at the NKI
AC_ARG_ENABLE(
	revision,
	[AS_HELP_STRING([--disable-revision], [Create a build number as revision])])

AS_IF([test "x$enable_revision" != "xno" ], [
	UPDATE_REVISION=1
])

AC_SUBST([UPDATE_REVISION], [$UPDATE_REVISION])

AX_FILESYSTEM

AC_CHECK_LIB([m], [sin])

AX_BOOST_BASE([1.71], [], [AC_MSG_ERROR([Sorry, your boost is not found or not up-to-date.])])
AX_BOOST_IOSTREAMS
AX_BOOST_DATE_TIME
AX_BOOST_REGEX

AC_CHECK_HEADER([endian.h], [
	HAVE_ENDIAN_H=1
	AC_DEFINE([HAVE_ENDIAN_H], [$HAVE_ENDIAN_H], [Define to 1 if you have the <endian.h> header file.])])
AC_CHECK_HEADER([sys/endian.h], [
	HAVE_SYS_ENDIAN_H=1
	AC_DEFINE([HAVE_SYS_ENDIAN_H], [$HAVE_SYS_ENDIAN_H], [Define to 1 if you have the <sys/endian.h> header file.])])

AS_IF([ test "x${HAVE_ENDIAN_H}" = "x" && test "x${HAVE_SYS_ENDIAN_H}" = "x" ], 
	[AC_MSG_WARN([The files <endian.h> and <sys/endian.h> are both missing, the code most likely will not build.])]
)

AX_LIBCIFPP
AX_LIBZEEP

AX_LIBCLIPPER

AC_MSG_CHECKING([clipper version])
AC_COMPILE_IFELSE(
	[read_test(clipper-test.cpp)],
	[AC_MSG_RESULT([ok])],
	[AC_MSG_ERROR([The version of clipper is not up to date])])

AX_NEWUOA

dnl gsl?

dnl check for GSL
AS_IF([test "x$GSL_LIBS" = "x" ],
	[
		if test -x "$PKG_CONFIG"
		then
			AX_PKG_CHECK_MODULES([GSL], [gsl], [], [], [AC_MSG_ERROR([the required package gsl-dev is not installed])])
		else
			AC_CHECK_HEADER(
				[gsl/gsl_version.h],
				[],
				[AC_MSG_ERROR([Can't find the gsl header, gsl/gsl_version.h.  Make sure that libgsl is installed.])])

			saved_LIBS="$LIBS"

			AX_CHECK_LIBRARY([GSLBLAS], [gsl/gsl_version.h], [gslcblas],
					[ LIBS="$LIBS -lgslcblas" ], [AC_MSG_ERROR([libgslcblas not found])])

			AX_CHECK_LIBRARY([GSL], [gsl/gsl_version.h], [gsl],
					[ GSL_LIBS="-lgsl -lgslcblas -lm" ],
					[AC_MSG_ERROR([libgsl not found])])
			
			LIBS="$saved_LIBS"
			HAVE_GSL=1
			AC_SUBST([HAVE_GSL], [1])
		fi
	]
)

dnl These are still needed outside the Debian environment
AX_CHECK_LIBRARY([LIBZ], [zlib.h], [z],
	[ LIBS="$LIBS -lz" ],
	[AC_MSG_ERROR([libz not found - compressed files not supported])])
AX_CHECK_LIBRARY([LIBBZ2], [bzlib.h], [bz2],
	[ LIBS="$LIBS -lbz2"],
	[AC_MSG_ERROR([libbz2 not found - compressed files not supported])])

AC_SUBST([LIBS], [$LIBS])

dnl Set output variables for the various directories
AC_ARG_VAR([DATA_LIB_DIR], [The location where to store auxiliary files])
AC_SUBST([DATA_LIB_DIR])
DATA_LIB_DIR=$datadir/libpdb-redo

dnl Process Makefile.in to create Makefile
AC_OUTPUT([GNUmakefile
           libpdb-redo.pc])
