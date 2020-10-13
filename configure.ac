AC_INIT([pdb-redo-tools], 1.0, [m.hekkelman@nki.nl])

dnl Switch to a decent C++ compiler, and check if it works.
AC_LANG(C++)
AX_CXX_COMPILE_STDCXX_17([noext])

AC_CONFIG_SRCDIR([src/pr-main.cpp])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_HEADERS([src/pdb-redo.h])
dnl AC_CONFIG_MACRO_DIR([config/m4])

AC_PREFIX_DEFAULT(/usr/local)

AC_PROG_INSTALL

AC_ARG_VAR([DEBUG], [Build a debug version of the application])

AX_CHECK_COMPILE_FLAG([-fstandalone-debug], , , [-Werror])

AC_PATH_PROG([PKG_CONFIG], [pkg-config])

AC_CHECK_LIB([m], [sin])
AX_PKG_CHECK_MODULES([GSL], [gsl], [], [], [AC_MSG_ERROR([the required package gsl-dev is not installed])])

AC_ARG_VAR([MRC], [Specify a location for the mrc executable])

dnl We really want to use mrc
if test "x$MRC" = "x"; then
	AC_PATH_PROG([MRC], [mrc])
fi

if test "x$MRC" = "x"; then
	AC_MSG_ERROR([mrc not found, the application will be built without resources])
else
	AC_DEFINE([USE_RSRC], [1])
fi

AX_PTHREAD

AC_CHECK_HEADER([filesystem], [], [AC_MSG_ERROR([The file <filesystem> is missing, perhaps you should install a more recent libstdc++ implementation.])])

AX_BOOST_BASE([1.65.1], [], [AC_MSG_ERROR([Could not find a recent version of boost])])
AX_BOOST_IOSTREAMS
AX_BOOST_THREAD
AX_BOOST_REGEX
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_DATE_TIME

AC_ARG_WITH([cif++],
	AS_HELP_STRING([--with-cif++=@<:@location@:>@],
		[Use the cif++ library as specified.]),
		[
			AS_IF([test -d ${withval}/include], [], [
				AC_MSG_ERROR(['${withval}'' is not a valid directory for --with-cif++])
			])
			dnl AC_SUBST([CIFPP_CFLAGS], ["-I ${withval}/include"])
			dnl AC_SUBST([CIFPP_LIBS], ["-L${withval}/.libs -lcif++"])

			CIFPP_CFLAGS="-I ${withval}/include"
			CIFPP_LIBS="-L${withval}/.libs -lcif++"
			CIFPP_RSRC="${withval}/rsrc"

			AC_SUBST([CIFPP_CFLAGS], [$CIFPP_CFLAGS])
			AC_SUBST([CIFPP_LIBS], [$CIFPP_LIBS])
		])

AC_ARG_WITH([newuoa],
	AS_HELP_STRING([--with-newuoa=@<:@location@:>@],
		[Use the newuoa library as specified.]),
		[
			CPPFLAGS="$CPPFLAGS -I ${withval}/include"
			CXXFLAGS="$CXXFLAGS -I ${withval}/include"
			LDFLAGS="$LDFLAGS -L${withval}/lib"
		])

AC_SUBST([CIFPP_RSRC], [$CIFPP_RSRC])

AS_IF([test "x$CIFPP_LIBS" = "x"], [
	AC_PATH_PROG([PKG_CONFIG], [pkg-config])
	if test -x "$PKG_CONFIG"
	then
		AX_PKG_CHECK_MODULES([CIFPP], [libcif++], [], [], [AC_MSG_ERROR([the required package libcif++ is not installed])])
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

			AX_CHECK_LIBRARY([CIFPP], [cif++/Config.hpp], [cif++],
				 [
					 LIBS="-lcif++ $LIBS"
				 ],
                 [AC_MSG_ERROR([libcif++ not found])])
	fi
])

dnl Check for libzeep

AC_ARG_WITH([zeep],
	AS_HELP_STRING([--with-zeep=@<:@location@:>@],
		[Use the libzeep library as specified.]),
		[
			AS_IF([test -d ${withval}/include], [], [
				AC_MSG_ERROR(['${withval}'' is not a valid directory for --with-zeep])
			])

			ZEEP_CFLAGS="-I ${withval}/include"
			ZEEP_LIBS="-L${withval}/.libs -lzeep-json"

			AC_SUBST([ZEEP_CFLAGS], [$ZEEP_CFLAGS])
			AC_SUBST([ZEEP_LIBS], [$ZEEP_LIBS])
		])

AS_IF([test "x$ZEEP_LIBS" = "x"], [
	if test -x "$PKG_CONFIG"
	then
		AX_PKG_CHECK_MODULES([ZEEP], [libzeep-json], [], [], [AC_MSG_ERROR([the required package libzeep is not installed])])
	else
		AC_CHECK_HEADER(
			[zeep/json/element.hpp],
			[],
			[AC_MSG_ERROR([
Can't find the libzeep header, zeep/json/element.hpp.  Make sure that libzeep
is installed, and either use the --with-zeep option or install
pkg-config.])])

			AX_CHECK_LIBRARY([ZEEP], [zeep/json/element.hpp], [zeep-json],
					[ LIBS="-lzeep-json $LIBS" ],
					[AC_MSG_ERROR([libzeep-json not found])])
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
	  ])


AC_ARG_WITH([clipper],
	AS_HELP_STRING([--with-clipper=@<:@location@:>@],
		[Use the clipper library as specified.
		@<:@location=$CCP4@:>@]),
		[
			CPPFLAGS="$CPPFLAGS -I ${withval}/include"
			CXXFLAGS="$CXXFLAGS -I ${withval}/include"
			LDFLAGS="$LDFLAGS -L${withval}/lib"
		])

AC_CHECK_HEADER(
	[clipper/clipper.h],
	[],
	[AC_MSG_ERROR([
Can't find the main clipper include file clipper/clipper.h. Please
specify either the CCP4 location of the location of the installed
clipper using --with-clipper=<path>
])])

AX_CHECK_LIBRARY([LIBCLIPPER], [clipper/clipper.h], [clipper-core],
	[
		LIBS="$LIBS -lclipper-ccp4 -lclipper-minimol -lclipper-mmdb -lclipper-cns -lclipper-phs -lclipper-contrib -lclipper-core -lrfftw -lfftw -lccp4c"
	],
	[AC_MSG_ERROR([libclipper not found - cannot continue])])

AC_CHECK_HEADER(
	[newuoa.h],
	[
		LIBS="$LIBS -lnewuoa"
	],
	[AC_MSG_ERROR([
Can't find the newuoa include file newuoa.h. Please install this
library and specify its location with --with-newuoa=<path>.
You can find newuoa at https://github.com/elsid/newuoa-cpp.git
])])

AX_CHECK_LIBRARY([LIBZ], [zlib.h], [z],
				 [ LIBS="$LIBS -lz" ],
                 [AC_MSG_ERROR([libz not found - compressed files not supported])])
AX_CHECK_LIBRARY([LIBBZ2], [bzlib.h], [bz2],
				 [ LIBS="$LIBS -lbz2"],
                 [AC_MSG_ERROR([libbz2 not found - compressed files not supported])])

dnl AC_MSG_CHECKING([clipper version])
dnl AC_COMPILE_IFELSE(
dnl 	[read_test(clipper-test.cpp)],
dnl 	[],
dnl 	[AC_MSG_ERROR([The version of clipper is not up to date])])
dnl 
dnl AC_CHECK_HEADER(
dnl 	[newuoa.h],
dnl 	[],
dnl 	[AC_MSG_ERROR([
dnl Can't find the newuoa include file newuoa.h. Please install this
dnl library and specify its location with --with-newuoa=<path>.
dnl You can find newuoa at https://github.com/elsid/newuoa-cpp.git
dnl ])])

dnl Process Makefile.in to create Makefile
AC_OUTPUT([GNUmakefile])
