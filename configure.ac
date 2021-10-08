
m4_define([idemuxCPP_version_major], [0])
m4_define([idemuxCPP_version_minor], [1])
m4_define([idemuxCPP_version_patch], [9])

my_url="my_url"
AC_PREREQ([2.59])
# 5 argument version only available with aclocal >= 2.64
AC_INIT( [idemuxcpp], 
         idemuxCPP_version_major.idemuxCPP_version_minor.idemuxCPP_version_patch, 
         [], 
         [idemuxcpp], [my_url] )

AC_CANONICAL_HOST
AC_CONFIG_AUX_DIR([.])
AC_CONFIG_SRCDIR([src/idemuxCPP.cpp])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_MACRO_DIR([m4])

# backwards compatibility for autoconf >= 2.64
# PACKAGE_URL should be the fifth argument of AC_INIT
m4_ifndef([AC_PACKAGE_URL],
	[m4_define([AC_PACKAGE_URL], [my_url])
	AC_DEFINE(PACKAGE_URL, "AC_PACKAGE_URL", [Package URL])
	AC_SUBST(PACKAGE_URL, AC_PACKAGE_URL)])

lt_enable_auto_import=""
case "$host_os" in 
	cygwin* | mingw* | cegcc*)
		AM_LDFLAGS="-Wl,--enable-auto-import $AM_LDFLAGS"
esac

# Checks for programs

# check for C++ compiler
# store current compiler flags to avoid default setup via AC_PROG_CXX and *_CC
OLD_CXXFLAGS=$CXXFLAGS
OLD_CFLAGS=$CFLAGS
# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

AC_CHECK_PROG([HAVE_PKG_CONFIG], [pkg-config], [yes],[no])

# reset compiler flags to initial flags
CXXFLAGS=$OLD_CXXFLAGS
CFLAGS=$OLD_CFLAGS

# automake initialisation (mandatory) and check for minimal automake API version
AM_INIT_AUTOMAKE([1.9])

# use the C++ compiler for the following checks
AC_LANG([C++])

# include openMP for multithreading.
dnl Check for OpenMP and prepare flags

dnl AX_OPENMP([],[])
AC_OPENMP
AC_SUBST(OPENMP_CXXFLAGS)
AC_SUBST(AM_CXXFLAGS,"$OPENMP_CXXFLAGS")

###############################################################################
###############################################################################

############  PARAMETERS  ########################################

###############################################################################
# SETUP FOR COMPILATION WITH 2011 ISO C++ STANDARD 
###############################################################################

AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])

###############################################################################
# DEBUG SUPPORT SETUP
###############################################################################

AC_MSG_CHECKING([whether to build with debug information])
debuger=no
AC_ARG_ENABLE([debug],
	[AS_HELP_STRING([--enable-debug],
	    [enable debug data generation (def=no)])],
  	[debuger="$enableval"])
AC_MSG_RESULT([$debuger])
if test x"$debuger" = x"yes"; then
	AC_DEFINE([DEBUG], [], [Run in DEBUG mode with additional assertions and debug output])
	AM_CXXFLAGS="$AM_CXXFLAGS -pg -g -Wno-uninitialized -O0 -pthread"
	AM_LDFLAGS="$AM_LDFLAGS -g -Wno-uninitialized -O0"
	#AM_CXXFLAGS="$AM_CXXFLAGS -g -O0 -Wall -Werror -Wno-uninitialized -Wno-deprecated"
else
	AC_DEFINE([NDEBUG], [], [Run in normal mode with minimal assertions])
	AM_CXXFLAGS="$AM_CXXFLAGS -O3 -fno-strict-aliasing -Wno-uninitialized -Wno-deprecated -pthread"
fi

############  CHECKS  ############################################

# check for gzlib, dl, boost etc.
AC_LIB_HAVE_LINKFLAGS([z], [], [#include <zlib.h>], [], [zlib is missing!])
if test "$HAVE_LIBZ" != "yes"; then
AC_MSG_FAILURE([ERROR: please install the zlib library!])
fi
AC_LIB_LINKFLAGS([dl])

# check for bamtools-dev
HAVE_LIBBAMTOOLS=no
AS_IF([test "$HAVE_PKG_CONFIG" = "yes"],[
  PKG_CHECK_MODULES([BAMTOOLS], [bamtools-1 >= 2.5.1], [
  AC_DEFINE([HAVE_LIBBAMTOOLS],[1])
  HAVE_LIBBAMTOOLS=yes
  ], [
       HAVE_LIBBAMTOOLS=no
       AC_DEFINE([HAVE_LIBBAMTOOLS],[0])
     ])
],[])

# add makefile vars
#AC_SUBST([HAVE_LIBBAMTOOLS],[1])
if test "$HAVE_LIBBAMTOOLS" = "yes"; then
  AC_SUBST([LIBBAMTOOLS])
# Warning: the file I need is bamtools/api/BamReader.h, however in the pkg_config it says that it is in ../include/
# the header references <api/BamReader.h> similar to headers referenced within BamReader.h which cannot be found!
# Thus, I add /include/bamtools/
  AC_SUBST([BAMTOOLS_CFLAGS],"${BAMTOOLS_CFLAGS}/bamtools/")
  AC_SUBST([BAMTOOLS_LIBS])
fi

AS_IF([test "$HAVE_PKG_CONFIG" != "yes" -o "$HAVE_LIBBAMTOOLS" != "yes"],[
 AC_CHECK_FILE("/usr/include/bamtools/api/BamReader.h", 
      [  AC_MSG_NOTICE([Found file in path /usr/include/bamtools/api. Test linking with updated CXXFLAGS...])
         CXXFLAGS="$CXXFLAGS -I/usr/include/bamtools/"
      ], 
      [ AC_MSG_NOTICE([Bamtools library not found in /usr/include/bamtools/api/BamReader.h! Let us hope it is in the include path])
 ])
 AC_LIB_HAVE_LINKFLAGS([bamtools], [z], [#include <bamtools/api/BamReader.h>], [], [Linking against bamtools library has failed!])
 if test "$HAVE_LIBBAMTOOLS" != "yes"; then
 AC_MSG_WARN([
********************************************************************************
   Please install the bamtools library!
   (https://github.com/pezmaster31/bamtools)
   If libbamtools-dev is already installed, you have to specify the path prefix 
   with --with-libbamtools-prefix or install it into a standard location!
 
   If you don't need bam files as input (fastq and fastq.gz should work), 
   then you can ignore this message.
********************************************************************************])
fi
],[])

AM_CONDITIONAL([HAVE_LIBBAMTOOLS], [test "$HAVE_LIBBAMTOOLS" != "no"])


###############################################################################
# Boost unit tests check
###############################################################################

AX_BOOST_BASE([1.55], [], AC_MSG_ERROR([Could not find a useful version of boost]))
AC_LIB_HAVE_LINKFLAGS([boost_system], [], [] , [], [Could not link the boost_system library!])
if test "$HAVE_LIBBOOST_SYSTEM" != "yes"; then
AC_MSG_FAILURE([ERROR: please install the boost_system library!])
fi
AC_LIB_HAVE_LINKFLAGS([boost_filesystem], [boost_system], [#include <boost/filesystem/path.hpp>] , [], [Could not link the boost_filesystem library!])
if test "$HAVE_LIBBOOST_FILESYSTEM" != "yes"; then
AC_MSG_FAILURE([ERROR: please install the boost_filesystem library!])
fi
AC_LIB_HAVE_LINKFLAGS([boost_iostreams], [], [] , [], [Could not link the boost_iostreams library!])
if test "$HAVE_LIBBOOST_IOSTREAMS" != "yes"; then
AC_MSG_FAILURE([ERROR: please install the boost_iostreams library!])
fi

AC_ARG_ENABLE([tests],
        [AS_HELP_STRING([--enable-tests],
            [enable boost unit tests (def=no)])],
        [tests="$enableval"])
AC_MSG_RESULT([$tests])
if test x"$tests" = x"yes"; then
  AC_LIB_HAVE_LINKFLAGS([boost_unit_test_framework], [], [#include <boost/test/unit_test.hpp>] , [using boost::unit_test::test_suite;
                                                         test_suite* test= BOOST_TEST_SUITE( "Unit test example 1" ); if (test == NULL) { return 1; } else { return 0; }], [Could not link the boost_test library!])
  if test "$HAVE_LIBBOOST_UNIT_TEST_FRAMEWORK" != "yes"; then
    AC_MSG_FAILURE([ERROR: please install the boost_test library!])
  else
    # lrt is requred by boost_unit_tests for clock_* in glibc < 2.1.7
    LIBS="$LIBS -lrt"
  fi
fi

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for header files.
AC_HEADER_STDC


###############################################################################
# help2man check
###############################################################################
AX_PROG_HELP2MAN([src/idemuxCPP])


# distribute additional compiler and linker flags
# --> set these variables instead of CXXFLAGS or LDFLAGS
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([LIBS])

# files to generate via autotools (.am or .in source files)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([doc/Makefile])
AC_CONFIG_FILES([misc/Makefile])
AC_CONFIG_FILES([tests/Makefile])

# generate the final Makefile etc.
AC_OUTPUT([packaging/idemuxcpp.spec packaging/PKGBUILD])

