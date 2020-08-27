
m4_define([idemuxCPP_version_major], [1])
m4_define([idemuxCPP_version_minor], [2])
m4_define([idemuxCPP_version_patch], [0])

my_url="my_url"
AC_PREREQ([2.59])
# 5 argument version only available with aclocal >= 2.64
AC_INIT( [idemuxCPP], 
         idemuxCPP_version_major.idemuxCPP_version_minor.idemuxCPP_version_patch, 
         [], 
         [idemuxCPP], [my_url] )

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
#AX_OPENMP([],[]) # not used

#AC_SUBST(OPENMP_CXXFLAGS)
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

# todo: ax_check for gzlib

AC_LIB_LINKFLAGS([z])
#AM_LDFLAGS="$AM_LDFLAGS -lz"

############  CHECKS  ############################################

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

# generate the final Makefile etc.
AC_OUTPUT([packaging/idemuxCPP.spec packaging/PKGBUILD])

