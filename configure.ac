#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(libcommon, 1.0, esler@uiuc.edu)
AM_INIT_AUTOMAKE(libcommon, 1.0)
AC_CONFIG_SRCDIR([Blitz.h])
AC_CONFIG_HEADER([config.h])

# Checks for libraries.
PKG_CHECK_MODULES(BLITZ, blitz, blitz_ok=yes, blitz_ok=no)
PKG_CHECK_MODULES(XML2, xml2, xml2_ok=yes, xml2_ok=no)
AC_CHECK_LIB([m], [sqrt])
AC_CHECK_LIB([hdf5], [H5open])
# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h sys/time.h unistd.h])
AC_SUBST(BLITZ_LIB)
AC_SUBST(BLITZ_CFLAGS)
AC_SUBST(XML2_LIB)
AC_SUBST(XML2_CFLAGS)


# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LIBTOOL

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_STRTOD
AC_CHECK_FUNCS([clock_gettime floor pow sqrt strtol])

AC_CONFIG_FILES([ Makefile IO/Makefile ])
AC_OUTPUT
