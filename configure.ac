#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59)
AC_INIT(pimc, 1.0, esler@uiuc.edu)
AM_INIT_AUTOMAKE(pimc, 1.0)
AC_CONFIG_SRCDIR([src/PathClass.h])
AC_CONFIG_HEADER([config.h])

AC_CANONICAL_HOST
if test "x$host" = "xpowerpc-apple-darwin7.7.0"; then
   CPPFLAGS="-DMAC $CPPFLAGS"
fi

PKG_CHECK_MODULES(COMMON, common, common_ok=yes, common_ok=no)

AC_CHECK_LIB([m], [sqrt])
AC_MSG_CHECKING([LIBS = $LIBS])
# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h sys/time.h unistd.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LIBTOOL

# Optionally compile visual
AC_ARG_ENABLE(pathvis, [  --enable-pathvis         compile path visualization GUI])
AM_CONDITIONAL([WANT_PATHVIS], [test "$enable_pathvis" = "yes"])
if test "x$enable_pathvis" = "xyes"; then
  PKG_CHECK_MODULES(GTKMM, gtkmm-2.4, gtkmm_ok=yes, gtkmm_ok=no)
  if test "x$gtkmm_ok" != "xyes"; then
     AC_MSG_ERROR([gtkmm libraries missing.  Please install.])
  fi
  PKG_CHECK_MODULES(GTKGLEXTMM, gtkglextmm-1.2, gtkglextmm_ok=yes, gtglextkmm_ok=no)
  if test "x$gtkglextmm_ok" != "xyes"; then
     AC_MSG_ERROR([gtkglextmm libraries missing.  Please install.])
  fi
  PKG_CHECK_MODULES(GTS, gts, gts_ok=yes, gts_ok=no)
  if test "x$gts_ok" != "xyes"; then
     AC_MSG_ERROR([GTS libraries missing.  Please install. (http://gts.sourceforge.net)])
  fi
  AC_CHECK_LIB([xvidcore],[decoder_create])
  AC_CHECK_LIB([revel],[main])
  AC_CHECK_LIB([glut],[glutSolidTeapot])
  AC_CHECK_LIB([gle],[gleSetNumSides])
fi


# Optionally set dimensionality
AC_ARG_WITH(ndim, [  --with-ndim=N            specify simulation dimensionality])
if (test "x$with_ndim" != "x"); then
   CPPFLAGS="-DNDIM=$with_ndim $CPPFLAGS"
else # Default to 3D
   CPPFLAGS="-DNDIM=3 $CPPFLAGS"
fi

# Checks if we want MPI support
AC_ARG_ENABLE(mpi,    [  --enable-mpi              enable MPI support (default=no)])

if test "$enable_mpi" = "yes"; then
    MPI_DEF="-DUSE_MPI"
else
    MPI_DEF=" "
fi
AC_SUBST(MPI_DEF)
AC_F77_WRAPPERS

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_STRTOD
AC_CHECK_FUNCS([clock_gettime floor pow sqrt strtol])

AC_SUBST(PKGDATADEF)

AC_CONFIG_FILES([ Makefile                  \
		  src/Makefile              \
                  src/Actions/Makefile      \
                  src/Moves/Makefile        \
                  src/Observables/Makefile  \
                  src/Visual/Makefile       \
                  doc/Makefile               ])
AC_MSG_RESULT([ COMMON_LIBS  = $COMMON_LIBS ])
AC_OUTPUT
