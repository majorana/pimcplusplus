#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59)
AC_INIT(pimc, 1.0, esler@uiuc.edu)
AM_INIT_AUTOMAKE(pimc, 1.0)
AC_CONFIG_SRCDIR([src/PathClass.h])
AC_CONFIG_HEADER([config.h])



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


AC_CANONICAL_HOST
AC_MSG_CHECKING([host cpu = ])
AC_MSG_RESULT([${host_cpu}])
case "${host_cpu}" in
   "powerpc")
      CPPFLAGS="-DMAC $CPPFLAGS"
      CXXFLAGS="-O3 -mcpu=970 -mtune=970 -mpowerpc64 -mpowerpc-gpopt -ffast-math -g"
      F77FLAGS="-O3 -mcpu=970 -mtune=970 -mpowerpc64 -mpowerpc-gpopt -ffast-math -g";;
   "i686")
      AC_MSG_CHECKING([ i686 cpu type])
      if grep -q XP /proc/cpuinfo; then
        AC_MSG_RESULT([Athlon XP])
        AX_CC_OPTION(TUNEOPT, [-mtune=athlon-xp], [-mtune=athlon-xp], [-mcpu=athlon-xp])
        AX_CC_OPTION(ARCHOPT, [-march=athlon-xp], [-march=athlon-xp], [-march=athlon])
      fi
      if grep -q "Pentium(R) M" /proc/cpuinfo; then
	AC_MSG_RESULT([Pentium M])
        AX_CC_OPTION(TUNEOPT, [-mtune=pentium-m], [-mtune=pentium-m], [-mcpu=pentium-m])
        AX_CC_OPTION(ARCHOPT, [-march=pentium-m], [-march=pentium-m], [-march=pentium-m])
      fi
      if grep -q "Intel(R) Xeon" /proc/cpuinfo; then
	AC_MSG_RESULT([Intel Xeon])
        AX_CC_OPTION(TUNEOPT, [-tpp7], [-tpp7], [-tpp7])
        AX_CC_OPTION(ARCHOPT, [-xW], [-xW], [-xW])
      fi
      CXXFLAGS="-O3 $TUNEOPT $ARCHOPT -ffast-math -g"
      F77FLAGS="-O3 -mtune=athlon -march=athlon -ffast-math -msse -m3dnow -g"
      if grep -q sse /proc/cpuinfo ; then
        CXXFLAGS="$CXXFLAGS -msse"
      fi
      if grep -q sse2 /proc/cpuinfo ; then
        CXXFLAGS="$CXXFLAGS -msse2"
      fi
      if grep -q sse3 /proc/cpuinfo ; then
        CXXFLAGS="$CXXFLAGS -msse3"
      fi
      if grep -q 3dnow /proc/cpuinfo ; then
        CXXFLAGS="$CXXFLAGS -m3dnow"
      fi;;
esac


# Optionally compile visual
AC_ARG_ENABLE(pathvis, [  --enable-pathvis        compile path visualization GUI])
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
#  PKG_CHECK_MODULES(GTS, gts, gts_ok=yes, gts_ok=no)
#  if test "x$gts_ok" != "xyes"; then
#     AC_MSG_ERROR([GTS libraries missing.  Please install. (http://gts.sourceforge.net)])
#  fi
  AC_CHECK_LIB([xvidcore],[decoder_create])
  AC_CHECK_LIB([revel],[main])
  AC_CHECK_LIB([glut],[glutSolidTeapot])
  AC_CHECK_LIB([gle],[gleSetNumSides])
fi

AC_ARG_ENABLE(profile,[  --enable-profile        instrument code with profiling information])
if test "x$enable_profile" = "xyes"; then
   CXXFLAGS="-pg $CXXFLAGS"
   F77FLAGS="-pg $F77FLAGS"
   LDFLAGS="-pg $LDFLAGS"
fi

AC_ARG_ENABLE(debug,[  --enable-debug          enable code for debugging checks])
if test "x$enable_debug" = "xyes"; then
   CXXFLAGS="-DDEBUG $CXXFLAGS"
   F77FLAGS="-DDEBUG $F77FLAGS"
fi

AC_ARG_ENABLE(open,[ --enable-open            enable code for open loops])
if test "x$enable_open" = "xyes"; then
   CXXFLAGS="-DOPEN_LOOPS $CXXFLAGS"
   F77FLAGS="-DOPEN_LOOPS $CXXFLAGS"
fi

AC_ARG_ENABLE(open,[ --enable-on-fermions  enable code for order-n fermions])
if test "x$enable_on_fermions" = "xyes"; then
   AC_CHECK_LIB([dctmp], [cg])
   CXXFLAGS="-DORDER_N_FERMIONS $CXXFLAGS"
fi

# Optionally set dimensionality
AC_ARG_WITH(ndim, [  --with-ndim=N           specify simulation dimensionality])
if (test "x$with_ndim" != "x"); then
   CPPFLAGS="-DNDIM=$with_ndim $CPPFLAGS"
else # Default to 3D
   CPPFLAGS="-DNDIM=3 $CPPFLAGS"
fi

# Checks if we want MPI support
AC_ARG_ENABLE(mpi,    [  --enable-mpi            enable MPI support (default=no)])

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
                  src/Langevin/Makefile     \
                  doc/Makefile               ])
AC_MSG_RESULT([ COMMON_LIBS  = $COMMON_LIBS ])
AC_OUTPUT
