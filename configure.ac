#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59)
AC_INIT(pimc++, 0.5, esler@uiuc.edu, pimc++)
AM_INIT_AUTOMAKE(pimc, 0.5)
AC_CONFIG_SRCDIR([src/PathClass.h])
AC_CONFIG_HEADER([config.h])

# Check which compilers we're using
AX_CXX_MAXOPT
AC_MSG_CHECKING([optimal compiler flags])
AC_MSG_RESULT([${CXXFLAGS}])

AC_CHECK_LIB([m], [sqrt])
AC_MSG_CHECKING([LIBS = $LIBS])
# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h sys/time.h unistd.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_F77_LIBRARY_LDFLAGS
AC_F77_MAIN
AC_PROG_LIBTOOL
PKG_PROG_PKG_CONFIG


AC_CANONICAL_HOST
AC_MSG_CHECKING([host cpu = ])
AC_MSG_RESULT([${host_cpu}])
case "${host_cpu}" in
   "powerpc")
      CPPFLAGS="-DMAC $CPPFLAGS"
      CXXFLAGS="$CXX_FLAGS -O3 -mcpu=970 -mtune=970 -mpowerpc64 -mpowerpc-gpopt -ffast-math -g"
      FFLAGS="-O3 -mcpu=970 -mtune=970 -mpowerpc64 -mpowerpc-gpopt -ffast-math -g";;
   "i686")
      AC_MSG_CHECKING([ i686 cpu type])
      if grep -q XP /proc/cpuinfo; then
        AC_MSG_RESULT([Athlon XP])
        AX_CC_OPTION(TUNEOPT, [-mtune=athlon-xp], [-mtune=athlon-xp], [-mcpu=athlon-xp])
        AX_CC_OPTION(ARCHOPT, [-march=athlon-xp], [-march=athlon-xp], [-march=athlon])
	FFLAG="-O3 -march=athlon -ffast-math -msse -m3dnow -g"
      fi
      if grep -q "Pentium(R) M" /proc/cpuinfo; then
	AC_MSG_RESULT([Pentium M])
#        AX_CC_OPTION(TUNEOPT, [-mtune=pentium-m], [-mtune=pentium-m], [-mcpu=pentium-m])
#        AX_CC_OPTION(ARCHOPT, [-march=pentium-m], [-march=pentium-m], [-march=pentium-m])
#	FFLAGS="-O3 $TUNEOPT $ARCHOPT -ffast-math -g"
        AX_CC_OPTION(TUNEOPT, [-mtune=pentium-m], [-mtune=pentium-m], [-mcpu=pentium4])
	AX_CC_OPTION(ARCHOPT, [-march=pentium-m], [-march=pentium-m], [-mcpu=pentium4])
        AX_F77_OPTION(FARCHOPT, [-march=pentium4], [-march=pentium4], [ ])
	FFLAGS="-O3 $FARCHOPT -ffast-math -g"
      fi
      if grep -q "Intel(R) Xeon" /proc/cpuinfo; then
	AC_MSG_RESULT([Intel Xeon])
        AX_CC_OPTION(TUNEOPT, [-tpp7], [-tpp7], [-tpp7])
        AX_CC_OPTION(ARCHOPT, [-xW], [-xW], [-xW])
	FFLAGS="-O3 $TUNEOPT $ARCHOPT -g"
      fi
      CXXFLAGS="$CXXFLAGS -O3 $TUNEOPT $ARCHOPT -ffast-math -g"
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
      
    "ia64")
       AC_MSG_CHECKING([ Itanium cpu type])
       if grep -q "Itanium 2" /proc/cpuinfo; then
         AC_MSG_RESULT([Intel Itanium 2])
         AX_CC_OPTION(TUNEOPT, [-mcpu=itanium2], [-mcpu=itanium2], [-mcpu=itanium2])
       fi
       CXXFLAGS="$CXXFLAGS -O3 $TUNEOPT -g"
       FFLAGS="-O3 $TUNEOPT -g";;
esac

CXXFLAGS="$CXXFLAGS -g -pthread"
FFLAGS="$FFLAGS -g"

# Checks if we want MPI support
AC_ARG_ENABLE(mpi,    [  --enable-mpi            enable MPI support (default=no)])
AM_CONDITIONAL([PARALLEL], [test "$enable_mpi" = "yes"])
# If so, use MPI version of libcommon.  Otherwise, use serial version
if test "x$enable_mpi" != "xyes"; then
  PKG_CHECK_MODULES(COMMON, common, common_ok=yes, common_ok=no)
  AC_MSG_RESULT([Linking against serial version of libcommon.])
else
  PKG_CHECK_MODULES(COMMON, mpicommon, common_ok=yes, common_ok=no)
  AC_MSG_RESULT([Linking against parallel version of libmpicommon.])
fi

if test "x$common_ok" != "xyes"; then
   AC_MSG_ERROR([Could not find common library.  Please install first.])
fi

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
  PKG_CHECK_MODULES(GTHREAD, gthread-2.0, gthread_ok=yes, gthread_ok=no)
  if test "x$gthread_ok" != "xyes"; then
     AC_MSG_ERROR([gthread libraries missing.  Please install.])
  fi
#  PKG_CHECK_MODULES(GTS, gts, gts_ok=yes, gts_ok=no)
#  if test "x$gts_ok" != "xyes"; then
#     AC_MSG_ERROR([GTS libraries missing.  Please install. (http://gts.sourceforge.net)])
#  fi
  AC_CHECK_LIB([xvidcore],[xvid_decore])
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  AC_CHECK_LIB([revel],[Revel_EncodeStart])
  AC_LANG_RESTORE
  AC_CHECK_LIB([glut],[glutSolidTeapot])
  AC_CHECK_LIB([gle],[gleSetNumSides])
fi

AC_ARG_ENABLE(profile,[  --enable-profile        instrument code with profiling information])
if test "x$enable_profile" = "xyes"; then
   CXXFLAGS="-pg $CXXFLAGS"
   FFLAGS="-pg $FFLAGS"
   LDFLAGS="-pg $LDFLAGS"
fi

AC_ARG_ENABLE(debug,[  --enable-debug          enable code for debugging checks])
if test "x$enable_debug" = "xyes"; then
   CXXFLAGS="-DDEBUG -DBZ_DEBUG $CXXFLAGS"
   FFLAGS="-DDEBUG $FFLAGS"
fi

AC_ARG_ENABLE(open, [  --enable-open           enable code for open loops])
if test "x$enable_open" = "xyes"; then
   CXXFLAGS="-DOPEN_LOOPS $CXXFLAGS"
   FFLAGS="-DOPEN_LOOPS $CXXFLAGS"
fi

# Optionally set dimensionality
AC_ARG_WITH(ndim, [  --with-ndim=N           specify simulation dimensionality])
if (test "x$with_ndim" != "x"); then
   CPPFLAGS="-DNDIM=$with_ndim $CPPFLAGS"
else # Default to 3D
   CPPFLAGS="-DNDIM=3 $CPPFLAGS"
fi



if test "$enable_mpi" = "yes"; then
    MPI_DEF="-DUSE_MPI"
else
    MPI_DEF=" "
fi
AC_SUBST(MPI_DEF)

AC_ARG_ENABLE(open, [  --enable-on-fermions    enable code for order-n fermions])
if test "x$enable_on_fermions" = "xyes"; then
   AC_LANG_SAVE
   AC_LANG_FORTRAN77
   AC_CHECK_LIB([dctmp], [cg], [HAVE_DCTMP="yes"], [HAVE_DCTMP="no"], [-lf2c -lm -lgsl -lgslcblas -llapack -lblas])
   if test "x$HAVE_DCTMP" = "xyes"; then
      ORDER_N_LIBS="-ldctmp"
   fi
   if test "x$HAVE_DCTMP" = "xno"; then
      AC_MSG_ERROR([  Cannot find libdctmp.  Aborting.])
   fi
   AC_LANG_RESTORE
   AC_SUBST(ORDER_N_LIBS)      
   ORDER_N_LIBS="-ldctmp"
   CXXFLAGS="-DORDER_N_FERMIONS $CXXFLAGS"
fi

# Checks if we want qmcPACK functionality for CEIMC
AC_ARG_WITH(qmc,    [  --with-qmcpack=/path/to/my/qmcpack/builddir	enable QMCPACK functionality])
if (test "x$with_qmcpack" != "x"); then
   	QMC_DEF="-DUSE_QMC -DAPPNAMESPACE=qmcplusplus -I$with_qmcpack/src -I$with_qmcpack/../src"
	 	QMC_LIB="$with_qmcpack/lib/libqmc.a $with_qmcpack/lib/libqmcdriver.a $with_qmcpack/lib/libqmcham.a $with_qmcpack/lib/libqmcwfs.a $with_qmcpack/lib/libqmcbase.a $with_qmcpack/lib/libqmcutil.a"
else # don't link with qmcpack
    QMC_DEF=" "
		QMC_LIB=" "
fi

AC_ARG_WITH(boost,    [  --with-boost=/path/to/my/boost	include boost libraries; needed for qmcpack])
if (test "x$with_boost" != "x"); then
	 	BOOST_LIB="-I$with_boost/include -L$with_boost/lib"
else # don't link with boost
		BOOST_LIB=" "
fi

AC_ARG_WITH(libxml,    [  --with-libxml=/path/to/my/libxml2  include libxml2 libraries; needed for qmcpack])
if (test "x$with_libxml" != "x"); then
	 	LIBXML_LIB="-I$with_libxml/include/libxml2 $with_libxml/lib/libxml2.so"
else # don't link with libxml
		LIBXML_LIB=" "
fi

AC_SUBST(QMC_DEF)
AC_SUBST(QMC_LIB)
AC_SUBST(BOOST_LIB)
AC_SUBST(LIBXML_LIB)

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
