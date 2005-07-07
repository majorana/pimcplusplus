#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(libcommon, 1.0, esler@uiuc.edu)
AM_INIT_AUTOMAKE(libcommon, 1.0)
AC_CONFIG_SRCDIR([Blitz.h])
AC_CONFIG_HEADER([config.h])

AC_CANONICAL_HOST
if test "x$host" = "xpowerpc-apple-darwin7.7.0"; then
   CPPFLAGS="-DMAC $CPPFLAGS"
fi

# Checks for libraries.
# this creates and AC_SUBST's vars BLITZ_CFLAGS and BLITZ_LIBS
PKG_CHECK_MODULES(BLITZ, blitz, blitz_ok=yes, blitz_ok=no)
# this creates and AC_SUBST's vars XML2_CFLAGS and XML2_LIBS
PKG_CHECK_MODULES(XML2, xml2, xml2_ok=yes, xml2_ok=no)
# this creates and AC_SUBST's vars FFTW3_CFLAGS and FFTW3_LIBS
PKG_CHECK_MODULES(FFTW3, fftw3, fftw3_ok=yes, fftw3_ok=no)
# this creates and AC_SUBST's vars SPRNG2_CFLAGS and SPRNG2_LIBS
PKG_CHECK_MODULES(SPRNG2, sprng2, sprng2_ok=yes, sprng2_ok=no)

AC_CHECK_LIB([m], [sqrt])
AC_CHECK_LIB([hdf5], [H5open],[hdf5_ok=yes], hdf5_ok=no)
if test "x$hdf5_ok" = "xyes"; then
   LIBS="-lhdf5 $LIBS"
fi
AC_CHECK_LIB([cblas],[cblas_dgemm], cblas_ok=yes, cblas_ok=no)
if test "x$cblas_okay" = "xyes"; then
   LIBS="-lcblas $LIBS"
fi
AC_CHECK_LIB([g2c],[do_lio])
AM_PATH_GSL(1.4,gsl_ok=yes,gsl_ok=no)
AC_CHECK_LIB([lapack], [dgesvd_])
AC_CHECK_LIB([blas],[daxp])
AC_MSG_CHECKING([LIBS = $LIBS])
# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h sys/time.h unistd.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LIBTOOL

# Checks if we want MPI support
AC_ARG_ENABLE(mpi,    [  --enable-mpi              enable MPI support (default=no)])
AC_ARG_WITH(mpi-inc,  [  --with-mpi-cflags=FLAGS   specify MPI compliler flags])
AC_ARG_WITH(mpi-libs, [  --with-mpi-libs=LIBS      specify MPI linker flags])
AC_ARG_WITH(mpi-run,  [  --with-mpi-run=cmd        specify MPI run command (default=mpirun)])

if test "$enable_mpi" = "yes"; then
    MPI_DEF="-DUSE_MPI"
    if test "X$with_mpi_libs" = "X"; then
        MPI_LIBS="-lmpi"
    else
        MPI_LIBS="$with_mpi_libs"
    fi
  
    if test "X$with_mpi_cflags" = "X"; then
        MPI_CFLAGS=""
    else
        MPI_CFLAGS="$with_mpi_cflags"
    fi

    if test "X$with_mpi_run" = "X"; then
        MPI_RUN="mpirun"
    else
        MPI_RUN="$with_mpi_run"
    fi
else
    MPI_DEF=" "
fi
AC_SUBST(MPI_DEF)
AC_SUBST(MPI_LIBS)
AC_SUBST(MPI_CFLAGS)
AC_SUBST(MPI_RUN)


AC_ARG_WITH(hdf5-cflags,[  --with-hdf5-cflags=FLAGS specify HDF5 linker flags])
AC_ARG_WITH(hdf5-libs, [  --with-hdf5-libs=LIBS    specify HDF5 linker flags])
if test "x$with_hdf5_cflags" != "x"; then
   HDF5_CFLAGS="$with_hdf5_cflags"
fi
if test "x$with_hdf5_libs" != "x"; then
   HDF5_LIBS="$with_hdf5_libs"
fi
AC_SUBST(HDF5_CFLAGS)
AC_SUBST(HDF5_LIBS)

if test "x$hdf5_ok" = "xno"; then
   AC_MSG_WARN([Didn't find libhdf5 in standard place.])
   if test "x$with_hdf5_libs" = "x"; then
       AC_MSG_ERROR([HDF5 library could not be found.  Please specify location with --with-hdf5-libs=])
   fi
fi

AC_ARG_WITH(cblas-cflags,[  --with-cblas-cflags=FLAGS specify CBLAS linker flags])
AC_ARG_WITH(cblas-libs, [  --with-cblas-libs=LIBS    specify CBLAS linker flags])
if test "x$with_cblas_cflags" != "x"; then
   CBLAS_CFLAGS="$with_cblas_cflags"
fi
if test "x$with_cblas_libs" != "x"; then
   CBLAS_LIBS="$with_cblas_libs"
fi
AC_SUBST(CBLAS_CFLAGS)
AC_SUBST(CBLAS_LIBS)
if test "x$cblas_ok" = "xno"; then
   AC_MSG_WARN([Didn't find libcblas in standard place.])
   if test "x$with_cblas_libs" = "x"; then
       AC_MSG_ERROR([CBLAS library could not be found.  Please specify location with --with-cblas-libs=])
   fi
fi

AC_ARG_WITH(xml2-cflags,[  --with-xml2-cflags=FLAGS specify XML2 linker flags])
AC_ARG_WITH(xml2-libs, [  --with-xml2-libs=LIBS    specify XML2 linker flags])
if test "x$with_xml2_cflags" != "x"; then
   XML2_CFLAGS="$with_xml2_cflags"
fi
if test "x$with_xml2_libs" != "x"; then
   XML2_LIBS="$with_xml2_libs"
fi
AC_SUBST(XML2_CFLAGS)
AC_SUBST(XML2_LIBS)
if test "x$xml2_ok" = "xno"; then
   AC_MSG_WARN([Didn't find libxml2 in standard place.])
   if test "x$with_xml2_libs" = "x"; then
       AC_MSG_ERROR([XML2 library could not be found.  Please specify location with --with-xml2-libs=])
   fi
fi



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

AC_CONFIG_FILES([ Makefile                  \
                  common.pc                 \
                  IO/Makefile               \
                  Splines/Makefile          \
                  SpecialFunctions/Makefile \
                  Random/Makefile           \
                  Fitting/Makefile          \
                  Atom/Makefile             \
                  DFT/Makefile              \
                  PH/Makefile               \
                  Integration/Makefile      \
                  MatrixOps/Makefile        \
                  MPI/Makefile              \
                  Distributed/Makefile      \
                  Ewald/Makefile            \
                  FFT/Makefile              \
                  PairAction/Makefile       \
                  PlaneWavePHDFT/Makefile ] )
AC_MSG_RESULT([ HDF5_LIBS = $HDF5_LIBS ])
AC_MSG_RESULT([ XML2_LIBS = $XML2_LIBS ])
AC_OUTPUT
