#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(libcommon, 1.0, esler@uiuc.edu)
AM_INIT_AUTOMAKE(libcommon, 1.0)
AC_CONFIG_SRCDIR([Blitz.h])
AC_CONFIG_HEADER([config.h])

# Checks for libraries.
# this creates and AC_SUBST's vars BLITZ_CFLAGS and BLITZ_LIBS
PKG_CHECK_MODULES(BLITZ, blitz, blitz_ok=yes, blitz_ok=no)
# this creates and AC_SUBST's vars XML2_CFLAGS and XML2_LIBS
PKG_CHECK_MODULES(XML2, xml2, xml2_ok=yes, xml2_ok=no)
# this creates and AC_SUBST's vars FFTW3_CFLAGS and FFTW3_LIBS
PKG_CHECK_MODULES(FFTW3, fftw3, fftw3_ok=yes, fftw3_ok=no)

AC_CHECK_LIB([m], [sqrt])
AC_CHECK_LIB([hdf5], [H5open])
AC_CHECK_LIB([g2c],[do_lio])
AM_PATH_GSL(1.4,gsl_ok=yes,gsl_ok=no)
AC_CHECK_LIB([lapack], [dgesvd_])
AC_CHECK_LIB([blas],[daxp])
# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h sys/time.h unistd.h])


# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LIBTOOL

# Checks if we want MPI support
AC_ARG_ENABLE(mpi,    [  --enable-mpi            enable MPI support (default=no)])
AC_ARG_WITH(mpi-inc,  [  --with-mpi-cflags=FLAGS specify MPI compliler flags])
AC_ARG_WITH(mpi-libs, [  --with-mpi-libs=LIBS    specify MPI linker flags])
AC_ARG_WITH(mpi-run,  [  --with-mpi-run=cmd      specify MPI run command (default=mpirun)])

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


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_STRTOD
AC_CHECK_FUNCS([clock_gettime floor pow sqrt strtol])

AC_CONFIG_FILES([ Makefile IO/Makefile      \
                  Splines/Makefile          \
                  SpecialFunctions/Makefile \
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
AC_OUTPUT
