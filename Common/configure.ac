#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.


##add m4_include([gsl.m4])
AC_PREREQ(2.59)
AC_INIT(libcommon, 0.5, esler@uiuc.edu)
AM_INIT_AUTOMAKE(libcommon, 0.5)
AC_CONFIG_SRCDIR([Blitz.h])
AC_CONFIG_HEADER([config.h])

# Check which compilers we're using
AX_CXX_MAXOPT
AX_F77_MAXOPT
AX_CC_MAXOPT

AC_MSG_CHECKING([optimal compiler flags])
AC_MSG_RESULT([${CXXFLAGS}])
#CXXFLAGS="$CXXFLAGS -pthread"
CFLAGS="$CFLAGS -std=c99 -g"


# Checks for libraries.
# this creates and AC_SUBST's vars FFTW3_CFLAGS and FFTW3_LIBS

# this creates and AC_SUBST's vars BLITZ_CFLAGS and BLITZ_LIBS
PKG_CHECK_MODULES(BLITZ, blitz, blitz_ok=yes, blitz_ok=no)
# this creates and AC_SUBST's vars XML2_CFLAGS and XML2_LIBS
# PKG_CHECK_MODULES(XML2, xml2, xml2_ok=yes, xml2_ok=no)
# this creates and AC_SUBST's vars SPRNG2_CFLAGS and SPRNG2_LIBS
PKG_CHECK_MODULES(SPRNG2, sprng2, sprng2_ok=yes, sprng2_ok=no)
# this creates and AC_SUBST's vars GSL_CFLAGS and GSL_LIBS
PKG_CHECK_MODULES(GSL, gsl, gsl_ok=yes, gsl_ok=no)

AC_CHECK_LIB([m], [sqrt])
AC_CHECK_LIB([hdf5], [H5open],[hdf5_ok=yes], hdf5_ok=no)
if test "x$hdf5_ok" = "xyes"; then
   LIBS="-lhdf5 $LIBS"
fi
#AC_CHECK_LIB([cblas],[cblas_dgemm], cblas_ok=yes, cblas_ok=no)
#if test "x$cblas_okay" = "xyes"; then
#   LIBS="-lcblas $LIBS"
#fi
#AC_CHECK_LIB([g2c],[do_lio],[],[],[frtbegin])
#AC_CHECK_LIB([gfortran],[do_lio])

#if test "x$gsl_ok" = "xno"; then
#   AM_PATH_GSL(1.4,gsl_ok=yes,gsl_ok=no)
#fi

AC_MSG_CHECKING([LIBS = $LIBS])
# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h sys/time.h unistd.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC_STDC
AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
AC_F77_MAIN
AC_PROG_LIBTOOL
AC_C_RESTRICT

# Optimization settings.
AC_CANONICAL_HOST
AC_MSG_CHECKING([host cpu = ])
AC_MSG_RESULT([${host_cpu}])
case "${host_cpu}" in
   "powerpc")
      CPPFLAGS="-DMAC $CPPFLAGS"
      CXXFLAGS="-O3 -mcpu=970 -mtune=970 -mpowerpc64 -mpowerpc-gpopt -ffast-math -g"
      FFLAGS="-O3 -mcpu=970 -mtune=970 -mpowerpc64 -mpowerpc-gpopt -ffast-math -g";;

   "i686")
      AC_MSG_CHECKING([ i686 cpu type])
      if grep -q XP /proc/cpuinfo; then
        AC_MSG_RESULT([Athlon XP])
        AX_CC_OPTION(TUNEOPT, [-mtune=athlon-xp], [-mtune=athlon-xp], [-mcpu=athlon-xp])
        AX_CC_OPTION(ARCHOPT, [-march=athlon-xp], [-march=athlon-xp], [-march=athlon])
	FFLAGS="-O3 -march=athlon -ffast-math -msse -m3dnow -g"
      fi
      if grep -q "Pentium(R) M" /proc/cpuinfo; then
	AC_MSG_RESULT([Pentium M])
#        AX_CC_OPTION(TUNEOPT, [-mcpu=pentium4], [-mcpu=pentium4], [-mcpu=pentium4])
#        AX_CC_OPTION(ARCHOPT, [-march=axB], [-maxB], [-axB])
#        AX_F77_OPTION(FARCHOPT, [-arch=SSE2], [-arch=SSE2], [-arch=SSE2])
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

CXXFLAGS="$CXXFLAGS -g"

# Check python stuff for IO wrappers
AC_ARG_ENABLE(pthread, AC_HELP_STRING([--enable-pthread], [compile with -pthread]))
if test "x$enable_pthread" = "xyes"; then
   CXXFLAGS="$CXXFLAGS -pthread"
   PTHREAD_FLAG="-pthread"
else
   PTHREAD_FLAG=""
fi
AC_SUBST(PTHREAD_FLAG)

AC_ARG_ENABLE(python,      [  --enable-python         enable IO python wrappers (default=yes)], ,enable_python=yes)
AM_CONDITIONAL([WANT_PYTHON_WRAPPERS], [test "$enable_python" = "yes"])
if test "x$enable_python" = "xyes"; then
  AM_PATH_PYTHON(2.3)
  AM_CHECK_PYTHON_HEADERS(,[AC_MSG_ERROR(could not find Python headers)])
  AX_PYTHON_RUN([[import numpy ; print numpy.get_include()]])
  if test $ax_python_cc -ne 0
  then
    AC_MSG_ERROR([cannot find numpy header files])
  fi
  AC_MSG_RESULT([$ax_python_output])
  AC_SUBST([NUMPY_INC], [-I$ax_python_output])
fi


# Check if we want to use the fftw2mkl wrappers
AC_ARG_ENABLE(fftw2mkl,    [  --enable-fftw2mkl       enable MKL FFTW wrappers (default=no)])
if test "x$enable_fftw2mkl" = "xyes"; then
   FFT_DEF="-DUSE_FFTW2MKL"
   PKG_CHECK_MODULES(FFTW3, fftw2mkl, fftw3_ok=yes, fftw3_ok=no)	   
else
   FFT_DEF=" "
   PKG_CHECK_MODULES(FFTW3, fftw3, fftw3_ok=yes, fftw3_ok=no)
fi

# Check if we want to use single-precision FFTs for SIMD spead improvements
AC_ARG_ENABLE(fft-single,    [  --enable-fft-single      enable single-precision FFTs (default=no)])
if test "$enable_fft_single" = "yes"; then
    FFT_DEF="$FFT_DEF -DFFT_USE_SINGLE"
    if test "$enable-fftw2mkl" != "yes"; then
      FFTW3_LIBS="$FFTW3_LIBS -lfftw3f"
    fi
else
    FFT_DEF="$FFT_DEF "
fi
AC_SUBST(FFT_DEF)

# Check if we have specified FFTW location
AC_ARG_WITH(fftw3-cflags,[  --with-fftw3-cflags=FLAGS specify FFTW3 linker flags])
AC_ARG_WITH(fftw3-libs, [  --with-fftw3-libs=LIBS    specify FFTW3 linker flags])
if test "x$with_fftw3_cflags" != "x"; then
   FFTW3_CFLAGS="$with_fftw3_cflags"
fi
if test "x$with_fftw3_libs" != "x"; then
   FFTW3_LIBS="$with_fftw3_libs"
fi
AC_SUBST(FFTW3_CFLAGS)
AC_SUBST(FFTW3_LIBS)

if test "x$fftw3_ok" = "xno"; then
   AC_MSG_WARN([Didn't find libfftw3 in standard place.])
   if test "x$with_fftw3_libs" = "x"; then
       AC_MSG_ERROR([FFTW3 library could not be found.  Please specify location with --with-fftw3-libs=])
   fi
fi


# Check if we want static-only binaries
AC_ARG_ENABLE(all-static, AC_HELP_STRING([--enable-all-static], [build static binaries (default=no)]))
if test "$enable_all_static" = "yes"; then
  ALL_STATIC="-all-static"
else
  ALL_STATIC=""
fi
AC_SUBST(ALL_STATIC)

# Checks if we want MPI support
AC_ARG_ENABLE(mpi,    [  --enable-mpi             enable MPI support (default=no)])
AM_CONDITIONAL([PARALLEL], [test "$enable_mpi" = "yes"])
AC_ARG_WITH(mpi-inc,  [  --with-mpi-cflags=FLAGS   specify MPI compliler flags])
AC_ARG_WITH(mpi-libs, [  --with-mpi-libs=LIBS      specify MPI linker flags])
AC_ARG_WITH(mpi-run,  [  --with-mpi-run=cmd        specify MPI run command (default=mpirun)])

if test "$enable_mpi" = "yes"; then
    MPI_DEF="-DUSE_MPI"
    if test "X$with_mpi_libs" = "X"; then
        MPI_LIBS=""
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

AC_ARG_ENABLE(debug,[  --enable-debug          enable code for debugging checks])
if test "x$enable_debug" = "xyes"; then
   CXXFLAGS="-DDEBUG -DBZ_DEBUG $CXXFLAGS"
   FFLAGS="-DDEBUG $FFLAGS"
fi



AC_ARG_ENABLE(profile,[  --enable-profile        instrument code with profiling information])
if test "x$enable_profile" = "xyes"; then
   CXXFLAGS="-pg $CXXFLAGS"
   FFLAGS="-pg $FFLAGS"
   LDFLAGS="-pg $LDFLAGS"
fi


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



#AC_ARG_WITH(cblas-cflags,[  --with-cblas-cflags=FLAGS specify CBLAS linker flags])
#AC_ARG_WITH(cblas-libs, [  --with-cblas-libs=LIBS    specify CBLAS linker flags])
#if test "x$with_cblas_cflags" != "x"; then
#   CBLAS_CFLAGS="$with_cblas_cflags"
#fi
#if test "x$with_cblas_libs" != "x"; then
#   CBLAS_LIBS="$with_cblas_libs"
#fi
#AC_SUBST(CBLAS_CFLAGS)
#AC_SUBST(CBLAS_LIBS)
#if test "x$cblas_ok" = "xno"; then
#   AC_MSG_WARN([Didn't find libcblas in standard place.])
#   if test "x$with_cblas_libs" = "x"; then
#       AC_MSG_ERROR([CBLAS library could not be found.  Please specify location with --with-cblas-libs=])
#   fi
#fi



AC_ARG_WITH(blas-cflags, [  --with-blas-cflags=CFLAGS  specify BLAS compile flags ])
if test "x$with_blas_cflags" != "x"; then
   BLAS_CFLAGS="$with_blas_cflags"
else
   AC_CHECK_LIB([blas],[dgemv_])
fi
AC_SUBST(BLAS_CFLAGS)

AC_ARG_WITH(blas-libs,   [  --with-blas-libs=LIBS      specify BLAS linker flags ])
if test "x$with_blas_libs" != "x"; then
   BLAS_LIBS="$with_blas_libs"
else
   AC_CHECK_LIB([blas],[dgemv_])
fi
AC_SUBST(BLAS_LIBS)

AC_ARG_ENABLE(cblas,      [  --enable-cblas            use     CBLAS instead of BLAS (default=no)])
AC_ARG_ENABLE(mkl_cblas,  [  --enable-mkl-cblas        use MKL CBLAS instead of BLAS (default=no)])
if test "$enable_mkl_cblas" = "yes"; then
  CBLAS_DEF="-DUSE_MKL_CBLAS -DUSE_CBLAS"
else
  CBLAS_DEF=" "
fi
if test "$enable_cblas" = "yes"; then
  CBLAS_DEF="-DUSE_CBLAS"
fi
AC_SUBST(CBLAS_DEF)



AC_ARG_WITH(lapack-cflags,[  --with-lapack-cflags=FLAGS specify LAPACK linker flags])
AC_ARG_WITH(lapack-libs, [  --with-lapack-libs=LIBS    specify LAPACK linker flags])
if test "x$with_lapack_cflags" != "x"; then
   LAPACK_CFLAGS="$with_lapack_cflags"
fi
if test "x$with_lapack_libs" != "x"; then
   LAPACK_LIBS="$with_lapack_libs"
else
   AC_CHECK_LIB([lapack], [dgesvd_])
fi
AC_SUBST(LAPACK_CFLAGS)
AC_SUBST(LAPACK_LIBS)
if test "x$lapack_ok" = "xno"; then
   AC_MSG_WARN([Didn't find liblapack in standard place.])
   if test "x$with_lapack_libs" = "x"; then
       AC_MSG_ERROR([LAPACK library could not be found.  Please specify location with --with-lapack-libs=])
   fi
fi

#AC_ARG_WITH(xml2-cflags,[  --with-xml2-cflags=FLAGS specify XML2 linker flags])
#AC_ARG_WITH(xml2-libs, [  --with-xml2-libs=LIBS    specify XML2 linker flags])
#if test "x$with_xml2_cflags" != "x"; then
#   XML2_CFLAGS="$with_xml2_cflags"
#fi
#if test "x$with_xml2_libs" != "x"; then
#   XML2_LIBS="$with_xml2_libs"
#fi
#AC_SUBST(XML2_CFLAGS)
#AC_SUBST(XML2_LIBS)
#if test "x$xml2_ok" = "xno"; then
#   AC_MSG_WARN([Didn't find libxml2 in standard place.])
#   if test "x$with_xml2_libs" = "x"; then
#       AC_MSG_ERROR([XML2 library could not be found.  Please specify location with --with-xml2-libs=])
#   fi
#fi

AC_ARG_WITH(gsl-cflags,[  --with-gsl-cflags=FLAGS specify GSL linker flags])
AC_ARG_WITH(gsl-libs,  [  --with-gsl-libs=LIBS    specify GSL linker flags])
if test "x$with_gsl_cflags" != "x"; then
   GSL_CFLAGS="$with_gsl_cflags"
fi
if test "x$with_gsl_libs" != "x"; then
   GSL_LIBS="$with_gsl_libs"
fi
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)
if test "x$gsl_ok" = "xno"; then
   AC_MSG_WARN([Didn't find libgsl in standard place.])
   if test "x$with_gsl_libs" = "x"; then
       AC_MSG_ERROR([GSL library could not be found.  Please specify location with --with-gsl-libs=])
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
AC_CHECK_FUNCS([clock_gettime floor pow sqrt strtol sincos])

AC_CONFIG_FILES([ Makefile                  \
                  common.pc                 \
                  mpicommon.pc              \
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
AC_MSG_RESULT([ HDF5_LIBS   = $HDF5_LIBS ])
#AC_MSG_RESULT([ XML2_LIBS   = $XML2_LIBS ])
AC_MSG_RESULT([ GSL_LIBS    = $GSL_LIBS  ])
AC_MSG_RESULT([ BLAS_LIBS   = $BLAS_LIBS ])
AC_MSG_RESULT([ LAPACK_LIBS = $LAPACK_LIBS ])
AC_MSG_RESULT([ FFTW3_LIBS  = $FFTW3_LIBS ])
AC_MSG_RESULT([ FFT_DEF     = $FFT_DEF ])
#AC_MSG_RESULT([ CBLAS_LIBS  = $CBLAS_LIBS ])
AC_OUTPUT
