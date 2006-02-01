check_PROGRAMS = TestPW2 BandStructure MPIBandStructure

TestPW2_SOURCES = TestPW2.cc
BandStructure_SOURCES = BandStructure.cc BandStructure.h
MPIBandStructure_SOURCES = MPIBandStructure.cc MPIBandStructure.h

LDADD = libpw.la                 ../IO/libIO.la                           \
        ../FFT/libfft.la       ../PH/libPH.la                             \
        ../MatrixOps/libmatops.la ../Splines/libspline.la                 \
        ../Fitting/libfitting.la ../Integration/libquad.la                \
        ../MPI/libcomm.la                                                 \
        $(XML2_LIBS) $(HDF5_LIBS) $(BLITZ_LIBS) $(G2C_LIBS) $(FFTW3_LIBS) \
        $(LAPACK_LIBS) $(BLAS_LIBS) $(GSL_LIBS) $(FLIBS) $(MPI_LIBS)

AM_CXXFLAGS = $(BLITZ_CFLAGS) $(XML2_CFLAGS)  $(FFTW3_CFLAGS)  \
              $(HDF5_CFLAGS)  $(CBLAS_CFLAGS) $(MPI_DEF) $(MPI_CFLAGS)

noinst_LTLIBRARIES = libpw.la
libpw_la_SOURCES = ConjGrad2.h        \
                   ConjGrad2.cc       \
                   ConjGradMPI.h      \
                   ConjGradMPI.cc     \
                   Coulomb.h          \
                   Coulomb.cc         \
                   CoulombFFT.h       \
                   CoulombFFT.cc      \
                   FFTBox.h           \
                   FFTBox.cc          \
                   GVecs.h            \
                   GVecs.cc           \
                   HamiltonianBase.h  \
                   HamiltonianBase.cc \
                   Hamiltonians.h     \
                   Hamiltonians.cc    \
                   LocalPotFFT.h      \
                   LocalPotFFT.cc     \
                   PHPot.h            \
                   PHPot.cc           \
                   PHPotFFT.h         \
                   PHPotFFT.cc        \
                   PlaneWaves.h       \
                   PlaneWaves.cc      \
                   PlaneWavesMPI.h    \
                   PlaneWavesMPI.cc   \
                   VectorOps.h

#TestXML_SOURCES = TestXML.cc
#TestXML_LDADD = libciao.la $(XML2_LIBS) $(BLITZ_LIBS)

#TestHDF5_SOURCES = TestHDF5.cc
#TestHDF5_LDADD = libciao.la $(XML2_LIBS) $(BLITZ_LIBS)

includedir=$(prefix)/include/Common/PlaneWavePHDFT


#if WITH_F77
#     libsprng_la_SOURCES += common/fwrap_mpi.c
#endif
libpw_la_CXXFLAGS = $(BLITZ_CFLAGS) $(XML2_CFLAGS)  $(FFTW3_CFLAGS) \
                    $(HDF5_CFLAGS)  $(CBLAS_CFLAGS) $(AM_CXXFLAGS)
libpw_la_LIBADD   =

#if HAVE_GMP
#     libsprng_la_SOURCES += pmlcg/pmlcg.h pmlcg/pmlcg.c pmlcg/info.h pmlcg/basic.h 
#endif


include_HEADERS = ConjGrad2.h FFTBox.h GVecs.h Hamiltonian2.h PlaneWaves.h \
                  VectorOps.h PlaneWavesMPI.h ConjGradMPI.h Hamiltonians.h \
                  Coulomb.h CoulombFFT.h HamiltonianBase.h LocalPotFFT.h   \
                  PHPot.h PHPotFFT.h

