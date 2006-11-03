check_PROGRAMS =  TestFastCubic TestFastMultiCubic TestBicubic TestComplexMultiTricubicSpline \
                  TestGrid TestMultiTricubicSpline TestMyTricubic TestPeriodic TestQuintic    \
                  TestCubicBspline TestTricubicBspline TestBlips

TestFastCubic_SOURCES = TestFastCubic.cc
TestFastMultiCubic_SOURCES = TestFastMultiCubic.cc
TestBicubic_SOURCES = TestBicubic.cc 
TestComplexMultiTricubicSpline_SOURCES = TestComplexMultiTricubicSpline.cc
TestGrid_SOURCES = TestGrid.cc
TestMultiTricubicSpline_SOURCES = TestMultiTricubicSpline.cc
TestMyTricubic_SOURCES = TestMyTricubic.cc
TestPeriodic_SOURCES = TestPeriodic.cc
TestQuintic_SOURCES = TestQuintic.cc
TestCubicBspline_SOURCES = TestCubicBspline.cc
TestTricubicBspline_SOURCES = TestTricubicBspline.cc
TestBlips_SOURCES = TestBlips.cc

LDADD = ../IO/libIO.la ../PlaneWavePHDFT/libpw.la ../FFT/libfft.la libspline.la  $(HDF5_LIBS) $(FLIBS) $(XML2_LIBS) $(BLITZ_LIBS) $(FFTW3_LIBS) $(LAPACK_LIBS)
AM_CXXFLAGS = $(BLITZ_CFLAGS) $(XML2_CFLAGS) $(HDF5_CFLAGS) $(FFTW3_CFLAGS)

noinst_LTLIBRARIES = libspline.la
libspline_la_SOURCES = BicubicSpline.h                \
                       BicubicSpline.cc               \
                       CubicBspline.cc                \
                       CubicBspline.h                 \
                       CubicSpline.h                  \
                       CubicSpline.cc                 \
                       FastCubicSpline.cc             \
                       FastCubicSpline.h              \
                       FastMultiCubicSpline.cc        \
                       FastMultiCubicSpline.h         \
                       Grid.h                         \
                       Grid.cc                        \
                       ComplexMultiTricubicSpline.h   \
                       ComplexMultiTricubicSpline.cc  \
                       MultiTricubicSpline.h          \
                       MultiTricubicSpline.cc         \
                       MyTricubicSpline.h             \
                       MyTricubicSpline.cc            \
                       PeriodicSpline.h               \
                       PeriodicSpline.cc              \
                       QuinticSpline.h                \
                       QuinticSpline.cc               \
                       QuinticSplines.f               \
                       TricubicSpline.h               \
                       TricubicBspline.cc             \
                       TricubicBspline.h              \
                       r3spline.f                     \
                       z3spline.f                     \
                       bcspeval.f                     \
                       bcspline.f                     \
                       evtricub.f                     \
                       herm3ev.f                      \
                       ibc_ck.f                       \
                       mktricubw.f                    \
                       splinck.f                      \
                       tcspeval.f                     \
                       tcspline.f                     \
                       v_spline.f                     \
                       zonfind.f



#TestXML_SOURCES = TestXML.cc
#TestXML_LDADD = libciao.la $(XML2_LIBS) $(BLITZ_LIBS)

#TestHDF5_SOURCES = TestHDF5.cc
#TestHDF5_LDADD = libciao.la $(XML2_LIBS) $(BLITZ_LIBS)

includedir=$(prefix)/include/Common/Splines


#if WITH_F77
#     libsprng_la_SOURCES += common/fwrap_mpi.c
#endif
libspline_la_CXXFLAGS = $(BLITZ_CFLAGS) $(XML2_CFLAGS) $(AM_CXXFLAGS)
libspline_la_LIBADD   =

#if HAVE_GMP
#     libsprng_la_SOURCES += pmlcg/pmlcg.h pmlcg/pmlcg.c pmlcg/info.h pmlcg/basic.h 
#endif


include_HEADERS = BicubicSpline.h ComplexMultiTricubicSpline.h CubicSpline.h Grid.h MultiTricubicSpline.h MyTricubicSpline.h PeriodicSpline.h QuinticSpline.h TricubicSpline.h FastCubicSpline.h

EXTRA_DIST = CubicSplines.tex
