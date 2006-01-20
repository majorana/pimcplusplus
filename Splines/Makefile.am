check_PROGRAMS =  TestFastCubic TestBicubic TestComplexMultiTricubicSpline TestGrid TestMultiTricubicSpline TestMyTricubic TestPeriodic TestQuintic

TestFastCubic_SOURCES = TestFastCubic.cc
TestBicubic_SOURCES = TestBicubic.cc 
TestComplexMultiTricubicSpline_SOURCES = TestComplexMultiTricubicSpline.cc
TestGrid_SOURCES = TestGrid.cc
TestMultiTricubicSpline_SOURCES = TestMultiTricubicSpline.cc
TestMyTricubic_SOURCES = TestMyTricubic.cc
TestPeriodic_SOURCES = TestPeriodic.cc
TestQuintic_SOURCES = TestQuintic.cc

LDADD = ../IO/libIO.la libspline.la $(XML2_LIBS) $(HDF5_LIBS) $(BLITZ_LIBS) $(FLIBS) 
AM_CXXFLAGS = $(BLITZ_CFLAGS) $(XML2_CFLAGS) $(HDF5_CFLAGS)

noinst_LTLIBRARIES = libspline.la
libspline_la_SOURCES = BicubicSpline.h                \
                       BicubicSpline.cc               \
                       CubicSpline.h                  \
                       CubicSpline.cc                 \
                       FastCubicSpline.cc             \
                       FastCubicSpline.h              \
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
