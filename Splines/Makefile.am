check_PROGRAMS = TestBicubic TestComplexMultiTricubicSpline TestGrid TestMultiTricubicSpline TestMyTricubic TestPeriodic TestQuintic

TestBicubic_SOURCES = TestBicubic.cc 
TestComplexMultiTricubicSpline_SOURCES = TestComplexMultiTricubicSpline.cc
TestGrid_SOURCES = TestGrid.cc
TestMultiTricubicSpline_SOURCES = TestMultiTricubicSpline.cc
TestMyTricubic_SOURCES = TestMyTricubic.cc
TestPeriodic_SOURCES = TestPeriodic.cc
TestQuintic_SOURCES = TestQuintic.cc


LDADD = ../IO/libciao.la libspline.la $(XML2_LIBS) $(HDF5_LIBS)
CXXFLAGS = $(BLITZ_CFLAGS) $(XML2_CFLAGS) $(AM_CXXFLAGS)

lib_LTLIBRARIES = libspline.la
libspline_la_SOURCES = BicubicSpline.h                \
                       BicubicSpline.cc               \
                       CubicSpline.h                  \
                       CubicSpline.cc                 \
                       Grid.h                         \
                       Grid.cc                        \
                       ComplexMulitTricubicSpline.h   \
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
                       r3spline.f                     \
                       z3spline.f          

#TestXML_SOURCES = TestXML.cc
#TestXML_LDADD = libciao.la $(XML2_LIBS) $(BLITZ_LIBS)

#TestHDF5_SOURCES = TestHDF5.cc
#TestHDF5_LDADD = libciao.la $(XML2_LIBS) $(BLITZ_LIBS)

includedir=$(prefix)/include/Common/Splines


#if WITH_F77
#     libsprng_la_SOURCES += common/fwrap_mpi.c
#endif
libspline_la_CXXFLAGS = $(BLITZ_CFLAGS) $(XML2_CFLAGS)
libspline_la_LIBADD   =

#if HAVE_GMP
#     libsprng_la_SOURCES += pmlcg/pmlcg.h pmlcg/pmlcg.c pmlcg/info.h pmlcg/basic.h 
#endif


include_HEADERS = BicubicSpline.h ComplexMultiTricubicSpline.h CubicSpline.h Grid.h MultiTricubicSpline.h MyTricubicSpline.h PeriodicSpline.h QuinticSpline.h 

