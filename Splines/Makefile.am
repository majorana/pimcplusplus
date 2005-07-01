check_PROGRAMS = TestBicubic TestComplexMultiTricubicSpline TestComplexMultiTricubicSpline2 TestGrid TestMultiTricubicSpline TestMultiTricubicSpline2 TestMultiTricubicSpline2 TestMyTricubic TestPeriodic TestQuintic

lib_LTLIBRARIES = libspline.la
libspline_la_SOURCES = BicubicSpline.h                \
                       BicubicSpline.cc               \
                       CubicSpline.h                  \
                       CubicSpline.cc                 \
                       Grid.h                         \
                       Grid.cc                        \
                       ComplexMulitTricubicSpline.h   \
                       ComplexMultiTricubicSpline.cc  \
                       ComplexMulitTricubicSpline2.h  \
                       ComplexMultiTricubicSpline2.cc \
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


include_HEADERS = BicubicSpline.h ComplexMultiTricubicSpline.h ComplexMultiTricubicSpline2.h CubicSpline.h Grid.h MultiTricubicSpline.h MultiTricubicSpline2.h MultiTricubicSpline3.h MyTricubicSpline.h PeriodicSpline.h QuinticSpline.h 

