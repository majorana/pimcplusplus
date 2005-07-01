check_PROGRAMS = TestIntegrate TestHermite

TestIntegrate_SOURCES = TestIntegrate.cc
TestHermite_SOURCES = TestHermite.cc 

LDADD = libquad.la ../IO/libciao.la  $(XML2_LIBS) $(BLITZ_LIBS)
CXXFLAGS = $(BLITZ_CFLAGS) $(XML2_CFLAGS) $(AM_CXXFLAGS)

lib_LTLIBRARIES = libquad.la
libquad_la_SOURCES = GKIntegration.h    \
                     GKIntegration.cc   \
                     HermiteQuad.h      \
                     HermiteQuad.cc     \
                     Integrate.h        \
                     Integrate.cc       \
                     RungeKutta.h       \
                     Standard.h         

libquad_la_CXXFLAGS = $(XML2_CFLAGS) $(HDF5_CLFAGS) $(AM_CXXFLAGS)

includedir=$(prefix)/include/Common/Integration

include_HEADERS = GKIntegration.h HermiteQuad.h Integrate.h RungeKutta.h Standard.h
