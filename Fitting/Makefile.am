check_PROGRAMS = TestFitting TestEOS

LDADD = libfitting.la ../MatrixOps/libmatops.la $(LAPACK_LIBS) $(BLITZ_LIBS)

TestFitting_SOURCES = TestFitting.cc
TestEOS_SOURCES     = TestEOS.cc

AM_CXXFLAGS = $(BLITZ_CFLAGS)

noinst_LTLIBRARIES = libfitting.la

libfitting_la_SOURCES = Fitting.cc      \
                        Fitting.h       \
                        NonlinearFit.h 

libfitting_la_CXXFLAGS = $(BLITZ_CFLAGS) $(AM_CXXFLAGS)

includedir=$(prefix)/include/Common/Fitting

