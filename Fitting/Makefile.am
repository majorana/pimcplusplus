check_PROGRAMS = TestFitting

LDADD = libfitting.la ../MatrixOps/libmatops.la $(LAPACK_LIBS) 

TestFitting_SOURCES = TestFitting.cc

noinst_LTLIBRARIES = libfitting.la

libfitting_la_SOURCES = Fitting.cc    \
                        Fitting.h     

libfitting_la_CXXFLAGS = $(AM_CXXFLAGS)

includedir=$(prefix)/include/Common/Fitting

