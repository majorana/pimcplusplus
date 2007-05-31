check_PROGRAMS = TestFitting 
bin_PROGRAMS   = eosfit

LDADD = libfitting.la ../MatrixOps/libmatops.la ../IO/libIO.la $(LAPACK_LIBS) $(BLITZ_LIBS) $(HDF5_LIBS)

TestFitting_SOURCES = TestFitting.cc
eosfit_SOURCES      = TestEOS.cc

AM_CXXFLAGS = $(BLITZ_CFLAGS) $(HDF5_CFLAGS) $(XML2_CFLAGS)

noinst_LTLIBRARIES = libfitting.la

libfitting_la_SOURCES = Fitting.cc      \
                        Fitting.h       \
                        NonlinearFit.h 

libfitting_la_CXXFLAGS = $(BLITZ_CFLAGS) $(AM_CXXFLAGS)
include_HEADERS = Fitting.h NonlinearFit.h
includedir=$(prefix)/include/Common/Fitting

