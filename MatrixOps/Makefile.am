check_PROGRAMS = Test

Test_SOURCES = Test.cc 

LDADD = libmatops.la $(LAPACK_LIBS) $(BLAS_LIBS)
CXXFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(AM_CXXFLAGS)

lib_LTLIBRARIES = libmatops.la
libmatops_la_SOURCES = MatrixOps.h MatrixOps.cc

includedir=$(prefix)/include/Common/MatrixOps

libmatops_la_CXXFLAGS = $(BLITZ_CFLAGS) $(XML2_CFLAGS)

include_HEADERS = MatrixOps.h

