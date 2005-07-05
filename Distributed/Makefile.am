check_PROGRAMS = DistributedTest

DistributedTest_SOURCES = TestDistributedMat.cc
LDADD = ../MPI/libcomm.la libdistrib.la $(BLITZ_LIBS) $(MPI_LIBS)
CXXFLAGS = $(BLITZ_CFLAGS) $(MPI_CFLAGS) $(MPI_DEF) $(AM_CXXFLAGS)

noinst_LTLIBRARIES = libdistrib.la
libdistrib_la_SOURCES = DistributedMat.cc
libdistrib_la_CXXFLAGS = $(BLITZ_CFLAGS)  $(MPI_CFLAGS) $(MPI_DEF) $(AM_CXXFLAGS)

includedir = $(prefix)/Common/Distributed

include_HEADERS = DistributedMat.h

