check_PROGRAMS = TestRandom

TestRandom_SOURCES = TestRandom.cc

CXXFLAGS = $(MPI_CFLAGS) $(MPI_DEF) $(BLITZ_CFLAGS) $(SPRNG2_CFLAGS) $(AM_CXXFLAGS) 
LDADD    = ../MPI/libcomm.la $(SPRNG2_LIBS) $(MPI_LIBS)

includedir=$(prefix)/include/Common/Random

include_HEADERS = Random.h
