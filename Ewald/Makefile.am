check_PROGRAMS = TestOptimizedBreakup

TestOptimizedBreakup_SOURCES = TestOptimizedBreakup.cc
LDADD =  libbreakup.la ../Integration/libquad.la ../MatrixOps/libmatops.la ../MPI/libcomm.la $(MPI_LIBS) $(BLITZ_LIBS) $(LAPACK_LIBS)
AM_CXXFLAGS = $(BLITZ_CFLAGS) $(MPI_CFLAGS) $(MPI_DEF)

noinst_LTLIBRARIES = libbreakup.la
libbreakup_la_SOURCES = OptimizedBreakup.cc
libbreakup_la_CXXFLAGS = $(BLITZ_CFLAGS)  $(MPI_CFLAGS) $(MPI_DEF) $(AM_CXXFLAGS)

includedir = $(prefix)/Common/Ewald/

include_HEADERS = OptimizedBreakup.h

