AM_LDFLAGS = $(ALL_STATIC)

check_PROGRAMS = TestFitting 
bin_PROGRAMS   = LiFit VinetFit VinetFit2

LDADD = libfitting.la ../MatrixOps/libmatops.la ../IO/libIO.la $(LAPACK_LIBS) $(BLITZ_LIBS) $(HDF5_LIBS) $(GSL_LIBS) $(SPRNG2_LIBS) $(FLIBS)

TestFitting_SOURCES = TestFitting.cc
LiFit_SOURCES       = TestEOS.cc LiEOS.h
VinetFit_SOURCES    = VinetFit.cc VinetEOS.h
VinetFit2_SOURCES   = VinetFit2.cc PhononFreeEnergy.h PhononFreeEnergy.cc \
                      DebyeModel.h DebyeModel.cc VinetEOS2.h RamanFreq.h RamanFreq.cc


AM_CXXFLAGS = $(BLITZ_CFLAGS) $(HDF5_CFLAGS) $(XML2_CFLAGS) $(GSL_CFLAGS)

noinst_LTLIBRARIES = libfitting.la

libfitting_la_SOURCES = Fitting.cc      \
                        Fitting.h       \
                        NonlinearFit.h 

#libfitting_la_CXXFLAGS = $(BLITZ_CFLAGS) $(AM_CXXFLAGS)
include_HEADERS = Fitting.h NonlinearFit.h
includedir=$(prefix)/include/Common/Fitting

