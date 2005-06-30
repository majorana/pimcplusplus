check_PROGRAMS = TestXML TestHDF5

TestXML_SOURCES = TestXML.cc
TestXML_LIBS = ciao

TestHDF5_SOURCES = TestHDF5.cc
TestHDF5_LIBS = ciao

lib_LTLIBRARIES = libciao.la
libciao_la_SOURCES = FileExpand.h        \
                     InputOutput.h       \
                     InputOutput.cc      \
                     InputOutputBase.h   \
                     InputOutputASCII.h  \
                     InputOutputASCII.cc \
                     InputOutputHDF5.h   \
                     InputOutputHDF5.cc  \
                     InputOutputXML.h    \
                     InputOutputXML.cc   
#if WITH_F77
#     libsprng_la_SOURCES += common/fwrap_mpi.c
#endif

libciao_la_CFLAGS = $(MPI_DEF) $(ENABLE_F77) $(MPI_CFLAGS) 
libciao_la_LIBADD = $(GMP_LIBS)

#if HAVE_GMP
#     libsprng_la_SOURCES += pmlcg/pmlcg.h pmlcg/pmlcg.c pmlcg/info.h pmlcg/basic.h 
#endif


include_HEADERS = InputOutput.h InputOutputASCII.h InputOutputBase.h InputOutputHDF5.h InputOutputXML.h

