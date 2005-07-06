check_PROGRAMS = TestXML TestHDF5

LDADD = libciao.la $(XML2_LIBS) $(HDF5_LIBS) $(BLITZ_LIBS) $(AM_LDFLAGS)
AM_CXXFLAGS = $(XML2_CFLAGS) $(HDF5_CFLAGS) $(BLITZ_CFLAGS)

TestXML_SOURCES = TestXML.cc
#TestXML_LDADD = libciao.la $(XML2_LIBS) $(BLITZ_LIBS)
#TestXML_CXXFLAGS = $(BLITZ_CFLAGS) $(XML2_CFLAGS) $(AM_CXXFLAGS)

TestHDF5_SOURCES = TestHDF5.cc
#TestHDF5_LDADD = libciao.la $(XML2_LIBS) $(BLITZ_LIBS)
#TestHDF5_CXXFLAGS = $(BLITZ_CFLAGS) $(XML2_CFLAGS) $(AM_CXXFLAGS)

includedir=$(prefix)/include/Common/IO

noinst_LTLIBRARIES = libciao.la
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
libciao_la_CXXFLAGS = $(BLITZ_CFLAGS) $(XML2_CFLAGS) $(HDF5_CFLAGS)

#if HAVE_GMP
#     libsprng_la_SOURCES += pmlcg/pmlcg.h pmlcg/pmlcg.c pmlcg/info.h pmlcg/basic.h 
#endif


include_HEADERS = InputOutput.h InputOutputASCII.h InputOutputBase.h InputOutputHDF5.h InputOutputXML.h

