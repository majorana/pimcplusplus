check_PROGRAMS = TestDFTAtom TestRadialWF

TestDFTAtom_SOURCES  = TestDFTAtom.cc
TestDFTAtom_LDADD    = ../IO/libciao.la ../Splines/libspline.la libatom.la $(XML2_LIBS) $(BLITZ_LIBS)

TestRadialWF_SOURCES = TestRadialWF.cc
TestRadialWF_LDADD   = ../IO/libciao.la ../Splines/libspline.la libatom.la $(XML2_LIBS) $(BLITZ_LIBS)

includedir=$(prefix)/include/Common/Atom

lib_LTLIBRARIES = libatom.la
libatom_la_SOURCES = AtomBase.h          \
                     DFTAtom.h           \
                     DFTAtom.cc          \
                     NewAtom.h           \
                     RadialWF.h          \
                     RadialWF.cc

#if WITH_F77
#     libsprng_la_SOURCES += common/fwrap_mpi.c
#endif
libatom_la_CXXFLAGS = $(BLITZ_CFLAGS) $(XML2_CFLAGS)
libatom_la_LIBADD = $(GMP_LIBS) 

#if HAVE_GMP
#     libsprng_la_SOURCES += pmlcg/pmlcg.h pmlcg/pmlcg.c pmlcg/info.h pmlcg/basic.h 
#endif


include_HEADERS = AtomBase.h DFTAtom.h RadialWF.h


