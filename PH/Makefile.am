check_PROGRAMS = TestkSpacePH

TestkSpacePH_SOURCES = TestkSpacePH.cc

LDADD = libPH.la ../IO/libciao.la   ../MatrixOps/libmatops.la ../Fitting/libfitting.la ../Integration/libquad.la ../Splines/libspline.la $(XML2_LIBS) $(HDF5_LIBS) $(BLITZ_LIBS) $(G2C_LIBS) $(GSL_LIBS) 
AM_CXXFLAGS = $(BLITZ_CFLAGS) $(XML2_CFLAGS) $(GSL_CFLAGS)

noinst_LTLIBRARIES = libPH.la
libPH_la_SOURCES = CoreTransform.h    \
                   CoreTransform.cc   \
                   CoulombPot.h       \
                   CoulombPot.cc      \
                   HeAzizPot.h        \
                   HeAzizPot.cc       \
                   Potential.h        \
                   Potential.cc       \
                   PotentialBase.h    \
                   QuinticPH.h        \
                   QuinticPH.cc       \
                   ScreenedPot.h      \
                   ScreenedPot.cc     \
                   SplinePot.h        \
                   SplinePot.cc       \
                   kSpacePH.h         \
                   kSpacePH.cc      


includedir=$(prefix)/include/Common/PH


#if WITH_F77
#     libsprng_la_SOURCES += common/fwrap_mpi.c
#endif
libPH_la_CXXFLAGS = $(BLITZ_CFLAGS) $(XML2_CFLAGS)
libPH_la_LIBADD   =

#if HAVE_GMP
#     libsprng_la_SOURCES += pmlcg/pmlcg.h pmlcg/pmlcg.c pmlcg/info.h pmlcg/basic.h 
#endif


include_HEADERS = CoreTransform.h CoulombPot.h HeAzizPot.h Potential.h PotentialBase.h \
                  QuinticPH.h ScreenedPot.h SplinePot.h kSpacePH.h 

