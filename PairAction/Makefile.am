noinst_LTLIBRARIES = libpairact.la
libpairact_la_SOURCES = DavidPAClass.cc    \
                        DavidPAClass.h     \
                        PADipoleFit.cc     \
                        PADipoleFit.h      \
                        PAFit.h            \
                        PAFitBase.h        \
                        PAMonopoleFit.cc   \
                        PAMonopoleFit.h    \
                        PATripoleFit.cc    \
                        PATripoleFit.h     \
                        PAclassicalFit.cc  \
                        PAclassicalFit.h   \
                        PAcoulombBCFit.cc  \
                        PAcoulombBCFit.h   \
                        PAtricubicFit.cc   \
                        PAtricubicFit.h    \
                        PAtricubicFit2.cc  \
                        PAtricubicFit2.h   \
                        PAzeroFit.cc       \
                        PAzeroFit.h        \
                        Particle.h


includedir=$(prefix)/include/Common/PairAction


#if WITH_F77
#     libsprng_la_SOURCES += common/fwrap_mpi.c
#endif
libpairact_la_CXXFLAGS = $(BLITZ_CFLAGS) $(XML2_CFLAGS) $(AM_CXXFLAGS)
libpairact_la_LIBADD   =

#if HAVE_GMP
#     libsprng_la_SOURCES += pmlcg/pmlcg.h pmlcg/pmlcg.c pmlcg/info.h pmlcg/basic.h 
#endif


include_HEADERS = DavidPAClass.h    PAMonopoleFit.h   PAtricubicFit2.h \
                  PADipoleFit.h     PATripoleFit.h    PAzeroFit.h      \
                  PAFit.h           PAclassicalFit.h  Particle.h       \
                  PAFitBase.h       PAcoulombBCFit.h  PAtricubicFit.h


