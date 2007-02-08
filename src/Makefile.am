if WANT_PATHVIS
   MAYBE_PATHVIS = Visual
endif

#QMCPACK_ROOT = /home/gergely/qmcPlusPlusNew
#QMCPACK_HOME = $(QMCPACK_ROOT)/build

SUBDIRS = Actions Moves Observables Langevin $(MAYBE_PATHVIS)

if WANT_PIMC_DEV
  bin_PROGRAMS = pimc++ pimcDev FreeParticles VacancyLoc vmc
else
  bin_PROGRAMS = pimc++ FreeParticles VacancyLoc vmc
endif
  


AM_CXXFLAGS = $(COMMON_CFLAGS) $(MPI_DEF) $(QMC_DEF) $(LIBXML_LIB)

pimc___CXXFLAGS = $(AM_CXXFLAGS)
pimcDev_CXXFLAGS = -DBUILD_DEV $(AM_CXXFLAGS)
FreeParticles_CXXFLAGS = $(AM_CXXFLAGS)
VacancyLoc_CXXFLAGS = $(AM_CXXFLAGS)
vmc_CXXFLAGS = $(AM_CXXFLAGS)

SHARED_SOURCES =          \
    det_calc_uekt.h       \
    det_calc_uekt.c       \
    Common.h              \
    EventClass.h          \
    EventClass.cc         \
    JosephsonJunctionHelper.cc \
    GridClass.h           \
    GridClass.cc          \
    LoopClass.cc          \
    LoopClass.h           \
    MirroredClass.cc      \
    MirroredClass.h       \
    PIMCClass.cc          \
    PIMCClass.h           \
    PathDataClass.cc      \
    PathDataClass.h       \
    QMCWrapper.cc         \
    QMCWrapper.h          \
    RunInfoClass.h        \
    SpeciesClass.cc       \
    SpeciesClass.h        \
    WriteData.h           \
    WriteData.cc          

pimc___SOURCES =          \
    $(SHARED_SOURCES)     \
    InitPaths.cc          \
    PathClass.cc          \
    PathClass.h           \
    WormHelper.cc         \
    Main.cc               

pimcDev_SOURCES =         \
    $(SHARED_SOURCES)     \
    InitPathsDev.cc       \
    PathClassDev.cc       \
    PathClassDev.h        \
    WormHelperDev.cc      \
    Main.cc               

FreeParticles_SOURCES =  FreeParticles.cc
FreeParticles_LDADD   = $(COMMON_LIBS)

VacancyLoc_SOURCES = \
                      Observables/lsapr.f           \
                      VacancyLoc.cc  
VacancyLoc_LDADD = $(COMMON_LIBS)

vmc_SOURCES = VMC6.cc 
vmc_LDADD = $(COMMON_LIBS) $(LAPACK_LIBS)

pimc___LDADD = Moves/libmoves.la Actions/libactions.la \
               Observables/libobservables.la $(COMMON_LIBS) $(QMC_LIB) $(ORDER_N_LIBS) $(BOOST_LIB) $(LIBXML_LIB)

pimcDev_LDADD = Moves/libmovesDev.la Actions/libactionsDev.la \
               Observables/libobservablesDev.la $(COMMON_LIBS) $(QMC_LIB) $(ORDER_N_LIBS) $(BOOST_LIB) $(LIBXML_LIB)
