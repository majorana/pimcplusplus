if WANT_PATHVIS
   MAYBE_PATHVIS = Visual
endif

#QMCPACK_ROOT = /home/gergely/qmcPlusPlusNew
#QMCPACK_HOME = $(QMCPACK_ROOT)/build

SUBDIRS = Actions Moves Observables Langevin $(MAYBE_PATHVIS)
bin_PROGRAMS = pimc++ FreeParticles VacancyLoc vmc


AM_CXXFLAGS = $(COMMON_CFLAGS) $(MPI_DEF) $(QMC_DEF) 


SHARED_SOURCES =          \
    det_calc_uekt.h       \
    det_calc_uekt.c       \
    Common.h              \
    EventClass.h          \
    EventClass.cc         \
    InitPaths.cc          \
    JosephsonJunctionHelper.cc \
    GridClass.h           \
    GridClass.cc          \
    LoopClass.cc          \
    LoopClass.h           \
    MirroredClass.cc      \
    MirroredClass.h       \
    PIMCClass.cc          \
    PIMCClass.h           \
    PathClass.cc          \
    PathClass.h           \
    PathDataClass.cc      \
    PathDataClass.h       \
    QMCWrapper.cc         \
    QMCWrapper.h          \
    RunInfoClass.h        \
    SpeciesClass.cc       \
    SpeciesClass.h        \
    WormHelper.cc         \
    WriteData.h           \
    WriteData.cc          

pimc___SOURCES =          \
    $(SHARED_SOURCES)     \
    Main.cc               

FreeParticles_SOURCES =  FreeParticles.cc
FreeParticles_LDADD   = $(COMMON_LIBS)

VacancyLoc_SOURCES = \
                      Observables/lsapr.f           \
                      VacancyLoc.cc  
VacancyLoc_LDADD = $(COMMON_LIBS)

vmc_SOURCES = vmc.cc 
vmc_LDADD = $(COMMON_LIBS)

#INCLUDES = -I$(QMCPACK_HOME)/src -I$(QMCPACK_ROOT)/src 

pimc___LDADD = Moves/libmoves.la Actions/libactions.la \
               Observables/libobservables.la $(COMMON_LIBS) $(QMC_LIB) $(ORDER_N_LIBS) $(BOOST_LIB)
#							 $(QMCPACK_HOME)/lib/libqmc.a \
#							 $(QMCPACK_HOME)/lib/libqmcutil.a \
#							 $(QMCPACK_HOME)/lib/libqmcbase.a
