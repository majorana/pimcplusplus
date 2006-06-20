if WANT_PATHVIS
   MAYBE_PATHVIS = Visual
endif

#QMCPACK_ROOT = /home/gergely/qmcPlusPlusNew
#QMCPACK_HOME = $(QMCPACK_ROOT)/build

SUBDIRS = Actions Moves Observables Langevin $(MAYBE_PATHVIS)
bin_PROGRAMS = pimc++


AM_CXXFLAGS = $(COMMON_CFLAGS) $(MPI_DEF) $(QMC_DEF)

SHARED_SOURCES =          \
    Common.h              \
    EventClass.h          \
    InitPaths.cc          \
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
    RunInfoClass.h        \
    SpeciesClass.cc       \
    SpeciesClass.h        \
    WriteData.h           \
    WriteData.cc

pimc___SOURCES =          \
    $(SHARED_SOURCES)     \
    Main.cc               

#INCLUDES = -I$(QMCPACK_HOME)/src -I$(QMCPACK_ROOT)/src 

pimc___LDADD = Moves/libmoves.la Actions/libactions.la \
               Observables/libobservables.la $(COMMON_LIBS) $(QMC_LIB)
#							 $(QMCPACK_HOME)/lib/libqmc.a \
#							 $(QMCPACK_HOME)/lib/libqmcutil.a \
#							 $(QMCPACK_HOME)/lib/libqmcbase.a
