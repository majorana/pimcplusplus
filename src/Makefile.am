if WANT_PATHVIS
   MAYBE_PATHVIS = Visual
endif

SUBDIRS = Actions Moves Observables $(MAYBE_PATHVIS)
bin_PROGRAMS = pimc++ FreeParticles


AM_CXXFLAGS = $(COMMON_CFLAGS) $(MPI_DEF)

pimc___SOURCES =          \
    EventClass.h          \
    Main.cc               \
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
    WrapClass.cc          \
    WrapClass.h
pimc___LDADD = Moves/libmoves.la Actions/libactions.la \
               Observables/libobservables.la $(COMMON_LIBS)



FreeParticles_SOURCES =  FreeParticles.cc 
FreeParticles_LDADD = $(COMMON_LIBS)
