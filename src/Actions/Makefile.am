if WANT_PIMC_DEV
  noinst_LTLIBRARIES = libactionsDev.la libactions.la
else 
  noinst_LTLIBRARIES = libactions.la
endif

libactions_la_CXXFLAGS = $(COMMON_CFLAGS) $(QMC_DEF) $(BOOST_LIB) $(LIBXML_LIB) $(MPI_DEF)
libactions_la_SOURCES =            \
    ActionBase.cc                  \
    ActionBase.h                   \
    ActionsClass.cc                \
   ../det_calc_uekt.h              \
   ../det_calc_uekt.c              \
    ActionsClass.h                 \
    DavidLongRangeClass.cc         \
    DavidLongRangeClass.h          \
    FixedPhaseActionClass.cc       \
    FixedPhaseActionClass.h        \
    FreeNodalActionClass.cc        \
    FreeNodalActionClass.h         \
    GroundStateNodalActionClass.cc \
    GroundStateNodalActionClass.h  \
    Hermele.h                      \
    DualHermele.h                  \
    Hermele.cc                     \
    DualHermele.cc                 \
    Josephson.h                    \
    Josephson.cc                   \
    KineticClass.cc                \
    KineticClass.h                 \
    KineticSphereClass.cc          \
    KineticSphereClass.h           \
    LongRangeClass.cc              \
    LongRangeClass.h               \
    LongRangePotClass.cc           \
    LongRangePotClass.h            \
    LongRangeRPAClass.cc           \
    LongRangeRPAClass.h            \
		MoleculeInteractionsClass.cc   \
		MoleculeInteractionsClass.h    \
    Mu.h                           \
    Mu.cc                          \
    NodalActionClass.cc            \
    NodalActionClass.h             \
    OpenLoopImportance.cc          \
    OpenLoopImportance.h           \
    PairFixedPhase.h               \
    PairFixedPhase.cc              \
    popen2.h                       \
    popen2.cc                      \
		QBoxAction.cc                  \
		QBoxAction.h                   \
    QMCSamplingClass.h             \
    QMCSamplingClass.cc            \
    ShortRangeApproximateClass.cc  \
    ShortRangeApproximateClass.h   \
    ShortRangeClass.cc             \
    ShortRangeClass.h              \
    ShortRangeOnClass.h            \
    ShortRangeOnClass.cc           \
    ShortRangePotClass.cc          \
    ShortRangePotClass.h           \
    ShortRangePrimitive.h          \
    ShortRangePrimitive.cc         \
		ST2WaterClass.cc               \
		ST2WaterClass.h                \
    StructureReject.cc             \
    Tether.h                       \
    Tether.cc                      \
		TIP5PWaterClass.cc             \
		TIP5PWaterClass.h              \
    VariationalPI.h                \
    VariationalPI.cc               \
    TruncatedInverse.h             \
    TruncatedInverse.cc            \
    VariationalPI.h                \
    VariationalPI.cc               

libactionsDev_la_CXXFLAGS = -DBUILD_DEV $(COMMON_CFLAGS) $(QMC_DEF) $(BOOST_LIB) $(LIBXML_LIB) $(MPI_DEF)
libactionsDev_la_SOURCES =         \
    ActionBase.cc                  \
    ActionBase.h                   \
    ActionsClass.cc                \
   ../det_calc_uekt.h              \
   ../det_calc_uekt.c              \
    ActionsClass.h                 \
    DavidLongRangeClass.cc         \
    DavidLongRangeClass.h          \
    FixedPhaseActionClass.cc       \
    FixedPhaseActionClass.h        \
    FreeNodalActionClass.cc        \
    FreeNodalActionClass.h         \
    GroundStateNodalActionClass.cc \
    GroundStateNodalActionClass.h  \
    Hermele.h                      \
    DualHermele.h                  \
    Hermele.cc                     \
    DualHermele.cc                 \
    Josephson.h                    \
    Josephson.cc                   \
    KineticClass.cc                \
    KineticClass.h                 \
    KineticSphereClass.cc          \
    KineticSphereClass.h           \
    LongRangeClass.cc              \
    LongRangeClass.h               \
    LongRangePotClass.cc           \
    LongRangePotClass.h            \
    LongRangeRPAClass.cc           \
    LongRangeRPAClass.h            \
		MoleculeInteractionsClass.cc   \
		MoleculeInteractionsClass.h    \
    Mu.h                           \
    Mu.cc                          \
    NodalActionClass.cc            \
    NodalActionClass.h             \
    OpenLoopImportance.cc          \
    OpenLoopImportance.h           \
    PairFixedPhase.h               \
    PairFixedPhase.cc              \
    popen2.h                       \
    popen2.cc                      \
		QBoxAction.cc                  \
		QBoxAction.h                   \
    QMCSamplingClass.h             \
    QMCSamplingClass.cc            \
    ShortRangeApproximateClass.cc  \
    ShortRangeApproximateClass.h   \
    ShortRangeClass.cc             \
    ShortRangeClass.h              \
    ShortRangeOnClass.h            \
    ShortRangeOnClass.cc           \
    ShortRangePotClass.cc          \
    ShortRangePotClass.h           \
    ShortRangePrimitive.cc         \
    ShortRangePotClass.cc          \
		ST2WaterClass.cc               \
		ST2WaterClass.h                \
    StructureReject.cc             \
    Tether.h                       \
    Tether.cc                      \
		TIP5PWaterClass.cc             \
		TIP5PWaterClass.h              \
    VariationalPI.h                \
    VariationalPI.cc               \
    TruncatedInverse.h             \
    TruncatedInverse.cc            \
    VariationalPI.h                \
    VariationalPI.cc               



EXTRA_DIST = FixedPhase.tex
