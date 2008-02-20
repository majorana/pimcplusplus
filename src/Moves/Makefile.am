#check_PROGRAMS = TestSpaceWarp

#TestSpaceWarp_SOURCES = TestSpaceWarp.cc

LDADD = $(COMMON_LIBS) libmoves.la

AM_CXXFLAGS = $(COMMON_FLAGS)

  noinst_LTLIBRARIES = libmoves.la

libmoves_la_CXXFLAGS = $(COMMON_CFLAGS) $(QMC_DEF) $(LIBXML_LIB) $(MPI_DEF)
libmoves_la_SOURCES  =           \
    AVBMove.cc                   \
    AVBMove.h                    \
    BisectionBlock.cc            \
    BisectionBlock.h             \
    BisectionSphereBlock.cc      \
    BisectionSphereBlock.h       \
    BisectionStage.cc            \
    BisectionStage.h             \
    BisectionJosephsonStage.cc   \
    BisectionJosephsonStage.h    \
    BisectionSphereStage.cc      \
    BisectionSphereStage.h       \
    CenterDropletMove.h          \
    CenterDropletMove.cc         \
    CenterofMassMove.h           \
    CenterofMassMove.cc          \
    CorrelatedBisectionBlock.cc  \
    CorrelatedBisectionBlock.h   \
    CoupledPermuteStage.cc       \
    CoupledPermuteStage.h        \
    CouplingMove.cc              \
    CouplingMove.h               \
    CouplingStage.h              \
    CouplingStage.cc             \
    DisplaceMove.cc              \
    DisplaceMove.h               \
    DisplaceMoveFast.cc              \
    DisplaceMoveFast.h               \
    EmptyStage.h                 \
    EmptyStage.cc                \
    EndStage.cc                  \
    EndStage.h                   \
    IonDisplaceMove.h            \
    IonDisplaceMove.cc           \
    IonMoveStage.h               \
    IonMoveStage.cc              \
    IonMoveManager.h             \
    IonMoveManager.cc            \
    GlobalJosephsonMove.h        \
    GlobalJosephsonMove.cc       \
    HermeleFourierMove.h         \
    HermeleFourierMove.cc        \
    LangevinMove.cc              \
    LangevinMove.h               \
    LeviFlightStage.h            \
    LeviFlightStage.cc           \
    MetaMoves.cc                 \
    MetaMoves.h                  \
    MoleculeBias.cc              \
    MoleculeBias.h               \
    MoleculeMove.cc              \
    MoleculeMove.h               \
    MoleculeMoveBase.cc          \
    MoleculeMoveBase.h           \
    MoleculeMoveManager.cc       \
    MoleculeMoveManager.h        \
    MoveBase.cc                  \
    MoveBase.h                   \
    MoveClass.h                  \
    MoveUtils.cc                 \
    MoveUtils.h                  \
    MultiStage.cc                \
    MultiStage.h                 \
    NoPermuteClass.h             \
    NoPermuteClass.cc            \
    OpenStage.h                  \
    OpenStage.cc                 \
    OpenEndMove.cc               \
    OpenEndMove.h                \
    OpenStage.h                  \
    OpenStage.cc                 \
    PermuteStage.h               \
    PermuteStage.cc              \
    PermuteTableOnClass.cc       \
    PermuteTableOnClass.h        \
    PermuteTableClass.cc         \
    PermuteTableClass.h          \
    PreSampling.cc               \
    PreSampling.h                \
    RandomPermClass.cc           \
    RandomPermClass.h            \
    RandomPermMove.h             \
    RefSliceMove.cc              \
    RefSliceMove.h               \
    SpaceWarp.cc                 \
    SpaceWarp.h                  \
    SwapMove.h                   \
    SwapMove.cc                  \
    StructureRejectStage.cc      \
    StructureRejectStage.h       \
    StageClass.h                 \
    StageClass.cc                \
    TablePermuteStage.cc         \
    TablePermuteStage.h          \
    ClusterMove.cc               \
    ClusterMove.h                \
    VariationalDisplaceMove.h    \
    VariationalDisplaceMove.cc   \
		SPS.cc                       \
		SPS.h                        \
    WormStage.h                  \
    WormStage.cc                 \
    WormGrow.h                   \
    WormGrow.cc                  \
    WormPermuteStage.h           \
    WormPermuteStage.cc   
