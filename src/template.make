PSPLINELIB = -L$(PWD)/Common/Splines/fortran -lpspline

ifeq ($(HOSTTYPE),alpha)
include /usr/users/0/kesler/lib/Make.include
    MPILIB = -lmpi -lelan
    LIBS = $(MPILIB) $(LAPACKLIB) $(BLITZLIB) $(SPRNGLIB) $(GSLLIB) \
         $(G2CLIB) $(HDF5LIB) $(XMLLIB) -lm 
    INCL = $(BLITZINC) $(SPRNGINC) $(GSLINC) $(HDF5INC) $(XMLINC)
    MAKE = gmake
    CC = g++
    LD = g++ #-static
    F77 = g77
    EXTRADEFS = -DNOCUSERID
    CCFLAGS = -c -g  -Wno-deprecated  #-pg 
endif
ifeq ($(HOSTTYPE),rs6000)
    include /users/uiuc/ux455254/lib/Make.include
    LIBS = $(MPILIB) $(LAPACKLIB) $(BLITZLIB) $(SPRNGLIB) $(GSLLIB) \
           $(G2CLIB) $(HDF5LIB) $(XMLLIB) $(MATHLIB) -lm -lz -lf
    INCL = $(BLITZINC) $(SPRNGINC) $(GSLINC) $(HDF5INC) $(XMLINC)
    MAKE = gmake
    CC = mpCC
    LD = mpCC -LOOK_IT_UP #-qstaticinline #-static
    F77 = f77
    EXTRADEFS = -DNOCUSERID -DNOUNDERSCORE
    CCFLAGS = -c -g 
endif
ifeq ($(HOSTTYPE),i386-linux)
    include /home/common/Codes/Make.include	
    LIBS = $(BLITZLIB) $(SPRNGLIB) $(GSLLIB) $(G2CLIB) $(LAPACKLIB) \
           $(G2CLIB) $(HDF5LIB) $(XMLLIB) -lm 
    INCL = $(BLITZINC) $(SPRNGINC) $(GSLINC) $(HDF5INC) $(XMLINC)
    CC = mpiCC
    LD = mpiCC  -Bstatic 
    CCFLAGS = -c -g  -Wno-deprecated  #-pg 
endif

MAKECC = g++

# Gets the subversion revision number
VER = \"`svn info | grep Revision | sed -e 's/Revision: //'`\"





DEFS = $(EXTRADEFS) -DVERSION=$(VER) -DTHREE_D -DNO_COUT -O3 #-DUSE_MPI#-DDEBUG -DBZ_DEBUG # -DPARALLEL -DUSE_MPI # -DDEBUG -DBZ_DEBUG #-ffast-math#  -DDEBUG -DBZ_DEBUG  # -DUSE_MPI #  DPARALLEL  # -DDEBUG -DBZ_DEBUG  -g #-DUSE_MPI 

PIMCobjs =                            \
  Main.o                              \
  Observables/ObservableEnergy.o      \
  Observables/StructureFactor.o       \
  Moves/BisectionClass.o              \
  Moves/NoPermuteClass.o              \
  Moves/PermuteStageClass.o           \
  Moves/BisectionBlock.o              \
  Moves/RefSliceMove.o                \
  Moves/BisectionStageClass.o         \
  PIMCClass.o                         \
  Moves/MetaMoves.o 	              \
  Moves/BlockMove.o                   \
  Observables/ObservableBase.o        \
  Observables/ObservableCorrelation.o \
  Observables/PathDump.o              \
  Observables/WindingNumber.o         \
  SpeciesClass.o                      \
  Common.o                            \
  Moves/PermuteTableClass.o	      \
  Moves/TablePermuteStageClass.o      \
  Moves/RandomPermClass.o             \
  Moves/OpenBisectionMoveClass.o      \
  Moves/BisectionMoveClass.o          \
  Moves/MoveBase.o                    \
  Actions/ActionBase.o                \
  Actions/ShortRangeClass.o           \
  Actions/LongRangeClass.o            \
  Actions/ShortRangePotClass.o        \
  Actions/LongRangePotClass.o         \
  Actions/LongRangeRPAClass.o         \
  Actions/KineticClass.o              \
  Actions/ActionsClass.o              \
  Actions/NodalActionClass.o          \
  Moves/MultiStage.o                  \
  ActionClass.o                       \
  LongRangeRPA.o                      \
  PathDataClass.o                     \
  CommunicatorClass.o                 \
  PathClass.o                         \
  WrapClass.o			      \
  NodalAction.o                       \
  Common/Splines/CubicSpline.o        \
  Common/Splines/MyTricubicSpline.o   \
  Common/Splines/Grid.o               \
  Common/Splines/QuinticSpline.o      \
  Common/Splines/QuinticSplines.o     \
  Common/MPI/Communication.o	      \
  Common/IO/InputOutput.o             \
  Common/IO/InputOutputHDF5.o         \
  Common/IO/InputOutputASCII.o        \
  Common/IO/InputOutputXML.o          \
  Common/PairAction/PAcoulombBCFit.o  \
  Common/PairAction/PACoulombFit.o    \
  Common/PairAction/PADipoleFit.o     \
  Common/PairAction/PATripoleFit.o    \
  Common/PairAction/PAclassicalFit.o  \
  Common/PairAction/PAzeroFit.o       \
  Common/Splines/BicubicSpline.o      \
  Common/PH/Potential.o               \
  Common/PH/QuinticPH.o               \
  Common/PH/CoulombPot.o              \
  Common/PH/ScreenedPot.o             \
  Common/PH/SplinePot.o               \
  Common/PH/HeAzizPot.o               \
  Common/PairAction/DavidPAClass.o    \
  Common/Ewald/OptimizedBreakup.o     \
  Common/MatrixOps/MatrixOps.o        \
  Common/Integration/GKIntegration.o  \
  Common/Fitting/Fitting.o            \
  MirroredClass.o
#  Common/PairAction/PAcoulombFit.o   \
#  Common/PairAction/PAszFit.o        \
#  Common/PairAction/PAsFit.o         \
#  Common/PairAction/PAtricubicFit.o  \


TestPermobjs =                        \
  TestPermutation.o                   \
  Moves/BisectionClass.o              \
  Moves/PermuteStageClass.o           \
  Moves/BisectionBlock.o              \
  Moves/RefSliceMove.o                \
  Moves/NoPermuteClass.o              \
  Moves/BlockMove.o                   \
  Moves/MetaMoves.o                   \
  Moves/TablePermuteStageClass.o      \
  Moves/BisectionStageClass.o         \
  Actions/ActionBase.o                \
  Actions/ShortRangeClass.o           \
  Actions/ActionsClass.o              \
  Actions/LongRangeClass.o            \
  Actions/LongRangeRPAClass.o         \
  Actions/ShortRangePotClass.o        \
  Actions/LongRangePotClass.o         \
  Actions/KineticClass.o              \
  Actions/NodalActionClass.o          \
  Moves/MultiStage.o                  \
  PIMCClass.o                         \
  Observables/ObservableBase.o        \
  Observables/ObservableEnergy.o      \
  Observables/ObservableCorrelation.o \
  Observables/StructureFactor.o       \
  Observables/PathDump.o              \
  Observables/WindingNumber.o         \
  SpeciesClass.o                      \
  Common.o                            \
  Moves/PermuteTableClass.o	      \
  Moves/RandomPermClass.o             \
  Moves/BisectionMoveClass.o          \
  Moves/OpenBisectionMoveClass.o      \
  Moves/MoveBase.o                    \
  ActionClass.o                       \
  LongRangeRPA.o                      \
  PathDataClass.o                     \
  CommunicatorClass.o                 \
  PathClass.o                         \
  WrapClass.o			      \
  NodalAction.o                       \
  Common/Splines/CubicSpline.o        \
  Common/Splines/MyTricubicSpline.o   \
  Common/Splines/Grid.o               \
  Common/Splines/QuinticSpline.o      \
  Common/Splines/QuinticSplines.o     \
  Common/MPI/Communication.o	      \
  Common/IO/InputOutput.o             \
  Common/IO/InputOutputHDF5.o         \
  Common/IO/InputOutputASCII.o        \
  Common/IO/InputOutputXML.o          \
  Common/PairAction/PAcoulombBCFit.o  \
  Common/PairAction/PACoulombFit.o    \
  Common/PairAction/PADipoleFit.o     \
  Common/PairAction/PATripoleFit.o    \
  Common/PairAction/PAclassicalFit.o  \
  Common/PairAction/PAzeroFit.o       \
  Common/Splines/BicubicSpline.o      \
  Common/PH/Potential.o               \
  Common/PH/QuinticPH.o               \
  Common/PH/CoulombPot.o              \
  Common/PH/ScreenedPot.o             \
  Common/PH/SplinePot.o               \
  Common/PH/HeAzizPot.o               \
  Common/PairAction/DavidPAClass.o    \
  Common/Ewald/OptimizedBreakup.o     \
  Common/MatrixOps/MatrixOps.o        \
  Common/Integration/GKIntegration.o  \
  Common/Fitting/Fitting.o            \
  MirroredClass.o  
#  Common/PairAction/PAcoulombFit.o   \
#  Common/PairAction/PAszFit.o        \
#  Common/PairAction/PAsFit.o         \
#  Common/PairAction/PAtricubicFit.o  \


TestEwaldobjs =                       \
  TestEwald.o                         \
  PathClass.o                         \
  MirroredClass.o                     \
  SpeciesClass.o                      \
  ActionClass.o                       \
  NodalAction.o                       \
  PathDataClass.o                     \
  LongRangeRPA.o                      \
  Actions/ActionBase.o                \
  Actions/ShortRangeClass.o           \
  Actions/LongRangeClass.o            \
  Actions/LongRangeRPAClass.o         \
  Actions/ShortRangePotClass.o        \
  Actions/LongRangePotClass.o         \
  Actions/KineticClass.o              \
  Actions/NodalActionClass.o          \
  Common/MPI/Communication.o	      \
  Common/IO/InputOutput.o             \
  Common/IO/InputOutputHDF5.o         \
  Common/IO/InputOutputASCII.o        \
  Common/IO/InputOutputXML.o          \
  Common/PairAction/PAcoulombBCFit.o  \
  Common/PairAction/PACoulombFit.o    \
  Common/PairAction/PADipoleFit.o     \
  Common/PairAction/PATripoleFit.o    \
  Common/PairAction/PAclassicalFit.o  \
  Common/PairAction/PAzeroFit.o       \
  Common/PairAction/DavidPAClass.o    \
  Common/Splines/BicubicSpline.o      \
  Common/PH/Potential.o               \
  Common/PH/QuinticPH.o               \
  Common/PH/CoulombPot.o              \
  Common/PH/ScreenedPot.o             \
  Common/PH/SplinePot.o               \
  Common/PH/HeAzizPot.o               \
  Common/Splines/CubicSpline.o        \
  Common/Splines/MyTricubicSpline.o   \
  Common/Splines/Grid.o               \
  Common/Splines/QuinticSpline.o      \
  Common/Ewald/OptimizedBreakup.o     \
  Common/MatrixOps/MatrixOps.o        \
  Common/Integration/GKIntegration.o  \
  Common/Fitting/Fitting.o            \
  Common/Splines/QuinticSplines.o    


FreeParticleObjs =                   \
  FreeParticles.o                    \
  Common/IO/InputOutput.o            \
  Common/IO/InputOutputHDF5.o        \
  Common/IO/InputOutputASCII.o       \
  Common/IO/InputOutputXML.o         \
  Common/MPI/Communication.o


PASS_DEFS = "CC=${CC}" "LD=${LD}" "CCFLAGS=${CCFLAGS}" "DEFS=${DEFS}" "INCL=${INCL}" "LIBS=${LIBS}" "F77=${F77}"\
  "MAKECC=${MAKECC}"

MAKE_ALL = $(MAKE) all $(PASS_DEFS)
MAKE_NEWMAKE = $(MAKE) -f template.make newmake $(PASS_DEFS)


all:    pimc++ TestPerm TestEwald FreeParticles

pimc++: Common_obj observables moves actions Tests $(PIMCobjs)
	$(LD) -o $@ $(PIMCobjs) $(LIBS) $(PSPLINELIB)

TestPerm: Common_obj Tests $(TestPermobjs)
	$(LD) -o $@ $(TestPermobjs) $(LIBS) $(PSPLINELIB)

TestEwald: Common_obj Tests $(TestEwaldobjs)
	$(LD) -o $@ $(TestEwaldobjs) $(LIBS) $(PSPLINELIB)

FreeParticles: Common_obj $(FreeParticleObjs)
	$(LD) -o $@ $(FreeParticleObjs) $(LIBS) $(PSPLINELIB)

Common_obj:
	cd Common; ${MAKE_ALL}

observables:
	cd Observables; ${MAKE_ALL}

moves:
	cd Moves; ${MAKE_ALL}

actions:
	cd Actions; ${MAKE_ALL}

Common_clean:
	cd Common; ${MAKE} clean

CodeTests:    
	cd Tests; ${MAKE_ALL}



TestHDF5:	Common_obj TestHDF5.o Common/IO/InputOutput.o Common/IO/InputOutputHDF5.o Common/IO/InputOutputASCII.o
	$(LD) -o $@ TestHDF5.o Common/IO/InputOutput.o Common/IO/InputOutputHDF5.o Common/IO/InputOutputASCII.o $(LIBS)

TestASCII:	Common_obj TestASCII.o Common/IO/InputOutput.o Common/IO/InputOutputASCII.o Common/IO/InputOutputHDF5.o
	$(LD) -o $@ TestASCII.o Common/IO/InputOutput.o Common/IO/InputOutputASCII.o Common/IO/InputOutputHDF5.o $(LIBS)

#TestSubarrays: 	$(TestSubarrayObjs)
#	pushd ..; make; pushd
#	$(LD) -o $@ $(TestSubarrayObjs) $(LIBS)

clean:	Common_clean
	rm *.o

.cc.o: $(HEADERS)
	$(CC) $(CCFLAGS) $(DEFS) $(INCL) $<
.f.o:
	g77 -c $<



SOURCES =  myprog.cc SpeciesClass.cc Common.cc ActionClass.cc PathDataClass.cc  CommunicatorClass.cc PathClass.cc TestSubarrays.cc  WrapClass.cc TestHDF5.cc TestASCII.cc  Main.cc PIMCClass.cc TestPermutation.cc MirroredClass.cc TestEwald.cc LongRangeRPA.cc NodalAction.cc FreeParticles.cc


newmake: Common_newmake Tests_newmake Observables_newmake Moves_newmake Actions_newmake
	$(MAKE) -f template.make Makefile FRC=force_rebuild

Common_newmake:
	cd Common; $(MAKE_NEWMAKE)

Tests_newmake:
	cd Tests; $(MAKE_NEWMAKE)

Observables_newmake:
	cd Observables; ${MAKE_NEWMAKE}

Moves_newmake:
	cd Moves; ${MAKE_NEWMAKE}

Actions_newmake:
	cd Actions; ${MAKE_NEWMAKE}

Makefile:	$(FRC)
	rm -f $@
	cp template.make $@
	echo 'Automatically generated dependency list:' >> $@
	$(MAKECC) $(CCFLAGS) $(INCL) -MM $(SOURCES) >> $@
	chmod -w $@


force_rebuild:



