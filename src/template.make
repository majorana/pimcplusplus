include /home/common/Codes/Make.include

PSPLINELIB = -L$(PWD)/Common/Splines/fortran -lpspline
LIBS = $(BLITZLIB) $(SPRNGLIB) $(GSLLIB) $(G2CLIB) $(LAPACKLIB) \
       $(G2CLIB) $(HDF5LIB) $(XMLLIB) -lm 
INCL = $(BLITZINC) $(SPRNGINC) $(GSLINC) $(HDF5INC) $(XMLINC)

CCFLAGS = -c -g  -Wno-deprecated  #-pg
CC = mpiCC
LD = mpiCC  -Bstatic 
DEFS = -DTHREE_D -DNO_COUT -O3 #-DDEBUG -DBZ_DEBUG #-ffast-math#  -DDEBUG -DBZ_DEBUG  # -DUSE_MPI #  DPARALLEL  # -DDEBUG -DBZ_DEBUG  -g #-DUSE_MPI 

PIMCobjs =                           \
  Main.o                             \
  ObservableEnergy.o		     \
  BisectionClass.o                   \
  PIMCClass.o                        \
  MetaMoves.o 			     \
  BlockMove.o                        \
  ObservableClass.o                  \
  Common/Splines/CubicSpline.o       \
  Common/Splines/MyTricubicSpline.o  \
  Common/Splines/Grid.o              \
  Common/Splines/QuinticSpline.o     \
  Common/Splines/QuinticSplines.o    \
  SpeciesClass.o                     \
  Common.o                           \
  PermuteTableClass.o		     \
  BisectionMoveClass.o               \
  MoveClass.o                        \
  ActionClass.o                      \
  PathDataClass.o                    \
  CommunicatorClass.o                \
  PathClass.o                        \
  MirroredArrayClass.o               \
  WrapClass.o			     \
  Common/MPI/Communication.o	     \
  Common/IO/InputOutput.o            \
  Common/IO/InputOutputHDF5.o        \
  Common/IO/InputFile.o              \
  Common/IO/InputOutputASCII.o       \
  Common/IO/InputOutputXML.o         \
  Common/PairAction/PAcoulombBCFit.o \
  Common/PairAction/PAclassicalFit.o \
  Common/PairAction/PAzeroFit.o      \
  Common/Splines/BicubicSpline.o     \
  Common/PH/Potential.o              \
  Common/PH/QuinticPH.o              \
  Common/PH/CoulombPot.o             \
  Common/PH/ScreenedPot.o            \
  Common/PH/SplinePot.o              \
  Common/PH/HeAzizPot.o              \
  Common/PairAction/DavidPAClass.o   \
  Common/Ewald/OptimizedBreakup.o    \
  Common/MatrixOps/MatrixOps.o       \
  Common/Integration/GKIntegration.o \
  MirroredClass.o
#  Common/PairAction/PAcoulombFit.o   \
#  Common/PairAction/PAszFit.o        \
#  Common/PairAction/PAsFit.o         \
#  Common/PairAction/PAtricubicFit.o  \


TestPermobjs =                       \
  TestPermutation.o                  \
  BisectionClass.o                   \
  BlockMove.o                        \
  MetaMoves.o                        \
  PIMCClass.o                        \
  ObservableClass.o                  \
  ObservableEnergy.o                 \
  Common/Splines/CubicSpline.o       \
  Common/Splines/MyTricubicSpline.o  \
  Common/Splines/Grid.o              \
  Common/Splines/QuinticSpline.o     \
  Common/Splines/QuinticSplines.o    \
  SpeciesClass.o                     \
  Common.o                           \
  PermuteTableClass.o		     \
  BisectionMoveClass.o               \
  MoveClass.o                        \
  ActionClass.o                      \
  PathDataClass.o                    \
  CommunicatorClass.o                \
  PathClass.o                        \
  MirroredArrayClass.o               \
  WrapClass.o			     \
  Common/MPI/Communication.o	     \
  Common/IO/InputOutput.o            \
  Common/IO/InputOutputHDF5.o        \
  Common/IO/InputFile.o              \
  Common/IO/InputOutputASCII.o       \
  Common/IO/InputOutputXML.o         \
  Common/PairAction/PAcoulombBCFit.o \
  Common/PairAction/PAclassicalFit.o \
  Common/PairAction/PAzeroFit.o      \
  Common/Splines/BicubicSpline.o     \
  Common/PH/Potential.o              \
  Common/PH/QuinticPH.o              \
  Common/PH/CoulombPot.o             \
  Common/PH/ScreenedPot.o            \
  Common/PH/SplinePot.o              \
  Common/PH/HeAzizPot.o              \
  Common/PairAction/DavidPAClass.o   \
  Common/Ewald/OptimizedBreakup.o    \
  Common/MatrixOps/MatrixOps.o       \
  Common/Integration/GKIntegration.o \
  MirroredClass.o  
#  Common/PairAction/PAcoulombFit.o   \
#  Common/PairAction/PAszFit.o        \
#  Common/PairAction/PAsFit.o         \
#  Common/PairAction/PAtricubicFit.o  \


TestEwaldobjs =                      \
  TestEwald.o                        \
  PathClass.o                        \
  MirroredClass.o                    \
  SpeciesClass.o                     \
  ActionClass.o                      \
  PathDataClass.o                    \
  Common/MPI/Communication.o	     \
  Common/IO/InputOutput.o            \
  Common/IO/InputOutputHDF5.o        \
  Common/IO/InputFile.o              \
  Common/IO/InputOutputASCII.o       \
  Common/IO/InputOutputXML.o         \
  Common/PairAction/PAcoulombBCFit.o \
  Common/PairAction/PAclassicalFit.o \
  Common/PairAction/PAzeroFit.o      \
  Common/PairAction/DavidPAClass.o   \
  Common/Splines/BicubicSpline.o     \
  Common/PH/Potential.o              \
  Common/PH/QuinticPH.o              \
  Common/PH/CoulombPot.o             \
  Common/PH/ScreenedPot.o            \
  Common/PH/SplinePot.o              \
  Common/PH/HeAzizPot.o              \
  Common/Splines/CubicSpline.o       \
  Common/Splines/MyTricubicSpline.o  \
  Common/Splines/Grid.o              \
  Common/Splines/QuinticSpline.o     \
  Common/Ewald/OptimizedBreakup.o    \
  Common/MatrixOps/MatrixOps.o       \
  Common/Integration/GKIntegration.o \
  Common/Splines/QuinticSplines.o    


PASS_DEFS = "CC=${CC}" "LD=${LD}" "CCFLAGS=${CCFLAGS}" "DEFS=${DEFS}" "INCL=${INCL}" "LIBS=${LIBS}"

MAKE_ALL = $(MAKE) all $(PASS_DEFS)
MAKE_NEWMAKE = $(MAKE) -f template.make newmake $(PASS_DEFS)


all:   pimc++ TestPerm TestEwald
	
pimc++: Common_obj Tests $(PIMCobjs)
	$(LD) -o $@ $(PIMCobjs) $(LIBS) $(PSPLINELIB)

TestPerm: Common_obj Tests $(TestPermobjs)
	$(LD) -o $@ $(TestPermobjs) $(LIBS) $(PSPLINELIB)

TestEwald: Common_obj Tests $(TestEwaldobjs)
	$(LD) -o $@ $(TestEwaldobjs) $(LIBS) $(PSPLINELIB)

Common_obj:
	cd Common; ${MAKE_ALL}

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


SOURCES = ObservableClass.cc ObservableEnergy.cc myprog.cc SpeciesClass.cc Common.cc BisectionMoveClass.cc MoveClass.cc ActionClass.cc PathDataClass.cc  MirroredArrayClass.cc CommunicatorClass.cc PathClass.cc TestSubarrays.cc DistanceTablePBCClass.cc DistanceTableFreeClass.cc DistanceTableClass.cc  WrapClass.cc TestHDF5.cc TestASCII.cc PermuteTableClass.cc Main.cc PIMCClass.cc TestPermutation.cc BisectionClass.cc  MetaMoves.cc BlockMove.cc MirroredClass.cc TestEwald.cc

newmake: Common_newmake Tests_newmake
	make -f template.make Makefile FRC=force_rebuild

Common_newmake:
	cd Common; $(MAKE_NEWMAKE)

Tests_newmake:
	cd Tests; $(MAKE_NEWMAKE)

Makefile:	$(FRC)
	rm -f $@
	cp template.make $@
	echo 'Automatically generated dependency list:' >> $@
	$(CC) $(CCFLAGS) $(INCL) -M $(SOURCES) >> $@
	chmod -w $@


force_rebuild:



