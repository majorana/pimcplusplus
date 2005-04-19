PSPLINELIB = -L$(PWD)/Common/Splines/fortran -lpspline

ifeq ($(HOSTTYPE),alpha)
    include /usr/users/0/kesler/lib/Make.include
    MPILIB = -lmpi -lelan
    LIBS =  $(LAPACKLIB) $(BLITZLIB) $(SPRNGLIB) $(GSLLIB) \
         $(G2CLIB) $(HDF5LIB) $(XMLLIB) $(MPILIB) -lm 
    INCL = $(BLITZINC) $(SPRNGINC) $(GSLINC) $(HDF5INC) $(XMLINC)
    MAKE = gmake
    CC = g++
    LD = g++ #-static
    F77 = g77
    EXTRADEFS = -DNOCUSERID
    CCFLAGS = -c -g  -Wno-deprecated  -O3 #-pg 
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

ifeq ($(HOSTTYPE),powermac)
   include /turing/home/esler/lib/Make.include
   LIBS = $(LAPACKLIB) $(BLITZLIB) $(SPRNGLIB) $(GSLLIB) \
          $(HDF5LIB) $(XMLLIB) $(PSPLINELIB) $(FORTRANLIB) \
          $(FFTW3LIB) $(CBLASLIB) -lm
   INCL = $(BLITZINC) $(SPRNGINC) $(GSLINC) $(HDF5INC) $(XMLINC) \
          $(FFTW3INC) $(CBLASINC)
   CC = mpiCC
   LD = mpiCC
   F77 = f77
   CCFLAGS = -c -g  -Wno-long-double
   EXTRADEFS = -DNOUNDERSCORE -DNOCUSERID -DMAC
   MAKE = make
endif
ifeq ($(HOSTTYPE),i386-linux)
    ifeq ($(GROUP),tvi)
       include /u/ac/esler/lib/Make.include
       CC = cmpic++ -ccl icpc
       LD = cmpic++ -ccl icpc 
#       CC = icc
#       LD = icc
       F77 = ifort
       EXTRADEFS = -DUSE_MKL -w1 -wr654,1011
    else
       include /usr/lib/Make.include	
       CC = distcc mpiCC
       LD = mpiCC  -Bstatic 
       F77 = g77 
       EXTRADEFS = -Wno-deprecated -march=athlon -mcpu=athlon -ffast-math
    endif
    LIBS = $(BLITZLIB) $(SPRNGLIB) $(GSLLIB) $(G2CLIB) $(LAPACKLIB) \
           $(G2CLIB) $(HDF5LIB) $(XMLLIB) $(FFTW3LIB) $(CBLASLIB) \
           $(FORTLIB) -lm 
    INCL = $(BLITZINC) $(SPRNGINC) $(GSLINC) $(HDF5INC) $(XMLINC) \
           $(FFTW3INC) $(CBLASINC)
#    CC = mpiCC


    CCFLAGS = -c -g  #-pg 
endif

ifeq ($(GROUP),tvi)
   MAKECC = icc
else
  MAKECC = g++
endif

# Gets the subversion revision number
VER = \"`svn info | grep Revision | sed -e 's/Revision: //'`\"
COMMONVER = \"`svn info Common | grep Revision | sed -e 's/Revision: //'`\"


DEFS = $(EXTRADEFS) -DNDIM=3 -DVERSION=$(VER)  -DNO_COUT  -O3 -DUSE_MPI #-DBZ_DEBUG -DDEBUG ##-ffast-math#  -DDEBUG -DBZ_DEBUG  # -DUSE_MPI #  DPARALLEL  # -DDEBUG -DBZ_DEBUG  -g #-DUSE_MPI -DTHREE_D 


PIMCobjs =                                \
  Common.o                                \
  Main.o                                  \
  Observables/ObservableEnergy.o          \
  Observables/Time.o                      \
  Observables/DistanceToHead.o            \
  Observables/ObservableModifiedEnergy.o  \
  Observables/Weight.o                    \
  Observables/StructureFactor.o           \
  Moves/NoPermuteClass.o                  \
  Moves/EndStageClass.o                   \
  Moves/PermuteStageClass.o               \
  Moves/CoupledPermuteStageClass.o        \
  Moves/BisectionBlock.o                  \
  Moves/RefSliceMove.o                    \
  Moves/BisectionStageClass.o             \
  Moves/StructureReject.o                 \
  PIMCClass.o                             \
  Moves/MetaMoves.o 	                  \
  Observables/ObservableBase.o            \
  Observables/ObservableCorrelation.o     \
  Observables/PathDump.o                  \
  Observables/WindingNumber.o             \
  SpeciesClass.o                          \
  Moves/PermuteTableClass.o	          \
  Moves/TablePermuteStageClass.o          \
  Moves/RandomPermClass.o                 \
  Moves/DisplaceMove.o                    \
  Moves/OpenEndMove.o                     \
  Moves/MoveBase.o                        \
  Moves/WaterMove.o                       \
  Moves/WaterMoveRing.o                   \
  Actions/StructureReject.o               \
  Actions/ActionBase.o                    \
  Actions/ShortRangeClass.o               \
  Actions/OpenLoopImportance.o            \
  Actions/ShortRangeApproximateClass.o    \
  Actions/LongRangeClass.o                \
  Actions/ShortRangePotClass.o            \
  Actions/LongRangePotClass.o             \
  Actions/LongRangeRPAClass.o             \
  Actions/DavidLongRangeClass.o           \
  Actions/KineticClass.o                  \
  Actions/ActionsClass.o                  \
  Actions/NodalActionClass.o              \
  Actions/GroundStateNodalActionClass.o   \
  Actions/TIP5PWaterClass.o               \
  Actions/ST2WaterClass.o                 \
  Moves/MultiStage.o                      \
  LongRangeRPA.o                          \
  PathDataClass.o                         \
  CommunicatorClass.o                     \
  PathClass.o                             \
  WrapClass.o			          \
  Common/Splines/CubicSpline.o            \
  Common/Splines/MyTricubicSpline.o       \
  Common/Splines/Grid.o                   \
  Common/Splines/QuinticSpline.o          \
  Common/Splines/QuinticSplines.o         \
  Common/Splines/MultiTricubicSpline3.o   \
  Common/MPI/Communication.o	          \
  Common/IO/InputOutput.o                 \
  Common/IO/InputOutputHDF5.o             \
  Common/IO/InputOutputASCII.o            \
  Common/IO/InputOutputXML.o              \
  Common/PairAction/PAcoulombBCFit.o      \
  Common/PairAction/PAtricubicFit.o       \
  Common/PairAction/PACoulombFit.o        \
  Common/PairAction/PADipoleFit.o         \
  Common/PairAction/PATripoleFit.o        \
  Common/PairAction/PAclassicalFit.o      \
  Common/PairAction/PAzeroFit.o           \
  Common/Splines/BicubicSpline.o          \
  Common/PH/Potential.o                   \
  Common/PH/QuinticPH.o                   \
  Common/PH/CoulombPot.o                  \
  Common/PH/ScreenedPot.o                 \
  Common/PH/SplinePot.o                   \
  Common/PH/HeAzizPot.o                   \
  Common/PH/kSpacePH.o                    \
  Common/PairAction/DavidPAClass.o        \
  Common/Ewald/OptimizedBreakup.o         \
  Common/MatrixOps/MatrixOps.o            \
  Common/Integration/GKIntegration.o      \
  Common/Fitting/Fitting.o                \
  Common/PlaneWavePHDFT/PlaneWaves.o      \
  Common/PlaneWavePHDFT/ConjGrad2.o       \
  Common/PlaneWavePHDFT/FFTBox.o          \
  Common/PlaneWavePHDFT/GVecs.o           \
  Common/PlaneWavePHDFT/Hamiltonian2.o    \
  Common/FFT/FFT.o                        \
  MirroredClass.o
#  Common/PairAction/PAcoulombFit.o       \
#  Common/PairAction/PAszFit.o            \
#  Common/PairAction/PAsFit.o             \
#  Common/PairAction/PAtricubicFit.o      \
#  Moves/BisectionClass.o                 \
#  Moves/BlockMove.o                      \
#  Moves/OpenBisectionMoveClass.o         \
#  Moves/BisectionMoveClass.o             \
#  ActionClass.o                          \
#  NodalAction.o                          \


TestPermobjs =                            \
  Common.o                                \
  TestPermutation.o                       \
  Moves/PermuteStageClass.o               \
  Moves/EndStageClass.o                   \
  Moves/BisectionBlock.o                  \
  Moves/DisplaceMove.o                    \
  Moves/RefSliceMove.o                    \
  Moves/NoPermuteClass.o                  \
  Moves/MetaMoves.o                       \
  Moves/TablePermuteStageClass.o          \
  Moves/CoupledPermuteStageClass.o        \
  Moves/BisectionStageClass.o             \
  Moves/OpenEndMove.o                     \
  Actions/ActionBase.o                    \
  Actions/ShortRangeClass.o               \
  Actions/ShortRangeApproximateClass.o    \
  Actions/ActionsClass.o                  \
  Actions/LongRangeClass.o                \
  Actions/LongRangeRPAClass.o             \
  Actions/ShortRangePotClass.o            \
  Actions/LongRangePotClass.o             \
  Actions/DavidLongRangeClass.o           \
  Actions/KineticClass.o                  \
  Actions/NodalActionClass.o              \
  Actions/GroundStateNodalActionClass.o   \
  Moves/MultiStage.o                      \
  PIMCClass.o                             \
  Observables/ObservableBase.o            \
  Observables/ObservableEnergy.o          \
  Observables/ObservableModifiedEnergy.o  \
  Observables/Weight.o                    \
  Observables/ObservableCorrelation.o     \
  Observables/StructureFactor.o           \
  Observables/PathDump.o                  \
  Observables/WindingNumber.o             \
  SpeciesClass.o                          \
  Moves/PermuteTableClass.o	          \
  Moves/RandomPermClass.o                 \
  Moves/MoveBase.o                        \
  LongRangeRPA.o                          \
  PathDataClass.o                         \
  CommunicatorClass.o                     \
  PathClass.o                             \
  WrapClass.o			          \
  Common/Splines/CubicSpline.o            \
  Common/Splines/MyTricubicSpline.o       \
  Common/Splines/Grid.o                   \
  Common/Splines/QuinticSpline.o          \
  Common/Splines/QuinticSplines.o         \
  Common/MPI/Communication.o	          \
  Common/IO/InputOutput.o                 \
  Common/IO/InputOutputHDF5.o             \
  Common/IO/InputOutputASCII.o            \
  Common/IO/InputOutputXML.o              \
  Common/PairAction/PAcoulombBCFit.o      \
  Common/PairAction/PAtricubicFit.o       \
  Common/PairAction/PACoulombFit.o        \
  Common/PairAction/PADipoleFit.o         \
  Common/PairAction/PATripoleFit.o        \
  Common/PairAction/PAclassicalFit.o      \
  Common/PairAction/PAzeroFit.o           \
  Common/Splines/BicubicSpline.o          \
  Common/PH/Potential.o                   \
  Common/PH/QuinticPH.o                   \
  Common/PH/CoulombPot.o                  \
  Common/PH/ScreenedPot.o                 \
  Common/PH/SplinePot.o                   \
  Common/PH/HeAzizPot.o                   \
  Common/PH/kSpacePH.o                    \
  Common/PairAction/DavidPAClass.o        \
  Common/Ewald/OptimizedBreakup.o         \
  Common/MatrixOps/MatrixOps.o            \
  Common/Integration/GKIntegration.o      \
  Common/Fitting/Fitting.o                \
  Common/Splines/MultiTricubicSpline3.o   \
  Common/PlaneWavePHDFT/PlaneWaves.o      \
  Common/PlaneWavePHDFT/ConjGrad2.o       \
  Common/PlaneWavePHDFT/FFTBox.o          \
  Common/PlaneWavePHDFT/GVecs.o           \
  Common/PlaneWavePHDFT/Hamiltonian2.o    \
  Common/FFT/FFT.o                        \
  MirroredClass.o  
#  Moves/BlockMove.o                      \
#  Moves/BisectionMoveClass.o             \
#  Moves/OpenBisectionMoveClass.o         \
#  Moves/BisectionClass.o                 \
#  Common/PairAction/PAcoulombFit.o       \
#  Common/PairAction/PAszFit.o            \
#  Common/PairAction/PAsFit.o             \
#  Common/PairAction/PAtricubicFit.o      \
#  ActionClass.o                          \
#  NodalAction.o                          \



TestEwaldobjs =                           \
  TestEwald.o                             \
  PathClass.o                             \
  MirroredClass.o                         \
  SpeciesClass.o                          \
  PathDataClass.o                         \
  LongRangeRPA.o                          \
  Actions/ActionBase.o                    \
  Actions/ShortRangeClass.o               \
  Actions/ShortRangeApproximateClass.o    \
  Actions/LongRangeClass.o                \
  Actions/DavidLongRangeClass.o           \
  Actions/LongRangeRPAClass.o             \
  Actions/ShortRangePotClass.o            \
  Actions/LongRangePotClass.o             \
  Actions/KineticClass.o                  \
  Actions/NodalActionClass.o              \
  Common/MPI/Communication.o	          \
  Common/IO/InputOutput.o                 \
  Common/IO/InputOutputHDF5.o             \
  Common/IO/InputOutputASCII.o            \
  Common/IO/InputOutputXML.o              \
  Common/PairAction/PAtricubicFit.o       \
  Common/PairAction/PAcoulombBCFit.o      \
  Common/PairAction/PACoulombFit.o        \
  Common/PairAction/PADipoleFit.o         \
  Common/PairAction/PATripoleFit.o        \
  Common/PairAction/PAclassicalFit.o      \
  Common/PairAction/PAzeroFit.o           \
  Common/PairAction/DavidPAClass.o        \
  Common/Splines/BicubicSpline.o          \
  Common/PH/Potential.o                   \
  Common/PH/QuinticPH.o                   \
  Common/PH/CoulombPot.o                  \
  Common/PH/ScreenedPot.o                 \
  Common/PH/SplinePot.o                   \
  Common/PH/HeAzizPot.o                   \
  Common/Splines/CubicSpline.o            \
  Common/Splines/MyTricubicSpline.o       \
  Common/Splines/Grid.o                   \
  Common/Splines/QuinticSpline.o          \
  Common/Ewald/OptimizedBreakup.o         \
  Common/MatrixOps/MatrixOps.o            \
  Common/Integration/GKIntegration.o      \
  Common/Fitting/Fitting.o                \
  Common/Splines/MultiTricubicSpline3.o   \
  Common/PlaneWavePHDFT/ConjGrad2.o       \
  Common/PlaneWavePHDFT/FFTBox.o          \
  Common/PlaneWavePHDFT/GVecs.o           \
  Common/PlaneWavePHDFT/Hamiltonian2.o    \
  Common/PlaneWavePHDFT/PlaneWaves.o      \
  Common/FFT/FFT.o                        \
  Common/Splines/QuinticSplines.o    
#  ActionClass.o                          \
#  NodalAction.o                          \



FreeParticleObjs =                       \
  FreeParticles.o                        \
  Common/IO/InputOutput.o                \
  Common/IO/InputOutputHDF5.o            \
  Common/IO/InputOutputASCII.o           \
  Common/IO/InputOutputXML.o             \
  Common/MPI/Communication.o


PASS_DEFS = "CC=${CC}" "LD=${LD}" "CCFLAGS=${CCFLAGS}" "DEFS=${DEFS}" "INCL=${INCL}" "LIBS=${LIBS}" "F77=${F77}"\
  "MAKECC=${MAKECC}"

MAKE_ALL = $(MAKE) all $(PASS_DEFS)
MAKE_NEWMAKE = $(MAKE) -f template.make newmake $(PASS_DEFS)


all:    pimc++  FreeParticles # Visual_obj #TestPerm TestEwald 

pimc++: Common_obj Observables_obj Moves_obj Actions_obj Tests $(PIMCobjs)
	$(LD) -o $@ $(PIMCobjs) $(LIBS) $(PSPLINELIB)

TestPerm: Common_obj Tests $(TestPermobjs)
	$(LD) -o $@ $(TestPermobjs) $(LIBS) $(PSPLINELIB)

TestEwald: Common_obj Tests $(TestEwaldobjs)
	$(LD) -o $@ $(TestEwaldobjs) $(LIBS) $(PSPLINELIB)

FreeParticles: Common_obj $(FreeParticleObjs)
	$(LD) -o $@ $(FreeParticleObjs) $(LIBS) $(PSPLINELIB)

Common_obj:
	+cd Common; ${MAKE_ALL}

Visual_obj:
	+cd Visual; ${MAKE_ALL}

Observables_obj:
	+cd Observables; ${MAKE_ALL}

Moves_obj:
	+cd Moves; ${MAKE_ALL}

Actions_obj:
	+cd Actions; ${MAKE_ALL}

Common_clean:
	+cd Common; ${MAKE} clean

Visual_clean:
	+cd Visual; ${MAKE} clean

Tests_clean:
	+cd Tests; ${MAKE} clean

Actions_clean:
	+cd Actions; ${MAKE} clean

Moves_clean:
	+cd Moves; ${MAKE} clean

Observables_clean:
	+cd Observables; ${MAKE} clean

CodeTests:    
	+cd Tests; ${MAKE_ALL}



TestHDF5:	Common_obj TestHDF5.o Common/IO/InputOutput.o Common/IO/InputOutputHDF5.o Common/IO/InputOutputASCII.o
	$(LD) -o $@ TestHDF5.o Common/IO/InputOutput.o Common/IO/InputOutputHDF5.o Common/IO/InputOutputASCII.o $(LIBS)

TestASCII:	Common_obj TestASCII.o Common/IO/InputOutput.o Common/IO/InputOutputASCII.o Common/IO/InputOutputHDF5.o
	$(LD) -o $@ TestASCII.o Common/IO/InputOutput.o Common/IO/InputOutputASCII.o Common/IO/InputOutputHDF5.o $(LIBS)

#TestSubarrays: 	$(TestSubarrayObjs)
#	pushd ..; make; pushd
#	$(LD) -o $@ $(TestSubarrayObjs) $(LIBS)

clean:	Common_clean Tests_clean Actions_clean Moves_clean Observables_clean Visual_clean
	rm *.o

.cc.o: $(HEADERS)
	$(CC) $(CCFLAGS) $(DEFS) $(INCL) $<
.f.o:
	g77 -c $<



SOURCES =  Common.cc myprog.cc SpeciesClass.cc  ActionClass.cc PathDataClass.cc  CommunicatorClass.cc PathClass.cc TestSubarrays.cc  WrapClass.cc TestHDF5.cc TestASCII.cc  Main.cc PIMCClass.cc TestPermutation.cc MirroredClass.cc TestEwald.cc LongRangeRPA.cc NodalAction.cc FreeParticles.cc


newmake: Common_newmake Tests_newmake Observables_newmake Moves_newmake Actions_newmake 
	+$(MAKE) -f template.make Makefile FRC=force_rebuild

Common_newmake:
	+cd Common; $(MAKE_NEWMAKE)

Visual_newmake:
	+cd Visual; $(MAKE_NEWMAKE)

Tests_newmake:
	+cd Tests; $(MAKE_NEWMAKE)

Observables_newmake:
	+cd Observables; ${MAKE_NEWMAKE}

Moves_newmake:
	+cd Moves; ${MAKE_NEWMAKE}

Actions_newmake:
	+cd Actions; ${MAKE_NEWMAKE}

Makefile:	$(FRC)
	rm -f $@
	cp template.make $@
	echo 'Automatically generated dependency list:' >> $@
	$(MAKECC) $(CCFLAGS) $(DEFS) $(INCL) -MM $(SOURCES) >> $@
	chmod -w $@


force_rebuild:



