include /home/common/Codes/Make.include

LIBS = $(BLITZLIB) $(SPRNGLIB) $(GSLLIB) $(G2CLIB) $(LAPACKLIB) \
       $(G2CLIB) $(HDF5LIB) -lm 
INCL = $(BLITZINC) $(SPRNGINC) $(GSLINC) $(HDF5INC) 

CCFLAGS = -c -g  -Wno-deprecated  #-O3   #-pg
CC = mpicc
LD = mpicc  -Bstatic 
DEFS = -DNO_COUT -DUSE_MPI   -DDEBUG -DBZ_DEBUG  -g #-DUSE_MPI 

TestObjs =                      \
  ObservableClass.o             \
  Common/Splines/CubicSpline.o  \
  Common/Splines/Grid.o         \
  SpeciesClass.o                \
  Common.o                      \
  BisectionMoveClass.o          \
  MoveClass.o                   \
  ActionClass.o                 \
  PathDataClass.o               \
  CommunicatorClass.o           \
  PathClass.o                   \
  test.o                        \
  DistanceTablePBCClass.o       \
  DistanceTableFreeClass.o      \
  DistanceTableClass.o          \
  MirroredArrayClass.o          \
  Common/IO/InputOutput.o       \
  Common/IO/InputOutputHDF5.o   \
  Common/IO/InputFile.o         \
  Common/IO/InputOutputASCII.o 

MakeInputObjs =                 \
  Common/IO/InputOutput.o       \
  Common/IO/InputOutputHDF5.o   \
  Common/IO/InputOutputASCII.o  \
  makeInput.o

PASS_DEFS = "CC=${CC}" "LD=${LD}" "CCFLAGS=${CCFLAGS}" "DEFS=${DEFS}" "INCL=${INCL}" "LIBS=${LIBS}"

MAKE_ALL = $(MAKE) all $(PASS_DEFS)
MAKE_NEWMAKE = $(MAKE) -f template.make newmake $(PASS_DEFS)



Test: 	Common_obj $(TestObjs)
	pushd ..; make; pushd
	$(LD) -o $@ $(TestObjs) $(LIBS)

MakeInput: 	Common_obj $(MakeInputObjs)
	pushd ..; make; pushd
	$(LD) -o $@ $(MakeInputObjs) $(LIBS)

Common_obj:
	cd Common; ${MAKE_ALL}

Common_clean:
	cd Common; ${MAKE} clean

TestHDF5:	Common_obj TestHDF5.o Common/IO/InputOutput.o Common/IO/InputOutputHDF5.o
	$(LD) -o $@ TestHDF5.o Common/IO/InputOutput.o Common/IO/InputOutputHDF5.o $(LIBS)

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


SOURCES = ObservableClass.cc myprog.cc SpeciesClass.cc Common.cc BisectionMoveClass.cc MoveClass.cc ActionClass.cc PathDataClass.cc  MirroredArrayClass.cc CommunicatorClass.cc PathClass.cc test.cc TestSubarrays.cc DistanceTablePBCClass.cc DistanceTableFreeClass.cc DistanceTableClass.cc TestHDF5.cc TestASCII.cc

newmake: Common_newmake
	make -f template.make Makefile FRC=force_rebuild

Common_newmake:
	cd Common; $(MAKE_NEWMAKE)

Makefile:	$(FRC)
	rm -f $@
	cp template.make $@
	echo 'Automatically generated dependency list:' >> $@
	$(CC) $(CCFLAGS) $(INCL) -M $(SOURCES) >> $@
	chmod -w $@


force_rebuild:



