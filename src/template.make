include /home/common/Codes/Make.include

LIBS = $(BLITZLIB) $(SPRNGLIB) $(GSLLIB) $(G2CLIB) $(LAPACKLIB) $(G2CLIB) $(HDF5LIB) -lm #-lstdc++
INCL = $(BLITZINC) $(SPRNGINC) $(GSLINC) $(HDF5INC) 

CCFLAGS = -c -g  -Wno-deprecated -DBZ_DEBUG # -O3   #-pg
CC = mpiCC
LD = mpiCC  -Bstatic 
DEFS = -DNO_COUT -DUSE_MPI  -DDEBUG #-DBZ_DEBUG  -g #-DUSE_MPI 

TestObjs = ObservableClass.o CubicSpline.o Grid.o InputFile.o SpeciesClass.o Common.o BisectionMoveClass.o MoveClass.o ActionClass.o PathDataClass.o   CommunicatorClass.o PathClass.o test.o DistanceTablePBCClass.o DistanceTableFreeClass.o DistanceTableClass.o MirroredArrayClass.o InputOutput.o InputOutputHDF5.o

TestSubarrayObjs = TestSubarrays.o

Test: 	$(TestObjs)
	pushd ..; make; pushd
	$(LD) -o $@ $(TestObjs) $(LIBS)

TestHDF5:	TestHDF5.o InputOutput.o InputOutputHDF5.o
	$(LD) -o $@ TestHDF5.o InputOutput.o InputOutputHDF5.o $(LIBS)

TestSubarrays: 	$(TestSubarrayObjs)
	pushd ..; make; pushd
	$(LD) -o $@ $(TestSubarrayObjs) $(LIBS)

.cc.o: $(HEADERS)
	$(CC) $(CCFLAGS) $(DEFS) $(INCL) $<
.f.o:
	g77 -c $<


SOURCES = ObservableClass.cc CubicSpline.cc Grid.cc InputFile.cc myprog.cc SpeciesClass.cc Common.cc BisectionMoveClass.cc MoveClass.cc ActionClass.cc PathDataClass.cc  MirroredArrayClass.cc CommunicatorClass.cc PathClass.cc test.cc TestSubarrays.cc DistanceTablePBCClass.cc DistanceTableFreeClass.cc DistanceTableClass.cc InputOutput.cc InputOutputHDF5.cc TestHDF5.cc


newmake: 
	make -f template.make Makefile FRC=force_rebuild

Makefile:	$(FRC)
	rm -f $@
	cp template.make $@
	echo 'Automatically generated dependency list:' >> $@
	$(CC) $(CCFLAGS) $(INCL) -M $(SOURCES) >> $@
	chmod -w $@


force_rebuild:



