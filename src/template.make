include /home/common/Codes/Make.include

LIBS = $(BLITZLIB) $(SPRNGLIB) $(GSLLIB) $(G2CLIB) $(LAPACKLIB) $(G2CLIB) -lm #-lstdc++
INCL = $(BLITZINC) $(SPRNGINC) $(GSLINC) 

CCFLAGS = -c -g -O3 -DBZ_DEBUG -Wno-deprecated #-pg
CC = mpiCC
LD = mpiCC  -Bstatic 
DEFS = -DNO_COUT -DUSE_MPI -DBZ_DEBUG  -g #-DUSE_MPI 

TestObjs = ObservableClass.o CubicSpline.o Grid.o InputFile.o SpeciesClass.o Common.o BisectionMoveClass.o MoveClass.o ActionClass.o PathDataClass.o SpeciesArrayClass.o MirroredArrayClass.o CommunicatorClass.o PathClass.o test.o

Test: 	$(TestObjs)
	pushd ..; make; pushd
	$(LD) -o $@ $(TestObjs) $(LIBS)

.cc.o: $(HEADERS)
	$(CC) $(CCFLAGS) $(DEFS) $(INCL) $<
.f.o:
	g77 -c $<


SOURCES = ObservableClass.cc CubicSpline.cc Grid.cc InputFile.cc myprog.cc SpeciesClass.cc Common.cc BisectionMoveClass.cc MoveClass.cc ActionClass.cc PathDataClass.cc SpeciesArrayClass.cc MirroredArrayClass.cc CommunicatorClass.cc PathClass.cc test.cc


newmake: 
	make -f template.make Makefile FRC=force_rebuild

Makefile:	$(FRC)
	rm -f $@
	cp template.make $@
	echo 'Automatically generated dependency list:' >> $@
	$(CC) $(CCFLAGS) $(INCL) -M $(SOURCES) >> $@
	chmod -w $@


force_rebuild:



