include /home/common/Codes/Make.include

LIBS = $(BLITZLIB) $(SPRNGLIB) $(GSLLIB) $(G2CLIB) $(LAPACKLIB) $(G2CLIB) -lm #-lstdc++
INCL = $(BLITZINC) $(SPRNGINC) $(GSLINC) 

CCFLAGS = -c -g -DBZ_DEBUG #-pg
CC = mpiCC
LD = mpiCC  -Bstatic 
DEFS = -DNO_COUT -DUSE_MPI -DBZ_DEBUG -g #-DUSE_MPI 

TestObjs = CubicSpline.o Grid.o InputFile.o myprog.o IdenticleParticleClass.o Common.o BisectionMoveClass.o MoveClass.o ActionClass.o PathDataClass.o ArrayOfIdenticalParticlesClass.o MirroredArrayClass.o

Test: 	$(TestObjs)
	pushd ..; make; pushd
	$(LD) -o $@ $(TestObjs) $(LIBS)

.cc.o: $(HEADERS)
	$(CC) $(CCFLAGS) $(DEFS) $(INCL) $<
.f.o:
	g77 -c $<


SOURCES = CubicSpline.cc Grid.cc InputFile.cc myprog.cc IdenticleParticleClass.cc Common.cc BisectionMoveClass.cc MoveClass.cc ActionClass.cc PathDataClass.cc ArrayOfIdenticalParticlesClass.cc MirroredArrayClass.cc


newmake: 
	make -f template.make Makefile FRC=force_rebuild

Makefile:	$(FRC)
	rm -f $@
	cp template.make $@
	echo 'Automatically generated dependency list:' >> $@
	$(CC) $(CFLAGS) $(INCL) -M $(SOURCES) >> $@
	chmod -w $@


force_rebuild:



