SOURCES = CubicSpline.cc Grid.cc BicubicSpline.cc TestBicubic.cc TestGrid.cc TestTricubic.cc MyTricubicSpline.cc TestMyTricubic.cc QuinticSpline.cc TestQuintic.cc DyutimanTest.cc

IOobjs = ../IO/InputOutput.o ../IO/InputOutputHDF5.o ../IO/InputOutputASCII.o  ../IO/InputOutputXML.o

F77Objs = fortran/evtricub.o  fortran/herm3ev.o  fortran/mktricubw.o  fortran/tcspline.o fortran/ibc_ck.o fortran/splinck.o fortran/zonfind.o fortran/tcspeval.o fortran/v_spline.o fortran/bcspline.o fortran/bcspeval.o

Objs:	MyTricubicSpline.o CubicSpline.o QuinticSpline.o DyutimanTest.o

all:	FortranObjs TestBicubic TestGrid  TestMyTricubic Objs TestQuintic DyutimanTest

TestQuintic:	QuinticSpline.o QuinticSplines.o TestQuintic.o
	$(LD) -o TestQuintic QuinticSpline.o QuinticSplines.o TestQuintic.o $(LIBS)

TestBicubic:	CubicSpline.o Grid.o BicubicSpline.o TestBicubic.o 
	$(LD) -o TestBicubic CubicSpline.o Grid.o BicubicSpline.o TestBicubic.o $(IOobjs) $(LIBS)

TestMyTricubic:	CubicSpline.o Grid.o  TestMyTricubic.o MyTricubicSpline.o
	$(LD) -o TestMyTricubic CubicSpline.o Grid.o TestMyTricubic.o MyTricubicSpline.o $(IOobjs) -Lfortran -lpspline $(LIBS)


DyutimanTest:	CubicSpline.o Grid.o  TestMyTricubic.o MyTricubicSpline.o
	$(LD) -o DyutimanTest CubicSpline.o Grid.o DyutimanTest.o MyTricubicSpline.o $(IOobjs) -Lfortran -lpspline $(LIBS)


TestTricubic:	Grid.o TestTricubic.o FortranObjs TestMyTricubic
	$(LD) -o TestTricubic Grid.o TestTricubic.o $(F77Objs) $(IOobjs) $(LIBS)

TestGrid:	 Grid.o  TestGrid.o 
	$(LD) -o TestGrid Grid.o TestGrid.o $(IOobjs) $(LIBS)

FortranObjs:
	cd fortran; $(MAKE)


clean:
	rm -f *.o

.cc.o: 
	$(CC) $(CCFLAGS) $(DEFS) $(INCL) -o $*.o $< 
.f.o:
	$(F77) -c -O3 -o $*.o $<

newmake:
	$(MAKE) -f template.make Makefile FRC=force_rebuild

Makefile:	$(FRC)
	rm -f $@
	cp template.make $@
	echo 'Automatically generated dependency list:' >> $@
	$(MAKECC) $(CCFLAGS) $(INCL) -MM $(SOURCES) >> $@
	chmod -w $@


force_rebuild:
