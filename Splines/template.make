SOURCES = CubicSpline.cc Grid.cc BicubicSpline.cc TestBicubic.cc TestGrid.cc TestTricubic.cc

IOobjs = ../IO/InputOutput.o ../IO/InputOutputHDF5.o ../IO/InputOutputASCII.o ../IO/InputFile.o ../IO/InputOutputXML.o

F77Objs = fortran/evtricub.o  fortran/herm3ev.o  fortran/mktricubw.o  fortran/tcspline.o fortran/ibc_ck.o fortran/splinck.o fortran/zonfind.o fortran/tcspeval.o fortran/v_spline.o fortran/bcspline.o fortran/bcspeval.o

all:	FortranObjs TestBicubic TestGrid  TestTricubic 

TestBicubic:	CubicSpline.o Grid.o BicubicSpline.o TestBicubic.o 
	$(LD) -o TestBiCubic CubicSpline.o Grid.o BicubicSpline.o TestBicubic.o $(IOobjs) $(LIBS)

TestTricubic:	Grid.o TestTricubic.o FortranObjs
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
	g77 -c -o $*.o $<

newmake:
	$(MAKE) -f template.make Makefile FRC=force_rebuild

Makefile:	$(FRC)
	rm -f $@
	cp template.make $@
	echo 'Automatically generated dependency list:' >> $@
	$(CC) $(CCFLAGS) $(INCL) -M $(SOURCES) >> $@
	chmod -w $@


force_rebuild:
