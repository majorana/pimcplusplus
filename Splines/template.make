SOURCES = CubicSpline.cc Grid.cc BicubicSpline.cc TestBicubic.cc TestGrid.cc

IOobjs = ../IO/InputOutput.o ../IO/InputOutputHDF5.o ../IO/InputOutputASCII.o ../IO/InputFile.o ../IO/InputOutputXML.o

all:	TestBicubic TestGrid
TestBicubic:	CubicSpline.o Grid.o BicubicSpline.o TestBicubic.o 
	$(LD) -o TestBiCubic CubicSpline.o Grid.o BicubicSpline.o TestBicubic.o $(IOobjs) $(LIBS)

TestGrid:	 Grid.o  TestGrid.o 
	$(LD) -o TestGrid Grid.o TestGrid.o $(IOobjs) $(LIBS)

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
