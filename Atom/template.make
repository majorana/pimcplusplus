SOURCES = Atom.cc RadialWF.cc TestRadialWF.cc DFTAtom.cc TestDFTAtom.cc

all:	Atom.o RadialWF.o TestRadialWF DFTAtom.o TestDFTAtom

TestWFobjs = TestRadialWF.o RadialWF.o ../PH/CoulombPot.o ../Splines/CubicSpline.o 

TestDFTAtomobjs = RadialWF.o ../PH/CoulombPot.o ../Splines/CubicSpline.o TestDFTAtom.o DFTAtom.o ../PH/ScreenedPot.o ../DFT/Functionals.o ../DFT/ExCorr.o

TestRadialWF:  $(TestWFobjs)
	$(LD) -o TestRadialWF $(TestWFobjs) $(LIBS)

TestDFTAtom:  $(TestDFTAtomobjs)
	$(LD) -o TestDFTAtom $(TestDFTAtomobjs) $(LIBS)

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
