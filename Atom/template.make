SOURCES = Atom.cc RadialWF.cc TestRadialWF.cc NewAtom.cc

all:	Atom.o RadialWF.o TestRadialWF NewAtom.o

TestWFobjs = TestRadialWF.o RadialWF.o ../PH/CoulombPot.o ../Splines/CubicSpline.o 

TestRadialWF:  $(TestWFobjs)
	$(LD) -o TestRadialWF $(TestWFobjs) $(LIBS)

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
