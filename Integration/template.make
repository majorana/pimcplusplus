SOURCES = Integrate.cc GKIntegration.cc HermiteQuad.cc TestHermite.cc

all:	Integrate.o GKIntegration.o HermiteQuad.o TestHermite.o TestHermite

TestHermite: TestHermite.o HermiteQuad.o
	$(LD) -o TestHermite  HermiteQuad.o TestHermite.o $(IOobjs) $(LIBS)

clean:
	rm -f *.o


.cc.o: $(HEADERS)
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
