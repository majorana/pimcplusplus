SOURCES = EwaldBase.cc SimpleEwald.cc NaClTest.cc OptimizedBreakup.cc TestOptimizedBreakup.cc

all:	EwaldBase.o SimpleEwald.o OptimizedBreakup.o TestOptimizedBreakup.o \
        TestOptimizedBreakup #NaClTest

TestOptimizedBreakup: OptimizedBreakup.o TestOptimizedBreakup.o
	$(LD) -o TestOptimizedBreakup TestOptimizedBreakup.o OptimizedBreakup.o $(LIBS)

NaClTest: NaClTest.o EwaldBase.o SimpleEwald.o 
	$(LD) -o NaClTest NaClTest.o EwaldBase.o SimpleEwald.o $(IOobjs) $(LIBS)
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
