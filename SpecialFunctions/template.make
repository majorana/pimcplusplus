SOURCES = SpecialFunctions.cc HermitePoly.cc TestLegendre.cc

all:	SpecialFunctions.o HermitePoly.o TestLegendre.o
	$(LD) -o TestHermite HermitePoly.o $(LIBS)
	$(LD) -o TestLegendre TestLegendre.o SpecialFunctions.o $(LIBS)

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
