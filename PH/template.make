SOURCES = PH.cc Chebyshev.cc Potential.cc Atom.cc Cost.cc CoreTransform.cc \
	  ScreenedPot.cc QuinticPH.cc CoulombPot.cc SplinePot.cc HeAzizPot.cc

all:	PH.o Chebyshev.o Potential.o Atom.o Cost.o CoreTransform.o        \
	ScreenedPot.o QuinticPH.o CoulombPot.o SplinePot.o HeAzizPot.o 

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
	$(MAKECC) $(CCFLAGS) $(INCL) -MM $(SOURCES) >> $@
	chmod -w $@


force_rebuild:
