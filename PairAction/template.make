SOURCES = PAszFit.cc PAcoulombFit.cc PAcoulombBCFit.cc PAclassicalFit.cc PAsFit.cc PAtricubicFit.cc PAzeroFit.cc DavidPAClass.cc PACoulombFit.cc PADipoleFit.cc PATripoleFit.cc

#all:	PAszFit.o PAcoulombFit.o PAcoulombBCFit.o PAclassicalFit.o PAsFit.o PAtricubicFit.o PAzeroFit.o DavidPAClass.o
all:	 PAcoulombBCFit.o PAtricubicFit.o PAclassicalFit.o PAzeroFit.o DavidPAClass.o PACoulombFit.o PADipoleFit.o PATripoleFit.o

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
