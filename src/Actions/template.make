SOURCES = ShortRangeClass.cc  LongRangeRPAClass.cc LongRangeClass.cc ActionsClass.cc ActionBase.cc ShortRangePotClass.cc LongRangePotClass.cc KineticClass.cc

objs = ShortRangeClass.o LongRangeClass.o LongRangeRPAClass.o ActionsClass.o ActionBase.o ShortRangePotClass.o LongRangePotClass.o KineticClass.o


all: Actions	

Actions: ${objs}  



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
