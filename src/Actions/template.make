SOURCES = ShortRangeClass.cc  LongRangeRPAClass.cc LongRangeClass.cc ActionsClass.cc ActionBase.cc ShortRangePotClass.cc LongRangePotClass.cc KineticClass.cc NodalActionClass.cc DavidLongRangeClass.cc ShortRangeApproximateClass.cc OpenLoopImportance.cc TIP5PWaterClass.cc StructureReject.cc GroundStateNodalActionClass.cc ST2WaterClass.cc


objs = ShortRangeClass.o LongRangeClass.o LongRangeRPAClass.o ActionsClass.o ActionBase.o ShortRangePotClass.o LongRangePotClass.o KineticClass.o NodalActionClass.o DavidLongRangeClass.o ShortRangeApproximateClass.o OpenLoopImportance.o TIP5PWaterClass.o ST2WaterClass.o StructureReject.o GroundStateNodalActionClass.o

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
	$(MAKECC) $(CCFLAGS) $(INCL) -M $(SOURCES) >> $@
	chmod -w $@


force_rebuild:
