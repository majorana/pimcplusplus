SOURCES = ShiftTest.cc

myobjs = ShiftTest.o 
objs = $(myobjs) ../MirroredArrayClass.o ../CommunicatorClass.o
all:	$(objs) ShiftTest

ShiftTest:	$(myobjs) 
	$(LD) -o $@ $(objs)  $(LIBS)

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
