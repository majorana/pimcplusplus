SOURCES = InputOutput.cc InputOutputHDF5.cc InputFile.cc  InputOutputASCII.cc TestXML.cc

objs = InputOutput.o  InputOutputHDF5.o InputFile.o   InputOutputASCII.o InputOutputXML.o
all:	$(objs) TestXML



TestXML:	$(objs) TestXML.o
	$(LD) -o $@ $(objs) TestXML.o $(LIBS)

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
