SOURCES = PathVis.cc GLObject.cc PathObject.cc trackball.c ViewClass.cc BoxObject.cc Visual.cc SphereObject.cc BoxClass.cc SmoothClass.cc

GLINC = `pkg-config gtkglextmm-1.2 --cflags`
GLLIBS = `pkg-config gtkglextmm-1.2 --libs` -lgle -lglut

IOobjs = ../Common/IO/InputOutput.o ../Common/IO/InputOutputHDF5.o \
         ../Common/IO/InputOutputASCII.o  ../Common/IO/InputOutputXML.o

Objs =	PathVis.o GLObject.o PathObject.o trackball.o ViewClass.o BoxObject.o Visual.o SphereObject.o BoxClass.o SmoothClass.o

all:	PathVis

PathVis:  $(IOobjs) $(Objs)
	$(LD) -o PathVis $(Objs) $(IOobjs) $(LIBS) $(GLLIBS)

clean:
	rm -f *.o

.cc.o: 
	$(CC) $(CCFLAGS) $(DEFS) $(INCL) $(GLINC) -o $*.o $< 
.c.o: 
	gcc $(CCFLAGS) $(DEFS) $(INCL) $(GLINC) -o $*.o $< 
.f.o:
	$(F77) -c -O3 -o $*.o $<

newmake:
	$(MAKE) -f template.make Makefile FRC=force_rebuild

Makefile:	$(FRC)
	rm -f $@
	cp template.make $@
	echo 'Automatically generated dependency list:' >> $@
	$(MAKECC) $(CCFLAGS) $(INCL) -MM $(SOURCES) >> $@
	chmod -w $@


force_rebuild:
