SOURCES = ShiftTest.cc


permobjs =                            \
  PermutationTest.o 		      \
  ../ObservableClass.o                \
  ../Common/Splines/CubicSpline.o     \
  ../Common/Splines/Grid.o            \
  ../SpeciesClass.o                   \
  ../Common.o                         \
  ../BisectionMoveClass.o             \
  ../MoveClass.o                      \
  ../ActionClass.o                    \
  ../PathDataClass.o                  \
  ../CommunicatorClass.o              \
  ../PathClass.o                      \
  ../DistanceTablePBCClass.o          \
  ../DistanceTableFreeClass.o         \
  ../DistanceTableClass.o             \
  ../MirroredArrayClass.o             \
  ../Common/IO/InputOutput.o          \
  ../Common/IO/InputOutputHDF5.o      \
  ../Common/IO/InputFile.o            \
  ../Common/IO/InputOutputASCII.o     \
  ../Common/IO/InputOutputXML.o       \
  ../Common/PairAction/PAcoulombFit.o \
  ../Common/PairAction/PAszFit.o      \
  ../Common/PH/PH.o                   \
  ../Common/PH/Potential.o


myobjs = ShiftTest.o 
objs = $(myobjs) ../MirroredArrayClass.o ../CommunicatorClass.o
all:	$(objs) ShiftTest 

ShiftTest:	$(myobjs) 
	$(LD) -o $@ $(objs)  $(LIBS)

PermTest: $(permobjs) 
	$(LD) -o $@ $(permobjs) $(LIBS)
	
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
	$(MAKECC) $(CCFLAGS) $(INCL) $(DEFS) -M $(SOURCES) >> $@
	chmod -w $@


force_rebuild:
