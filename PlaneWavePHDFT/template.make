SOURCES = PlaneWavePHDFT.cc Hamiltonian.cc ConjGrad.cc TestPW.cc

Objs = Hamiltonian.o ConjGrad.o

TestObjs  = TestPW.o

PHObjs = ../PH/Potential.o ../PH/ScreenedPot.o ../PH/QuinticPH.o              \
         ../PH/CoulombPot.o ../PH/SplinePot.o ../PH/HeAzizPot.o               \
         ../PH/kSpacePH.o

IOObjs = ../IO/InputOutput.o ../IO/InputOutputASCII.o ../IO/InputOutputHDF5.o \
 ../IO/InputOutputXML.o

SplineObjs = ../Splines/BicubicSpline.o ../Splines/CubicSpline.o              \
             ../Splines/Grid.o ../Splines/MyTricubicSpline.o                  \
             ../Splines/QuinticSpline.o ../Splines/QuinticSplines.o 

MiscObjs =   ../Integration/GKIntegration.o ../Fitting/Fitting.o              \
             ../MatrixOps/MatrixOps.o
all:    $(Objs) TestPW

TestPW:  $(TestObjs) $(Objs) $(PHObjs) $(IOObjs) $(SplineObjs) $(MiscObjs)
	$(LD) -o $@ $(TestObjs) $(Objs) $(PHObjs) $(IOObjs) $(SplineObjs) \
        $(MiscObjs) $(LIBS)

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
