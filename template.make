PASS_DEFS = "CC=${CC}" "LD=${LD}" "CCFLAGS=${CCFLAGS}" \
"DEFS=${DEFS}" "INCL=${INCL}" "LIBS=${LIBS}"

MAKE_ALL = ${MAKE} all $(PASS_DEFS)
MAKE_NEWMAKE = ${MAKE} -f template.make newmake $(PASS_DEFS)
MAKE_CLEAN = ${MAKE} clean

all:	PH_obj Splines_obj Integration_obj IO_obj DFT_obj Random_obj \
	MPI_obj Optimize_obj SpecialFunctions_obj MatrixOps_obj \
	Fitting_obj


PH_obj:
	cd PH; $(MAKE_ALL) 
Splines_obj:
	cd Splines; $(MAKE_ALL)

Integration_obj:
	cd Integration; $(MAKE_ALL)

IO_obj:
	cd IO; $(MAKE_ALL)

DFT_obj:
	cd DFT; $(MAKE_ALL)

Random_obj:
	cd Random; $(MAKE_ALL)

MPI_obj:
	cd MPI; $(MAKE_ALL)

Optimize_obj:
	cd Optimize; $(MAKE_ALL)

SpecialFunctions_obj:
	cd SpecialFunctions; $(MAKE_ALL)

MatrixOps_obj:
	cd MatrixOps; $(MAKE_ALL)

Fitting_obj:
	cd Fitting; $(MAKE_ALL)


CLEANS = PH_clean Splines_clean Integration_clean IO_clean DFT_clean \
	 Random_clean MPI_clean Optimize_clean SpecialFunctions_clean \
	 MatrixOps_clean Fitting_clean


clean:	$(CLEANS)
	rm *.o


PH_clean:
	cd PH; $(MAKE_CLEAN)

Splines_clean:
	cd Splines; $(MAKE_CLEAN)

Integration_clean:
	cd Integration; $(MAKE_CLEAN)

IO_clean:
	cd IO; $(MAKE_CLEAN)

DFT_clean:
	cd DFT; $(MAKE_CLEAN)

Random_clean:
	cd Random; $(MAKE_CLEAN)

MPI_clean:
	cd MPI; $(MAKE_CLEAN)

Optimize_clean:
	cd Optimize; $(MAKE_CLEAN)

SpecialFunctions_clean:
	cd SpecialFunctions; $(MAKE_CLEAN)

MatrixOps_clean:
	cd MatrixOps; $(MAKE_CLEAN)

Fitting_clean:
	cd Fitting; $(MAKE_CLEAN)




PH_newmake:
	cd PH; $(MAKE_NEWMAKE)

Splines_newmake:
	cd Splines; $(MAKE_NEWMAKE)

Integration_newmake:
	cd Integration; $(MAKE_NEWMAKE)

IO_newmake:
	cd IO; $(MAKE_NEWMAKE)

DFT_newmake:
	cd DFT; $(MAKE_NEWMAKE)

Random_newmake:
	cd Random; $(MAKE_NEWMAKE)

MPI_newmake:
	cd MPI; $(MAKE_NEWMAKE)

Optimize_newmake:
	cd Optimize; $(MAKE_NEWMAKE)

SpecialFunctions_newmake:
	cd SpecialFunctions; $(MAKE_NEWMAKE)

MatrixOps_newmake:
	cd MatrixOps; $(MAKE_NEWMAKE)

Fitting_newmake:
	cd Fitting; $(MAKE_NEWMAKE)



NEW_MAKES = PH_newmake Splines_newmake Integration_newmake IO_newmake \
DFT_newmake Random_newmake MPI_newmake Optimize_newmake \
SpecialFunctions_newmake MatrixOps_newmake Fitting_newmake

SOURCES = `*.cc`

ifeq ($(strip $(SOURCES)),)
	TO_MAKE = $(CC) $(CCFLAGS) $(INCL) -M $(SOURCES) >> $@
else
	TO_MAKE = 
endif


newmake: $(NEW_MAKES)
	$(MAKE) -f template.make Makefile FRC=force_rebuild

Makefile:	$(FRC)
	rm -f $@
	cp template.make $@
	echo 'Automatically generated dependency list:' >> $@
	$(TO_MAKE)
	chmod -w $@


force_rebuild:
