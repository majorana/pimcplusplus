SOURCES = BisectionClass.cc BlockMove.cc OpenBisectionMoveClass.cc  MetaMoves.cc PermuteTableClass.cc BisectionMoveClass.cc MoveBase.cc RandomPermClass.cc  MultiStage.cc NoPermuteClass.cc BisectionBlock.cc PermuteStageClass.cc TablePermuteStageClass.cc BisectionStageClass.cc RefSliceMove.cc

objs = BisectionClass.o BlockMove.o OpenBisectionMoveClass.o MetaMoves.o PermuteTableClass.o BisectionMoveClass.o  MoveBase.o RandomPermClass.o  MultiStage.o NoPermuteClass.o BisectionBlock.o PermuteStageClass.o TablePermuteStageClass.o BisectionStageClass.o RefSliceMove.o


all:	Moves	

Moves:	${objs}  



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
