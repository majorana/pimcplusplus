all:	evtricub.o  herm3ev.o  mktricubw.o  tcspline.o ibc_ck.o splinck.o\
        zonfind.o tcspeval.o v_spline.o bcspline.o bcspeval.o

clean:
	rm -f *.o

.f.o:
	g77 -g -O3 -c -o $*.o $<
