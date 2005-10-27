fname = 'H2_2.0_2.6_beta25.0.h5';

wA   = hdf5read(fname, '/Moves_4/CorrelatedBisectionBlock_1/wA');
wB   = hdf5read(fname, '/Moves_4/CorrelatedBisectionBlock_1/wB');
wAEA = hdf5read(fname, '/Moves_4/CorrelatedBisectionBlock_1/wAEA');
wBEB = hdf5read(fname, '/Moves_4/CorrelatedBisectionBlock_1/wBEB');
wASA = hdf5read(fname, '/Moves_4/CorrelatedBisectionBlock_1/wASA');
wBSB = hdf5read(fname, '/Moves_4/CorrelatedBisectionBlock_1/wBSB');

l = length(wAEA);
[EA, EAerr] = stats(wAEA);
EA    = l*EA/sum(wA);
EAerr = l*EAerr/sum(wA);
[EB, EBerr] = stats(wBEB);
EB    = l*EB/sum(wB);
EBerr = l*EBerr/sum(wB);
[SA, SAerr] = stats(wASA);
SA    = l*SA/sum(wA);
SAerr = l*SAerr/sum(wA);
[SB, SBerr] = stats(wBSB);
SB    = l*SB/sum(wB);
SBerr = l*SBerr/sum(wB);

deltaE = EA-EB;
[dE, deltaEerr] = stats(wAEA./wA-wBEB./wB);
[dS, deltaSerr] = stats(wASA./wA-wBSB./wB);
%deltaEerr = l*deltaEerr/sum(wA+wB);

figure(1);
plot (wAEA./wA-wBEB./wB);
figure(2)
plot (wASA./wA-wBSB./wB);
sprintf ('EA = %1.5f +/- %1.5f', EA, EAerr)
sprintf ('EB = %1.5f +/- %1.5f', EB, EBerr)
%sprintf ('SA = %1.5f +/- %1.5f', SA, SAerr)
%sprintf ('SB = %1.5f +/- %1.5f', SB, SBerr)
sprintf ('delta E = %1.5f +/- %1.5f', EA-EB, deltaEerr)
%sprintf ('delta S = %1.5f +/- %1.5f', SA-SB, deltaSerr)
sprintf ('log(wA/wB)/25 = %1.5f', log(sum(wA)/sum(wB))/25)
[wAavg, wAerr] = stats(wA);
[wBavg, wBerr] = stats(wB);
sprintf ('wA = %1.5f +/- %1.5f', wAavg/(wAavg+wBavg), wAerr/(wAavg+wBavg))

%Eerror = zeros(1,500);
%for i=[1:500]
%  blockE = block(wA, wAEA, wB, wBEB,i);
%  [dE, Eerror(i)] = stats(blockE);
%end;
%plot (Eerror)
