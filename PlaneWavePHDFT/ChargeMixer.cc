#include "ChargeMixer.h"

KerkerMixerClass::KerkerMixerClass (FFTBox &fft) :
  ChargeMixerBase(fft), Lambda(2.0), HaveLastCharge(false)
{
  LastCharge.resize(FFT.GVecs.DeltaSize());
  NewCharge.resize(FFT.GVecs.DeltaSize());
}

void
KerkerMixerClass::SetLambda(double lambda)
{  Lambda = lambda; }

void
KerkerMixerClass::SetNewFraction(double newFrac)
{ NewFraction = newFrac; }

void
KerkerMixerClass::Reset()
{
  HaveLastCharge = false;
}

void
KerkerMixerClass::Mix (const Array<double,3> &newCharge, 
		       Array<double,3>       &mixedCharge_r,
		       zVec                  &mixedCharge_G)
{
  copy (newCharge, FFT.rBox);
  FFT.r2k();
  int nx, ny,nz;
  FFT.GVecs.GetFFTBoxSize(nx,ny,nz);
  double meshCellVol = FFT.GVecs.GetBoxVol()/(double)(nx*ny*nz);
  FFT.GetkVec(NewCharge);
  NewCharge *= meshCellVol;
  if (HaveLastCharge) 
    for (int i=0; i<NewCharge.size(); i++) {
      double G2 = dot(FFT.GVecs.DeltaG(i), FFT.GVecs.DeltaG(i));
      double fracNew = NewFraction * G2/(G2 + Lambda*Lambda);
      mixedCharge_G(i) = fracNew * NewCharge(i) + (1.0-fracNew)*LastCharge(i);
    }
  else
    mixedCharge_G = NewCharge;
  LastCharge = NewCharge;
  HaveLastCharge = true;
  // Now FFT mixed charge to real space
  FFT.PutkVec (mixedCharge_G);
  FFT.k2r();
  mixedCharge_r = real(FFT.rBox);
}

BroydenMixerClass::BroydenMixerClass (FFTBox &fft) :
  ChargeMixerBase(fft)
{
  // nothing more for now
}
