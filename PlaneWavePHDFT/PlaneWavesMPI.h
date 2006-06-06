#ifndef PLANE_WAVES_MPI_H
#define PLANE_WAVES_MPI_H

#include "ConjGradMPI.h"

class CommunicatorClass;

class MPISystemClass
{
protected:
  Array<Vec3,1> Rions;
  FFTBox FFT;
  HamiltonianClass H;
  Array<complex<double>,2> Bands, LastBands;
  ConjGradMPI CG;
  Vec3 Box;
  double kCut;
  int NumBands;
  Potential *PH;
  bool UseFFT, UseLDA;
  CommunicatorClass &BandComm, &kComm;
  bool MDExtrap, FirstTime;
public:
  GVecsClass GVecs;
  void Setup(Vec3 box, Vec3 k, double kcut, Potential &ph, bool useFFT=true);
  void Setup(Vec3 box, Vec3 k, double kcut, double z, bool useFFT=true);
  void SetIons (const Array<Vec3,1> &rions);
  inline Vec3 GetIonPos(int i) { return Rions(i); }
  void Setk (Vec3 k);
  void DiagonalizeH();
  inline double GetEnergy(int band) { return CG.Energies(band); }
  inline int GetNumBands() { return NumBands; }

  /// Gets the FFT box dimensions.
  inline void GetBoxDims(int &nx, int &ny, int &nz)
  { FFT.GetDims(nx, ny, nz); }
  /// This FFT's the desired band into real space.
  void SetRealSpaceBandNum(int num);
  inline complex<double> RealSpaceBand (int ix, int iy, int iz)
  { return FFT.rBox(ix, iy, iz); }
  

  void CalcChargeDensity(Array<double,3> &rho);
  void WriteXSFFile(string filename);

  MPISystemClass(int numBands, int numElecs,
		 CommunicatorClass &bandcomm, CommunicatorClass &kcomm,
		 bool useLDA=false, bool mdextrap=false) 
    : CG(H, numElecs, Bands, bandcomm, kcomm, FFT), FFT(GVecs), H(GVecs, FFT), 
      NumBands(numBands), BandComm(bandcomm), kComm(kcomm), 
      MDExtrap(mdextrap), FirstTime(true), UseLDA(useLDA)
  {

  }
};

#endif
