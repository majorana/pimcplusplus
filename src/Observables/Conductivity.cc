#include "Conductivity.h"

void
ConductivityClass::CalcCurrentT()
{
  PathClass &Path = PathData.Path;
  MyCurrent = 0.0;
  for (int si=0; si<Path.NumSpecies(); si++) {
    SpeciesClass &species = Path.Species(si);
    if (species.lambda != 0.0) {
      int first = species.FirstPtcl;
      int last = species.LastPtcl;
      for (int slice=0; slice<(Path.NumTimeSlices()-1); slice++) 
	for (int ptcl=first; ptcl<=last; ptcl++)
	  MyCurrent(slice) += Path.Velocity(slice, slice+1, ptcl);
    }
  }
  // Now that we've calculated this processors current, let us gather
  // to processor 0
  Path.Communicator.Gather(MyCurrent, TempCurrent, ProcNumLinks);
}

void
ConductivityClass::AccumulateSlow()
{
  CalcCurrentT();
  for (int sep=0; sep<M; sep++) {
    dVec TT = 0.0;
    for (int slice1=0; slice1<M; slice1++) {
      int slice2 = (slice1+sep);
      slice2 -= (slice2>=M) ? M : 0;
      for (int dim=0; dim<NDIM; dim++)
	TT[dim] += TempCurrent(slice1)[dim]*TempCurrent(slice2)[dim];
    }
    CorrSumTT(sep) = TT;
  }
}

inline double mag2(complex<double> z)
{
  return z.real()*z.real() + z.imag()*z.imag();
}

void
ConductivityClass::AccumulateFast()
{
  CalcCurrentT();
  double Minv = 1.0/(double)M;
  for (int dim=0; dim<NDIM; dim++) {
    for (int i=0; i<M; i++)
      FFTtemp.rBox(i) = TempCurrent(i)[dim];
    FFTtemp.r2k();
    for (int i=0; i<M; i++)
      FFTtemp.kBox(i) = Minv*mag2(FFTtemp.kBox(i));
    FFTtemp.k2r();
    for (int i=0; i<M; i++)
      CorrSumTT(i)[dim] += FFTtemp.rBox(i).real();
  }
}

void
ConductivityClass::Accumulate()
{
//   AccumulateSlow();
//   Array<dVec,1> temp;
//   temp.resize(CorrSumTT.size());
//   temp = CorrSumTT;
//   CorrSumTT = 0.0;
  AccumulateFast();
//   fprintf (stderr, "CorrSumTT = \n");
//   for (int i=0; i<CorrSumTT.size(); i++)
//     fprintf (stderr, "%1.12e %1.12e\n", temp(i)[0], CorrSumTT(i)[0]);
  NumSamples++;
}


void
ConductivityClass::Read(IOSectionClass &in)
{
  ObservableClass::Read(in);
  M = PathData.Path.TotalNumSlices;
  Current_T.resize  (M);
  CorrSumTT.resize  (M);
  TempCurrent.resize(M);
  MyCurrent.resize(PathData.Path.NumTimeSlices()-1);
  WriteArray.resize(M, NDIM);
  Current_T = 0.0;
  CorrSumTT = 0.0;
  MyCurrent = 0.0;
  ProcNumLinks.resize(PathData.Path.Communicator.NumProcs());
  FFTtemp.resize(PathData.Path.TotalNumSlices);
  for (int proc=0; proc<PathData.Path.Communicator.NumProcs(); proc++) {
    int slice1, slice2;
    PathData.Path.SliceRange(proc, slice1, slice2);
    ProcNumLinks = (slice2-slice1);
  }
}


void 
ConductivityClass::WriteBlock()
{
  double norm = 1.0/(double)NumSamples;
  for (int i=0; i<M; i++)
    for (int j=0; j<NDIM; j++)
      WriteArray(i,j) = norm*CorrSumTT(i)[j];
  KineticVar.Write(WriteArray);

  CorrSumTT = 0.0;
  NumSamples = 0;
}
