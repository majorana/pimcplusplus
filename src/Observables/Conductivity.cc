/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

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
    RealSumTT(sep) += TT;
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
  if (PathData.Path.Communicator.MyProc()==0) {
    double Minv = 1.0/(double)M;
    for (int dim=0; dim<NDIM; dim++) {
      for (int i=0; i<M; i++)
	FFTtemp.rBox(i) = TempCurrent(i)[dim];
      FFTtemp.r2k();
      for (int i=0; i<M; i++)
	FFTtemp.kBox(i) = Minv*mag2(FFTtemp.kBox(i));
      for (int i=0; i<M; i++)
	FreqSumTT(i)[dim] += FFTtemp.rBox(i).real();
      FFTtemp.k2r();
      for (int i=0; i<M; i++)
	RealSumTT(i)[dim] += FFTtemp.rBox(i).real();
    }
  }
  NumSamples++;
}

void
ConductivityClass::Accumulate()
{
//   AccumulateSlow();
//   Array<dVec,1> temp;
//   temp.resize(RealSumTT.size());
//   temp = RealSumTT;
//   RealSumTT = 0.0;
  AccumulateFast();
//   fprintf (stderr, "RealSumTT = \n");
//   for (int i=0; i<RealSumTT.size(); i++)
//     fprintf (stderr, "%1.12e %1.12e\n", temp(i)[0], RealSumTT(i)[0]);
}


void
ConductivityClass::Read(IOSectionClass &in)
{
  ObservableClass::Read(in);
  M = PathData.Path.TotalNumSlices;
  Current_T.resize  (M);
  RealSumTT.resize  (M);
  FreqSumTT.resize  (M);
  TempCurrent.resize(M);
  MyCurrent.resize(PathData.Path.NumTimeSlices()-1);
  WriteArray.resize(M, NDIM);
  Current_T = 0.0;
  RealSumTT = 0.0;
  MyCurrent = 0.0;
  RealSumTT = 0.0;
  FreqSumTT = 0.0;
  ProcNumLinks.resize(PathData.Path.Communicator.NumProcs());
  FFTtemp.resize(PathData.Path.TotalNumSlices);
  for (int proc=0; proc<PathData.Path.Communicator.NumProcs(); proc++) {
    int slice1, slice2;
    PathData.Path.SliceRange(proc, slice1, slice2);
    ProcNumLinks(proc) = (slice2-slice1);
  }
}


void 
ConductivityClass::WriteBlock()
{
  double norm = 1.0/(double)NumSamples;
  for (int i=0; i<M; i++)
    for (int j=0; j<NDIM; j++)
      WriteArray(i,j) = norm*RealSumTT(i)[j];
  RealKineticVar.Write(WriteArray);
  for (int i=0; i<M; i++)
    for (int j=0; j<NDIM; j++)
      WriteArray(i,j) = norm*FreqSumTT(i)[j];
  FreqKineticVar.Write(WriteArray);

  RealSumTT = 0.0;
  FreqSumTT = 0.0;
  NumSamples = 0;
}
