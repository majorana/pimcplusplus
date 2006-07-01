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

#include "PairFixedPhase.h"
#include "../PathDataClass.h"
#include <Common/MatrixOps/MatrixOps.h>

double 
PairFixedPhaseClass::d_dBeta(int slice1, int slice2, int level)
{
  //do nothing for now

}

/// Return essentiall 0 or infinity
double 
PairFixedPhaseClass::SingleAction (int startSlice, int endSlice,
				   const Array<int,1> &activeParticles, 
				   int level)
{
  int skip = 1<<level;
  int myStart, myEnd;
  Path.SliceRange(PathData.Path.Communicator.MyProc(), myStart, myEnd);
  int refSlice = Path.GetRefSlice() - myStart;

  double uNode=0.0;
  for (int slice=startSlice; slice <= endSlice; slice+=skip) {
    if ((slice != refSlice) && (slice != refSlice+Path.TotalNumSlices)){
      //      complex<double> detValue=Det(slice);
      complex<double> detValue=CalcDegenerateDet(slice);
      cerr<<detValue<<endl;
      if (detValue.real() < 0.0)
	uNode += 1.0e50;
    }
  }
  return uNode;
}


PairFixedPhaseClass::PairFixedPhaseClass (PathDataClass &pathData) :
  //  NodalActionClass (pathData), 
  ActionBaseClass (pathData), 
  PathData(pathData),
  Path (pathData.Path)
{
  //do nothing for now
}


complex<double>
PairFixedPhaseClass::Det (int slice, dVec kVec)
{


  ///Figure out the difference between the slice and the refslice

  int myStartSlice, myEndSlice;
  int myProc = PathData.Path.Communicator.MyProc();
  Path.SliceRange (myProc, myStartSlice, myEndSlice);
  int refSlice = Path.GetRefSlice()-myStartSlice;
  int sliceDiff = abs(slice-refSlice);
  sliceDiff = min (sliceDiff, Path.TotalNumSlices-sliceDiff);
  assert (sliceDiff <= Path.TotalNumSlices);
  assert (sliceDiff > 0);


  SpeciesClass &species1 = Path.Species(Species1Num);
  SpeciesClass &species2 = Path.Species(Species2Num);

  int firstS1 = species1.FirstPtcl;
  int lastS1 = species1.LastPtcl;

  int firstS2 = species2.FirstPtcl;
  int lastS2 = species2.LastPtcl;
  assert ((lastS1-firstS1)-(lastS2-firstS2)-2==0);

  int minParticleNum=lastS2-firstS2+1;
  // Fill up determinant matrix
  //We start by filling up the determinant matrix in the area where
  //there is an equal number of up and down electrons and so there is
  //an overlap of the spin.
  for (int detIndexA=0;detIndexA<minParticleNum+2;detIndexA++){
    for (int detIndexB=0;detIndexB<minParticleNum;detIndexB++){
      int ptclSpecies1=detIndexA;
      int ptclSpecies2=detIndexB+firstS2;
      dVec disp;
      double dist;
      Path.DistDisp(slice,ptclSpecies1,ptclSpecies2,dist,disp);
      double pairingDist2=dot(disp,disp);
      double action = 10.0*pairingDist2; //ActionSplines(sliceDiff)[dim](diff[dim]);
      DetMatrix(detIndexA, detIndexB) = exp(-action);
    }
  }


  for (int detIndexA=0;detIndexA<minParticleNum+2;detIndexA++){
    int detIndexB=minParticleNum;
    DetMatrix(detIndexA,detIndexB)=1.0;
  }


  for (int detIndexA=0;detIndexA<minParticleNum+2;detIndexA++){
    int detIndexB=minParticleNum+1;
    int ptclSpecies1=detIndexA;
    dVec s=PathData.Path(slice,ptclSpecies1);
    PathData.Path.PutInBox(s);
    double kdots=dot(kVec,s);
    complex<double> action(cos(kdots),-sin(kdots));
    DetMatrix(detIndexA,detIndexB)=action;
  }
  
  return Determinant (DetMatrix);
}



complex<double>
PairFixedPhaseClass::CalcDegenerateDet(int slice)
{

//   double k_x=2*M_PI/PathData.Path.GetBox()[0];
//   double k_y=2*M_PI/PathData.Path.GetBox()[1];
//   double k_z=2*M_PI/PathData.Path.GetBox()[2];
//   double totalDet=1.0;

//   dVec kVec;
//   kVec(0)=k_x;
//   kVec(1)=0;
//   kVec(2)=0;
//   totalDet*=DegenerateRefSliceDeterminates(0)*Det(slice,kVec);
//   kVec(0)=-k_x;
//   totalDet*=DegenerateRefSliceDeterminates(1)*Det(slice,kVec);

//   kVec(0)=0;
//   kVec(1)=k_y;
//   totalDet*=DegenerateRefSliceDeterminates(2)*Det(slice,kVec);
//   kVec(1)=-k_y;
//   totalDet*=DegenerateRefSliceDeterminates(3)*Det(slice,kVec);


//   kVec(1)=0;
//   kVec(2)=k_z;
//   totalDet*=DegenerateRefSliceDeterminates(4)*Det(slice,kVec);
//   kVec(2)=-k_z;
//   totalDet*=DegenerateRefSliceDeterminates(5)*Det(slice,kVec);
  
  
//   //  int sliceDiff = abs(slice-refSlice);
//   //  sliceDiff = min (sliceDiff, Path.TotalNumSlices-sliceDiff);
//   //  assert (sliceDiff <= Path.TotalNumSlices);
//   //  assert (sliceDiff > 0);
  
//   return totalDet;

}




void 
PairFixedPhaseClass::BuildRefDet ()
{

  double k_x=2*M_PI/PathData.Path.GetBox()[0];
  double k_y=2*M_PI/PathData.Path.GetBox()[1];
  double k_z=2*M_PI/PathData.Path.GetBox()[2];
  
  int refSliceOwner=PathData.Path.SliceOwner(Path.GetRefSlice());
  int myProc = PathData.Path.Communicator.MyProc();
  if (myProc!=refSliceOwner)
    return;
  int myStartSlice, myEndSlice;
  Path.SliceRange (myProc, myStartSlice, myEndSlice);
  int refSliceLocal = Path.GetRefSlice()-myStartSlice;
  dVec kVec;
  kVec(0)=k_x;
  kVec(1)=0;
  kVec(2)=0;
  DegenerateRefSliceDeterminates(0)=Det(refSliceLocal,kVec);
  kVec(0)=-k_x;
  DegenerateRefSliceDeterminates(1)=Det(refSliceLocal,kVec);

  kVec(0)=0;
  kVec(1)=k_y;
  DegenerateRefSliceDeterminates(2)=Det(refSliceLocal,kVec);
  kVec(1)=-k_y;
  DegenerateRefSliceDeterminates(3)=Det(refSliceLocal,kVec);


  kVec(1)=0;
  kVec(2)=k_z;
  DegenerateRefSliceDeterminates(4)=Det(refSliceLocal,kVec);
  kVec(2)=-k_z;
  DegenerateRefSliceDeterminates(5)=Det(refSliceLocal,kVec);
  
  
  //  int sliceDiff = abs(slice-refSlice);
  //  sliceDiff = min (sliceDiff, Path.TotalNumSlices-sliceDiff);
  //  assert (sliceDiff <= Path.TotalNumSlices);
  //  assert (sliceDiff > 0);
  


}


// void 
// BuildRefDet()
// {
//   for (int ptcl1=0;ptcl1<DistanceMatrix.extent(0);ptcl1++){
//     for (int ptcl2=0;ptcl2<DistanceMatrix.extent(1);ptcl2++){
//       TrialRefDet(ptcl1,ptcl2)=exp(-dot(DistanceMatrix(ptc1,ptcl2),DistanceMatrix(ptcl1,ptcl2))*sigma2Inverse);
//     }
//   }
  

//   for (int ptcl1=0;ptcl1<DistanceMatrix.extent(0);ptcl1++){
//     TrialRefDet(ptcl1,DistanceMatrix.extent(1)+1);
//   }
//   //compute k \dot r
//   ///exp(ik \dot r)
  
  




// }


///Lapack routine: cgdi zgedi

//The second species shoudl have 2 less particles then the first species
void 
PairFixedPhaseClass::Read (IOSectionClass &in)
{
  string species1String;
  string species2String;
  assert(in.ReadVar("Species1",species1String));
  assert(in.ReadVar("Species2",species2String));
  Species1Num=Path.SpeciesNum(species1String);
  Species2Num=Path.SpeciesNum(species2String);
  SpeciesClass &species1=Path.Species(Species1Num);
  SpeciesClass &species2=Path.Species(Species2Num);
  NumberAssymetry=
    (species1.LastPtcl-species1.FirstPtcl)-(species2.LastPtcl-species2.FirstPtcl);
  //////HACK!  assert (NumberAssymetry==2);
  int N = Path.Species(Species2Num).LastPtcl - 
    Path.Species(Species2Num).FirstPtcl+1;
  DetMatrix.resize(N+2,N+2);
  ///Will have to be different for different forms of assymetry
  DegenerateRefSliceDeterminates.resize(6);
}

string
PairFixedPhaseClass::GetName()
{
  return "PairFixedPhase";
}
