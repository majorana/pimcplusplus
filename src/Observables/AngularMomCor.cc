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

#include "AngularMomCor.h"

///////////////////////////////////////////////////////
//Angular Correlation Integral[f(t)*f(0),{t,0,beta}]///
//where f(t) = Sum[dot[R(t),V(t)] over ptcl & slices]//
///////////////////////////////////////////////////////




void AngularMomCor::WriteBlock()
{
  //  Array<double,1> CorSum(Correlation.size());
  //  Path.Communicator.Sum(Correlation, CorSum);

  Correlation=Correlation/TotalCounts;
  CorVar.Write(Correlation);
  CorVar.Flush();
  Correlation=0;
  TotalCounts=0;

}


void AngularMomCor::Read(IOSectionClass &in)
{  

  ObservableClass::Read(in);
  string speciesName;
  Species=-1;
  assert(in.ReadVar("Species1",speciesName));
  for (int spec=0;spec<PathData.NumSpecies();spec++){ //???what is Species ?
    if (PathData.Species(spec).Name==speciesName){
      Species=spec;
    }
  }
  if (PathData.Path.Communicator.MyProc()==0){
    IOSection.WriteVar("Type","CorrelationFunction");
    IOSection.WriteVar("Cumulative", false);
  }
  Correlation.resize(PathData.Path.NumTimeSlices()+1);
  Correlation=0;
  TotalCounts=0;
  /// Now write the one-time output variables
//   if (PathData.Path.Communicator.MyProc()==0)
//     WriteInfo();

}

void AngularMomCor::Accumulate()
{
#if NDIM==3
  TotalCounts++;
#ifdef BUILD_DEV
  PathClassDev &Path= PathData.Path;
#else
  PathClass &Path= PathData.Path;
#endif
  SpeciesClass &species=PathData.Path.Species(Species);
  dVec r;
  dVec angMom1,angMom2;
  dVec v;
  int sliceskip;
  for (int skip=0;skip<=PathData.NumTimeSlices();skip++){
    for (int ptcl=species.FirstPtcl;ptcl<=species.LastPtcl;ptcl++){
      for (int slice=0;slice<=PathData.NumTimeSlices()-1;slice++){

	r=Path(slice,ptcl);

	v=Path.Velocity(slice,(slice+1)%PathData.NumTimeSlices(),ptcl);
	//cerr<<"skip "<<skip<<" ptcl "<<ptcl<<" slice "<<slice<<endl;
	//cerr<<"r="<<r<<" v= "<<v<<endl;
	angMom1=cross(r,v);

	sliceskip=(slice+skip)%PathData.NumTimeSlices();
	r=Path(sliceskip,ptcl);
	v=Path.Velocity(sliceskip,(sliceskip+1)%PathData.NumTimeSlices(),ptcl);
	angMom2=cross(r,v);
	//cerr<<"later, r="<<r<<" v= "<<v<<endl;
	Correlation(skip)+=dot(angMom1,angMom2);

	//cerr<<"Corr "<<Correlation(skip)<<endl;
      }
    }
    
  }
#else
  cerr<<"AngularMomCor doesn't work in 3d"<<endl;
  assert(NDIM==3);
#endif
}
