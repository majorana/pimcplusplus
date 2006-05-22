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


#include "Josephson.h"
#include "../PathDataClass.h"

///This has to be called after pathdata knows how many
///particles it has
void JosephsonClass::Read(IOSectionClass& in)
{
  assert(in.ReadVar("alpha",alpha));
  TotalTime=0;
}

JosephsonClass::JosephsonClass(PathDataClass &pathData) : 
  ActionBaseClass (pathData)
{
}

double 
JosephsonClass::SingleAction (int slice1, int slice2,
			       const Array<int,1> &changedParticles,
			       int level)
{ 
  double levelTau=Path.tau* (1<<level);
  double totalU=0.0;
  int skip = 1<<level;
  //  cerr<<"My skip is "<<skip<<" and my level Tau is "<<levelTau<<endl;
  PathClass &Path=PathData.Path;
  //  for (int slice=0;slice<Path.NumTimeSlices()-1;slice+=skip){
  for (int slice=slice1;slice<slice2;slice+=skip){
    dVec diff=(Path(slice+skip,0)-Path(slice,0));
    double diff2=diff[0]*diff[0]/(levelTau)*1.0/16.0;
    totalU+=diff2;
  }
  //  cerr<<"My total u is "<<totalU;
  //  return totalU;


  //  cerr<<"My alpha is "<<alpha<<endl;
  //  for (int slice=0;slice<Path.NumTimeSlices();slice+=skip){
  for (int slice=slice1;slice<=slice2;slice+=skip){
    totalU -= cos(Path(slice,0)[0])*levelTau;
  }
  //  return totalU;
  double totalPotential=0.0;
  double  PiOverBeta=M_PI/(Path.tau*Path.NumTimeSlices());
  for (int sliceA=slice1;sliceA<=slice2;sliceA+=skip)
    for (int sliceB=0;sliceB<Path.NumTimeSlices();sliceB+=skip){
      if (sliceA!=sliceB){
	double diffPhi=Path(sliceA,0)[0]-Path(sliceB,0)[0];
	double diffTau=(sliceA-sliceB)*PathData.Path.tau;
	double denom=sin(PiOverBeta*diffTau);
	totalPotential+=2.0*alpha/(8.0*M_PI*M_PI)*(diffPhi*diffPhi)/(denom*denom)*levelTau*levelTau;
      }
    }
  for (int sliceA=slice1;sliceA<=slice2;sliceA+=skip)
    for (int sliceB=slice1;sliceB<=slice2;sliceB+=skip){
      if (sliceA!=sliceB){
	double diffPhi=Path(sliceA,0)[0]-Path(sliceB,0)[0];
	double diffTau=(sliceA-sliceB)*PathData.Path.tau;
	double denom=sin(PiOverBeta*diffTau);
	totalPotential-=alpha/(8.0*M_PI*M_PI)*(diffPhi*diffPhi)/(denom*denom)*levelTau*levelTau;
      }
    }
  
  double massTerm=0.0;
  for (int sliceA=slice1;sliceA<=slice2;sliceA+=skip)
    massTerm += 
      10*PathData.Path.tau*PathData.Path(sliceA,0)[0]*PathData.Path(sliceA,0)[0];
  
  totalU+=totalPotential*PiOverBeta*PiOverBeta+massTerm;
  return totalU;
}



double 
JosephsonClass::d_dBeta (int slice1, int slice2, int level)
{
  return 0.0;
}
