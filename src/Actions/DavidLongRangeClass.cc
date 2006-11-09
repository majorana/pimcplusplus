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

#include "DavidLongRangeClass.h"
#include "../PathDataClass.h"


void DavidLongRangeClass::Read(IOSectionClass &in)
{
  double myNum;
  cerr<<"The current size of the thing is "<<Path.kVecs.size()<<endl;
  uk.resize(Path.MagK.size());
  duk.resize(Path.MagK.size());
  for (int counter=0;counter<duk.size();counter++){
    duk(counter)=0.0;
  }
  string fileName;
  assert(in.ReadVar("LongRangeFile",fileName));
  ifstream infile;
  ///BUG: Currently hardcoded for actual file
  infile.open(fileName.c_str());
  cerr<<"Beginning now"<<endl;
  for (int lvl=0;lvl<2;lvl++)
    for (int isEnergy=0;isEnergy<3;isEnergy++)
      ///BUG: 
      ///Currently hard coded for 20. Ugly 
      for (int kVec=0;kVec<20;kVec++){
	infile>>myNum;
	if (lvl==0 && isEnergy==1){
	  uk(kVec)=myNum;
	}
	if ((lvl==0 && isEnergy==2) || (lvl==0 && isEnergy==0)){
	  duk(kVec)+=myNum;
	  cout<<"My energy is "<<myNum<<endl;
	}

      }
  cerr<<"Done"<<endl;
  infile.close();
  
}

inline double mag2 (const complex<double> &z)
{
  return (z.real()*z.real() + z.imag()*z.imag());
}


/// Calculates the long range part of the action using David's breakup.
/// The short range part must be supplied as a dm file without the long
/// range part in it.  It ignores active particles.
double 
DavidLongRangeClass::SingleAction (int slice1, int slice2, 
				   const Array<int,1> &activeParticles, 
				   int level)
{
  double total=0;
  double factor;
  for (int slice=slice1;slice<=slice2;slice++){
    if ((slice == slice1) || (slice==slice2))
      factor = 0.5;
    else
      factor = 1.0;
    for (int species=0; species<Path.NumSpecies(); species++) {
      Path.CalcRho_ks_Fast(slice,species);
      //      PairActionFitClass &pa = *PairMatrix(species,species);
      //      if (pa.IsLongRange()) {
      for (int ki=0; ki<Path.kVecs.size(); ki++) {
	double rhok2 = mag2(Path.Rho_k(slice,species,ki));
	double kmagnitude=sqrt(Path.kVecs(ki)[0]*Path.kVecs(ki)[0]+
			Path.kVecs(ki)[1]*Path.kVecs(ki)[1]);
	int kcounter=0;
	while (abs(kmagnitude-Path.MagK(kcounter))>1e-10)
	  kcounter++;
	assert(kcounter<Path.MagK.size());
	//	cerr<<rhok2<<" "<<Path.MagKint(kcounter)<<" "<<uk(Path.MagKint(kcounter))<<endl;
	total +=  factor*rhok2 * uk(Path.MagKint(kcounter));
	
      }
    }
  }
  //  cerr<<"The action total is "<<total<<endl;
  return total;

}

  ///Not really d_dbeta but total energy
double DavidLongRangeClass::d_dBeta (int slice1, int slice2,  int level)
{
  cerr<<"Calculating long range energy"<<endl;
  //  cerr<<"My level is "<<level<<endl;
  double total=0.0;
  double factor=1.0;
  for (int slice=slice1;slice<=slice2;slice++){
    double sliceTotal=0.0;
    //    if ((slice == slice1) || (slice==slice2))
    //      factor = 0.5;
    //    else
    //      factor = 1.0;
    //    cerr<<"Starting loop\n";
    for (int species=0; species<Path.NumSpecies(); species++) {
      Path.CalcRho_ks_Fast(slice,species);
      //      PairActionFitClass &pa = *PairMatrix(species,species);
      //      if (pa.IsLongRange()) {
      for (int ki=0; ki<Path.kVecs.size(); ki++) {
	double rhok2 = mag2(Path.Rho_k(slice,species,ki));
	double kmagnitude=sqrt(Path.kVecs(ki)[0]*Path.kVecs(ki)[0]+
			Path.kVecs(ki)[1]*Path.kVecs(ki)[1]);
	int kcounter=0;
	while (abs(kmagnitude-Path.MagK(kcounter))>1e-10)
	  kcounter++;
	assert(kcounter<Path.MagK.size());
	////	cerr<<"K counter is "<<kcounter<<" "<<Path.MagKint(kcounter)<<endl;
	//	cerr<<"My ki is "<<ki<<endl;
	//	cerr<<"The spot I'm acessing is "<<Path.MagKint(ki)<<endl;
	//	cerr<<"The value of this spot is "<<uk(Path.MagKint(ki))<<endl;
	sliceTotal +=  factor*rhok2 * duk(Path.MagKint(kcounter));
	
      }
    }
    //    cerr<<"Ending loop";
    total += sliceTotal;
    cerr<<"My slice total is "<<slice<<" "<<sliceTotal<<endl;
    
  }
  //  cerr<<"I am being called"<<endl;
  //  cerr<<"The energy total is "<<total<<endl;

  return total;


}

DavidLongRangeClass::DavidLongRangeClass (PathDataClass &pathData) :
  ActionBaseClass (pathData)
{
  //Do  nothing for now
}


string
DavidLongRangeClass::GetName()
{
  return "DavidsLongRange";
}
