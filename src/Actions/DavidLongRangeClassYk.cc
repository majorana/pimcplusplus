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

#include <Common/MPI/Communication.h>
#include "DavidLongRangeClassYk.h"
#include "../PathDataClass.h"

bool fequals(double a,double b,double tol)
{
  return abs(a-b)<tol;
}

bool vecEquals(dVec &a, dVec &b,double tol)
{
  bool equals=true;
  for (int dim=0;dim<NDIM;dim++)
    equals= equals && fequals(a[dim],b[dim],tol);
  return equals;
}

void DavidLongRangeClassYk::BuildRPA_SingleSpecies()
{
  int speciesNum=0;
  double vol=1.0;
  for (int dim=0;dim<NDIM;dim++)
    vol*=Path.GetBox()[dim];
  DavidPAClass &pa(*(DavidPAClass*)PairArray(speciesNum));
  int ncomps=Path.Species(speciesNum).NumParticles;
  double lambda=Path.Species(speciesNum).lambda;
  cerr<<"lambda is "<<lambda<<endl;
  cerr<<"ncomps is "<<ncomps<<endl;
  cerr<<"Vol is "<<vol<<endl;
  for (int i=0;i<pa.kVals.size();i++){
    double k=pa.kVals(i);
    double arg=1.0+ncomps*2.0*pa.uk_long(i)/(lambda*k*k*vol);
    double theta=0.5*k*k*lambda*Path.tau;
    double q;
    double s;
    double tn;
    if (arg<0){
      q=sqrt(-arg);
      tn=tan(q*theta);
      s=q*(1-q*tn)/(q+tn);
    }
    else if (arg==0)
      s=1.0/(1.0+theta);
    else {
      q=sqrt(arg);
      tn=tanh(theta*q);
      s=q*(1.0+q*tn)/(q+tn);
    }
    for (int j=0;j<Path.kVecs.size();j++){
      if (fequals(sqrt(blitz::dot(Path.kVecs(j),Path.kVecs(j))),k,1e-10)){
	uk(j)=(-1.0+s)/ncomps;
	duk(j)=pa.uk_long(i)/(vol)-lambda*k*k*uk(j)*((1.0+0.5*ncomps*uk(j)));
	Vlong_k(j)=pa.uk_long(i)/(vol);
	cerr<<"VLONG IS "<<Vlong_k(j)<<endl;
      }
    }
  }
  cerr<<"I have built the rpa"<<endl;
  for (int i=0;i<uk.size();i++)
    cerr<<"KVecs: "<<sqrt(blitz::dot(Path.kVecs(i),Path.kVecs(i)))<<" "<<uk(i)<<" "<<duk(i)<<endl;
}


void DavidLongRangeClassYk::ReadYk()
{ 
  assert(PairArray.size()==1);
  DavidPAClass &pa(*((DavidPAClass*)PairArray(0)));
  assert(pa.LongRangeDim==NDIM);
  for (int dim=0;dim<NDIM;dim++)
    assert(pa.LongRangeBox(dim)==Path.GetBox()[dim]);
  assert(pa.LongRangeMass1==pa.LongRangeMass2);
  assert(pa.LongRangeMass1==Path.Species(0).lambda);
  uk.resize(Path.kVecs.size());
  duk.resize(Path.kVecs.size());
  Vlong_k.resize(Path.kVecs.size());
  BuildRPA_SingleSpecies();  

//   uk=-999;
//   for (int i=0;i<Path.kVals.size();i++){
//     for (int j=0;j<pa.kVals.extent(0);j++){
//       dVec temp;
//       for (int dim=0;dim<NDIM;dim++)
// 	temp[dim]=pa.kVecs(j,dim);
//       if (vecEquals(temp,Path.kVecs(i),1e-10)){
// 	uk(i)=pa.uk_long(j);
//       }
//     }
//   }
//   for (int i=0;i<uk.size();i++)
//     if (uk(i)==-999){
//       cerr<<"There is a missing value of uk"<<endl;
//       assert(1==2);
//     }
}

void DavidLongRangeClassYk::Read(IOSectionClass &in)
{
  assert(1==2);
  double myNum;
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
  cerr<<" of "<<fileName.c_str()<<endl;
  string isRank;
  infile >> isRank;
  assert(isRank=="RANK");
  int is5; infile >> is5; assert(is5==5);
  int numkVec; infile >> numkVec;
  int is1; infile >>is1; assert(is1==1); infile >>is1; assert(is1==1);
  int is3; infile >> is3; assert(is3==3);
  int numLvl; infile >> numLvl; 
  infile >> isRank;
  assert(isRank=="BEGIN");
  infile >> isRank;
  assert(isRank=="k-space");
  infile >>isRank;
  assert(isRank=="action");
		    
  cerr<<"Loading David Long Range: "<<numLvl<<" "<<numkVec<<endl;
  //  for (int lvl=0;lvl<2;lvl++)
  for (int lvl=0;lvl<numLvl;lvl++)
    for (int isEnergy=0;isEnergy<3;isEnergy++)
      ///BUG: 
      ///Currently hard coded for 20. Ugly 
      //      for (int kVec=0;kVec<20;kVec++){
      for (int kVec=0;kVec<numkVec;kVec++){
	infile>>myNum;
	cerr<<"My num is "<<myNum<<endl;
	if (lvl==0 && isEnergy==1){
	  uk(kVec)=myNum;
	}
	if ((lvl==0 && isEnergy==2) || (lvl==0 && isEnergy==0)){
	  duk(kVec)+=myNum;
	  cout<<"My energy is "<<myNum<<endl;
	}

      }
  //  cerr<<"Done"<<endl;
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
DavidLongRangeClassYk::SingleAction (int slice1, int slice2, 
				   const Array<int,1> &activeParticles, 
				   int level)

{
  int species=0;
  if (GetMode() == NEWMODE)
    Path.UpdateRho_ks(slice1, slice2, activeParticles, level);

  double total=0;
  double factor;
  for (int ki=0; ki<Path.kVecs.size(); ki++) {
    for (int slice=slice1;slice<=slice2;slice++){
      if ((slice == slice1) || (slice==slice2))
	factor = 0.5;
      else
	factor = 1.0;
      double rhok2 = mag2(Path.Rho_k(slice,species,ki));
      total +=  factor*rhok2 * uk(ki);
    }
  }

  return total;
    
}

  ///Not really d_dbeta but total energy
double DavidLongRangeClassYk::d_dBeta (int slice1, int slice2,  int level)
{
  double total=0.0;
  double factor=1.0;
  for (int slice=slice1;slice<=slice2;slice++){
    double sliceTotal=0.0;
    if ((slice == slice1) || (slice==slice2))
      factor = 0.5;
    else
      factor = 1.0;
    for (int species=0; species<Path.NumSpecies(); species++) {
      Path.CalcRho_ks_Fast(slice,species);
      for (int ki=0; ki<Path.kVecs.size(); ki++) {
	double rhok2 = mag2(Path.Rho_k(slice,species,ki));
	sliceTotal +=  factor*rhok2 * duk(ki);
      }
    }
    total += sliceTotal;
  }
  return total;

}


  ///Not really d_dbeta but total energy
double DavidLongRangeClassYk::V (int slice1, int slice2,  int level)
{
  double total=0.0;
  double factor=1.0;
  for (int slice=slice1;slice<=slice2;slice++){
    double sliceTotal=0.0;
    if ((slice == slice1) || (slice==slice2))
      factor = 0.5;
    else
      factor = 1.0;
    for (int species=0; species<Path.NumSpecies(); species++) {
      Path.CalcRho_ks_Fast(slice,species);
      for (int ki=0; ki<Path.kVecs.size(); ki++) {
	double rhok2 = mag2(Path.Rho_k(slice,species,ki));
	sliceTotal +=  factor*rhok2 * Vlong_k(ki);
      }
    }
    total += sliceTotal;
  }
  return total;

}



DavidLongRangeClassYk::DavidLongRangeClassYk(PathDataClass &pathData,
					     Array<PairActionFitClass* ,2> &pairMatrix,
					     Array<PairActionFitClass*, 1> &pairArray) :
  ActionBaseClass (pathData), PairMatrix(pairMatrix), PairArray(pairArray)
{

}



string
DavidLongRangeClassYk::GetName()
{
  return "DavidsLongRange";
}
