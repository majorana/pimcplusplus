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

#include "BisectionStage.h"



// double BisectionStageClass::LogSampleProb(int slice1, int slice2, 
// 					  Array<int,1> particles)
// {
//   PathClass &Path = PathData.Path;
//   dVec rpp;
//   int skip = 1<<(BisectionLevel+1);
//   double logSampleProb=0.0;

//   double levelTau = 0.5*PathData.Action.tau*skip;
//   for (int ptclIndex=0;ptclIndex<particles.size();ptclIndex++){
//     int ptcl = particles(ptclIndex);
//     double lambda=Path.ParticleSpecies(ptcl).lambda;
//     double sigma2=(1.0*lambda*levelTau);
//     double sigma=sqrt(sigma2);
//     double prefactorOfSampleProb=0.0;//-NDIM/2.0*log(2*M_PI*sigma2);
//     for (int slice=slice1;slice<slice2;slice+=skip){
//       dVec r = Path(slice,ptcl);
//       dVec rdiff = Path.Velocity(slice, slice+skip, ptcl);

//       dVec rpp=Path(slice+(skip>>1),ptcl);
      
//       // We've ignored boundary conditions here (well we think this is
//       // fixed but we're not sure)  
//       dVec rbar=r + 0.5*rdiff;
//       dVec Delta= rpp - rbar;
//       Path.PutInBox(Delta);
      
//       double GaussProd=1.0;
//       for (int dim=0; dim<NDIM; dim++) {
// 	double GaussSum = 0.0;
// 	int NumImage = 1;
// 	for (int image=-NumImage; image <= NumImage; image++) {
// 	  double dist = Delta[dim]+(double)image*Path.GetBox()[dim];
// 	  GaussSum += exp(-0.5*dist*dist/sigma2);
// 	}
// 	GaussProd *= GaussSum;
//       }
//       logSampleProb += prefactorOfSampleProb + log(GaussProd);

//     }
//   }
//   return logSampleProb;
// }

// double BisectionStageClass::SamplePaths(int startSlice, int endSlice, 
// 					Array<int,1> particles) 
// {
//   dVec rpp;
//   int skip = 1<<(BisectionLevel+1);
//   double logNewSampleProb=0.0;
//   PathClass &Path = PathData.Path;
//   ActionClass &Action = PathData.Action;
//   double levelTau = 0.5*PathData.Action.tau*skip;
//   for (int ptclIndex=0;ptclIndex<particles.size();ptclIndex++){
//     int  ptcl=particles(ptclIndex);
//     double lambda=PathData.Path.ParticleSpecies(ptcl).lambda;
//     double sigma2=(1.0*lambda*levelTau);
//     double sigma=sqrt(sigma2);
//     double prefactorOfSampleProb=0.0;//-NDIM/2.0*log(2*M_PI*sigma2);
//     for (int slice=startSlice;slice<endSlice;slice+=skip){
//       dVec r = Path(slice,ptcl);
//       //\\     dVec rp= Path(slice+skip,ptcl);
//       dVec rpp; //\\=Path(slice+(skip>>1),ptcl);
//       //      Path.PutInBox(r);
//       //      Path.PutInBox(rp);
//       //      Path.PutInBox(rpp);
//       dVec rdiff=Path.Velocity(slice, slice+skip, ptcl);
//       dVec rbar = r + 0.5*rdiff;
//       dVec newDelta;
//       Path.Random.LocalGaussianVec(sigma,newDelta);
//       PathData.Path.PutInBox(newDelta);
//       double GaussProd=1.0;
//       for (int dim=0; dim<NDIM; dim++) {
// 	  double GaussSum = 0.0;
// 	  int NumImage = 1;
// 	  for (int image=-NumImage; image <= NumImage; image++) {
// 	    double dist = newDelta[dim]+(double)image*Path.GetBox()[dim];
// 	    GaussSum += exp(-0.5*dist*dist/sigma2);
// 	  }
// 	  GaussProd *= GaussSum;
//       }
//       logNewSampleProb += prefactorOfSampleProb + log(GaussProd);
//       rpp=rbar+newDelta;

//       ///Here we've stored the new position in the path
//       Path.SetPos(slice+(skip>>1),ptcl,rpp);
//     }
//   }
//   //  cerr<<"My logNewSampleProb is "<<logNewSampleProb<<endl;
//   return logNewSampleProb;
// }

void BisectionStageClass::WriteRatio()
{ 
  AcceptRatioVar.Write((double)NumAccepted/(double)NumAttempted);
  AcceptRatioVar.Flush();
}

void BisectionStageClass::Accept()
{
  //do nothing for now
  
}

void BisectionStageClass::Reject()
{
  //do nothing for now

}

double BisectionStageClass::Sample(int &slice1,int &slice2,
				   Array<int,1> &activeParticles)
{
  PathClass &Path = PathData.Path;
  int skip = 1<<(BisectionLevel+1);
  double levelTau = 0.5*PathData.Path.tau*skip;

  int numImages = PathData.Actions.NumImages;
  int oldSlice2=slice2;
  slice2=slice1+(1<<TotalLevels);
  double logSampleProb=0.0;
  double logOldSampleProb=0.0;
  dVec firstDelta;
  for (int ptclIndex=0;ptclIndex<activeParticles.size();ptclIndex++){
    int ptcl=activeParticles(ptclIndex);
    double lambda=PathData.Path.ParticleSpecies(ptcl).lambda;
    double sigma2=(1.0*lambda*levelTau);
    double sigma=sqrt(sigma2);
    double prefactorOfSampleProb=0.0; //-NDIM/2.0*log(2*M_PI*sigma2);
    for (int slice=slice1;slice<slice2;slice+=skip){
      SetMode(OLDMODE);
      dVec rOld=Path(slice,ptcl);
      dVec rdiffOld=Path.Velocity(slice,slice+skip,ptcl);
      dVec rbarOld=rOld+ 0.5*rdiffOld;
      dVec rppOld=Path(slice+(skip>>1),ptcl);
      dVec DeltaOld=rppOld-rbarOld;
      Path.PutInBox(DeltaOld);
      
      SetMode(NEWMODE);
      dVec r=Path(slice,ptcl);
      dVec rdiff=Path.Velocity(slice,slice+skip,ptcl);
  
      dVec rbar=r+ 0.5*rdiff;
      ///Here we've stored the new position in the path
      dVec  Delta;
      Path.Random.LocalGaussianVec(sigma,Delta);
      PathData.Path.PutInBox(Delta);
//       if (ptclIndex==0)
//       	firstDelta=Delta;
//       if (ptclIndex==1){
//       	double length=sqrt(dot(Delta,Delta));
//       	double firstLength=sqrt(dot(firstDelta,firstDelta));
//       	Delta=-1.0*firstDelta*length/firstLength;
//       }
      double GaussProd=1.0;
      double GaussProdOld=1.0;
      for (int dim=0; dim<NDIM; dim++) {
	double GaussSum = 0.0;
	double GaussSumOld =0.0;
	for (int image=-numImages; image <= numImages; image++) {
	  double dist = Delta[dim]+(double)image*Path.GetBox()[dim];
	  double distOld=DeltaOld[dim]+(double)image*Path.GetBox()[dim];
	  GaussSum += exp(-0.5*dist*dist/sigma2);
	  GaussSumOld += exp(-0.5*distOld*distOld/sigma2);
	}
	GaussProd *= GaussSum;
	GaussProdOld *= GaussSumOld;
      }
      logOldSampleProb += prefactorOfSampleProb + log(GaussProdOld);
      logSampleProb += prefactorOfSampleProb  + log(GaussProd);
      
      dVec rpp=rbar+Delta;
      ////REALLY BAD HACK!
      ///      dVec rpp=rbar+100*Delta;
      ///Here we've stored the new position in the path
      Path.SetPos(slice+(skip>>1),ptcl,rpp);
    }
  }
//   SetMode (NEWMODE);
//   SamplePaths (slice1, slice2, activeParticles);
//   SetMode(OLDMODE);
//   double oldSample = LogSampleProb (slice1, slice2, activeParticles);
//   SetMode(NEWMODE);
//   double newSample = LogSampleProb (slice1, slice2, activeParticles);
  //cerr << "oldSample = " << oldSample << endl;
  //cerr << "oldSampleProb = " << logOldSampleProb << endl;
  //cerr << "newSample = " << newSample << endl;
  //cerr << "newSampleProb = " << logSampleProb << endl;
  //  if (activeParticles.size()>1)
  //    cerr<<"Value is "<<BisectionLevel<<": "<<exp(-logSampleProb+logOldSampleProb)<<endl;
  slice2=oldSlice2;
  return exp(-logSampleProb+logOldSampleProb);
  //return (exp (-newSample + oldSample));
}

