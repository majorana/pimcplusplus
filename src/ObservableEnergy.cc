#include "ObservableEnergy.h"


// Fix to include final link between link M and 0
void TotalEnergyClass::Accumulate()
{
  TimesCalled++;
  if (TimesCalled % DumpFreq==0){
    WriteBlock();
  }

  if ((TimesCalled % Freq)!=0){
    return;
  }
  //Move the join to the end so we don't have to worry about permutations
  PathData.MoveJoin(PathData.NumTimeSlices()-1);
  // Loop over all links
  int numPtcls = PathData.NumParticles();
  int numLinks = PathData.NumTimeSlices()-1; 
  //  cerr<<"My numLinks is "<<numLinks<<endl;
  double tau = PathData.Action.tau;
  // Add constant part.  Note: we should really check the number of
  // dimensions. 
  double sum = 0.0;
  double vSum=0.0;
  double sSum=0.0;
  double fSum=0.0;
  double prefact=0.0;
  int NumImage=1;
  for (int ptcl=0; ptcl<numPtcls; ptcl++)
    if (PathData.Path.ParticleSpecies(ptcl).lambda != 0.0){
      sum += 1.5/tau * (double)numLinks;
      sSum+= 1.5/tau * (double)numLinks;
    }
  for (int link=0; link<numLinks; link++) {
    for (int ptcl1=0; ptcl1<numPtcls; ptcl1++) {
      // Do free-particle part
      int species1 = PathData.Path.ParticleSpeciesNum(ptcl1);
      double lambda = PathData.Path.ParticleSpecies(ptcl1).lambda;
      if (lambda != 0.0) {
	double FourLambdaTauInv=1.0/(4.0*PathData.Path.Species(species1).lambda*tau);
	dVec vel;
	vel = PathData.Path.Velocity(link, link+1, ptcl1);
	double Z = 1.0;
	dVec GaussSum=0.0;
	for (int dim=0; dim<NDIM; dim++) {
	  for (int image=-NumImage; image<=NumImage; image++) {
	    double dist = vel[dim]+(double)image*PathData.Path.GetBox()[dim];
	    GaussSum[dim] += exp(-dist*dist*FourLambdaTauInv);
	  }
	  Z *= GaussSum[dim];
	}
	dVec numSum=0.0;
	for (int dim=0;dim<NDIM;dim++){
	  for (int image=-NumImage;image<=NumImage;image++){
	    double dist = vel[dim]+(double)image*PathData.Path.GetBox()[dim];
	    numSum[dim] += 
	      (-dist*dist*FourLambdaTauInv/tau)*exp(-dist*dist*FourLambdaTauInv);
	  }
	}
	double scalarnumSum=0.0;
	for (int dim=0;dim<NDIM;dim++){
	  dVec numProd=1.0;
	  for (int dim2=0;dim2<NDIM;dim2++){
	    if (dim2!=dim){
	      numProd[dim] *= GaussSum[dim2];
	    }
	    else {
	      numProd[dim] *=  numSum[dim2];
	    }
	    
	  }
	  scalarnumSum += numProd[dim];
	}
	sum += scalarnumSum/Z; //NOT HACK!!!!
	sSum +=scalarnumSum/Z;
	//	sum += log(scalarnumSum/Z);
	//	sum -= log(Z);
      }
           
     
      for (int ptcl2=0; ptcl2<ptcl1; ptcl2++) {
	dVec r, rp;
	double rmag, rpmag;
	PathData.Path.DistDisp(link, link+1, ptcl1, ptcl2,
					 rmag, rpmag, r, rp); 

	double s2 = dot(r-rp, r-rp);
	double q = 0.5*(rmag+rpmag);
	double z = (rmag-rpmag);
	double dU;
	double dV;
	int PairIndex = 
	  PathData.Action.PairMatrix(species1, 
				     PathData.Path.ParticleSpeciesNum(ptcl2));
	dU=PathData.Action.PairActionVector(PairIndex)->dU(q, z, s2, 0);
	
	PairActionFitClass &PA=*PathData.Action.PairActionVector(PairIndex);
	sum += dU; // HACK!
	fSum += dU;
      }
    }
  }
  

  
  vSum=0.0;
  for (int ptcl1=0; ptcl1<numPtcls; ptcl1++){
    int species1=PathData.Path.ParticleSpeciesNum(ptcl1);
    for (int ptcl2=0;ptcl2<ptcl1;ptcl2++){
      int species2=PathData.Path.ParticleSpeciesNum(ptcl2);
      int PairIndex =PathData.Action.PairMatrix(species1,species2);
      for (int slice=0;slice<PathData.NumTimeSlices()-1;slice++){
	dVec r;
	double rmag;
	PathData.Path.DistDisp(slice, ptcl1, ptcl2,rmag, r); 
	double dV;
	dV=(PathData.Action.PairActionVector(PairIndex))->V(rmag);
	vSum +=dV;
      }
    }
  }
  
  
  ESum += sum; //HACK!
  VSum += vSum;
  SSum += sSum;
  FSum += fSum; 
  NumSamples++;

  /// CHECK code
  double Echeck = 0.0;
  double spring, dU = 0.0;
  spring = dU = 0.0;
  for (int slice=0; slice<numLinks; slice++) {
    double sp, du;
    PathData.Action.Energy(slice, 0, sp, du);
    dU += du;
    spring += sp;
  }
  Echeck = spring + dU;
  if (fabs(sum-Echeck) > 1.0e-10*max(1.0,fabs(sum))) {
    cerr << "sum   = " << sum << endl;
    cerr << "Echeck = " << Echeck << endl;
    cerr << "diff   = " << sum-Echeck << endl;
    cerr << "fSum   = " << fSum << endl;
    cerr << "dU     = " << dU << endl;
    cerr << "sSum   = " << sSum << endl;
    cerr << "spring = " << spring << endl << endl;
  }
}

void TotalEnergyClass::ShiftData (int NumTimeSlices)
{
  // Do nothing
}

void TotalEnergyClass::WriteBlock()
{
  double totSum;
  double totNumSamples;
  
  double myAvg = ESum/(double)NumSamples; //everybody should have the same number of samples for this to be happy
  double myVAvg= VSum/(double)NumSamples;
  double mySAvg= SSum/(double)NumSamples;
  double myFAvg= FSum/(double)NumSamples;
  double avg = PathData.Communicator.Sum(myAvg);
  double vavg =PathData.Communicator.Sum(myVAvg);
  double savg =PathData.Communicator.Sum(mySAvg);
  double favg =PathData.Communicator.Sum(myFAvg);
  avg  = avg/(double)PathData.Path.TotalNumSlices;
  vavg =vavg/(double)PathData.Path.TotalNumSlices;
  savg =savg/(double)PathData.Path.TotalNumSlices;
  favg =favg/(double)PathData.Path.TotalNumSlices;
  // Only processor 0 writes.
  if (PathData.Communicator.MyProc()==0) {
    cerr << "myAvg = " << myAvg << endl;
    cerr << "avg = " << avg << endl;
    cerr << "Pot avg = " << vavg << endl;
    cerr << "S avg = " << savg << endl;
    cerr << "U avg = " <<favg <<endl;
    if (FirstTime) {
      FirstTime = false;
      Array<double,1> dummy(1);
      dummy(0)=avg;
      IOSection.WriteVar ("TotalEnergy", dummy);
      dummy(0)=vavg;
      IOSection.WriteVar ("PotentialEnergy",dummy);
      dummy(0)=savg;
      IOSection.WriteVar ("SpringEnergy",dummy);
      dummy(0)=favg;
      IOSection.WriteVar ("DBetaEnergy",dummy);
      IOVar = IOSection.GetVarPtr("TotalEnergy");
      IOVVar= IOSection.GetVarPtr("PotentialEnergy");
      IOSVar= IOSection.GetVarPtr("SpringEnergy");
      IOUVar= IOSection.GetVarPtr("DBetaEnergy");
    }
    else {
      IOVar->Append(avg);
      IOVVar->Append(vavg);
      IOSVar->Append(savg);
      IOUVar->Append(favg);
      IOSection.FlushFile();
    }
  }
  ESum = 0.0;
  VSum = 0.0;
  SSum = 0.0;
  FSum=0.0;
  NumSamples = 0;
}

