#include "ObservableClass.h"

// Fix to include final link between link M and 0
void TotalEnergyClass::Accumulate()
{
  // Loop over all links
  int numPtcls = PathData.NumParticles();
  int numLinks = PathData.NumTimeSlices()-1;
  double tau = PathData.Action.tau;
  // Add constant part.  Note: we should really check the number of
  // dimensions. 
  double sum = 0.0;
  for (int ptcl=0; ptcl<numPtcls; ptcl++)
    if (PathData.Path.ParticleSpecies(ptcl).lambda != 0.0)
      sum += 1.5/tau * (double)numLinks;
  for (int link=0; link<numLinks; link++) {
    for (int ptcl1=0; ptcl1<numPtcls; ptcl1++) {
      // Do free-particle part
      int species1 = PathData.Path.ParticleSpeciesNum(ptcl1);
      double lambda = PathData.Path.ParticleSpecies(ptcl1).lambda;
      if (lambda != 0.0) {
	dVec vel = PathData.DistanceTable->Velocity(link, link+1, ptcl1);
	sum -= dot(vel,vel)/(4.0*lambda*tau*tau);
      }
      for (int ptcl2=0; ptcl2<ptcl1; ptcl2++) {
	dVec r, rp;
	double rmag, rpmag;
	PathData.DistanceTable->DistDisp(link, link+1, ptcl1, ptcl2,
					 rmag, rpmag, r, rp);
	double s2 = dot(r-rp, r-rp);
	double q = 0.5*(rmag+rpmag);
	double z = (rmag-rpmag);
	int PairIndex = 
	  PathData.Action.PairMatrix(species1, 
				     PathData.Path.ParticleSpeciesNum(ptcl2));
	double dU;

	dU=PathData.Action.PairActionVector(PairIndex)->dU(q, z, s2, 0);
	sum += dU;
      }
    }
  }
  ESum += sum;
  NumSamples++;
}

void TotalEnergyClass::ShiftData (int NumTimeSlices)
{
  // Do nothing
}

void TotalEnergyClass::WriteBlock()
{
  double totSum;
  double totNumSamples;

  double myAvg = ESum/(double)NumSamples;
  double avg = PathData.Communicator.Sum(myAvg);
  
  cerr << "myAvg = " << myAvg << endl;
  cerr << "avg = " << avg << endl;

  // Only processor 0 writes.
  if (PathData.Communicator.MyProc()==0) {
    if (FirstTime) {
      FirstTime = false;
      Array<double,1> dummy(1);
      dummy(0)=avg;
      IOSection.WriteVar ("TotalEnergy", dummy);
      IOVar = IOSection.GetVarPtr("TotalEnergy");
    }
    else
      IOVar->Append(avg);
  }
  ESum = 0.0;
  NumSamples = 0;
}







void PairCorrelationClass::WriteBlock()
{
  if (FirstTime){
    FirstTime=false;
    IOSection.NewSection("grid");
    grid.Write(IOSection);
    IOSection.CloseSection();
    Array<double,2> gofrArray(1,Histogram.size());
    for (int i=0; i<(grid.NumPoints-1); i++){
	double r1 = grid(i);
	double r2 = grid(i+1);
	double r = 0.5*(r1+r2);
	double vol = 4.0*M_PI/3 * (r2*r2*r2-r1*r1*r1);
	gofrArray(0,i) = (double) Histogram(i) / (vol*TotalCounts);
    }
    IOSection.WriteVar("gofr",gofrArray);
    IOVar = IOSection.GetVarPtr("gofr");
  }
  else {
    Array<double,1> gofrArray(Histogram.size());
    for (int i=0; i<(grid.NumPoints-1); i++){
	double r1 = grid(i);
	double r2 = grid(i+1);
	double r = 0.5*(r1+r2);
	double vol = 4.0*M_PI/3 * (r2*r2*r2-r1*r1*r1);
	gofrArray(i) = (double) Histogram(i) / (vol*TotalCounts);
    }
    IOVar->Append(gofrArray);
  }
  

		       

}


void PairCorrelationClass::Print()
{
  for (int i=0; i<(grid.NumPoints-1); i++)
    {
      double r1 = grid(i);
      double r2 = grid(i+1);
      double r = 0.5*(r1+r2);
      double vol = 4.0*M_PI/3 * (r2*r2*r2-r1*r1*r1);
      double gofr = (double) Histogram(i) / (vol*TotalCounts);
      fprintf (stderr, "%1.12e %1.12e\n", r, gofr);
    }
}



/// Fix me to accumulate data only between the two species I'm
/// interested in.
void PairCorrelationClass::Accumulate()
{
  Array<bool,1> doPtcl2(PathData.NumParticles());

  /// HACK HACK HACK
  for (int slice=0;slice<PathData.NumTimeSlices();slice++){
    for (int ptcl1=0;ptcl1<PathData.NumParticles();ptcl1++){
      for (int ptcl2=ptcl1+1;ptcl2<PathData.NumParticles();ptcl2++){
	dVec r1=PathData(slice,ptcl1);
	dVec r2=PathData(slice,ptcl2);
	
	dVec disp;
	double dist;
	PathData.DistanceTable->DistDisp(slice,ptcl1,ptcl2,dist,disp);
	//	cout<<"For periodic bc "<<disp<<" ";
	//PathData.Path.MinDistance(PathData(slice,ptcl1),
	//			  PathData(slice,ptcl2),dist,disp);
	//cout<<disp<<endl;
#ifdef OLDDEBUG
	dVec dispDummy=r2-r1;
	double distDummy=sqrt(dot(dispDummy,dispDummy));
	for (int i=0; i<NDIM; i++)
	  if (disp[i] != dispDummy[i])
	    cerr << "Bad bad evil inconsistency is DistTable.\n";
	if (dist != distDummy)
	  cerr << "Bad bad evil inconsistency is DistTable.\n";
#endif
	
	TotalCounts++;
	if (dist<grid.End){
	  int index=grid.ReverseMap(dist);
	  Histogram(index)++;
	}
	else {
	  //	  cerr<<"The distance is really "<<dist<<endl;
	}
      }
    }
  }
}




void PairCorrelationClass::Initialize()
{
  int numTimeSlices=PathData.Path.NumTimeSlices();
  grid.Init (0.0, 5.0, numTimeSlices);
  TotalCounts = 0;
  Histogram.resize(numTimeSlices-1);
  Histogram = 0;
}



