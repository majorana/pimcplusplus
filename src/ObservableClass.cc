#include "ObservableClass.h"
#include "Common/IO/InputOutput.h"


// // void TotalEnergyClass::Accumulate()
// // {///Currently only works for square boxes
// //   TimesCalled++;
// //   if (TimesCalled % DumpFreq==0){
// //     WriteBlock();
// //   }
// //   if ((TimesCalled % Freq)!=0){
// //     return;
// //   }
// //   int numPtcls=PathData.NumParticles();
// //   int numLinks=PathData.NumTimeSlices()-1;
// //   double tau=PathData.Action.tau;
// //   double sum=0.0;
// //   double prefact=pow(PathData.Path.Box[0],-3*numPtcls);
// //   for (int ptcl=0;ptcl<numPtcls;ptcl++){
// //       int species1 = PathData.Path.ParticleSpeciesNum(ptcl);
// //       double lambda = PathData.Path.ParticleSpecies(ptcl).lambda;
// //       for (int link=0;link<numLinks;link++){
// // 	for (int n=0;n<10;n++){
// // 	  dVec vel=PathData.Path(link+1,ptcl)-PathData.Path(link,ptcl);
// // 	  //>Velocity(link, link+1, ptcl);       
// // 	  double tempSum=1.0;
// // 	  for (int dim=0;dim<NDIM;dim++){
// // 	    double K_n=2*M_PI*n/PathData.Path.Box[dim];
// // 	    tempSum += exp(-tau*lambda*K_n*K_n)*cos(K_n*vel[dim]);
// // 	  }
// // 	  sum += tempSum;
// // 	}
// //       }
// //   }
// //   sum *= prefact;
// //   ESum += sum;
// //   NumSamples++;

// // }

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
  // Loop over all links
  int numPtcls = PathData.NumParticles();
  int numLinks = PathData.NumTimeSlices()-1;
  double tau = PathData.Action.tau;
  // Add constant part.  Note: we should really check the number of
  // dimensions. 
  double sum = 0.0;
  double prefact=0.0;
  int NumImage=4;
  for (int ptcl=0; ptcl<numPtcls; ptcl++)
    if (PathData.Path.ParticleSpecies(ptcl).lambda != 0.0)
      sum += 1.5/tau * (double)numLinks;
  for (int link=0; link<numLinks; link++) {
    for (int ptcl1=0; ptcl1<numPtcls; ptcl1++) {
      // Do free-particle part
      int species1 = PathData.Path.ParticleSpeciesNum(ptcl1);
      double lambda = PathData.Path.ParticleSpecies(ptcl1).lambda;
      if (lambda != 0.0) {
	double FourLambdaTauInv=1.0/(4.0*PathData.Path.Species(species1).lambda*tau);
	dVec vel;
	vel = PathData.DistanceTable->Velocity(link, link+1, ptcl1);
	double Z = 1.0;
	dVec GaussSum=0.0;
	for (int dim=0; dim<NDIM; dim++) {
	  for (int image=-NumImage; image<=NumImage; image++) {
	    double dist = vel[dim]+(double)image*PathData.Path.Box[dim];
	    GaussSum[dim] += exp(-dist*dist*FourLambdaTauInv);
	  }
	  Z *= GaussSum[dim];
	}
	dVec numSum=0.0;
	for (int dim=0;dim<NDIM;dim++){
	  for (int image=-NumImage;image<=NumImage;image++){
	    double dist = vel[dim]+(double)image*PathData.Path.Box[dim];
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
	sum += scalarnumSum/Z;
	
	//	sum += log(scalarnumSum/Z);
	//	sum -= log(Z);
      }
      
     
      
      for (int ptcl2=0; ptcl2<ptcl1; ptcl2++) {
	dVec r, rp;
	double rmag, rpmag;
	PathData.DistanceTable->DistDisp(link, link+1, ptcl1, ptcl2,
					 rmag, rpmag, r, rp);
	// 	dVec r1 = PathData(link,ptcl1);
	// 	dVec r2 = PathData(link,ptcl2);
	// 	dVec rp1 = PathData(link+1,ptcl1);
	// 	dVec rp2 = PathData(link+1,ptcl2);
	// 	r=r2-r1;
	// 	rp=rp2-rp1;
	// 	rmag=sqrt(dot(r,r));
	// 	rpmag=sqrt(dot(rp,rp));
	double s2 = dot(r-rp, r-rp);
	double q = 0.5*(rmag+rpmag);
	double z = (rmag-rpmag);
	double dU;
	int PairIndex = 
	  PathData.Action.PairMatrix(species1, 
				     PathData.Path.ParticleSpeciesNum(ptcl2));
	dU=PathData.Action.PairActionVector(PairIndex)->dU(q, z, s2, 0);
	PairActionFitClass &PA=*PathData.Action.PairActionVector(PairIndex);
	// 	cerr << "ptcl1 = " << ptcl1 << endl;
	// 	cerr << "ptcl2 = " << ptcl2 << endl;
	// 	cerr << "species1 = " << PathData.Path.ParticleSpecies(ptcl1).Name
	// 	     << endl;
	// 	cerr << "species2 = " << PathData.Path.ParticleSpecies(ptcl2).Name
	// 	     << endl;
	// 	cerr << "PA species1 = " << PA.Particle1.Name << endl;
	// 	cerr << "PA species2 = " << PA.Particle2.Name << endl;
	//       	if (((ptcl1==2) && (ptcl2==1)) || ((ptcl1==3) && (ptcl2==0)))
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

  double myAvg = ESum/(double)NumSamples; //everybody should have the same number of samples for this to be happy
  double avg = PathData.Communicator.Sum(myAvg);
  avg=avg/(double)PathData.Path.TotalNumSlices;


  // Only processor 0 writes.
  if (PathData.Communicator.MyProc()==0) {
    cerr << "myAvg = " << myAvg << endl;
    cerr << "avg = " << avg << endl;
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





void PairCorrelationClass::Read(IOSectionClass& IO)
{
  string type1;
  string type2;
  Species1=-1;
  Species2=-1;
  assert(IO.ReadVar("name",Name));
  assert(IO.ReadVar("type1",type1));
  assert(IO.ReadVar("type2",type2));
  assert(IO.ReadVar("freq",Freq));
  assert(IO.ReadVar("dumpFreq",DumpFreq));
  for (int spec=0;spec<PathData.NumSpecies();spec++){
    if (PathData.Species(spec).Name==type1){
      Species1=spec;
    }
    if (PathData.Species(spec).Name==type2){
      Species2=spec;
    }
  }
  assert(Species1!=-1);
  assert(Species2!=-1);
  assert(IO.OpenSection("Grid"));
  string gridType;
  double gridStart;
  double gridEnd;
  int numGridPoints;
  assert(IO.ReadVar("type",gridType));
  assert(gridType=="Linear");
  assert(IO.ReadVar("start",gridStart));
  assert(IO.ReadVar("end",gridEnd));
  assert(IO.ReadVar("NumPoints",numGridPoints));
  grid.Init(gridStart,gridEnd,numGridPoints);
  TotalCounts=0;
  Histogram.resize(numGridPoints);
  Histogram=0;
  IO.CloseSection();
    
  

}

void PairCorrelationClass::WriteBlock()
{
  Array<int,1> HistSum(Histogram.size());

  PathData.Communicator.Sum(Histogram, HistSum);

  if (PathData.Communicator.MyProc()==0) {
    if (FirstTime){
      FirstTime=false;
      IOSection.NewSection("grid");
      grid.Write(IOSection);
      IOSection.CloseSection();
      IOSection.WriteVar("Species1", PathData.Species(Species1).Name);
      IOSection.WriteVar("Species2", PathData.Species(Species2).Name);
      Array<double,2> gofrArray(1,HistSum.size());
      for (int i=0; i<(grid.NumPoints-1); i++){
	double r1 = grid(i);
	double r2 = grid(i+1);
	double r = 0.5*(r1+r2);
	double vol = 4.0*M_PI/3 * (r2*r2*r2-r1*r1*r1);
	gofrArray(0,i) = (double) HistSum(i) / (vol*TotalCounts);
      }
      IOSection.WriteVar("gofr",gofrArray);
      IOVar = IOSection.GetVarPtr("gofr");
    }
    else {
      Array<double,1> gofrArray(HistSum.size());
      for (int i=0; i<(grid.NumPoints-1); i++){
	double r1 = grid(i);
	double r2 = grid(i+1);
	double r = 0.5*(r1+r2);
	double vol = 4.0*M_PI/3 * (r2*r2*r2-r1*r1*r1);
	gofrArray(i) = (double) HistSum(i) / (vol*TotalCounts);
      }
      IOVar->Append(gofrArray);
    }
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
  SpeciesClass &species1=PathData.Path.Species(Species1);
  SpeciesClass &species2=PathData.Path.Species(Species2);

  TimesCalled++;

  if (TimesCalled % DumpFreq==0){
    WriteBlock();
  }
  if ((TimesCalled % Freq)!=0){
    return;
  }


  /// HACK HACK HACK
  if (Species1==Species2) {
    for (int slice=0;slice<PathData.NumTimeSlices();slice++) 
      for (int ptcl1=species1.FirstPtcl;ptcl1<=species1.LastPtcl;ptcl1++)
	for (int ptcl2=ptcl1+1;ptcl2<=species1.LastPtcl;ptcl2++){
	  dVec r1=PathData(slice,ptcl1);
	  dVec r2=PathData(slice,ptcl2);
	  
	  dVec disp;
	  double dist;
	  PathData.DistanceTable->DistDisp(slice,ptcl1,ptcl2,dist,disp);
	
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
	}
  }
  else {
    for (int slice=0;slice<PathData.NumTimeSlices();slice++) 
      for (int ptcl1=species1.FirstPtcl;ptcl1<=species1.LastPtcl;ptcl1++)
	for (int ptcl2=species2.FirstPtcl;ptcl2<=species2.LastPtcl;ptcl2++){
	  dVec r1=PathData(slice,ptcl1);
	  dVec r2=PathData(slice,ptcl2);
	  
	  dVec disp;
	  double dist;
	  PathData.DistanceTable->DistDisp(slice,ptcl1,ptcl2,dist,disp);
	  
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
	}
  }
}



void PairCorrelationClass::Initialize()
{
  int numTimeSlices=PathData.Path.NumTimeSlices();//<--What is that there for...it's not used
  grid.Init (0.0, 12.0, 100);
  TotalCounts = 0;
  Histogram.resize(100);
  Histogram = 0;
}


void PathDumpClass::Accumulate()
{//Do nothing!
}
void PathDumpClass::WriteBlock()
{
  int numPtcls = PathData.NumParticles();
  int numTimeSlices = PathData.NumTimeSlices();
  if (FirstTime){
    FirstTime=false;

    Array<string,1> speciesNames(numPtcls);
    for (int speciesIndex=0; speciesIndex<PathData.NumSpecies(); speciesIndex++) {
      SpeciesClass &species = PathData.Path.Species(speciesIndex);
      for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++)
      	speciesNames(ptcl)=species.Name;
    }
    IOSection.WriteVar("SpeciesNames", speciesNames);
    
    Array<double,4> pathArray(1,numPtcls,numTimeSlices,NDIM);
    for (int ptcl=0;ptcl<numPtcls;ptcl++){
      for (int slice=0;slice<numTimeSlices;slice++){
	for (int dim=0;dim<NDIM;dim++){
	  pathArray(0,ptcl,slice,dim)=PathData(slice,ptcl)[dim];
	}
      }
    }
    
    
    // Write the first path here
    IOSection.WriteVar("Path",pathArray);
    // Now get the pointer to it
    IOVar = IOSection.GetVarPtr("Path");
  }
  else {
    // Append the new path here.
    Array<double,3> pathArray(numPtcls,numTimeSlices,NDIM);
    for (int ptcl=0;ptcl<numPtcls;ptcl++)
      for (int slice=0;slice<numTimeSlices;slice++)
	for (int dim=0;dim<NDIM;dim++)
	  pathArray(ptcl,slice,dim)=PathData(slice,ptcl)[dim];
    IOVar->Append(pathArray);
  }
}
