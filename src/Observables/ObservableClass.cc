#include "ObservableClass.h"
#include "../Common/IO/InputOutput.h"



void ObservableClass::WriteInfo()
{
  IOSection.WriteVar("Description",Description);
}

void ObservableClass::Read(IOSectionClass &in)
{
  assert(in.ReadVar("Name",Name));
  if(!(in.ReadVar("Description",Description))){
    Description="No description available";
  }
  
}



void PairCorrelationClass::Read(IOSectionClass& in)
{
  
  ObservableClass::Read(in);
  string species1Name;
  string species2Name;
  Species1=-1;
  Species2=-1;
  assert(in.ReadVar("Species1",species1Name));
  assert(in.ReadVar("Species2",species2Name));
  assert(in.ReadVar("freq",Freq));
  assert(in.ReadVar("dumpFreq",DumpFreq));
  for (int spec=0;spec<PathData.NumSpecies();spec++){
    if (PathData.Species(spec).Name==species1Name){
      Species1=spec;
    }
    if (PathData.Species(spec).Name==species2Name){
      Species2=spec;
    }
  }
  assert(Species1!=-1);
  assert(Species2!=-1);
  assert(in.OpenSection("Grid"));
  string gridType;
  double gridStart;
  double gridEnd;
  int numGridPoints;
  assert(in.ReadVar("type",gridType));
  assert(gridType=="Linear");
  bool readStartGrid=in.ReadVar("start",gridStart);
  bool readEndGrid=in.ReadVar("end",gridEnd);
  if (!readStartGrid)
    gridStart=0.0;
  if (!readEndGrid){
    if (PathData.Path.IsPeriodic(0)){
      gridEnd=PathData.Path.GetBox()[0];
    }
    else {
      cerr<<"I don't know where you want me to end this grid"<<endl;
      assert(1==2);
    }
  }
  assert(in.ReadVar("NumPoints",numGridPoints));
  grid.Init(gridStart,gridEnd,numGridPoints);
  TotalCounts=0;
  Histogram.resize(numGridPoints-1);
  Histogram=0;
  in.CloseSection();
}



void nofrClass::Read(IOSectionClass& in)
{
  
  ObservableClass::Read(in);
  string species1Name;
  string species2Name;
  Species1=-1;
  Species2=-1;
  assert(in.ReadVar("Species1",species1Name));
  assert(in.ReadVar("Species2",species2Name));
  assert(in.ReadVar("freq",Freq));
  assert(in.ReadVar("dumpFreq",DumpFreq));
  for (int spec=0;spec<PathData.NumSpecies();spec++){
    if (PathData.Species(spec).Name==species1Name){
      Species1=spec;
    }
    if (PathData.Species(spec).Name==species2Name){
      Species2=spec;
    }
  }
  assert(Species1!=-1);
  assert(Species2!=-1);
  assert(in.OpenSection("Grid"));
  string gridType;
  double gridStart;
  double gridEnd;
  int numGridPoints;
  assert(in.ReadVar("type",gridType));
  assert(gridType=="Linear");
  assert(in.ReadVar("start",gridStart));
  assert(in.ReadVar("end",gridEnd));
  assert(in.ReadVar("NumPoints",numGridPoints));
  grid.Init(gridStart,gridEnd,numGridPoints);
  TotalCounts=0;
  Histogram.resize(numGridPoints-1);
  Histogram=0;
  in.CloseSection();
}


void PairCorrelationClass::WriteInfo()
{
  ObservableClass::WriteInfo();
  IOSection.NewSection("grid");
  grid.Write(IOSection);
  IOSection.CloseSection();

  int numBins = grid.NumPoints-1;
  Array<double,1> r(numBins);
  for (int i=0; i<numBins; i++) {
    double ra = grid(i);
    double rb = grid(i+1);
    r(i) = 0.75 * (rb*rb*rb*rb-ra*ra*ra*ra)/(rb*rb*rb-ra*ra*ra);
  }
  IOSection.WriteVar("x", r);
  IOSection.WriteVar("xlabel", "r");
  IOSection.WriteVar("ylabel", "g(r)");
  IOSection.WriteVar("Species1", PathData.Species(Species1).Name);
  IOSection.WriteVar("Species2", PathData.Species(Species2).Name);
  IOSection.WriteVar("Type","CorrelationFunction");
}



void nofrClass::WriteInfo()
{
  ObservableClass::WriteInfo();
  IOSection.NewSection("grid");
  grid.Write(IOSection);
  IOSection.CloseSection();

  int numBins = grid.NumPoints-1;
  Array<double,1> r(numBins);
  for (int i=0; i<numBins; i++) {
    double ra = grid(i);
    double rb = grid(i+1);
    r(i) = 0.75 * (rb*rb*rb*rb-ra*ra*ra*ra)/(rb*rb*rb-ra*ra*ra);
  }
  IOSection.WriteVar("x", r);
  IOSection.WriteVar("xlabel", "r");
  IOSection.WriteVar("ylabel", "g(r)");
  IOSection.WriteVar("Species1", PathData.Species(Species1).Name);
  IOSection.WriteVar("Species2", PathData.Species(Species2).Name);
  IOSection.WriteVar("Type","CorrelationFunction");
}
  

void PairCorrelationClass::WriteBlock()
{
  Array<int,1> HistSum(Histogram.size());
  double norm=0.0;
  int N1 = PathData.Species(Species1).NumParticles;
  int N2 = PathData.Species(Species2).NumParticles;
  norm = TotalCounts * N1*N2/PathData.Path.GetVol();

  if (Species1==Species2){//Normalizes things when species are same
    norm *= 0.5;
    //double N = (double)PathData.Species(Species1).NumParticles;
    //norm = N/(N-1.0)/PathData.Path.GetVol();
    //norm = (double)PathData.Species(Species1).NumParticles/(double)(PathData.Species(Species1).NumParticles-1)*1.0/PathData.Path.GetVol();
  }
  
  PathData.Communicator.Sum(Histogram, HistSum);
  if (PathData.Communicator.MyProc()==0) {
    if (FirstTime) {
      FirstTime=false;
      WriteInfo();
      Array<double,2> gofrArray(1,HistSum.size());
      for (int i=0; i<grid.NumPoints-1; i++){
	double r1 = grid(i);
	double r2 = (i<(grid.NumPoints-1)) ? grid(i+1):(2.0*grid(i)-grid(i-1));
	double r = 0.5*(r1+r2);
	double binVol = 4.0*M_PI/3 * (r2*r2*r2-r1*r1*r1);
	gofrArray(0,i) = (double) HistSum(i) / (binVol*norm);
      }
      IOSection.WriteVar("y",gofrArray);
      IOVar = IOSection.GetVarPtr("y");
    }
    else {
      Array<double,1> gofrArray(HistSum.size());
      for (int i=0; i<grid.NumPoints-1; i++){
	double r1 = grid(i);
	double r2 = (i<(grid.NumPoints-1)) ? grid(i+1):(2.0*grid(i)-grid(i-1));
	double r = 0.5*(r1+r2);
	double binVol = 4.0*M_PI/3 * (r2*r2*r2-r1*r1*r1);
	gofrArray(i) = (double) HistSum(i) / (binVol*norm);
      }
      IOVar->Append(gofrArray);
    }
  }
}




void nofrClass::WriteBlock()
{
  Array<double,1> HistSum(Histogram.size());
  double norm=0.0;
  int N1 = PathData.Species(Species1).NumParticles;
  int N2 = PathData.Species(Species2).NumParticles;
  norm = TotalCounts * N1*N2/PathData.Path.GetVol();

  if (Species1==Species2){//Normalizes things when species are same
    norm *= 0.5;
    //double N = (double)PathData.Species(Species1).NumParticles;
    //norm = N/(N-1.0)/PathData.Path.GetVol();
    //norm = (double)PathData.Species(Species1).NumParticles/(double)(PathData.Species(Species1).NumParticles-1)*1.0/PathData.Path.GetVol();
  }
  
  PathData.Communicator.Sum(Histogram, HistSum);
  if (PathData.Communicator.MyProc()==0) {
    if (FirstTime) {
      FirstTime=false;
      WriteInfo();
      Array<double,2> gofrArray(1,HistSum.size());
      for (int i=0; i<grid.NumPoints-1; i++){
	double r1 = grid(i);
	double r2 = (i<(grid.NumPoints-1)) ? grid(i+1):(2.0*grid(i)-grid(i-1));
	double r = 0.5*(r1+r2);
	double binVol = 4.0*M_PI/3 * (r2*r2*r2-r1*r1*r1);
	gofrArray(0,i) = (double) HistSum(i) / (norm);
      }
      IOSection.WriteVar("y",gofrArray);
      IOVar = IOSection.GetVarPtr("y");
    }
    else {
      Array<double,1> gofrArray(HistSum.size());
      for (int i=0; i<grid.NumPoints-1; i++){
	double r1 = grid(i);
	double r2 = (i<(grid.NumPoints-1)) ? grid(i+1):(2.0*grid(i)-grid(i-1));
	double r = 0.5*(r1+r2);
	double binVol = 4.0*M_PI/3 * (r2*r2*r2-r1*r1*r1);
	gofrArray(i) = (double) HistSum(i) / (norm);
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



void nofrClass::Print()
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
    for (int slice=0;slice<PathData.NumTimeSlices()-1;slice++) {
      TotalCounts++;
      for (int ptcl1=species1.FirstPtcl;ptcl1<=species1.LastPtcl;ptcl1++)
	for (int ptcl2=ptcl1+1;ptcl2<=species1.LastPtcl;ptcl2++){
	  
	  dVec disp;
	  double dist;
	  PathData.Path.DistDisp(slice,ptcl1,ptcl2,dist,disp);
	
	  #ifdef OLDDEBUG
	  dVec r1=PathData(slice,ptcl1);
	  dVec r2=PathData(slice,ptcl2);
	  dVec dispDummy=r2-r1;
	  double distDummy=sqrt(dot(dispDummy,dispDummy));
	  for (int i=0; i<NDIM; i++)
	    if (disp[i] != dispDummy[i]){
	      cerr << "Bad bad evil inconsistency is DispTable.\n";
	      cerr<<r1<<" "<<r2<<" "<<dispDummy<<" "<<disp<<endl;
	    }
	  if (dist != distDummy)
	    cerr << "Bad bad evil inconsistency is DistTable.\n";
	  #endif
	  if (dist<grid.End){
	    int index=grid.ReverseMap(dist);
	    Histogram(index)++;
	  } 
	}
    }
  }
  else {
    for (int slice=0;slice<PathData.NumTimeSlices()-1;slice++) {
      TotalCounts++;
      for (int ptcl1=species1.FirstPtcl;ptcl1<=species1.LastPtcl;ptcl1++)
	for (int ptcl2=species2.FirstPtcl;ptcl2<=species2.LastPtcl;ptcl2++){
	  
	  dVec disp;
	  double dist;
	  PathData.Path.DistDisp(slice,ptcl1,ptcl2,dist,disp);
	  
#ifdef OLDDEBUG
	  dVec r1=PathData(slice,ptcl1);
	  dVec r2=PathData(slice,ptcl2);
	  dVec dispDummy=r2-r1;
	  double distDummy=sqrt(dot(dispDummy,dispDummy));
	  for (int i=0; i<NDIM; i++)
	    if (disp[i] != dispDummy[i])
	      cerr << "Bad bad evil inconsistency in DistTable.\n";
	  if (dist != distDummy)
	    cerr << "Bad bad evil inconsistency in DistTable.\n";
#endif
	  if (dist<grid.End){
	    int index=grid.ReverseMap(dist);
	    Histogram(index)++;
	  }
	}
    }
  }
}




/// Fix me to accumulate data only between the two species I'm
/// interested in.
void nofrClass::Accumulate()
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


    
  int openLink=(int)(PathData.Path.OpenLink);
  int openPtcl=(int)(PathData.Path.OpenPtcl);
  //     if (openPtcl==-1){
  //       cerr<<"Something is broken"<<endl;
//       openPtcl=0;
//     }
  //    disp=PathData.Path.GetPosOpen(openLink+1,openPtcl)-
    //      PathData.Path.GetPosOpen(openLink,openPtcl);
    //    disp=PathData.Path.Velocity(openLink,openLink+1,openPtcl);
  dVec disp;
  double dist;
  PathData.Path.DistDisp(openLink,openPtcl,PathData.Path.NumParticles(),
			 dist,disp); //This is distance between head and tail!
  //  cerr<<PathData.Path(openLink,openPtcl)<<" "<<PathData.Path(openLink,PathData.Path.NumParticles())<<" "<<openLink<<" "<<openPtcl<<endl;
  
    //    cerr<<"Disp is "<<disp<<" "<<PathData.Path(0,0)<<" "
    //	<<PathData.Path(PathData.Path.NumTimeSlices()-1,0)<<endl;
    ////////    PathData.Path.PutInBox(disp);
    //    dist=sqrt(dot(disp,disp));
  if (dist<grid.End){
    int index=grid.ReverseMap(dist);
    //    Histogram(index)=Histogram(index)+(0.5)/(dist*dist)+(0.9*exp(-dist*dist)+0.1);
    Histogram(index)=Histogram(index)+1.0;
  }  
  TotalCounts++;  
  return; 


  /// HACK HACK HACK
  if (Species1==Species2) {
    for (int slice=0;slice<PathData.NumTimeSlices()-1;slice++) {
      TotalCounts++;
      for (int ptcl1=species1.FirstPtcl;ptcl1<=species1.LastPtcl;ptcl1++)
	for (int ptcl2=ptcl1+1;ptcl2<=species1.LastPtcl;ptcl2++){
	  
	  dVec disp;
	  double dist;
	  PathData.Path.DistDisp(slice,ptcl1,ptcl2,dist,disp);
	
	  #ifdef OLDDEBUG
	  dVec r1=PathData(slice,ptcl1);
	  dVec r2=PathData(slice,ptcl2);
	  dVec dispDummy=r2-r1;
	  double distDummy=sqrt(dot(dispDummy,dispDummy));
	  for (int i=0; i<NDIM; i++)
	    if (disp[i] != dispDummy[i]){
	      cerr << "Bad bad evil inconsistency is DispTable.\n";
	      cerr<<r1<<" "<<r2<<" "<<dispDummy<<" "<<disp<<endl;
	    }
	  if (dist != distDummy)
	    cerr << "Bad bad evil inconsistency is DistTable.\n";
	  #endif
	  if (dist<grid.End){
	    int index=grid.ReverseMap(dist);
	    Histogram(index)++;
	  } 
	}
    }
  }
  else {
    for (int slice=0;slice<PathData.NumTimeSlices()-1;slice++) {
      TotalCounts++;
      for (int ptcl1=species1.FirstPtcl;ptcl1<=species1.LastPtcl;ptcl1++)
	for (int ptcl2=species2.FirstPtcl;ptcl2<=species2.LastPtcl;ptcl2++){
	  
	  dVec disp;
	  double dist;
	  PathData.Path.DistDisp(slice,ptcl1,ptcl2,dist,disp);
	  
#ifdef OLDDEBUG
	  dVec r1=PathData(slice,ptcl1);
	  dVec r2=PathData(slice,ptcl2);
	  dVec dispDummy=r2-r1;
	  double distDummy=sqrt(dot(dispDummy,dispDummy));
	  for (int i=0; i<NDIM; i++)
	    if (disp[i] != dispDummy[i])
	      cerr << "Bad bad evil inconsistency in DistTable.\n";
	  if (dist != distDummy)
	    cerr << "Bad bad evil inconsistency in DistTable.\n";
#endif
	  if (dist<grid.End){
	    int index=grid.ReverseMap(dist);
	    Histogram(index)++;
	  }
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


void nofrClass::Initialize()
{
  int numTimeSlices=PathData.Path.NumTimeSlices();//<--What is that there for...it's not used
  grid.Init (0.0, 12.0, 100);
  TotalCounts = 0;
  Histogram.resize(100);
  Histogram = 0.0;
}


void PathDumpClass::Accumulate()
{
  TimesCalled++;
  if (TimesCalled % 5000==0){
    WriteBlock();
  }
}

void PathDumpClass::Read(IOSectionClass &in)
{
  ObservableClass::Read(in);
}

void PathDumpClass::WriteBlock()
{
  //Move the join to the end so we don't have to worry about permutations
  PathData.MoveJoin(PathData.NumTimeSlices()-1);

  int numPtcls = PathData.NumParticles();
  int numTimeSlices = PathData.NumTimeSlices();
  if (FirstTime){
    FirstTime=false;
    WriteInfo();
    Array<string,1> speciesNames(numPtcls);
    for (int speciesIndex=0;speciesIndex<PathData.NumSpecies();speciesIndex++){
      SpeciesClass &species = PathData.Path.Species(speciesIndex);
      for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++)
      	speciesNames(ptcl)=species.Name;
    }
    IOSection.WriteVar("SpeciesNames", speciesNames);
    IOSection.WriteVar("Type","Path");
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
    IOSection.FlushFile();
  }
}


///This currently tells you the winding number of all the species. It
///might make sense to fix this so it tells only the winding number of
///a specific species
void WindingNumberClass::Accumulate()
{
  TimesCalled++;

  if (TimesCalled % DumpFreq==0){
    WriteBlock();
  }
  if ((TimesCalled % Freq)!=0){
    return;
  }
  TotalDisp=0.0;
  int numLinks=PathData.Path.NumTimeSlices()-1;
  for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
    for (int slice=0;slice<numLinks;slice++) {
      NumSamples++;
      dVec disp;
      disp=PathData.Path.Velocity(slice,slice+1,ptcl);
      TotalDisp(ptcl) =TotalDisp(ptcl)+ disp;
    }
  }
  PathData.Communicator.Sum(TotalDisp,TempDisp);
  if (PathData.Communicator.MyProc()==0) {  
    for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
      for (int dim=0;dim<NDIM;dim++){
	TotalW2[dim] =TotalW2[dim]+TempDisp(ptcl)[dim]*TempDisp(ptcl)[dim];
      }
    }
  }
  //  for (dim=0;dim<NDIM;dim++){
  //    allPtclDisp[dim]=(allPtclDisp[dim]*allPtclDisp[dim]);
  //  }
  //  TotalW2 += allPtclDispl;
}

void WindingNumberClass::Read(IOSectionClass& in)
{
  ObservableClass::Read(in);
  assert(in.ReadVar("freq",Freq));
  assert(in.ReadVar("dumpFreq",DumpFreq));
}
void WindingNumberClass::WriteBlock()
{
//   double totSum;
//   double totNumSamples;
  
//   double myAvg = ESum/(double)NumSamples; //everybody should have the same number of samples for this to be happy
//   double avg = PathData.Communicator.Sum(myAvg);
//   double vavg =PathData.Communicator.Sum(myVAvg);
//   double savg =PathData.Communicator.Sum(mySAvg);
//   double favg =PathData.Communicator.Sum(myFAvg);
//   avg  = avg/(double)PathData.Path.TotalNumSlices;
//   vavg =vavg/(double)PathData.Path.TotalNumSlices;
//   savg =savg/(double)PathData.Path.TotalNumSlices;
//   favg =favg/(double)PathData.Path.TotalNumSlices;
//   // Only processor 0 writes.
  if (PathData.Communicator.MyProc()==0) {
//     cerr << "myAvg = " << myAvg << endl;
//     cerr << "avg = " << avg << endl;
//     cerr << "Pot avg = " << vavg << endl;
//     cerr << "S avg = " << savg << endl;
//     cerr << "U avg = " <<favg <<endl;
    if (FirstTime) {
      FirstTime = false;
      WriteInfo();
      IOSection.WriteVar("Type","Vector");
       Array<double,2> dummy(1,3);
       for (int dim=0;dim<NDIM;dim++)
	 dummy(0,dim)=TotalW2[dim];
       IOSection.WriteVar ("WindingNumber", dummy);
//       dummy(0)=vavg;
//       IOSection.WriteVar ("PotentialEnergy",dummy);
//       dummy(0)=savg;
//       IOSection.WriteVar ("SpringEnergy",dummy);
//       dummy(0)=favg;
//       IOSection.WriteVar ("DBetaEnergy",dummy);
       IOVar = IOSection.GetVarPtr("WindingNumber");
//       IOVVar= IOSection.GetVarPtr("PotentialEnergy");
//       IOSVar= IOSection.GetVarPtr("SpringEnergy");
//       IOUVar= IOSection.GetVarPtr("DBetaEnergy");
    }
    else {
      Array<double,1> dummy(3);
      for (int dim=0;dim<NDIM;dim++)
	dummy(dim)=TotalW2[dim];
      IOVar->Append(dummy);
//       IOVVar->Append(vavg);
//       IOSVar->Append(savg);
//       IOUVar->Append(favg);
      IOSection.FlushFile();
    }
  }
//   ESum = 0.0;
//   VSum = 0.0;
//   SSum = 0.0;
//   FSum=0.0;
  NumSamples = 0;
}

