#include "ObservableCorrelation.h"

////////////////////////////////////////
///Pair Correlation Class           ///
///////////////////////////////////////

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





////////////////////////////////////////
///N(r) Correlation Class            ///
///////////////////////////////////////



///Initializes the class. The grid should always be overwritten with
///some other size. 
void nofrClass::Initialize()
{
  grid.Init (0.0, 12.0, 100);
  TotalCounts = 0;
  Histogram.resize(100);
  Histogram = 0.0;
}


///Reads in this classes input.
void nofrClass::Read(IOSectionClass& in)
{
  
  ObservableClass::Read(in);
  assert(in.ReadVar("freq",Freq));
  assert(in.ReadVar("dumpFreq",DumpFreq));

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

///Writes the data relevant for this classes output including its
///name, axis to plot, etc.
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
  IOSection.WriteVar("ylabel", "n(r)");
  IOSection.WriteVar("Type","CorrelationFunction");
}
  

///Writes a block of histogram data.  Does not compensate for volume
///effects or for importance sampling. This has to be done after the
///fact currently
void nofrClass::WriteBlock()
{
  Array<double,1> HistSum(Histogram.size());
  double norm=0.0;
  norm = TotalCounts/PathData.Path.GetVol();
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
	//	double binVol = 4.0*M_PI/3 * (r2*r2*r2-r1*r1*r1);
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



///Bins the distance between the head and the tail of the open
///loop. This does not compensaite for volume effects or importance sampling
void nofrClass::Accumulate()
{
  TimesCalled++;
  if (TimesCalled % DumpFreq==0){
    WriteBlock();
  }
  if ((TimesCalled % Freq)!=0){
    return;
  }

  dVec disp;
  double dist;
  int openLink=(int)(PathData.Path.OpenLink);
  int openPtcl=(int)(PathData.Path.OpenPtcl);
  PathData.Path.DistDisp(openLink,openPtcl,PathData.Path.NumParticles(),
			 dist,disp); //This is distance between head and tail!

  if (dist<grid.End){
    int index=grid.ReverseMap(dist);
    //    Histogram(index)=Histogram(index)+(0.5)/(dist*dist)+(0.9*exp(-dist*dist)+0.1);
    Histogram(index)=Histogram(index)+1.0;
  }  
  TotalCounts++;  
  return; 
}
