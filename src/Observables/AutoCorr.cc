#include "AutoCorr.h"

int NumSlots;
int WaitToFill;
int LastTotal;
dVec v;

////////////////////////////////////////
///Autocorrelation Class           ///
///////////////////////////////////////

void AutoCorrClass::Read(IOSectionClass& in)
{
  ObservableClass::Read(in);
  assert(in.ReadVar("numSlots",NumSlots));
  assert(in.OpenSection("Grid"));
  string gridType;
  double gridStart;
  double gridEnd;
  int numGridPoints = NumSlots;
  WaitToFill = NumSlots;
  LastTotal = WaitToFill;
  assert(in.ReadVar("type",gridType));
  assert(gridType=="Linear");
//  bool readStartGrid=in.ReadVar("start",gridStart);
//  bool readEndGrid=in.ReadVar("end",gridEnd);
//  if (!readStartGrid)
    gridStart=0.0;
//  if (!readEndGrid){
 //   if (PathData.Path.GetPeriodic()[0]){
    //gridEnd=PathData.Path.GetBox()[0];
    gridEnd=NumSlots*Freq;
//    }
//    else {
//      cerr<<"I don't know where you want me to end this grid"<<endl;
 //     assert(1==2);
 //   }
//  }
  //assert(in.ReadVar("NumPoints",numGridPoints));
  grid.Init(gridStart,gridEnd,numGridPoints);
  //TotalCounts=0;
  Histogram.resize(NumSlots);
cerr << "Histogram size is " << Histogram.size() << endl;
  DipoleBin.resize(PathData.Path.numMol*(PathData.Path.NumTimeSlices()-1),NumSlots);
cerr << "DipoleBin size: " << DipoleBin.size() << endl;
  cerr << "slices " << PathData.NumTimeSlices()-1 << " and molecules " << PathData.Path.numMol << endl;
cerr << "NumSlots is " << NumSlots << endl;
  in.CloseSection();
}



void AutoCorrClass::WriteInfo()
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
  IOSection.WriteVar("Type","Autocorrelation");
}



void AutoCorrClass::WriteBlock()
{
  int countNeg = 0;
  PathClass &Path = PathData.Path;
  double norm=(double)(TotalCounts - LastTotal);
  cerr << "normalizing by " << norm << "; LastTotal is " << LastTotal << endl;
  LastTotal = TotalCounts;

  if (Path.Communicator.MyProc()==0) 
    if (FirstTime) {
      FirstTime=false;
      WriteInfo();
    }
  Array<double,1> gofrArray(Histogram.size());
  for (int i=0; i<grid.NumPoints; i++){
    gofrArray(i) =  Histogram(i)/norm;// / (binVol*norm);
    if (gofrArray(i) < 0.0) {
      countNeg++;
    }
  }
  HistVar.Write(gofrArray);
  cerr << countNeg << " negative entries." << endl;
}


void AutoCorrClass::Print()
{
  for (int i=0; i<(grid.NumPoints-1); i++)
    {
      double r1 = grid(i);
      double r2 = grid(i+1);
      double r = 0.5*(r1+r2);
      double vol = 4.0*M_PI/3 * (r2*r2*r2-r1*r1*r1);
      double gofr = Histogram(i)/(TotalCounts - WaitToFill);// / (vol*TotalCounts);
      fprintf (stderr, "%1.12e %1.12e\n", r, gofr);
    }
}


int AutoCorrClass::MapIndex(int slice, int molecule){
  int numMol = PathData.Path.numMol;
  return (molecule + slice*numMol); 
}

void AutoCorrClass::Advance(int& index,int limit){
  index ++;
  if (index > limit)
    index = 0;
}

dVec AutoCorrClass::MeasureDipole(int slice,int molecule){
  Array <int,1> activeParticles(5);
  for (int a = 0; a < 5; a++)
    activeParticles(a) = molecule + a*PathData.Path.numMol;
  dVec O = PathData.Path(slice,activeParticles(0));
  dVec P1 = PathData.Path(slice,activeParticles(3));
  dVec P2 = PathData.Path(slice,activeParticles(4));
  P1 -= O;
  P2 -= O;
  P1 = PathData.Actions.TIP5PWater.Normalize(P1);
  P2 = PathData.Actions.TIP5PWater.Normalize(P2);
  dVec pvec = PathData.Actions.TIP5PWater.Normalize(PathData.Actions.TIP5PWater.GetBisector(P1,P2));
  return pvec;
}

int AutoCorrClass::Locate(int i, int t, int limit){
  int place = i + t;
  if (place > limit)
    place -= (limit + 1);
  if (place < 0)
    place += (limit + 1);
  return place;
}

double AutoCorrClass::CalcDotProd(dVec v1, dVec v2){
  double total = 0.0;
  for (int i = 0; i < 3; i++){
    total += v1(i)*v2(i);
  }
/*  if (total < 0){
    cerr << "dotproduct returning " << total << endl;
    cerr << "from v1 " << v1 << endl;
    cerr << "and v2 " << v2 << endl;
  }*/
  return total;
}

double AutoCorrClass::CalcAutoCorr(int index, int t, int limit){
  int count = 0;
  int count2 = 0;
  double k = 0.0;
//cerr << "calculating autocorrelation between " << index << " and " << Locate(index,t,limit) << endl;;
  double norm = (PathData.NumTimeSlices()-1)*PathData.Path.numMol;
  for (int slice=0;slice<PathData.NumTimeSlices()-1;slice++) {
    for (int mol=0;mol<PathData.Path.numMol;mol++){
      dVec mu1 = DipoleBin(MapIndex(slice,mol),index);
      dVec mu2 = DipoleBin(MapIndex(slice,mol),Locate(index,t,limit));
//cerr << "     slice " << slice << " and mol " << mol << "; mu1 is " << mu1 << " and mu2 is " << mu2;
      double dot = CalcDotProd(mu1,mu2);
      if (dot == 1)
        count++;
      else if (dot < -0.2)
        count2++;
      k += dot;
    }
  }
//cerr << ": k " << k << " /norm " << norm <<  " is " << k/norm << endl;
  //cerr << "I counted " << count << " unchanged particles out of " << norm <<  endl;
  //cerr << "I counted " << count2 << " particles with dot < -0.2" << endl;
  return k/norm;
}

dVec AutoCorrClass::Rotate(dVec coord, double theta){
  double x0 = coord(0);
  double y0 = coord(1);
  double z0 = coord(2);
  double c = cos(theta);
  double s = sin(theta);
  double x = x0*c - y0*s;
  double y = y0*c + x0*s;
  dVec newcoord;
  newcoord(0) = x;
  newcoord(1) = y;
  newcoord(2) = 0;
  return newcoord;
}

// THIS FUNCTION IS ONLY FOR TESTING!!!
/* SET UP TO GENERATE RANDOM ORIENTATIONS -- UNCORRELATED
void AutoCorrClass::Accumulate()
{
  TimesCalled++;

  // measure and catalog dipole moments at now, then calculate autocorrelation 
  if ((TimesCalled % Freq)==0){
    TotalCounts++;
    double theta = M_PI/(NumSlots);
    theta = (0.5-PathData.Path.Random.Local())*(2*M_PI/28);
    v = Rotate(v,theta);
    // loop over slices
    for (int slice=0;slice<PathData.NumTimeSlices()-1;slice++) {
      /// loop over molecules 
      for (int mol=0;mol<PathData.Path.numMol;mol++){
        DipoleBin(MapIndex(slice,mol),now) = v;
      }
    }
    // calculate autocorrelations
    for (int t = 0;t < NumSlots; t++){
      double c = CalcAutoCorr(now,-t,NumSlots);
      if (t<=Histogram.size()){
//cerr << "adding " << c << " to histogram " << Histogram(t);
        Histogram(t) += c;
//cerr << " for " << Histogram(t) << endl;
      }
      else
        cerr << "array size error: " << t << " > " << Histogram.size() << endl;
    }
//cerr << "Accumulate: after entry is " << DipoleBin(50,now) << endl;
    Advance(now,NumSlots);
//cerr << "column updated to " << now << ". ";

    // Write to file
    if (((TotalCounts % DumpFreq) == 0) && (TotalCounts > WaitToFill)){
      cerr << TimesCalled << ", " << TotalCounts << ": Writing Autocorrelation" << endl;
      WriteBlock();
      Histogram = 0;
      //TotalCounts = 0;
    }
  }

  // do nothing this time
  else{
    return;
  }

}*/

// ORIGINAL
//
void AutoCorrClass::Accumulate()
{

  // measure and catalog dipole moments at now, then calculate autocorrelation 

    TotalCounts++;
    // loop over slices
    for (int slice=0;slice<PathData.NumTimeSlices()-1;slice++) {
      /// loop over molecules 
      for (int mol=0;mol<PathData.Path.numMol;mol++){
        DipoleBin(MapIndex(slice,mol),now) = MeasureDipole(slice,mol);
      }
    }
    // calculate autocorrelations
    if (TotalCounts > WaitToFill){
      for (int t = 0;t < NumSlots; t++){
        double c = CalcAutoCorr(now,-t,NumSlots-1);
        if (t<Histogram.size()){
          Histogram(t) += c;
//cerr << TotalCounts << ": adding " << c << " to bin " << t << endl;
        }
        else
          cerr << "array size error: " << t << " > " << Histogram.size() << endl;
      }
//cerr << "Accumulate: after entry is " << DipoleBin(50,now) << endl;
//cerr << "column updated to " << now << ". ";

      // Write to file
      if ((TotalCounts % DumpFreq) == 0){
        cerr << TimesCalled << ", " << TotalCounts << ": Writing Autocorrelation" << endl;
        WriteBlock();
        Histogram = 0;
        //TotalCounts = 0;
cerr << "cleared it out.  now is " << now << endl;
      }
    }
    Advance(now,NumSlots-1);



}//

void AutoCorrClass::Initialize()
{
  TotalCounts = 0;
  TimesCalled=0;
  now = 0;
  Histogram=0;
  v(0) = 1;
  v(1) = 0;
  v(2) = 0;
}
