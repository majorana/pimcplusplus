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

#include "AutoCorr.h"
#include "../Moves/MoveUtils.h"

//dVec v;

////////////////////////////////////////
///Autocorrelation Class           ///
///////////////////////////////////////

void AutoCorrClass::Read(IOSectionClass& in)
{
  ObservableClass::Read(in);
  assert(in.ReadVar("Frequency", Frequency));
  assert(in.ReadVar("dumpFrequency", dumpFrequency));
  assert(in.ReadVar("numSlots",NumSlots));
	if(!in.ReadVar("Species",dipoleSpecies))
		dipoleSpecies = "p";
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
    gridEnd=NumSlots*Frequency;
//    }
//    else {
//      cerr<<"I don't know where you want me to end this grid"<<endl;
 //     assert(1==2);
 //   }
//  }
  //assert(in.ReadVar("NumPoints",numGridPoints));
  grid.Init(gridStart,gridEnd,numGridPoints);
  //TotalCounts=0;
  OneOneHistogram.resize(NumSlots);
  OneNetHistogram.resize(NumSlots);
  NetNetHistogram.resize(NumSlots);
	cerr << "OneOneHistogram size is " << OneOneHistogram.size() << endl;
  DipoleBin.resize(PathData.Path.numMol*(PathData.Path.NumTimeSlices()-1),NumSlots);
	cerr << "DipoleBin size: " << DipoleBin.size() << endl;
  NetDipoleBin.resize((PathData.Path.NumTimeSlices()-1),NumSlots);
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
  IOSection.WriteVar("xlabel", "MC timestep");
  IOSection.WriteVar("ylabel", "Autocorrelation");
  IOSection.WriteVar("Type","CorrelationFunction");
  IOSection.WriteVar("Cumulative", false);
}



void AutoCorrClass::WriteBlock()
{
}

void AutoCorrClass::LocalWriteBlock()
{
  PathClass &Path = PathData.Path;
  double norm=(double)(TotalCounts - LastTotal);
	double sliceNorm = PathData.NumTimeSlices()-1;
  cerr << "TotalCounts is " << TotalCounts << "; normalizing by " << norm << " and sliceNorm is " << sliceNorm << "; LastTotal is " << LastTotal << endl;
  LastTotal = TotalCounts;

  if (Path.Communicator.MyProc()==0) 
    if (FirstTime) {
      FirstTime=false;
      WriteInfo();
    }
  Array<double,1> OneOnegofrArray(OneOneHistogram.size());
  for (int i=0; i<grid.NumPoints; i++){
    OneOnegofrArray(i) =  OneOneHistogram(i)/norm;// / (binVol*norm);
  }
  Array<double,1> OneNetgofrArray(OneNetHistogram.size());
  for (int i=0; i<grid.NumPoints; i++){
    OneNetgofrArray(i) =  OneNetHistogram(i)/norm;// / (binVol*norm);
  }
  Array<double,1> NetNetgofrArray(NetNetHistogram.size());
  for (int i=0; i<grid.NumPoints; i++){
    NetNetgofrArray(i) =  NetNetHistogram(i)/norm;// / (binVol*norm);
  }

//cerr << "going to write OneOne " << OneOnegofrArray << endl;
  OneOneHistVar.Write(OneOnegofrArray);
//cerr << "going to write OneNet " << OneNetgofrArray << endl;
  OneNetHistVar.Write(OneNetgofrArray);
//cerr << "going to write NetNet " << NetNetgofrArray << endl;
  NetNetHistVar.Write(NetNetgofrArray);

	TotalNetMu *= 1.0/(norm*sliceNorm);
	TotalNetMuVar.Write(TotalNetMu);
}


void AutoCorrClass::Print()
{
  for (int i=0; i<(grid.NumPoints-1); i++)
    {
      double r1 = grid(i);
      double r2 = grid(i+1);
      double r = 0.5*(r1+r2);
      double vol = 4.0*M_PI/3 * (r2*r2*r2-r1*r1*r1);
      double gofr = OneOneHistogram(i)/(TotalCounts - WaitToFill);// / (vol*TotalCounts);
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

// new one
dVec AutoCorrClass::MeasureDipole(int slice,int molecule){
	vector<int> dipoleIndex(0);
  for (int a = 0; a < PathData.Path.MolMembers(molecule).size(); a++){
    int ptcl = PathData.Path.MolMembers(molecule)(a);
		if(PathData.Path.ParticleSpeciesNum(ptcl) == PathData.Path.SpeciesNum(dipoleSpecies)){
			dipoleIndex.push_back(ptcl);
		}
	}
	assert(dipoleIndex.size() == 2);
	//cerr << "autocorr loaded molecule " << molecule << " and protons " << dipoleIndex << endl;
  dVec O = PathData.Path(slice,molecule);
  dVec P1 = PathData.Path(slice,dipoleIndex[0]);
  dVec P2 = PathData.Path(slice,dipoleIndex[1]);
  P1 -= O;
  P2 -= O;
  //P1 = Normalize(P1);
  //P2 = Normalize(P2);
  dVec pvec = Normalize(GetBisector(P1,P2));
  return pvec;
}

// old one; try to make it a little more robust
//dVec AutoCorrClass::MeasureDipole(int slice,int molecule){
//  Array <int,1> activeParticles(5);
//  for (int a = 0; a < 5; a++)
//    activeParticles(a) = molecule + a*PathData.Path.numMol;
//  dVec O = PathData.Path(slice,activeParticles(0));
//  dVec P1 = PathData.Path(slice,activeParticles(3));
//  dVec P2 = PathData.Path(slice,activeParticles(4));
//  P1 -= O;
//  P2 -= O;
//  P1 = Normalize(P1);
//  P2 = Normalize(P2);
//  dVec pvec = Normalize(GetBisector(P1,P2));
//  return pvec;
//}

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

void
AutoCorrClass::CalcAutoCorr(int index, int t, int limit, double& SingleSingle, double& SingleNet, double& NetNet){
	SingleSingle = 0.0;
	SingleNet = 0.0;
	NetNet = 0.0;
  double k = 0.0;
	//cerr << "calculating autocorrelation between " << index << " and " << Locate(index,t,limit) << endl;;
  double norm = (PathData.NumTimeSlices()-1)*PathData.Path.numMol;
  for (int slice=0;slice<PathData.NumTimeSlices()-1;slice++) {
		dVec muNet1 = NetDipoleBin(slice, index);
		dVec muNet2 = NetDipoleBin(slice, Locate(index,t,limit));
		double NetNetDot = CalcDotProd(muNet1, muNet2);
		NetNet += NetNetDot;
    for (int mol=0;mol<PathData.Path.numMol;mol++){
      dVec mu1 = DipoleBin(MapIndex(slice,mol),index);
      dVec mu2 = DipoleBin(MapIndex(slice,mol),Locate(index,t,limit));
//cerr << "     slice " << slice << " and mol " << mol << "; mu1 is " << mu1 << " and mu2 is " << mu2;
      double dot = CalcDotProd(mu1,mu2);
			double SingleNetDot = CalcDotProd(mu2,muNet2);
      SingleSingle += dot;
			SingleNet += SingleNetDot;
    }
  }
	SingleSingle *= 1.0/norm;
	SingleNet *= 1.0/norm;
	NetNet *= 1.0/(PathData.NumTimeSlices() - 1);
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

void AutoCorrClass::Accumulate()
{
	TimesCalled++;
  // measure and catalog dipole moments at now, then calculate autocorrelation 
  TotalCounts++;
  // loop over slices
	//cerr << "In Autocorr " << TimesCalled << " " << TotalCounts << endl;
  for (int slice=0;slice<PathData.NumTimeSlices()-1;slice++) {
    /// loop over molecules 
		dVec netMu(0.0, 0.0, 0.0);
    for (int mol=0;mol<PathData.Path.numMol;mol++){
			int index = MapIndex(slice,mol); 
      DipoleBin(index,now) = MeasureDipole(slice,mol);
			netMu += DipoleBin(index,now);
    }
		TotalNetMu += sqrt(CalcDotProd(netMu, netMu));
		netMu = Normalize(netMu);
		//cerr << "Adding " << netMu << " at " << slice << ", " << now << endl;
		NetDipoleBin(slice,now) = netMu;
  }
  // calculate autocorrelations
  if (TotalCounts > WaitToFill){
    for (int t=0; t<NumSlots; t++){
			double OneOneCorr, OneNetCorr, NetNetCorr;
      CalcAutoCorr(now,-t,NumSlots-1,OneOneCorr,OneNetCorr,NetNetCorr);
      assert(t<OneOneHistogram.size());
      OneOneHistogram(t) += OneOneCorr;
      OneNetHistogram(t) += OneNetCorr;
      NetNetHistogram(t) += NetNetCorr;
			//cerr << TotalCounts << ": adding " << OneOneCorr << " to bin " << t << endl;
    }
		//cerr << "Accumulate: after entry is " << DipoleBin(50,now) << endl;
		//cerr << "  and Net is " << NetDipoleBin(50,now) << endl;
		//cerr << "  column updated to " << now << ". ";

    // Write to file
		//cerr << "TotalCounts " << TotalCounts << "%" << dumpFrequency << " is " << TotalCounts%dumpFrequency << endl;
    if ((TotalCounts % dumpFrequency) == 0){
      cerr << TimesCalled << ", " << TotalCounts << ": Writing Autocorrelation" << endl;
      LocalWriteBlock();
      OneOneHistogram = 0;
      OneNetHistogram = 0;
      NetNetHistogram = 0;
			TotalNetMu = 0.0;
    }
  }
  Advance(now,NumSlots-1);
}

void AutoCorrClass::Initialize()
{
  TotalCounts = 0;
  TimesCalled=0;
  now = 0;
  OneOneHistogram=0;
  OneNetHistogram=0;
  NetNetHistogram=0;
	TotalNetMu = 0.0;
  //v(0) = 1;
  //v(1) = 0;
  //v(2) = 0;
}
