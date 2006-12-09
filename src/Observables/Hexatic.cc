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

#include "Hexatic.h"

double unit2angle(double x,double y)
{
  double angle=atan(y/x);
  if (x<0)
    angle=angle+M_PI;
  return angle;
}


complex<double>
HexaticClass::OrderParamater(int slice,int ptcl)
{
  complex<double> op=0.0;
  for (int nearPtcl=0;nearPtcl<PathData.Path.NumParticles();
       nearPtcl++){
    double r12dist;
    dVec r12disp;
    PathData.Path.DistDisp(slice,ptcl,nearPtcl,r12dist,r12disp);
    if (r12dist<DistCutoff){
      r12disp=r12disp/r12dist;
      if (dot(r12disp,r12disp)-1.0<0.001)
	cerr<<dot(r12disp,r12disp);
      assert(dot(r12disp,r12disp)-1.0<0.001);
      double theta_12=unit2angle(r12disp(0),r12disp(1));
      op=op+complex<double>(cos(theta_12*q),sin(theta_12*q));
    }
  }
  return op;
}


void
HexaticClass::Accumulate()
{
  
  //  PathData.MoveJoin(PathData.NumTimeSlices()-1);
  for (int slice=0;slice<PathData.Path.NumTimeSlices();slice++){
    for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
      ParticleOrder(ptcl)=OrderParamater(slice,ptcl);
    }
    for (int ptcl1=0;ptcl1<PathData.Path.NumParticles();ptcl1++){
      for (int ptcl2=0;ptcl2<PathData.Path.NumParticles();ptcl2++){
	double dist;
	dVec disp;
	PathData.Path.DistDisp(slice,ptcl1,ptcl2,dist,disp);
	if (dist<grid.End){
	  int index=grid.ReverseMap(dist);
	  Histogram(index)++;
	}
      }
    }
  }
  NumSamples++;
}

void
HexaticClass::WriteBlock()
{
  PathClass &Path = PathData.Path;
  Path.Communicator.Sum(Histogram,HistSum);
  double norm=1.0/((double)NumSamples*Path.TotalNumSlices);
  for (int i=0;i<HistDouble.size();i++)
    HistDouble(i)=(double)HistSum(i)*norm;
  HexaticVar.Write(HistDouble);
  HexaticVar.Flush();
  Histogram=0;
  NumSamples=0;
}

void 
HexaticClass::ReadGrid(IOSectionClass &in)
{
  assert(in.OpenSection("Grid"));
  string gridType;
  double gridStart;
  double gridEnd;
  int numGridPoints;
  assert(in.ReadVar("Type",gridType));
  assert(gridType=="Linear");
  bool readStartGrid=in.ReadVar("start",gridStart);
  bool readEndGrid=in.ReadVar("end",gridEnd);
  if (!readStartGrid)
    gridStart=0.0;
  if (!readEndGrid){
    if (PathData.Path.GetPeriodic()[0]){
      gridEnd=PathData.Path.GetBox()[0];
    }
    else {
      cerr<<"I don't know where you want me to end this grid"<<endl;
      assert(1==2);
    }
  }
  assert(in.ReadVar("NumPoints",numGridPoints));
  grid.Init(gridStart,gridEnd,numGridPoints);
  NumSamples=0;
  Histogram.resize(numGridPoints-1);
  HistSum.resize(Histogram.size());
  HistDouble.resize(Histogram.size());
  Histogram=0;
  in.CloseSection();
}

void
HexaticClass::Read (IOSectionClass &in)
{
  NumSamples=0;
  ParticleOrder.resize(PathData.Path.NumParticles());
  ///It's probably important that the grid is the same grid that is in
  ///the pair correlation function. Not sure how to authenticate this.
  ReadGrid(in);
  ObservableClass::Read(in);
}
