#include "VacancyLocation.h"


// Fix to include final link between link M and 0
void VacancyLocClass::Accumulate()
{
  TimesCalled++;
  if (TimesCalled % DumpFreq==0)
    WriteBlock();

  if ((TimesCalled % Freq)!=0){
    return;
  }
  dVec displaceAmount;
  double distanceAmount;
  for (int slice=0;slice<PathData.NumTimeSlices();slice++){
    TempLoc=0;
    for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
      double closestAmount=5*dot(PathData.Path.GetBox(),PathData.Path.GetBox());
      int closestLoc=0;
      for (int counter=0;counter<FixedLoc.size();counter++){
	displaceAmount=PathData.Path(slice,ptcl)-FixedLoc(counter);
	PathData.Path.PutInBox(displaceAmount);
	distanceAmount=dot(displaceAmount,displaceAmount);
	if (distanceAmount<closestAmount){
	  closestAmount=distanceAmount;
	  closestLoc=counter;
	}
      }
      Loc(closestLoc)++;
      TempLoc(closestLoc)++;
      R2Dist+=closestAmount;
    }
    for (int counter=0;counter<TempLoc.size();counter++)
      for (int counter2=0;counter2<TempLoc.size();counter2++)
	if (TempLoc(counter)==0 && TempLoc(counter2)==0){
	  dVec disp=FixedLoc(counter)-FixedLoc(counter2);
	  PathData.Path.PutInBox(disp);
	  double dist=sqrt(dot(disp,disp));
	  if (dist<Grid.End){
	    int index=Grid.ReverseMap(dist);
	    Histogram(index)++;
	  }
	}
  }
  NumSamples++;
}

void VacancyLocClass::WriteBlock()
{
  double norm = 1.0/((double)NumSamples);
  double R2Avg=R2Dist*norm;
  //  cerr<<"About to write the double"<<endl;
  //  R2Var.Write(R2Avg);
  //  cerr<<"Written teh double"<<endl;

  Array<int,1> locSum(Loc.size());
  Array<double,1> locWrite(Loc.size());
  locSum=0;
  PathData.Path.Communicator.Sum(Loc,locSum);
  for (int counter=0;counter<locSum.size();counter++){
    locWrite(counter)=locSum(counter)*norm;
  }
  VacancyLocVar.Write(locWrite);
  Array<double,1> histogramWrite(Histogram.size());
  for (int counter=0;counter<histogramWrite.size();counter++)
    histogramWrite(counter)=Histogram(counter)*norm;
  HistogramVar.Write(histogramWrite);
  HistogramVar.Flush();
  Loc=0;
  TempLoc=0;
  NumSamples = 0;
  R2Dist=0.0;
  //  cerr<<"I'm done with that"<<endl;
}

void VacancyLocClass::Read(IOSectionClass &in)
{  


  int numFixedPoints;
  ObservableClass::Read(in);
  assert(in.ReadVar("freq",Freq));
  assert(in.ReadVar("dumpFreq",DumpFreq));
  assert(in.ReadVar("NumFixedPoints",numFixedPoints));
  Loc.resize(numFixedPoints);
  Loc=0;
  TempLoc.resize(numFixedPoints);
  TempLoc=0;
  FixedLoc.resize(numFixedPoints);
  Array<double,2> positions;
  assert(in.ReadVar("LocationsToCompare",positions));

  ///Verify you used the right number of points to compare against
  assert(positions.extent(0)==Loc.size());
  assert(positions.extent(1)==NDIM);
  dVec pos;

  for (int loc=0;loc<FixedLoc.size(); loc++){
    for (int dim=0; dim<NDIM; dim++)
      pos(dim) = positions(loc,dim);
    FixedLoc(loc) = pos;
  }      
  

  if (PathData.Path.Communicator.MyProc()==0){
    WriteInfo();
    IOSection.WriteVar("Type","Scalar");
  }
  Array<int,1> toDivide(Histogram.size());
  for (int counter=0;counter<FixedLoc.size();counter++){
    for (int counter2=0;counter<FixedLoc.size();counter2++){
      disp=PathData.Path(slice,ptcl)-FixedLoc(counter);
      PathData.Path.PutInBox(disp);
      double dist=sqrt(dot(disp,disp));
      if (dist<Grid.End){
	int index=Grid.ReverseMap(dist);
	toDivide(index)++;
      }
    }
    IOSection.WriteVar("Multiplicity",toDivide);
    
    
  
}



void VacancyLocClass::WriteInfo()
{


}
