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
  for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
    for (int slice=0;slice<PathData.NumTimeSlices();slice++){
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
      R2Dist+=closestAmount;
    }
  }

  NumSamples++;
  cerr<<"Current thing to write is "<<R2Dist/NumSamples<<endl;
}

void VacancyLocClass::WriteBlock()
{
  cerr<<"I'm writing my block"<<endl;
  double norm = 1.0/((double)NumSamples);
  double R2Avg=R2Dist*norm;
  cerr<<"About to write the double"<<endl;
  R2Var.Write(R2Avg);
  cerr<<"Written teh double"<<endl;

  Array<int,1> locSum(Loc.size());
  Array<double,1> locWrite(Loc.size());
  locSum=0;
  PathData.Path.Communicator.Sum(Loc,locSum);
  for (int counter=0;counter<locSum.size();counter++){
    locWrite(counter)=locSum(counter)*norm;
  }
  cerr<<"About to wirt eh array"<<endl;
  VacancyLocVar.Write(locWrite);

  cerr<<"Written the array"<<endl;
  Loc=0;
  NumSamples = 0;
  R2Dist=0.0;
  cerr<<"I'm done with that"<<endl;
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
}



void VacancyLocClass::WriteInfo()
{


}
