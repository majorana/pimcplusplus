#include "DistanceToHead.h"


// Fix to include final link between link M and 0
void HeadLocClass::Accumulate()
{
  TimesCalled++;
  if (TimesCalled % DumpFreq==0)
    WriteBlock();

  if ((TimesCalled % Freq)!=0){
    return;
  }
  
  int closestHeadLoc=0;
  int closestTailLoc=0;
  ///Multiplication by 5 just to make sure closestHead and closestTail
  ///are larger then you could ever have the distance to the head or
  ///tail from any fixed location
  double closestHead=5*dot(PathData.Path.GetBox(),PathData.Path.GetBox());
  double closestTail=5*dot(PathData.Path.GetBox(),PathData.Path.GetBox());
  dVec headDisp,tailDisp;
  double headDist2,tailDist2; 
  int openSlice=PathData.Path.OpenLink;
  int headPtcl=PathData.Path.OpenPtcl;
  ///In non-open loops mode this particle wouldn't exist.  In open
  ///loop mood it stores the location of the tail.
  int tailPtcl=PathData.Path.NumTimeSlices(); 				   
  dVec headLoc=PathData.Path(openSlice,headPtcl);
  dVec tailLoc=PathData.Path(openSlice,tailPtcl);
  for (int counter=0;counter<FixedLoc.size();counter++){
    headDisp=FixedLoc(counter)-headLoc;
    tailDisp=FixedLoc(counter)-tailLoc;
    PathData.Path.PutInBox(headDisp);
    PathData.Path.PutInBox(tailDisp);
    headDist2=dot(headDisp,headDisp);
    tailDist2=dot(tailDisp,tailDisp);
    if (headDist2<closestHead){
      closestHead=headDist2;
      closestHeadLoc=counter;
    }
    if (tailDist2<closestTail){
      closestTail=tailDist2;
      closestTailLoc=counter;
    }
  }
  HeadLoc(closestHeadLoc)++;
  TailLoc(closestTailLoc)++;
  NumSamples++;
}

void HeadLocClass::WriteBlock()
{
  double norm = 1.0/((double)NumSamples);
  Array<int,1> headLocSum(HeadLoc.size());
  Array<int,1> tailLocSum(TailLoc.size());
  PathData.Path.Communicator.Sum(HeadLoc,headLocSum);
  PathData.Path.Communicator.Sum(TailLoc,tailLocSum);
  Array<double,1> headLocDouble(HeadLoc.size());
  Array<double,1> tailLocDouble(TailLoc.size());
  for (int counter=0;counter<headLocSum.size();counter++){
    headLocDouble(counter)=headLocSum(counter)*norm;
    tailLocDouble(counter)=tailLocSum(counter)*norm;
  }
  HeadLocVar.Write(headLocDouble);
  TailLocVar.Write(tailLocDouble);
  //  NumSamples = 0;
}

void HeadLocClass::Read(IOSectionClass &in)
{  
  int numFixedPoints;
  ObservableClass::Read(in);
  assert(in.ReadVar("freq",Freq));
  assert(in.ReadVar("dumpFreq",DumpFreq));
  assert(in.ReadVar("NumFixedPoints",numFixedPoints));
  HeadLoc.resize(numFixedPoints);
  TailLoc.resize(numFixedPoints);
  FixedLoc.resize(numFixedPoints);
  Array<double,2> positions;
  assert(in.ReadVar("LocationsToCompare",positions));
  ///Verify you used the right number of points to compare against
  assert(positions.extent(0)==HeadLoc.size());
  assert(positions.extent(1)==NDIM);
  dVec pos;
  for (int loc=0;loc<=FixedLoc.size(); loc++){
    for (int dim=0; dim<NDIM; dim++)
      pos(dim) = positions(loc,dim);
    FixedLoc(loc) = pos;
  }      
  
  if (PathData.Path.Communicator.MyProc()==0){
    WriteInfo();
    IOSection.WriteVar("Type","Scalar");
  }
}



void HeadLocClass::WriteInfo()
{


}
