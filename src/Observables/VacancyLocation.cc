#include "VacancyLocation.h"



void VacancyLocClass::PrintNearbySites()
{

  for (int site=0;site<FixedLoc.size();site++){
    list<int>::iterator neighborIter;
    cerr<<"Site: "<<site<<endl;
    for (neighborIter=Neighbors[site].begin();
	 neighborIter!=Neighbors[site].end();
	 neighborIter++){
      cerr<<(*neighborIter)<<", ";
    }
    cerr<<endl;
  }
}


void VacancyLocClass::TabulateNearbySites()
{
  double minDist=10*dot(PathData.Path.GetBox(),PathData.Path.GetBox());
  for (int siteA=0;siteA<FixedLoc.size();siteA++){
    for (int siteB=0;siteB<FixedLoc.size();siteB++){
      double currDist=
	PathData.Path.MinImageDistance(FixedLoc(siteA),FixedLoc(siteB));
      if (currDist<minDist && abs(currDist)>0.01)
	minDist=currDist;
    }
  }
  double epsilon=0.1;
  for (int siteA=0;siteA<FixedLoc.size();siteA++){
    for (int siteB=0;siteB<FixedLoc.size();siteB++){
      double currDist=
	PathData.Path.MinImageDistance(FixedLoc(siteA),FixedLoc(siteB));
      if (currDist<minDist+epsilon && abs(currDist)>0.01)
	Neighbors[siteA].push_back(siteB);
    }
  }
}

bool VacancyLocClass::NeighborsVacancyFree(int site)
{
  list<int>::iterator neighborIter;
  for (neighborIter=Neighbors[site].begin();
       neighborIter!=Neighbors[site].end();
       neighborIter++){
    if (TempLoc(*neighborIter)==0)
      return false;
  }
  return true;
}

bool VacancyLocClass::NeighborsNotDoublyOccupied(int site)
{
  list<int>::iterator neighborIter;
  for (neighborIter=Neighbors[site].begin();
       neighborIter!=Neighbors[site].end();
       neighborIter++){
    if (TempLoc(*neighborIter)>1)
      return false;
  }
  return true;
}

// Fix to include final link between link M and 0
void VacancyLocClass::Accumulate()
{
  dVec displaceAmount;
  double distanceAmount;
  SiteEmptyAtSomeTimeSlice=0;
  for (int slice=0;slice<PathData.NumTimeSlices()-1;slice++){
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
      //Commented out Mar 1      Loc(closestLoc)++;
      TempLoc(closestLoc)++;
      R2Dist+=closestAmount;
    }
    
    for (int counter=0;counter<TempLoc.size();counter++)
      if (TempLoc(counter)==0 &&  NeighborsNotDoublyOccupied(counter))
	SiteEmptyAtSomeTimeSlice(counter)=1;
    for (int counter=0;counter<TempLoc.size();counter++)
      if (TempLoc(counter)==0 && NeighborsNotDoublyOccupied(counter)){
	Loc(counter)++;
	cerr<<"My temp Loc is "<<TempLoc<<endl;
      }
    //    bool doPair=true;
    //    for (int counter=0;counter<TempLoc.size();counter++)
    //      if (TempLoc(counter)>1){
    //	doPair=false;
    //      }
    //    doPair=false;  //HACK! HACK! HACK!
    //    if (doPair){
    //      for (int counter=0;counter<TempLoc.size();counter++)
    //	for (int counter2=0;counter2<TempLoc.size();counter2++)
    //	  if (TempLoc(counter)==0 && TempLoc(counter2)==0){
    //	    dVec disp=FixedLoc(counter)-FixedLoc(counter2);
    //	    PathData.Path.PutInBox(disp);
    //	    double dist=sqrt(dot(disp,disp));
    //	    if (dist<Grid.End){
    //	      int index=Grid.ReverseMap(dist);
    //	      Histogram(index)++;
    //	    }
    //	  }
    //    }
    //  }






    if (PathData.Path.NumSpecies()>1){ //do vacancy-helium3 correlation
      for (int counter=0;counter<TempLoc.size();counter++)
	if (TempLoc(counter)==0 && NeighborsNotDoublyOccupied(counter)){
	  dVec disp=FixedLoc(counter)-PathData.Path(slice,PathData.Path.Species(1).FirstPtcl);
	  PathData.Path.PutInBox(disp);
	  double dist=sqrt(dot(disp,disp));
	  if (dist<Grid.End){
	    int index=Grid.ReverseMap(dist);
	    Histogram(index)++;
	  }
	}
    

    }
    else { //do vacancy-vacancy correlation
  ///TRying again 
    for (int counter=0;counter<TempLoc.size();counter++)
      for (int counter2=0;counter2<TempLoc.size();counter2++)
	if (TempLoc(counter)==0 && TempLoc(counter2)==0 &&
	    NeighborsNotDoublyOccupied(counter) && 
	    NeighborsNotDoublyOccupied(counter2)){
	  dVec disp=FixedLoc(counter)-FixedLoc(counter2);
	  PathData.Path.PutInBox(disp);
	  double dist=sqrt(dot(disp,disp));
	  if (dist<Grid.End){
	    int index=Grid.ReverseMap(dist);
	    Histogram(index)++;
	  }
	}   
    }

  }
  for (int counter=0;counter<TempLoc.size();counter++)
    if (SiteEmptyAtSomeTimeSlice(counter)==1)
      NumEmptyLatticeSites=NumEmptyLatticeSites+1;
  


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
  NumEmptyLatticeSitesVar.Write(NumEmptyLatticeSites/(double)NumSamples);
  Array<double,1> histogramWrite(Histogram.size());
  for (int counter=0;counter<histogramWrite.size();counter++)
    histogramWrite(counter)=Histogram(counter)*norm;
  HistogramVar.Write(histogramWrite);
  HistogramVar.Flush();
  Loc=0;
  TempLoc=0;
  NumSamples = 0;
  R2Dist=0.0;
  Histogram=0;
  NumEmptyLatticeSites=0.0;
  //  cerr<<"I'm done with that"<<endl;
}

void VacancyLocClass::Read(IOSectionClass &in)
{  


  int numFixedPoints;
  ObservableClass::Read(in);
  assert(in.ReadVar("NumFixedPoints",numFixedPoints));
  Loc.resize(numFixedPoints);
  Loc=0;
  TempLoc.resize(numFixedPoints);
  SiteEmptyAtSomeTimeSlice.resize(numFixedPoints);
  TempLoc=0;
  FixedLoc.resize(numFixedPoints);
  Neighbors.resize(numFixedPoints);
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
    for (int counter2=0;counter2<FixedLoc.size();counter2++){
      dVec disp=FixedLoc(counter)-FixedLoc(counter2);
      PathData.Path.PutInBox(disp);
      double dist=sqrt(dot(disp,disp));
      if (dist<Grid.End){
	int index=Grid.ReverseMap(dist);
	toDivide(index)++;
      }
    }
  }
  IOSection.WriteVar("Multiplicity",toDivide);
  TabulateNearbySites();
  PrintNearbySites();
  
  
}



void VacancyLocClass::WriteInfo()
{


}
