#include "ObservableClass.h"

void PairCorrelation::Accumulate()
{
  Array<bool,1> doPtcl2(PathData.NumParticles());

  for (int slice=0;slice<PathData.NumTimeSlices();slice++){
    for (int ptcl1=0;ptcl1<PathData.NumParticles();ptcl1++){
      for (int ptcl2=ptcl1+1;ptcl2<PathData.NumParticles();ptcl2++){
	dVec r1=PathData(slice,ptcl1);
	dVec r2=PathData(slice,ptcl2);
	dVec diff=r1-r2;
	dVec disp;
	//	PathData.DistanceTable->UpdateAll();
	double distDummy;
	PathData.DistanceTable->DistDisp(slice,ptcl1,ptcl2,distDummy,disp);
	cout<<"The two things are: "<<disp<<diff<<endl;
	double dist=sqrt(dot(diff,diff));
	TotalCounts++;
	if (dist<grid.End){
	  int index=grid.ReverseMap(dist);
	  Histogram(index)++;
	}
      }
    }
  }
}


	    //  int NumPtcl1 = PathData.SpeciesArray(Species1).NumParticles();
	    //  int NumPtcl2 = PathData.SpeciesArray(Species2).NumParticles();
//   bool DifferentSpecies = (Species1 != Species2);
  
  
//   if (DifferentSpecies)
//     for (int Slice=0; Slice<PathData.NumTimeSlices(); Slice++) {
//       for (int Ptcl1=0; Ptcl1<NumPtcl1; Ptcl1++) {
// 	dVec r1 =PathData.SpeciesArray(Species1).Path(Ptcl1, Slice);
// 	for (int Ptcl2=0; Ptcl2<NumPtcl2; Ptcl2++) {
// 	  dVec r2=PathData.SpeciesArray(Species2).Path(Ptcl2,Slice);
// 	  dVec diff = r1-r2;
// 	  double dist = sqrt(dot(diff,diff));

// 	  if (dist < grid.End)
// 	    {
// 	      int index = grid.ReverseMap(dist);
// 	      Histogram(index)++;
// 	    }
// 	  TotalCounts++;
// 	}
//       }
//     }
//   else
//     for (int Slice=0; Slice<PathData.NumTimeSlices(); Slice++) {
//       for (int Ptcl1=0; Ptcl1<NumPtcl1; Ptcl1++) {
// 	dVec r1 = PathData.SpeciesArray(Species1).Path(Ptcl1,Slice);
// 	for (int Ptcl2=Ptcl1+1; Ptcl2<NumPtcl2; Ptcl2++)  {
// 	  dVec r2=PathData.SpeciesArray(Species2).Path(Ptcl2,Slice);
// 	  dVec diff = r1-r2;
// 	  double dist = sqrt(dot(diff,diff));
	  
// 	  if (dist < grid.End)
// 	    {
// 	      int index = grid.ReverseMap(dist);
// 	      Histogram(index)++;
// 	    }
// 	  TotalCounts++;
// 	}
//       }
//     }
//}



void PairCorrelation::Initialize()
{
  grid.Init (0.0, 5.0, 50);
  TotalCounts = 0;
  Histogram.resize(49);
  Histogram = 0;
}


void PairCorrelation::Print()
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
