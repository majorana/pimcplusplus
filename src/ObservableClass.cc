#include "ObservableClass.h"

void PairCorrelation::Accumulate()
{
  int NumPtcl1 = PathData->IdenticalParticleArray(Species1).NumParticles;
  int NumPtcl2 = PathData->IdenticalParticleArray(Species2).NumParticles;
  bool DifferentSpecies = (Species1 != Species2);
  
  
  if (DifferentSpecies)
    for (int Slice=0; Slice<PathData->NumTimeSlices; Slice++) {
      for (int Ptcl1=0; Ptcl1<NumPtcl1; Ptcl1++) {
	dVec r1 =PathData->IdenticalParticleArray(Species1).Path(Ptcl1, Slice);
	for (int Ptcl2=0; Ptcl2<NumPtcl2; Ptcl2++) {
	  dVec r2=PathData->IdenticalParticleArray(Species2).Path(Ptcl2,Slice);
	  dVec diff = r1-r2;
	  double dist = sqrt(dot(diff,diff));

	  if (dist < grid.End)
	    {
	      int index = grid.ReverseMap(dist);
	      Histogram(index)++;
	    }
	  TotalCounts++;
	}
      }
    }
  else
    for (int Slice=0; Slice<PathData->NumTimeSlices; Slice++) {
      for (int Ptcl1=0; Ptcl1<NumPtcl1; Ptcl1++) {
	dVec r1 = PathData->IdenticalParticleArray(Species1).Path(Ptcl1,Slice);
	for (int Ptcl2=Ptcl1+1; Ptcl2<NumPtcl2; Ptcl2++)  {
	  dVec r2=PathData->IdenticalParticleArray(Species2).Path(Ptcl2,Slice);
	  dVec diff = r1-r2;
	  double dist = sqrt(dot(diff,diff));
	  
	  if (dist < grid.End)
	    {
	      int index = grid.ReverseMap(dist);
	      Histogram(index)++;
	    }
	  TotalCounts++;
	}
      }
    }
}



void PairCorrelation::Initialize()
{
  grid.Init (0.0, 5.0, 100);
  TotalCounts = 0;
  Histogram.resize(99);
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
