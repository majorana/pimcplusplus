#include "DistanceTableFreeClass.h"



void DistanceTableFreeClass::DistDispTest(int timeSlice, int ptcl1, int ptcl2, 
		       double &dist, dVec &disp)
{
  dVec p1=Path(timeSlice,ptcl1);
  dVec p2=Path(timeSlice,ptcl2);
  
    


  disp=p2-p1;
  dist=sqrt(dot(disp,disp));
  
}

void DistanceTableFreeClass::DistDispTest(int timeSliceA, int timeSliceB,
		       int ptcl1, int ptcl2, double &distA, double &distB,
		       dVec &dispA, dVec &dispB)
{
  dVec p1A=Path(timeSliceA,ptcl1);
  dVec p1B=Path(timeSliceB,ptcl1);
  dVec p2A=Path(timeSliceA,ptcl2);
  dVec p2B=Path(timeSliceB,ptcl2); 
    

  
  dispA=p2A-p1A;
  distA=sqrt(dot(dispA,dispA));
  dispB=p2B-p1B;
  distB=sqrt(dot(dispB,dispB));
  
  
}


void DistanceTableFreeClass::UpdateAll(int timeSlice)
{
  int index=0;

  
  dVec disp;
  double dist;


  ///Loop over all the particles that have NDIM dimensions
  for (int ptcl1=0; ptcl1<FirstLowerDimPtcl; ptcl1++) {
    for (int ptcl2=0; ptcl2 <=ptcl1; ptcl2++)
      {
	Displacement(timeSlice,ptcl1,ptcl2, disp, dist);
	DispTable.Set(timeSlice,index,disp);
        DistTable.Set(timeSlice,index,dist);
	index++;
      }
  }
  ///Loop over all the lower dimensional particles
  for (int species1=FirstLowerDimSpecies;
       species1<Path.NumSpecies(); species1++) {
    for (int ptcl1=Path.Species(species1).FirstPtcl;
	 ptcl1 <= Path.Species(species1).LastPtcl; ptcl1++) {
      for (int species2=0;
	   species2<=species1; species2++){
	dVec vectorMask = LowerDimVectorMasks(species1, species2);
	for (int ptcl2=Path.Species(species2).FirstPtcl;
	     ptcl2 <= Path.Species(species2).LastPtcl; ptcl2++) {
	  Displacement(timeSlice, ptcl1, ptcl2, disp, dist, 
		       vectorMask);
	  DispTable.Set(timeSlice,index, disp);
	  DistTable.Set(timeSlice,index,dist);
	  double dist = sqrt(dot(disp,disp));
	  DistTable.Set(timeSlice,index,dist);
	  index++;
	}
      }
    }
  }
}
  


//   /// Now zero out components of displacement due to inactive
//   /// dimensions of lower dimensional particles, eg. a cylinder in
//   /// 3-space.
//   for (int specCntr1=0; specCntr1<Path.NumSpecies; specCntr1++) {
//     SpeciesClass &species1=Path.Species(speciesCntr1);
//     if (species1.NumDim<NDIM){
//       int NotMask = 3;
//       for (int dim=0; dim<NDIM; dim++) {
// 	if (!species1.DimensionActive(dim)){
// 	  int Mask = ~NotMask;
// 	  for (int ptcl=0;ptcl<Path.NumParticles();ptcl++){
// 	    int index=ArrayIndex(ptcl1,ptcl2);
// 	    // This zeros out the image number corresponding to the
// 	    // dimension dim;
// 	    ImageNumTable(timeSlice,index) &= Mask;
// 	    DispTable(timeSlice,index)[dim] = 0.0;	
// 	  }
// 	}
// 	NotMask <<=2;
//       }
//     }
//   }
//   for (int ptcl1=0;ptcl1<Path.NumParticles;ptcl1++){
//     for (int ptcl2=0;ptcl2<=ptcl1;ptcl2++){
//       dVec& disp= DispTable(timeSlice,index);
//       DistanceTable(timeSlice,index)=sqrt(dot(disp,disp));
//       index++;
//     }
//   }
// }

void DistanceTableFreeClass::UpdateAll()
{
  for (int timeSlice=0;timeSlice<Path.NumTimeSlices();timeSlice++){
    UpdateAll(timeSlice);
  }
}


void DistanceTableFreeClass::Update(int timeSlice, 
				    const Array<int,1> &ptclArray)
{
  Path.DoPtcl=true;
  dVec disp;
  double dist;

  int index;
  double sign;

  for (int ptcl1Index=0;ptcl1Index<ptclArray.size();ptcl1Index++){
    int ptcl1=ptclArray(ptcl1Index);
    Path.DoPtcl(ptcl1)=false;
    int species1=Path.ParticleSpeciesNum(ptcl1);
    if (Path.Species(species1).NumDim<NDIM){
      for (int species2=0;species2<Path.NumSpecies();species2++){

	dVec vectorMask = LowerDimVectorMasks(species1, species2);
	for (int ptcl2=Path.Species(species2).FirstPtcl; 
	     ptcl2 <=  Path.Species(species2).LastPtcl; ptcl2++)
	  if (Path.DoPtcl(ptcl2)){
	    Displacement(timeSlice,ptcl1,ptcl2,disp,dist,
			 vectorMask);
	    ArrayIndex(ptcl1,ptcl2,index,sign);
	    disp=disp*sign;
	    DispTable.Set(timeSlice,index,disp);
	    DistTable.Set(timeSlice,index,dist);

	  }
      }
    }
    else {

      for (int ptcl2=0;ptcl2<FirstLowerDimPtcl;ptcl2++){
	if (Path.DoPtcl(ptcl2)){
	  Displacement(timeSlice,ptcl1,ptcl2,disp,dist);
	  ArrayIndex(ptcl1,ptcl2,index,sign);
	  disp=disp*sign;

	  DispTable.Set(timeSlice,index,disp);
	  DistTable.Set(timeSlice,index,dist);

	}
      }
      for (int species2=FirstLowerDimSpecies;species2<Path.NumSpecies();
	   species2++){

	dVec vectorMask = LowerDimVectorMasks(species1, species2);
	for (int ptcl2=Path.Species(species2).FirstPtcl; 
	     ptcl2 <=  Path.Species(species2).LastPtcl; ptcl2++)
	  if (Path.DoPtcl(ptcl2)){
	    Displacement(timeSlice,ptcl1,ptcl2,disp,dist,
			 vectorMask);
	    ArrayIndex(ptcl1,ptcl2,index,sign);
	    disp=disp*sign;

	    DispTable.Set(timeSlice,index,disp);
	    DistTable.Set(timeSlice,index,dist);

	  }
      }
    }
  }
}   


dVec DistanceTableFreeClass::Velocity(int timeSliceA, int timeSliceB, int ptcl)
{
  return (Path(timeSliceB,ptcl)-Path(timeSliceA,ptcl));
}



void DistanceTableFreeClass::PutInBox (dVec &r)
{
}
