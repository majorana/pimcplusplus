#include "DistanceTablePBCClass.h"

void DistanceTablePBCClass::UpdateAll(int timeSlice)
{
  int index=0;

  dVec disp;
  double dist;
  int imageNum;

  ///Loop over all the particles that have NDIM dimensions
  for (int ptcl1=0; ptcl1<FirstLowerDimPtcl; ptcl1++) {
    for (int ptcl2=0; ptcl2 <=ptcl1; ptcl2++)
      {
	Displacement(timeSlice,ptcl1,ptcl2, disp, dist, imageNum);
	DispTable.Set(timeSlice,index,disp);
	ImageNumTable.Set(timeSlice,index, imageNum);
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
	dVecInt imageMask = LowerDimImageMasks(species1,species2);
	dVec vectorMask = LowerDimVectorMasks(species1, species2);
	for (int ptcl2=Path.Species(species2).FirstPtcl;
	     ptcl2 <= Path.Species(species2).LastPtcl; ptcl2++) {
	  Displacement(timeSlice, ptcl1, ptcl2, disp, dist, imageNum,
		       vectorMask, imageMask);
	  DispTable.Set(timeSlice,index, disp);
	  ImageNumTable.Set(timeSlice,index,imageNum);
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

void DistanceTablePBCClass::UpdateAll()
{
  for (int timeSlice=0;timeSlice<Path.NumTimeSlices();timeSlice++){
    UpdateAll(timeSlice);
  }
}


void DistanceTablePBCClass::Update(int timeSlice,
				   const Array<int,1> &ptclArray)
{
  Path.DoPtcl=true;
  dVec disp;
  double dist;
  int imageNum;
  int index;
  double sign;

  for (int ptcl1Index=0;ptcl1Index<ptclArray.size();ptcl1Index++){
    int ptcl1=ptclArray(ptcl1Index);
    Path.DoPtcl(ptcl1)=false;
    int species1=Path.ParticleSpeciesNum(ptcl1);
    if (Path.Species(species1).NumDim<NDIM){
      for (int species2=0;species2<Path.NumSpecies();species2++){
	dVecInt imageMask = LowerDimImageMasks(species1,species2);
	dVec vectorMask = LowerDimVectorMasks(species1, species2);
	for (int ptcl2=Path.Species(species2).FirstPtcl; 
	     ptcl2 <=  Path.Species(species2).LastPtcl; ptcl2++)
	  if (Path.DoPtcl(ptcl2)){
	    Displacement(timeSlice,ptcl1,ptcl2,disp,dist,
			 imageNum,vectorMask,imageMask);
	    ArrayIndex(ptcl1,ptcl2,index,sign);
	    disp=disp*sign;
	    if (sign<0){
	      imageNum=NegateImageNum(imageNum);
	    }
	    DispTable.Set(timeSlice,index,disp);
	    DistTable.Set(timeSlice,index,dist);
	    ImageNumTable.Set(timeSlice,index,imageNum);
	  }
      }
    }
    else {

      for (int ptcl2=0;ptcl2<FirstLowerDimPtcl;ptcl2++){
	if (Path.DoPtcl(ptcl2)){
	  Displacement(timeSlice,ptcl1,ptcl2,disp,dist,imageNum);
	  ArrayIndex(ptcl1,ptcl2,index,sign);
	  disp=disp*sign;
	  if (sign<0){
	    imageNum=NegateImageNum(imageNum);
	  }
	  DispTable.Set(timeSlice,index,disp);
	  DistTable.Set(timeSlice,index,dist);
	  ImageNumTable.Set(timeSlice,index,imageNum);
	}
      }
      for (int species2=FirstLowerDimSpecies;species2<Path.NumSpecies();
	   species2++){
	dVecInt imageMask = LowerDimImageMasks(species1,species2);
	dVec vectorMask = LowerDimVectorMasks(species1, species2);
	for (int ptcl2=Path.Species(species2).FirstPtcl; 
	     ptcl2 <=  Path.Species(species2).LastPtcl; ptcl2++)
	  if (Path.DoPtcl(ptcl2)){
	    Displacement(timeSlice,ptcl1,ptcl2,disp,dist,
			 imageNum,vectorMask,imageMask);
	    ArrayIndex(ptcl1,ptcl2,index,sign);
	    disp=disp*sign;
	    if (sign<0){
	      imageNum=NegateImageNum(imageNum);
	    }
	    DispTable.Set(timeSlice,index,disp);
	    DistTable.Set(timeSlice,index,dist);
	    ImageNumTable.Set(timeSlice,index,imageNum);
	  }
      }
    }
  }
}   

