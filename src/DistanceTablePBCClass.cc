#include "DistanceTablePBCClass.h"

void DistanceTablePBCClass::DistDispTest(int timeSlice, int ptcl1, int ptcl2, 
		       double &dist, dVec &disp)
{
  dVec p1=Path(timeSlice,ptcl1);
  dVec p2=Path(timeSlice,ptcl2);
  
    
  dVec normalDisplace=p2-p1;
  for (int dim=0;dim<NDIM;dim++){
    while (normalDisplace(dim)>Path.Box(dim)/2){
      normalDisplace(dim)=normalDisplace(dim)-Path.Box(dim)/2;
    }
    while (normalDisplace(dim)< -Path.Box(dim)/2){
      normalDisplace(dim)=normalDisplace(dim)+Path.Box(dim)/2;
    }
    
  }
  disp=normalDisplace;
  dist=sqrt(dot(disp,disp));
  
}

void DistanceTablePBCClass::DistDispTest(int timeSliceA, int timeSliceB,
		       int ptcl1, int ptcl2, double &distA, double &distB,
		       dVec &dispA, dVec &dispB)
{
  dVec p1A=Path(timeSliceA,ptcl1);
  dVec p1B=Path(timeSliceB,ptcl1);
  dVec p2A=Path(timeSliceA,ptcl2);
  dVec p2B=Path(timeSliceB,ptcl2); 
    
  dVec normalDisplaceA=p2A-p1A;
  dVec normalDisplaceB=p2B-p1B;
  
  for (int dim=0;dim<NDIM;dim++){
    while (normalDisplaceA(dim)>Path.Box(dim)/2){
      normalDisplaceA(dim)=normalDisplaceA(dim)-Path.Box(dim)/2;
      normalDisplaceB(dim)=normalDisplaceB(dim)-Path.Box(dim)/2;
    }
    while (normalDisplaceA(dim)< -Path.Box(dim)/2){
      normalDisplaceA(dim)=normalDisplaceA(dim)+Path.Box(dim)/2;
      normalDisplaceB(dim)=normalDisplaceB(dim)+Path.Box(dim)/2;
    }
    
  }
  dispA=normalDisplaceA;
  distA=sqrt(dot(dispA,dispA));
  dispB=normalDisplaceB;
  distB=sqrt(dot(dispB,dispB));
  
}


void DistanceTablePBCClass::UpdateAll(int timeSlice)
{
  /// First, force all particles back into the box.
  for (int ptcl=0; ptcl < Path.NumParticles(); ptcl++)
    {
      dVec v = Path(timeSlice,ptcl);
      for (int dim=0; dim<NDIM;dim++)
	{
	  while (v[dim] > (0.5*Path.Box[dim]))
	    v[dim] -= Path.Box[dim];
	  while (v[dim] < -(0.5*Path.Box[dim]))
	    v[dim] += Path.Box[dim];
	}
      Path.SetPos(timeSlice, ptcl, v);
    } 

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
  /// First, force all particles that have moved back into the box.
  for (int ptclIndex=0; ptclIndex<ptclArray.size(); ptclIndex++)
    {
      int ptcl = ptclArray(ptclIndex);
      dVec v = Path(timeSlice,ptcl);
      for (int dim=0; dim<NDIM;dim++)
	{
	  while (v[dim] > (0.5*Path.Box[dim]))
	    v[dim] -= Path.Box[dim];
	  while (v[dim] < -(0.5*Path.Box[dim]))
	    v[dim] += Path.Box[dim];
	}
      Path.SetPos(timeSlice, ptcl, v);
    }
	    

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


dVec DistanceTablePBCClass::Velocity(int timeSliceA, int timeSliceB, int ptcl)
{
  dVec vel = Path(timeSliceB,ptcl) - Path(timeSliceA,ptcl);
   for (int dim=0; dim<NDIM; dim++)
     {
       while (vel[dim] > (0.5*Path.Box[dim]))
  	vel[dim] -= Path.Box[dim];
       while (vel[dim] < -(0.5*Path.Box[dim]))
  	vel[dim] += Path.Box[dim];
     }
  return (vel);
}

void DistanceTablePBCClass::PutInBox (dVec &r)
{
  for (int dim=0; dim<NDIM; dim++)
    {
      while (r[dim] > (0.5*Path.Box[dim]))
	r -= Path.Box[dim];
      while (r[dim] < -(0.5*Path.Box[dim]))
	r += Path.Box[dim];
    }
}
