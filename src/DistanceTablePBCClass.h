#ifndef DISTANCE_TABLE_PBC_CLASS_H
#define DISTANCE_TABLE_PBC_CLASS_H


class DistanceTablePBCClass : public DistanceTableClass
{
private:
  inline void Displacement (int timeslice, int ptcl1, int ptcl2,
			    dVec &disp, int imageNum);
public:
  void Update (int timeSlice, int ptcl);
  void UpdateAll();
  void UpdateAll(int timeSlice);
  /// Constructor
  DistanceTablePBCClass (PathClass &myPath) : DistanceTable(myPath)
  { /* Currently DistanceTable constructor does everything */ }
};


/// Determine the particle displacement
inline void DistanceTableClassPBC::Displacement(int timeSlice,
						int ptcl1,int ptcl2, 
						dVec &disp,
						dVecInt &imageNum)
{
  dVecInt image;
  disp=Path(timeSlice,ptcl2)-Path(timeSlice,ptcl1);
  // Determine the image number of ptcl2 w.r.t ptcl1
  for (int i=0; i<NDIM; i++) {
    image[i] = -(disp[i]<-0.5*Path.Box[i]) + (disp[i]>0.5*Path.Box[i]);
    disp[i] += image[i]*Path.Box[i];
  }
  imageNum = ImageNum(image);
  disp = disp + ImageVectors(imageNum);
}




void DistanceTablePBClass::UpdateAll(int timeSlice)
{
  int index=0;
  for (int ptcl1=0;ptcl1<Path.NumParticles;ptcl1++){
    for (int ptcl2=0;ptcl2<=ptcl1;ptcl2++){
      Displacement(timeSlice,ptcl1,ptcl2, DisplaceTable(timeSlice,index),
		   ImageNumTable(timeSlice, index));
      index++;
    }
  }
  /// Now zero out components of displacement due to inactive
  /// dimensions of lower dimensional particles, eg. a cylinder in
  /// 3-space.
  for (int specCntr1=0; specCntr1<Path.NumSpecies; specCntr1++) {
    SpeciesClass &species1=Path.Species(speciesCntr1);
    if (species1.NumDim<NDIM){
      int NotMask = 3;
      for (int dim=0; dim<NDIM; dim++) {
	if (!species1.DimensionActive(dim)){
	  int Mask = ~NotMask;
	  for (int ptcl=0;ptcl<Path.NumParticles();ptcl++){
	    int index=ArrayIndex(ptcl1,ptcl2);
	    // This zeros out the image number corresponding to the
	    // dimension dim;
	    ImageNumTable(timeSlice,index) &= Mask;
	    DispTable(timeSlice,index)[dim] = 0.0;	
	  }
	}
	NotMask <<=2;
      }
    }
  }
  for (int ptcl1=0;ptcl1<Path.NumParticles;ptcl1++){
    for (int ptcl2=0;ptcl2<=ptcl1;ptcl2++){
      int index = ArrayIndex(ptcl1,ptcl2);
      dVec& disp= DispTable(timeSlice,index);
      DistanceTable(timeSlice,index)=sqrt(dot(disp,disp));
      index++;
    }
  }
}

void DistanceTablePBClass::UpdateAll()
{
  for (int timeSlice=0;timeSlice<Path.NumTimeSlices();timeSlice++){
    updateAll(timeSlice);
  }
}


void DistanceTablePBCClass::Update(int timeSlice,int ptcl1)
{
  SpeciesClass &Species1=Path.Species(Path.ParticleSpeciesNum(ptcl1));
  
  for (int ptcl2; ptcl2<Path.NumParticles; ptcl2++){
    int index;
    double sign;
    ArrayIndex(ptcl1, ptcl2, index, sign);
    Displacement(timeSlice,ptcl1,ptcl2, DispTable(timeSlice,index),
		 ImageNumTable(timeSlice,index));
    DispTable(timeSlice,index) *=sign;
    if (sign<0.0)
      ImageNumTable(timeSlice,index) = 
	NegateImageNum(ImageNumTable(timeSlice,index); 
  }
  for (int speciesNum=0; speciesNum<Path.NumSpecies; speciesNum++) {
    SpeciesClass &Species = Path.Species(speciesNum);
    if (Species.NumDim < NDIM || Species1.NumDim<NDIM) {
      for (int dim=0; dim<NDIM; dim++) 
	// Zero out inactive dimensions
	if (!Species.DimensionActive(dim) || 
	    !Species1.DimensionActive(dim))
	  for (int ptcl2=Species.StartPtcl; ptcl2<Species.EndPtcl; ptcl2++) {
	    int index = ArrayIndex(ptcl1,ptcl2);
	    DispTable(timeSlice,index)[dim] = 0.0;
	  }
    } 
  }
  for (int ptcl2=0; ptcl2<Path.NumParticles; ptcl2++) {
    int index = ArrayIndex(ptcl1, ptcl2);
    dVec &disp = DispTable(TimeSlice,Index);
    DistTable(timeSlice,index) = sqrt(dot(disp,disp));
  }
}



#endif
