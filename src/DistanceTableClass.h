#ifndef DISTANCE_TABLE_CLASS_H
#define DISTANCE_TABLE_CLASS_H

#include "PathClass.h"

///Particle 1 is always in the middle of the box
///Particle 1 > Particle 2
///Particle 1 is in the middle of the box
///If Particle 1 < Particle 2 then you have to swap the sign of the
///displacement
///The imageNum is the box number of particle 2 if particle 1 is in the 
///middle of the box. 
///The displacement is always defined as particle 2 - particle 1
class DistanceTableClass
{
protected:
  ///Pointer to the SpeciesArray
  PathClass &Path;
  /// Stores the distances between all ptcls at all time slices
  MirroredArrayClass<double> DistTable; ///<(timeslice, ptcl x ptcl)
  /// Stores the displacements between all ptcls at all time slices
  MirroredArrayClass<dVec> DispTable; ///<(timeslice, ptcl x ptcl)  
  ///Table of Image numbers (stored 0-63 for 3d)
  MirroredArrayClass<int> ImageNumTable; ///<(timeslice, ptcl x ptcl)
  ///Stores what vectors the image numbers corresponds to
  Array<dVec,1> ImageVectors;
  Array<dVecInt,2> LowerDimImageMasks;
  Array<dVec,2> LowerDimVectorMasks;
  int FirstLowerDimSpecies;
  int FirstLowerDimPtcl;
  inline void ArrayIndex(int ptcl1, int ptcl2, 
			 int &index, double &sign) const;
public:
  virtual void Update (int timeSlice, Array<int,1> ptclArray) = 0;
  virtual void UpdateAll() = 0;
  virtual void UpdateAll(int timeSlice) = 0;
  inline void DistDisp(int timeSlice, int ptcl1, int ptcl2, 
		       double &distance, dVec &displacement);
  inline void DistDisp(int timeSliceA, int timeSliceB,
		       int ptcl1, int ptcl2, double &distA, double &distB,
		       dVec &dispA, dVec &dispB);
  inline DistanceTableClass (PathClass &myPath);
  inline void AcceptCopy(int startTimeSlice, int endTimeSlice,
			 const Array<int,1> &activeParticles);
  
  inline void RejectCopy(int startTimeSlice, int endTimeSlice,
			 const Array<int,1> &activeParticles);
};

inline void DistanceTableClass::AcceptCopy(int startTimeSlice, 
					   int endTimeSlice, 
					   const Array<int,1> &activeParticles)
{


  for (int slice=startTimeSlice;slice<=endTimeSlice;slice++){
    for (int ptcl1Index=0;ptcl1Index<activeParticles.size();ptcl1Index++){
      int ptcl1=activeParticles(ptcl1Index);
      for (int ptcl2=0;ptcl2<Path.NumParticles();ptcl2++){
	double dummySign;
	int index;
	ArrayIndex(ptcl1,ptcl2,index,dummySign);
	DistTable.AcceptCopy(slice,index);
	DispTable.AcceptCopy(slice,index);
	ImageNumTable.AcceptCopy(slice,index);
	
      }
    }
  }
}

inline void DistanceTableClass::RejectCopy(int startTimeSlice, 
					   int endTimeSlice,
					   const Array<int,1> &activeParticles)
{


  for (int slice=startTimeSlice;slice<=endTimeSlice;slice++){
    for (int ptcl1Index=0;ptcl1Index<activeParticles.size();ptcl1Index++){
      int ptcl1=activeParticles(ptcl1Index);
      for (int ptcl2=0;ptcl2<Path.NumParticles();ptcl2++){
	double dummySign;
	int index;
	ArrayIndex(ptcl1,ptcl2,index,dummySign);
	DistTable.RejectCopy(slice,index);
	DispTable.RejectCopy(slice,index);
	ImageNumTable.RejectCopy(slice,index);
	
      }
    }
  }
}
      

      
  


/// Return the image number correspoding to the negative of the
/// image.  Basically maps (i,j,k) into (-i,-j,-k) if
/// {i,j,k} are all int he range [-1,1].
inline int NegateImageNum(int imageNum)
{
  int mask = ~((~0)<<(2*NDIM));
  int sub = 0x55555555;  // binary: 01010101010101
  return (((~imageNum)-sub)&mask);
}



/// Takes a dVecInt, whose elements have the range [-1,1] and maps
/// it into a unique integer.  Works for dimensions < 16.
inline int ImageNum(dVecInt image)
{
  int imageNum;
  image +=1;
  if (NDIM==1)
    imageNum = image[0];
  else if (NDIM==2)
    imageNum = image[0] + (image[1] << 2);
  else if (NDIM==3)
    imageNum = image[0] + (image[1]<<2) + (image[2]<<4);
  return (imageNum);
}

/// Performs the reverse mapping of the previous function
inline dVecInt Image(int imageNum)
{
  dVecInt image;
  if (NDIM==1)
    image[0] = imageNum;
  else if (NDIM==2)
    {
      image[0]=(imageNum&3);
      imageNum >>=2;
      image[1]=(imageNum&3);
    }
  else if (NDIM==3)
    {
      image[0]=(imageNum&3);
      imageNum >>=2;
      image[1]=(imageNum&3);
      imageNum >>=2;
      image[2]=(imageNum&3);
    }
  image -= 1;
  return image;
}  



inline void DistanceTableClass::DistDisp(int timeSlice,int ptcl1, int ptcl2,
		     double &distance, dVec &displacement)
{
  int index;
  double sign;
  ArrayIndex(ptcl1,ptcl2,index,sign);
  distance=DistTable(timeSlice,index);
  displacement=sign*DispTable(timeSlice,index);
}

inline void DistanceTableClass::DistDisp(int timeSliceA, int timeSliceB,
		     int ptcl1, int ptcl2, double &distA, double &distB,
		     dVec &dispA, dVec &dispB)
{
  DistDisp(timeSliceA,ptcl1,ptcl2,distA,dispA);
  int index;
  double sign;
  ArrayIndex(ptcl1,ptcl2,index,sign);
  distA=DistTable(timeSliceA,index);
  dispA=sign*DispTable(timeSliceA,index);
  int imageNumA=ImageNumTable(timeSliceA,index);
  int imageNumB=ImageNumTable(timeSliceB,index);
  if (imageNumA!=imageNumB){
    dispB=Path(timeSliceB,ptcl1)-Path(timeSliceB,ptcl2)+
      ImageVectors(imageNumA)*sign;
    distB=sqrt(dot(dispB,dispB));
  }
  else {
    dispB=sign*DispTable(timeSliceB,index);
    distB=DistTable(timeSliceB,index); 
  }




}



inline DistanceTableClass::DistanceTableClass (PathClass &myPath) :
  Path(myPath)
{
  int size=(Path.NumParticles()*(Path.NumParticles()+1))/2;
  DistTable.Resize(Path.NumTimeSlices(),size);
  DispTable.Resize(Path.NumTimeSlices(),size);
  ImageNumTable.Resize(Path.NumTimeSlices(),size);
  		   
  ///Initialize the ImageVectors
  int NumVectors = 1;
  for (int i=0; i<NDIM; i++)
    NumVectors <<= 2;
  ImageVectors.resize (NumVectors);
  
  for (int i=0; i<NumVectors; i++)
    {
      dVecInt image = Image(i);
      for (int dim=0; dim<NDIM; dim++)
	ImageVectors(dim) = (double)image[dim]*Path.Box[dim];
    }
  // Now initialize the Masks for lower dimensionality;
  int NumSpecies = Path.NumSpecies();
  LowerDimImageMasks.resize(NumSpecies,NumSpecies);
  LowerDimVectorMasks.resize(NumSpecies,NumSpecies);
  FirstLowerDimSpecies = NumSpecies;
  for (int species1=0; species1<NumSpecies; species1++) {
    dVecInt imageMask;
    dVec vectorMask;
    imageMask = 1;
    vectorMask = 1.0;
    for (int dim=0; dim<NDIM; dim++)
      if (!Path.Species(species1).DimensionActive[dim]) {
	if (FirstLowerDimSpecies > species1)
	  FirstLowerDimSpecies = species1;
	imageMask[dim] = 0;
	vectorMask[dim] = 0.0;
      }
    for (int species2=0; species2<NumSpecies; species2++)
      {
	for (int dim=0; dim<NDIM; dim++)
	  if (!Path.Species(species2).DimensionActive[dim]) {
	    imageMask[dim] = 0;
	    vectorMask[dim] = 0.0;
	  }
	LowerDimImageMasks(species1,species2) = imageMask;
	LowerDimVectorMasks(species1,species2) = vectorMask;
      }
  }  
  if (FirstLowerDimSpecies<Path.NumSpecies()){
    FirstLowerDimPtcl=Path.Species(FirstLowerDimSpecies).FirstPtcl;
  }
  else {
    FirstLowerDimPtcl=Path.NumParticles();
  }

}


/// Maps particle x particle to a location in our distance or
/// displacement table. 
inline void DistanceTableClass::ArrayIndex(int ptcl1, int ptcl2, 
					   int &index, double &sign) const
{
  // Make sure ptclA > ptclB
  sign = 1.0;
  if (ptcl1<ptcl2)
    {
      int temp = ptcl1;
      ptcl1 = ptcl2;
      ptcl2 = temp;
      sign = -1.0;
    }
  //ptcl1*(ptcl1+1)/2+ptcl2;  
  index = ((ptcl1*(ptcl1+1))>>1)+ptcl2;   
}


#endif
