#ifndef DISTANCE_TABLE_BASE_CLASS_H
#define DISTANCE_TABLE_BASE_CLASS_H

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
  MirroredSymmetricMatrixClass<double> DistTable; ///<(timeslice, ptcl x ptcl)
  /// Stores the displacements between all ptcls at all time slices
  /// (timeslice, ptcl x ptcl)  
  MirroredAntiSymmetricMatrixClass<dVec> DispTable; 
  ///Table of Image numbers (stored 0-63 for 3d)
  ///(timeslice, ptcl x ptcl)
  MirroredAntiSymmetricMatrixClass<ImageNumClass> ImageNumTable; 
  ///Stores what vectors the image numbers corresponds to
  Array<dVec,1> ImageVectors;
  Array<dVecInt,2> LowerDimImageMasks;
  Array<dVec,2> LowerDimVectorMasks;
  int FirstLowerDimSpecies;
  int FirstLowerDimPtcl;
  inline void ArrayIndex(int ptcl1, int ptcl2, 
			 int &index, double &sign) const;
public:
  virtual void Update (int timeSlice, const Array<int,1> &ptclArray) = 0;
  virtual void UpdateAll() = 0;
  virtual void UpdateAll(int timeSlice) = 0;
  virtual void PutInBox(dVec &r)=0;
  void ShiftData(int numTimeSlicesToShift,PIMCCommunicatorClass& Communicator);
  inline void DistDisp(int timeSlice, int ptcl1, int ptcl2, 
		       double &distance, dVec &displacement);
  inline void DistDisp(int timeSliceA, int timeSliceB,
		       int ptcl1, int ptcl2, double &distA, double &distB,
		       dVec &dispA, dVec &dispB);

  ///Slow! Debugging purposes only!
  virtual void DistDispTest(int timeSlice, int ptcl1, int ptcl2, 
		       double &distance, dVec &displacement)=0;
  ///Slow! Debugging purposes only!
  virtual void DistDispTest (int timeSliceA, int timeSliceB,
			    int ptcl1, int ptcl2, double &distA, double &distB,
			    dVec &dispA, dVec &dispB) = 0;


  virtual dVec Velocity(int timeSliceA, int timeSliceB, int ptcl) = 0;
  inline DistanceTableClass (PathClass &myPath);
  inline void AcceptCopy(int startTimeSlice, int endTimeSlice,
			 const Array<int,1> &activeParticles);
  
  inline void RejectCopy(int startTimeSlice, int endTimeSlice,
			 const Array<int,1> &activeParticles);
  inline void MoveJoin (int oldJoin, int newJoin);



};

inline void DistanceTableClass::MoveJoin (int oldJoin, int newJoin)
{
  DistTable.MoveJoin(Path.Permutation, oldJoin, newJoin);
  DispTable.MoveJoin(Path.Permutation, oldJoin, newJoin);
  ImageNumTable.MoveJoin(Path.Permutation, oldJoin, newJoin);
//   SetMode(OLDMODE);
//   UpdateAll();
//   SetMode(NEWMODE);
//   UpdateAll();
}


inline void DistanceTableClass::AcceptCopy(int startTimeSlice, 
					   int endTimeSlice, 
					   const Array<int,1> &activeParticles)
{
  DistTable.AcceptCopy(startTimeSlice, endTimeSlice, activeParticles);
  DispTable.AcceptCopy(startTimeSlice, endTimeSlice, activeParticles);
  ImageNumTable.AcceptCopy(startTimeSlice, endTimeSlice, activeParticles);

}

inline void DistanceTableClass::RejectCopy(int startTimeSlice, 
					   int endTimeSlice,
					   const Array<int,1> &activeParticles)
{

  DistTable.RejectCopy(startTimeSlice, endTimeSlice, activeParticles);
  DispTable.RejectCopy(startTimeSlice, endTimeSlice, activeParticles);
  ImageNumTable.RejectCopy(startTimeSlice, endTimeSlice, activeParticles);



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
inline dVecInt Image(ImageNumClass myImageNum)
{
  int imageNum = myImageNum.ImageNum;
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
  
#ifdef DEBUG
  dVec dispTest;
  double distTest;
  DistDispTest(timeSlice,ptcl1,ptcl2,distTest,dispTest);
  if (!(dispTest==displacement)){
    cerr<<"Broken: "<<dispTest<<" "<<displacement<<endl;
  }
  if (!(distTest==distance)){
    cerr<<"Broken dist: "<<dispTest<<" "<<displacement<<endl;
    cerr<<"Broken dist: "<<distTest<<" "<<distance<<endl;
    cerr<<Path(timeSlice,ptcl1)<<" "<<Path(timeSlice,ptcl2);
  }
  assert(dispTest==displacement);
  assert(distTest==distance);
    
#endif
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
    dispB=Path(timeSliceB,ptcl2)-Path(timeSliceB,ptcl1)-
      ImageVectors(imageNumA)*sign;
    distB=sqrt(dot(dispB,dispB));
  }
  else {
    dispB=sign*DispTable(timeSliceB,index);
    distB=DistTable(timeSliceB,index); 
  }


#ifdef DEBUG
  dVec dispATest, dispBTest;
  double distATest, distBTest;

  DistDispTest(timeSliceA,timeSliceB,ptcl1,ptcl2,distATest,distBTest,
	       dispATest,dispBTest);
  //  cerr<<dispBTest<<" "<<dispB<<" "<<endl;
  if (!(dispBTest==dispB)){
    cerr<<dispBTest<<" "<<dispB<<" "<<endl;
    cerr<<ptcl1<<" "<<ptcl2<<" "<<timeSliceA<<" "<<timeSliceB<<" "<<endl;
    cerr<<dispATest<<" "<<dispA<<" "<<endl;
    cerr<<imageNumA<<" "<<imageNumB<<endl;
    cerr<<sign<<endl;
    cerr<<"IM: "<<ImageVectors(imageNumA)<<" "<<ImageVectors(imageNumB)<<endl;
  }
  assert(dispATest==dispA);
  assert(dispBTest==dispB);
  assert(distATest==distA);
  assert(distBTest==distB);

#endif

}



inline DistanceTableClass::DistanceTableClass (PathClass &myPath) :
  Path(myPath)
{
  int size=(Path.NumParticles()*(Path.NumParticles()+1))/2;
  DistTable.Resize(Path.NumTimeSlices(),Path.NumParticles());
  DispTable.Resize(Path.NumTimeSlices(),Path.NumParticles());
  ImageNumTable.Resize(Path.NumTimeSlices(),Path.NumParticles());
  		   
  ///Initialize the ImageVectors
  int NumVectors = 1;
  for (int i=0; i<NDIM; i++)
    NumVectors <<= 2;
  ImageVectors.resize (NumVectors);
  
  for (int i=0; i<NumVectors; i++)
    {
      dVecInt image = Image(i);
      for (int dim=0; dim<NDIM; dim++)
	ImageVectors(i)[dim] = (double)image[dim]*Path.Box[dim];
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
