#ifndef DISTANCE_TABLE_PBC_CLASS_H
#define DISTANCE_TABLE_PBC_CLASS_H

#include "DistanceTableClass.h"

class DistanceTablePBCClass : public DistanceTableClass
{
private:
  inline void Displacement (int timeslice, int ptcl1, int ptcl2,
			    dVec &disp, double &dist, int &imageNum);
  inline void Displacement (int timeSlice, int ptcl1, int ptcl2,
			    dVec &disp, double &dist, int &imageNum,
			    dVec vecMask, dVecInt imageMask);
public:
  void Update (int timeSlice, const Array<int,1> &ptclArray);
  void UpdateAll();
  void UpdateAll(int timeSlice);
  /// Constructor
  DistanceTablePBCClass (PathClass &myPath) : DistanceTableClass(myPath)
  { /* Currently DistanceTable constructor does everything */ }
};


/// Determine the particle displacement
inline void DistanceTablePBCClass::Displacement(int timeSlice,
						int ptcl1,int ptcl2, 
						dVec &disp,
						double &dist,
						int &imageNum)
{
  dVecInt image;
  disp=Path(timeSlice,ptcl2)-Path(timeSlice,ptcl1);
  // Determine the image number of ptcl2 w.r.t ptcl1
  for (int i=0; i<NDIM; i++) {
    image[i] = -(disp[i]<-0.5*Path.Box[i]) + (disp[i]>0.5*Path.Box[i]);
    //    disp[i] += image[i]*Path.Box[i];
  }
  imageNum = ImageNum(image);
  cerr << "ImageNum = " << imageNum << endl;
  cerr << "ImageVectors(imageNum) = " << ImageVectors(imageNum) << endl;
  disp = disp + ImageVectors(imageNum);
  dist = sqrt(dot(disp,disp));
}



///Calculates the displacement of two particles ensuring that they 
///both have the same image
inline void DistanceTablePBCClass::Displacement(int timeSlice,
						int ptcl1,int ptcl2, 
						dVec &disp,
						double &dist,
						int &imageNum,
						dVec vecMask,
						dVecInt imageMask)
{
  dVecInt image;

  // Determine the image number of ptcl2 w.r.t ptcl1
  for (int i=0; i<NDIM; i++) {
    disp[i] = vecMask[i]*(Path(timeSlice,ptcl2)[i] -
			  Path(timeSlice,ptcl1)[i]);
    image[i] = (-(disp[i]<-0.5*Path.Box[i])+(disp[i]>0.5*Path.Box[i]))
      *imageMask[i];
  }
  imageNum = ImageNum(image);
  disp = disp + ImageVectors(imageNum);
  dist = sqrt(dot(disp,disp));
}






#endif
