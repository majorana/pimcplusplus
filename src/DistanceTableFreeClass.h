#ifndef DISTANCE_TABLE_FREE_CLASS_H
#define DISTANCE_TABLE_FREE_CLASS_H






#include "DistanceTableClass.h"

class DistanceTableFreeClass : public DistanceTableClass
{
private:
  inline void Displacement (int timeslice, int ptcl1, int ptcl2,
			    dVec &disp, double &dist);
  inline void Displacement (int timeSlice, int ptcl1, int ptcl2,
			    dVec &disp, double &dist, dVec vecMask);
public:
  void Update (int timeSlice, Array<int,1> ptclArray);
  void UpdateAll();
  void UpdateAll(int timeSlice);
  /// Constructor
  DistanceTableFreeClass (PathClass &myPath) : DistanceTableClass(myPath)
  { /* Currently DistanceTable constructor does everything */ }
};


/// Determine the particle displacement
inline void DistanceTableFreeClass::Displacement(int timeSlice,
						 int ptcl1,int ptcl2, 
						 dVec &disp,
						 double &dist)
{
  disp=Path(timeSlice,ptcl2)-Path(timeSlice,ptcl1);
  dist = sqrt(dot(disp,disp));
}



///Calculates the displacement of two particles ensuring that they 
///both have the same image
inline void DistanceTableFreeClass::Displacement(int timeSlice,
						 int ptcl1,int ptcl2, 
						 dVec &disp,
						 double &dist,
						 dVec vecMask)
{

  double dist2;
  for (int i=0; i<NDIM; i++) {
    disp[i] = vecMask[i]*(Path(timeSlice,ptcl2)[i] -
			  Path(timeSlice,ptcl1)[i]);
    dist2+=disp[i]*disp[i];
  }

  
  dist = sqrt(dist2);
}




























// class DistanceTableFreeClass : public DistanceTableClass
// {
// public:
//   void Update (int timeSlice, int ptcl);
//   void UpdateAll();
//   void UpdateAll(int timeSlice);
//   /// Constructor
//   DistanceTableFreeClass (PathClass &myPath) : DistanceTableClass(myPath)
//   { /* Currently DistanceTable constructor does everything */ }
// };



// void DistanceTableFreeClass::UpdateAll(int timeSlice)
// {
//   int index=0;
//   for (int ptcl1=0;ptcl1<Path.NumParticles();ptcl1++){
//     for (int ptcl2=0;ptcl2<=ptcl1;ptcl2++){
//       DispTable.Set(timeSlice,index,	
// 		    Path(timeSlice,ptcl2) - Path(timeSlice,ptcl1));  
//       ImageNumTable.Set(timeSlice,index,0);
//       index++;
//     }
//   }
//   /// Now zero out components of displacement due to inactive
//   /// dimensions of lower dimensional particles, eg. a cylinder in
//   /// 3-space.
//   for (int specCntr1=0; specCntr1<Path.NumSpecies(); specCntr1++) {
//     SpeciesClass &species1=Path.Species(specCntr1);
//     if (species1.NumDim<NDIM){
//       for (int dim=0; dim<NDIM; dim++) {
// 	if (!species1.DimensionActive(dim)){
// 	  for (int ptcl=0;ptcl<Path.NumParticles();ptcl++){
// 	    int index=ArrayIndex(ptcl1,ptcl2);
// 	    // This zeros out the image number corresponding to the
// 	    // dimension dim;
// 	    dVec tempdVec=DispTable(timeSlice,index);
// 	    tempdVec[dim]=0.0;
// 	    DispTable.Set(timeSlice,index,tempdVec);
// 	  }
// 	}
//       }
//     }
//   }
//   index = 0;
//   for (int ptcl1=0;ptcl1<Path.NumParticles();ptcl1++){
//     for (int ptcl2=0;ptcl2<=ptcl1;ptcl2++){
//       dVec& disp= DispTable(timeSlice,index);
//       DistanceTable(timeSlice,index)=sqrt(dot(disp,disp));
//       index++;
//     }
//   }
// }



// void DistanceTableFreeClass::UpdateAll()
// {
//   for (int timeSlice=0;timeSlice<Path.NumTimeSlices();timeSlice++){
//     UpdateAll(timeSlice);
//   }
// }


// void DistanceTableFreeClass::Update(int timeSlice,int ptcl1)
// {
//   SpeciesClass &Species1=Path.Species(Path.ParticleSpeciesNum(ptcl1));
  
//   for (int ptcl2; ptcl2<Path.NumParticles(); ptcl2++){
//     int index;
//     double sign;
//     ArrayIndex(ptcl1, ptcl2, index, sign);
//     DispTable(timeSlice,index) = Path(timeSlice,ptcl2)-Path(timeSlice,ptcl1);
//     DispTable(timeSlice,index) *=sign;
//   }
//   for (int speciesNum=0; speciesNum<Path.NumSpecies; speciesNum++) {
//     SpeciesClass &Species = Path.Species(speciesNum);
//     if (Species.NumDim < NDIM || Species1.NumDim<NDIM) {
//       for (int dim=0; dim<NDIM; dim++) 
// 	// Zero out inactive dimensions
// 	if (!Species.DimensionActive(dim) || 
// 	    !Species1.DimensionActive(dim))
// 	  for (int ptcl2=Species.StartPtcl; ptcl2<Species.EndPtcl; ptcl2++) {
// 	    int index = ArrayIndex(ptcl1,ptcl2);
// 	    DispTable(timeSlice,index)[dim] = 0.0;
// 	  }
//     } 
//   }
//   for (int ptcl2=0; ptcl2<Path.NumParticles(); ptcl2++) {
//     int index = ArrayIndex(ptcl1, ptcl2);
//     dVec &disp = DispTable(timeSlice,index);
//     DistTable(timeSlice,index) = sqrt(dot(disp,disp));
//   }
// }



#endif
