#ifndef ARRAYOFIDENTICALPARTICLECLASS_H
#define ARRAYOFIDENTICALPARTICLECLASS_H

#include "SpeciesClass.h"

/// This is a wrapper class that holds an array of pointers to the species.
/// We have wrapped this array because we need it to be an array of SpeciesClass* 
class SpeciesArrayClass
{
  /// This is the actual array that holds the SpeciesClass objects.
  Array<SpeciesClass*,1> SpeciesArray; 
 public:
  /// Calculates the distance and displacement betweeen
  ///two particles and their respective timeslices
  inline void Distance(ParticleID particle1, ParticleID particle2,int timeSlice1, int timeSlice2,
	   double &dist,dVec &disp );
  
  /// Returns the size of the array
  inline int Size(){
    return SpeciesArray.size();
  }
  /// Resizes the array
  inline void Resize(int newsize)
    {
      //      cerr<<"The Array of SpeciesClass is being resized to ";
      //      cerr<<newsize<<endl;
      SpeciesArray.resize(newsize);
    }
  /// Sets one of the Species Objects in the Species Array
  inline void Set(int i, SpeciesClass &IDptcls)
    {
      SpeciesArray(i) = &IDptcls;
    }
  /// Returns a pointer to the Species Object in the Species Array
  inline SpeciesClass& operator()(int i){
    return (*(SpeciesArray(i)));
  }
  /// Sets the value of the Ptcl x Slice in the (Path Object) of the Species Object of the Species Arrya
  inline void SetPos(int Species, int Ptcl, int Slice, const dVec &pos){
    SpeciesArray(Species)->Path.SetPos(Ptcl, Slice, pos);
  }
  /// Returns  the value of the Ptcl x Slice in the (Path Object) of the Species Object of the Species Array
  inline dVec operator()(int Species, int Ptcl, int Slice) const {
    return (SpeciesArray(Species)->Path(Ptcl,Slice));
  }



};

  /// Calculates the distance and displacement betweeen
  ///two particles and their respective timeslices
inline void Distance(ParticleID particle1, ParticleID particle2,int timeSlice1, int timeSlice2,
	   double &dist,dVec &disp )
{
  


}


#endif

