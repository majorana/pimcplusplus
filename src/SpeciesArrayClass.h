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
  /// Returns the size of the array
  inline int size(){
    return SpeciesArray.size();
  }
  /// Resizes the array
  inline void resize(int newsize)
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



#endif

