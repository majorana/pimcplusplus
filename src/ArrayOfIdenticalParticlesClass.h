#ifndef ARRAYOFIDENTICALPARTICLECLASS_H
#define ARRAYOFIDENTICALPARTICLECLASS_H

#include "IdenticalParticlesClass.h"

class ArrayOfIdenticalParticlesClass
{
  Array<IdenticalParticlesClass*,1> IdenticalParticlesArray;
 public:
  inline int size(){
    return IdenticalParticlesArray.size();
  }
  inline void resize(int newsize)
    {
      //      cerr<<"The Array of IdenticalParticleClass is being resized to ";
      //      cerr<<newsize<<endl;
      IdenticalParticlesArray.resize(newsize);
    }
  inline void Set(int i, IdenticalParticlesClass &IDptcls)
    {
      IdenticalParticlesArray(i) = &IDptcls;
    }
  inline IdenticalParticlesClass& operator()(int i){
    return (*(IdenticalParticlesArray(i)));
  }
  inline void SetPos(int Species, int Ptcl, int Slice, dVec &pos){
    IdenticalParticlesArray(Species)->Path.SetPos(Ptcl, Slice, pos);
  }
  inline dVec operator()(int Species, int Ptcl, int Slice) const {
    return (IdenticalParticlesArray(Species)->Path(Ptcl,Slice));
  }


  //  inline IdenticalParticlesClass operator()(int i) const {
  //    return (*(IdenticlesParticleArray(i)));
  //  }

};



#endif

