#ifndef ARRAYOFIDENTICALPARTICLECLASS_H
#define ARRAYOFIDENTICALPARTICLECLASS_H

#include "IdenticleParticleClass.h"

class ArrayOfIdenticalParticlesClass
{
  Array<IdenticalParticlesClass*,1> IdenticalParticlesArray;
 public:
  inline int size(){
    return IdenticalParticlesArray.size();
  }
  inline IdenticalParticlesClass& operator()(int i){
    return (*(IdenticalParticlesArray(i)));
  }
  //  inline IdenticalParticlesClass operator()(int i) const {
  //    return (*(IdenticlesParticleArray(i)));
  //  }

};



#endif

