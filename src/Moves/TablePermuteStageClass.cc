#include "PermuteStageClass.h"


double TablePermuteStageClass::Sample(int &slice1,int &slice2,
				      Array<int,1> &changedParticles)
{
  //do nothing for now
}


bool TablePermuteStageClass::Attempt (int &slice1, int &slice2,
				      Array<int,1> &changedParticles, double &prevActionChange)
{
  if (changedParticles(0) == -1) {
    // Constrct a new permutation
    
    // Now return true always
    return true;
  }
  else {
    // Calculate real transition probability ratio
    
    // Accept or reject based on this ratio
    return true;
  }
}

void TablePermuteStageClass::InitBlock()
{


}
