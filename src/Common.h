#ifndef COMMON_H
#define COMMON_H

#define SIMPLE_SPRNG

/// Number of dimensions that dVec uses.
#define NDIM 3

#include <sprng.h>
#include "Blitz.h"


/// These are the different mode types for the MirroredArrayClass
typedef enum {OLDMODE, NEWMODE, BOTHMODE} ModeType;

///ParticleID=(species,particle number)

typedef TinyVector<int,2> ParticleID;

/// These are the global variables to be used to decide what part of
/// the mirrored array we are writing to and reading from  
extern int Write1;
extern int Write2; 

/// Returns the current monte-carlo time step
int GetCurrentTimeStamp();

/// Produces a guassian random vector with radius that has variance
/// sigma 
dVec GaussianRandomVec(double sigma);

/// Changes the mode the entire code is running in.
void SetMode(ModeType);

class ImageNumClass
{
public:
  int ImageNum;
  inline ImageNumClass operator-() const
  {
    ImageNumClass minusNum;
    int mask = ~((~0)<<(2*NDIM));
    int sub = 0x55555555;  // binary: 01010101010101
    minusNum.ImageNum= (((~ImageNum)-sub)&mask);
    return minusNum;
  }
  inline operator int() 
  {
    return (ImageNum);
  }
  inline ImageNumClass(int i)
  {
    ImageNum = i;
  }
  inline ImageNumClass()
  { }
};


#endif
