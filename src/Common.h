#ifndef COMMON_H
#define COMMON_H

#define SIMPLE_SPRNG


#define NDIM 3
#define OLD 0
#define NEW 1
#include <sprng.h>
#include "Blitz.h"
typedef enum {MOVEMODE, OBSERVABLEMODE} ModeType;

///ParticleID=(species,particle number)

typedef TinyVector<int,2> ParticleID;

extern int Write1;
extern int Write2; //These are the variables for writing 


int GetCurrentTimeStamp();
dVec GuassianRandomVec(double sigma);
dVec dVecSubtract(dVec a,dVec b); //a-b
void setMode(int);
double distSqrd(dVec a,dVec b); 
#endif
