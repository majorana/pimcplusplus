#ifndef COMMON_H
#define COMMON_H

#define NDIM 3
#define OLD 0
#define NEW 1
#include "Blitz.h"
typedef enum {MOVEMODE, OBSERVABLEMODE} ModeType;

///ParticleID=(species,particle number)

typedef TinyVector<int,2> ParticleID;

int GetCurrentTimeStamp();
dVec GuassianRandomVec(double sigma);
dVec dVecSubtract(dVec a,dVec b); //a-b
void setMode(int);
#endif
