#ifndef COMMON_H
#define COMMON_H

#include "Blitz.h"
typedef enum {MOVEMODE, OBSERVABLEMODE} ModeType;

///ParticleID=(species,particle number)

typedef TinyVector<int,2> ParticleID;

int GetCurrentTimeStamp();

#endif
