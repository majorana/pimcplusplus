#ifndef MOVE_UTILITIES_H
#define MOVE_UTILITIES_H

#include "MoveBase.h"

// A number of functions used by several 
// Move and Action classes
// -jg

// return the dot product of two vectors
double dotprod(dVec vec1, dVec vec2);

// return the dot product of two vectors,
// normalizing by a given value
// i.e. gives the dot product of two unit vectors
double dotprod(dVec vec1, dVec vec2, double mag);

// returns the magnitude of a vector
double Mag(dVec v);

// normalizes a vector
// i.e. generates a unit vector
dVec Normalize(dVec v);

// Scales a vector by a give scalar
dVec Scale(dVec v, double scale);

// Takes a vector u and a unit vector R
// and returns the vector aligned that is
// the component of u parallel to R, and
// perp, the component of u perpendicular to R
void Strip(dVec R, dVec u,dVec &aligned, dVec &perp);

#endif
