#include "MoveUtils.h"

double dotprod(dVec vec1, dVec vec2) {
  double total = 0;
  for(int i = 0; i<3; i++){
    total += vec1[i]*vec2[i];
  }
  return total;
}

double dotprod(dVec vec1, dVec vec2, double mag) {
  double norm = 1.0/mag;
  return dotprod(vec1, vec2)*norm;
}

double Mag(dVec v) {
  double mag = sqrt(dotprod(v,v));
  return mag;
}

dVec Normalize(dVec v) {
  double mag = Mag(v);
  dVec norm = v/mag;
  return norm;
}

dVec Scale(dVec v, double scale) {
  double mag = Mag(v);
  dVec norm = v/mag;
  norm *= scale;
  return norm;
}

void Strip(dVec R, dVec u,dVec &aligned, dVec &perp){
  aligned = Scale(R,dotprod(R,u));
  perp = u - aligned;
}
