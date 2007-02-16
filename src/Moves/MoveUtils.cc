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

dVec crossprod(dVec v1, dVec v2)
{
  dVec cross;
  cross[0] = v1[1]*v2[2] - v1[2]*v2[1];
  cross[1] = -v1[0]*v2[2] + v1[2]*v2[0];
  cross[2] = v1[0]*v2[1] - v1[1]*v2[0];
  return cross;
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

dVec Renormalize(dVec v, double scale) {
  double mag = Mag(v);
  dVec norm = v/mag;
  norm *= scale;
  return norm;
}

dVec Scale(dVec v, double scale) {
  cerr << "CALLING SCALE, NOW DEPRECATED.  PLEASE REPLACE WITH RENORMALIZE IN MOVEUTILS.  IT'S THE SAME THING BUT THE NAME IS CLEARER" << endl;
  double mag = Mag(v);
  dVec norm = v/mag;
  norm *= scale;
  return norm;
}

void Strip(dVec R, dVec u,dVec &aligned, dVec &perp){
  //aligned = Scale(R,dotprod(R,u));
  aligned = Renormalize(R,dotprod(R,u));
  perp = u - aligned;
}

double GetAngle(dVec v1, dVec v2)
{
  double mag = Mag(v1);
  mag *= Mag(v2);
  double dot = dotprod(v1,v2,mag);
  double angle = acos(dot);
  //cerr << "Computed angle " << angle << " in deg " << angle*180/M_PI << endl;
//  if (dot > 1.0){
//    cerr << "OH CRAP: DOT PRODUCT IS " << dot << " between " << v1 << " and " << v2 << "; I used mag " << mag << " and I'm going to return " << angle << endl;
//  }
  if (dot-1 < 0.0001 && dot-1 > 0){
    cerr << "correcting angle" << endl;
    angle = 0.0;
  }
  return angle;
}

dVec GetBisector(dVec v1, dVec v2)
{
  dVec bisector = v1 + v2;
  return bisector;
}
