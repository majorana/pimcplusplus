#ifndef BISECTION_STAGE_CLASS_H
#define BISECTION_STAGE_CLASS_H

#include "MultiStage.h"
#include "../Observables/ObservableBase.h"

class BisectionStageSphereClass : public LocalStageClass
{
public:
  void WriteRatio();
  ObservableDouble AcceptRatioVar;
  double Sample(int &slice1,int &slice2, 
		Array<int,1> &activeParticles);
  ///Really needs to go someplace else, but for the moment we will keep it here
  double SphereRadius; 
  void Accept();
  void Reject();
  void CartesianToSpherical(dVec &r, double &theta,double &phi);
  TinyVector<double,3> SphericalToCartesian(double &a,double &theta, double &phi);
  void ProjectOntoSphere(dVec &r, double a);
  void RotateAroundX(dVec &vec, double theta);
  void RotateAroundY(dVec &vec, double theta);
  void RotateAroundZ(dVec &vec, double theta);
  double angleInBetween(dVec &r1, dVec &r2);
  BisectionStageSphereClass(PathDataClass &pathData, int level,
		      IOSectionClass outSection) : 
    LocalStageClass(pathData,outSection),
    AcceptRatioVar("AcceptRatio",OutSection,pathData.Path.Communicator) 
  { 
    //do nothing for now
    BisectionLevel = level;
    SphereRadius=31.0;
  }
};

#endif
