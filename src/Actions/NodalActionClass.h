#ifndef NODAL_ACTION_CLASS_H
#define NODAL_ACTION_CLASS_H

#include "ActionBase.h"
#include <Common/Splines/CubicSpline.h>

class PathClass;

typedef enum { FREE_PARTICLE, GROUND_STATE, GROUND_STATE_FP } NodeType;

class NodalActionClass : public ActionBaseClass
{
public:
  virtual bool IsPositive (int slice) = 0;
  //  virtual double Det(int slice)       = 0;
  //  virtual Array<double,2> GetMatrix (int slice=0) = 0;
  virtual void AcceptCopy (int slice1, int slice2);
  virtual void RejectCopy (int slice1, int slice2);
  virtual void Init();
  virtual bool IsGroundState() = 0;
  virtual NodeType Type() = 0;
  virtual void Setk (Vec3 kVec);
  virtual void WriteInfo(IOSectionClass &out) = 0;
  virtual void Update();
  NodalActionClass (PathDataClass &pathData) :
    ActionBaseClass (pathData)
  {

  }
};


#endif
