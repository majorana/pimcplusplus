
#ifndef MULTI_STAGE_H
#define MULTI_STAGE_H


#include <list>
#include "MoveBase.h"
#include "../Actions/ActionBase.h"

class StageClass
{
protected:
  PathDataClass &PathData;
public:
  int NumAccepted, NumAttempted;
  int BisectionLevel;
  list<ActionBaseClass*> Actions;
  ///The first stage will set the slices and activeParticles
  ///This returns transition probability ratio T(new->old)/T(old->new)
  virtual double Sample (int &slice1,int &slice2,
			 Array<int,1> &activeParticles) = 0; 
  virtual void Read (IOSectionClass &in);
  virtual void Accept();
  virtual void Reject();
  inline double StageAction(int startSlice,int endSlice,
			    const Array<int,1> &changedParticles);
  inline double GlobalStageAction (const Array<int,1> &changeParticles);
  inline double AcceptRatio () 
  { return (double)NumAccepted / (double) NumAttempted; }

  virtual bool Attempt(int &slice1, int &slice2, 
		       Array<int,1> &activeParticles,
		       double &prevActionChange) = 0;

  StageClass(PathDataClass &pathData) :
    PathData(pathData), NumAccepted(0), NumAttempted(0),
    BisectionLevel(0)
  {
    // Do nothing for now
  }

};

class CommonStageClass : public StageClass
{
public:
  bool Attempt(int &slice1, int &slice2,
	       Array<int,1> &activeParticles,
	       double &prevActionChange);
  CommonStageClass(PathDataClass &pathData) :
    StageClass(pathData)
  {
    // Do nothing for now
  }
	       
};

class LocalStageClass : public StageClass
{
public:
  bool Attempt(int &slice1, int &slice2,
	       Array<int,1> &activeParticles,
	       double &prevActionChange);
  LocalStageClass(PathDataClass &pathData) :
    StageClass(pathData)
  {
    // Do nothing for now
  }
	       
};


inline double StageClass::StageAction(int startSlice,int endSlice,
				      const Array<int,1> &changedParticles)
{
  double TotalAction=0.0;
  list<ActionBaseClass*>::iterator actionIter=Actions.begin();
  while (actionIter!=Actions.end()){
    TotalAction += 
      ((*actionIter)->Action(startSlice, endSlice, changedParticles,
			     BisectionLevel));
    actionIter++;
  }
  return TotalAction;
}

inline double 
StageClass::GlobalStageAction (const Array<int,1> &changedParticles)
{
  int slice1 = 0;
  int slice2 = PathData.Path.NumTimeSlices()-1;
  double localAction = StageAction (slice1, slice2, changedParticles);
  double globalAction = PathData.Path.Communicator.AllSum (localAction);

  return globalAction;
}


class MultiStageClass : public ParticleMoveClass
{
protected:
  list<StageClass*> Stages;
  int NumSteps;
  int Slice1,Slice2;
public:
  void Read(IOSectionClass &io);
  void Accept();
  void Reject();
  virtual void WriteRatio()=0;

  ///Why was this MakeMove()=0 and virtual?
  void MakeMove();
  MultiStageClass(PathDataClass &pathData, IOSectionClass &outSection) : 
    ParticleMoveClass(pathData,outSection) 
  {
    //do nothing for now
  }
};

// class MultiStageLocalClass : public MultiStageClass
// {
// public:
//   void MakeMove();
//   MultiStageLocalClass(PathDataClass &pathData, IOSectionClass &outSection) :
//     MultiStageClass(pathData,outSection)
//   {
//     //do nothing for now
//   }
// };


// class MultiStageCommonClass : public MultiStageClass
// {
// public:
//   void MakeMove();
//   MultiStageCommonClass(PathDataClass &pathData, IOSectionClass &outSection) :
//     MultiStageClass(pathData,outSection)
//   {
//     //do nothing for now
//   }
// };


#endif
