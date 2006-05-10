#ifndef LEVI_FLIGHT_STAGE_CLASS_H
#define LEVI_FLIGHT_STAGE_CLASS_H

#include "MultiStage.h"
#include "../Observables/ObservableVar.h"


///The levi flight stage does the levi flight on a single
//particle specified by flight particle instead of activeParticles.
//This is intentional to allow the move to work on a section of the 
//path that's different from the one the action sees changed.
class LeviFlightStageClass : public LocalStageClass
{
public:
  void WriteRatio();
  double Sample(int &slice1,int &slice2, 
		Array<int,1> &activeParticles);
  void Accept();
  void Reject();
  void Read(IOSectionClass &in);
  int TotalLevels;
  int FlightPtcl;
  int EndSlice;
  int useStart;
  int useEnd;
  bool Attempt(int &slice1, int &slice2, 
	       Array<int,1> &activeParticles,
	       double &prevActionChange);

  int amountToOpen;

  double LeviFlight(int slice1,int slice2, int ptcl, double lambda);
  LeviFlightStageClass(PathDataClass &pathData, 
		      IOSectionClass outSection) : 
    LocalStageClass(pathData,outSection)
  { 
    //do nothing for now
    

  }
};

#endif
