//#include "Common/Blitz.h"
#include "../PathDataClass.h"
#ifndef RANDOM_PERM_CLASS_H
#define RANDOM_PERM_CLASS_H


class RandomPermClass
{
 private:
/*   double GetActionEstimate(int i, int j, double fourLambdaBetaInv, */
/* 			   Array<double,2> speciesData); */
  Array<double,2> ActionPairs;




 public:
  void GetActionPairs();
  void Apply();
  void Read(IOSectionClass &inSection);
  PathDataClass &PathData;
  RandomPermClass(PathDataClass &myPathData);
  int Slice1,Slice2;
  int SpeciesNum;
  Array<int,1> PermArray; //Must resize perm array
  void GeneratePerm();
  Array<int,1> CurrentParticles;
};




/* inline double RandomPermClass::GetActionEstimate(int i,int j, */
/* 						 double fourLambdaBetaInv, */
/* 						 Array<double,2> speciesData) */
/* { */
/*   dVec disp_ij=speciesData(Slice2,j)-speciesData(Slice1,i); */
/*   PathData.Path.PutInBox(disp_ij); */
/*   double dist_ij=dot(disp_ij,disp_ij); */
/*   return exp(-dist_ij*fourLambdaBetaInv); */
/* } */





#endif
