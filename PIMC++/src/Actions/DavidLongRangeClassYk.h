/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#ifndef DAVID_LONG_RANGE_CLASS_YK_H
#define DAVID_LONG_RANGE_CLASS_YK_H

#include "ActionBase.h"
#include <Common/PairAction/PAFit.h>




/// The LongRangeClass is an action class responsible for the long
/// wavelength components of the action that are summed in k-space.
/// This class reads in and uses an optimized breakup from a file that
/// David supplies. 

class DavidLongRangeClassYk : public ActionBaseClass
{
protected:
  Array<PairActionFitClass*,2> &PairMatrix;
  Array<PairActionFitClass*,1> &PairArray;
 //  LinearGrid LongGrid;


  /// This calculates the quantity 
  /// \f$ X_k \equiv -\frac{4\pi}{\Omega k} \int_{r_c}^\infty dr \, r \sin(kr) V(r).\f$
  //  double CalcXk (int paIndex, int level, double k, double rc, JobType type);
  // This must be called after all of the OptimizedBreakup_x's

  //  int Level, ki;

  ///The values of the action will be uk*\rho_k*\rho_{-k}
  //Array<double,1> uk;
  //Array<double,1> duk;
public:
  Array<double,2> uk;
  Array<double,2> duk;


  //Array<double,1> Vlong_k;
  Array<double,2> Vlong_k;
  // void Init(IOSectionClass &in);
  Array<double,1> yk_zero;
  void Read (IOSectionClass &in);
  void ReadYk();
  void BuildRPA_SingleSpecies();
  void Build_MultipleSpecies();
  double V(int slice1,int slice2,int level);

  double SingleAction (int slice1, int slice2, 
		       const Array<int,1> &activeParticles, int level);
  double d_dBeta (int slice1, int slice2,  int level);
  string GetName();
  DavidLongRangeClassYk(PathDataClass &pathData,
			Array<PairActionFitClass* ,2> &pairMatrix,
			Array<PairActionFitClass*, 1> &pairArray);

};

#endif
