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

#ifndef LONG_RANGE_POT_CLASS
#define LONG_RANGE_POT_CLASS

#include "ActionBase.h"
#include <Common/PairAction/PAFit.h>

class LongRangePotClass : public PotentialBaseClass
{
private:
  Array<PairActionFitClass*,2> &PairMatrix;
public:
  bool UseBackground;
  double V (int slice);

  LongRangePotClass (PathDataClass &pathData,
		     Array<PairActionFitClass*,2> &pairMatrix);
};



#endif
