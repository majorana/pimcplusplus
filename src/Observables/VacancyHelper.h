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

#ifndef VACANCY_HELPER_H
#define VACANCY_HELPER_H

#include "../PathDataClass.h"
#include "ObservableBase.h"
#include <Common/Blitz.h>
#include "../Common.h"
#include <Common/IO/IO.h>

class VacancyHelperClass
{

private:
  ///This is the set of locations you should compare against to decide
  ///the location of the head and the tail

  PathDataClass &PathData;
public:
  Array<dVec,1> FixedLoc;
  int NumVacancies;
  void Read(IOSectionClass &in);
  void FindVacancy(int slice);
  Array<int,1> VacancyArray;
  
  VacancyHelperClass(PathDataClass &myPathData) :
    PathData(myPathData)
  {

  }
};


#endif 
