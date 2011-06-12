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

#ifndef CONDUCTIVITY_H
#define CONDUCTIVITY_H

#include "ObservableBase.h"
#include <Common/FFT/FFT.h>

/// This class computes the DC conductivity for a system containing
/// electrons and ions.  This ions are assumed to be classical, and
/// hence do not contribute the the conductivity through kinetic
/// terms.  They do, however, contributed through the force they exert
/// on the electrons.  We distinguish electron and ions based on
/// whether lambda for the species is 0 or not.
class ConductivityClass : public ObservableClass
{
private:
  /// The following variables are used to store the sums of the
  /// current-current correlation functions.  They are indexed by the
  /// imaginary time separation, \beta_2 - \beta_1.
  /// The kinetic-kinetic contribution
  Array<dVec,1> RealSumTT;
  Array<dVec,1> FreqSumTT;
  /// The kinetic-potential contribution
  Array<dVec,1> CorrSumTV;
  /// The potential-potential contribution
  Array<dVec,1> CorrSumVV;
  Array<dVec,1> Current_T, Current_Short, Current_Long;
  Array<dVec,1> MyCurrent, TempCurrent;
  Array<int,1> ProcNumLinks;
  /// This FFT is used to do fast convolutions.
  FFT1D FFTtemp;
  void CalcCurrentT();
  void CalcCurrentShort();
  void CalcCurrentLong();
  void AccumulateFast();
  void AccumulateSlow();
  ObservableVecDouble2 RealKineticVar;
  ObservableVecDouble2 FreqKineticVar;
  Array<double,2> WriteArray;
  int M, NumSamples;
public:
  void Accumulate();
  void WriteBlock();
  void Read (IOSectionClass &in);
  ConductivityClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection),
      RealKineticVar("RealKinetic", IOSection, myPathData.Path.Communicator),
      FreqKineticVar("FreqKinetic", IOSection, myPathData.Path.Communicator),
      NumSamples(0)
  {
    // nothing for now
  }

};
#endif
