#ifndef PA_FIT_H
#define PA_FIT_H

#include "PAszFit.h"
#include "PAcoulombFit.h"
#include "PAcoulombBCFit.h"
#include "PAclassicalFit.h"

inline PairActionFitClass *ReadPAFit (IOSectionClass &in, 
				      double smallestBeta, int numLevels)
{
  assert (in.OpenSection("Fits"));
  string type;
  assert (in.ReadVar("Type", type));
  in.CloseSection (); // "Fits"
  PairActionFitClass *fit;
  if (type == "szfit")
    fit = new PAszFitClass;
  else if (type == "coulombfit")
    fit = new PAcoulombFitClass;
  else if (type == "coulombBCfit")
    fit = new PAcoulombBCFitClass;
  else if (type == "classical")
    fit = new PAclassicalFitClass;
  else {
    cerr << "Unrecognize pair action fit type \"" 
	 << type << "\".  Exitting.\n";
    exit(1);
  }
  fit->Read(in, smallestBeta, numLevels);
  return (fit);
}


#endif
