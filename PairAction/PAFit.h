#ifndef PA_FIT_H
#define PA_FIT_H

#include "PAszFit.h"
#include "PAcoulombFit.h"

inline PairActionFitClass *ReadFit (IOSectionClass &in)
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
  else {
    cerr << "Unrecognize pair action fit type \"" 
	 << type << "\".  Exitting.\n";
    exit(1);
  }
  return (fit);
}


#endif
