#ifndef NEW_ATOM_H
#define NEW_ATOM_H

#include "DFTAtom.h"

inline Atom* ReadAtom (IOSectionClass &in)
{
  Atom *atom;
  string type;
  assert (in.ReadVar ("Type", type));
  if (type == "DFT")
    atom = new DFTAtom;
  else {
    cerr << "Unrecognized atom type \"" << type 
	 << "\" in ReadAtom.  Exitting.\n";
    exit(1);
  }
  atom->Read(in);
  return (atom);
}

#endif
