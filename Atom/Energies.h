#ifndef ENERGIES_H
#define ENERGIES_H

#include "Atom.h"


scalar KEderiv(scalar r, scalar E, void *WFptr);
scalar PEderiv(scalar r, scalar E, void *Atomptr);
scalar HEderiv(scalar r, scalar E, void *Atomptr);
scalar XCEderiv(scalar r, scalar E, void *Atomptr);
scalar NUCderiv(scalar r, scalar E, void *Atomptr);

#endif
