#include "../Blitz.h"

void ExchangePotential (scalar nup, scalar ndown,
			scalar &Vup, scalar &Vdown);
void CorrelationPotential(scalar  nup, scalar ndown,
			  scalar &Vup, scalar &Vdown);
void FortranExCorr(scalar  nup, scalar  ndown,
		   scalar &Vup, scalar &Vdown);

double FortranXCE (scalar nup, scalar ndown);
