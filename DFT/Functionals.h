#ifndef FUNCTIONALS_H
#define FUNCTIONALS_H

void ExchangePotential (double nup, double ndown,
			double &Vup, double &Vdown);
void CorrelationPotential(double  nup, double ndown,
			  double &Vup, double &Vdown);
void FortranExCorr(double  nup, double  ndown,
		   double &Vup, double &Vdown);

double FortranXCE (double nup, double ndown);

#endif
