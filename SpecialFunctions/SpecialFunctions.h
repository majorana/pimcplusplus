#ifndef SPECFUNC_H
#define SPECFUNC_H
#include "../../Common/Blitz.h"

double ModBesselScaled(int l, double x);
double RegModBesselFracScaled(double nu, double x);
double RegModBesselFrac(double nu, double x);
double Legendre (int l, double CosTheta);
double LogRegModBesselScaled(int l, double x);

double il(int l, double z);
double dil_dz (int l, double z);
double d2il_dz2 (int l, double z);

double il_scaled(int l, double z);
double log_il_scaled(int l, double z);
double dil_dz_scaled (int l, double z);
double d2il_dz2_scaled (int l, double z);

void Test_il_Derivs();
void Test_il_series_large();
void Test_il_series_small();
#endif
