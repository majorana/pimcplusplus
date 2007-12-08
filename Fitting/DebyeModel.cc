#include "DebyeModel.h"
#include "Fitting.h"

double
DebyeModel::Chi2 (double theta, 
		  Array<double,1> &Fvals,
		  Array<double,1> &Tvals)
{
  SetTheta(theta);
  double chi2 = 0.0;
  for (int i=0; i<Fvals.size(); i++) {
    double diff = F(Tvals(i)) - Fvals(i);
    chi2 += diff *diff;
  }
  return chi2;
}

double
DebyeModel::MinTheta(Array<double,1> &Fvals, Array<double,1> &Tvals,
		      double startTheta, double endTheta)
{
  double minTheta = startTheta;
  double minChi2 = Chi2(minTheta, Fvals, Tvals);
  
  double dTheta = (endTheta - startTheta)/999.0;
  for (double theta = startTheta; theta <= (1.00000001*endTheta); 
       theta += dTheta) {
    double chi2 = Chi2(theta, Fvals, Tvals);
    // cerr << "Theta = " << theta << "  Chi2 = " << chi2 << endl;
    if (chi2 < minChi2) {
      minChi2 = chi2;
      minTheta = theta;
    }
  }
  // cerr << "minChi2 = " << minChi2 << endl;
  return minTheta;
}

double
DebyeModel::OptTheta (Array<double,1> &Fvals, Array<double,1> &Tvals)
{
  double range = 7998.0;
  double mid = 4000.0;

  for (int i=0; i< 5; i++) {
    double start = mid - 0.5 * range;
    double end   = mid + 0.5 * range;
    mid = MinTheta (Fvals, Tvals, start, end);
    range = 0.1 * range;
  }
  return mid;
}

void
DebyeFreeEnergy::AddModel (double Vval, vector<double> &Fvec,
			   vector<double> &Tvec)
{
  DebyeModel debye;
  debye.SetN(2);
  Array<double,1> F(Fvec.size()), T(Tvec.size());
  for (int i=0; i<Fvec.size(); i++) {
    F(i) = Fvec[i] - Fvec[0];
    T(i) = Tvec[i];
  }
  double theta = debye.OptTheta (F, T);
  V.push_back(Vval);
  F0.push_back(Fvec[0]);
  Theta.push_back(theta);
}


void
DebyeFreeEnergy::FitTheta_V (int numCoefs)
{
  int N = V.size();
  Array<double,2> basis(N, numCoefs);
  for (int i=0; i<N; i++) {
    basis(i,0) = 1.0;
    for (int j=1; j<numCoefs; j++)
      basis(i,j) = V[i]*basis(i,j-1);
  }
  ThetaCoefs.resize(numCoefs);
  Array<double,1> sigma(N), theta(N), errors(numCoefs);
  for (int i=0; i<N; i++) {
    sigma(i) = 1.0e-7;
    theta(i) = Theta[i];
  }
  LinFitLU(theta, sigma, basis, ThetaCoefs, errors);
}


double
DebyeFreeEnergy::F(double Vval, double T)
{
  Debye.SetTheta (Theta_V(Vval));
  return Debye.F(T);
}


double
DebyeFreeEnergy::P(double Vval, double T)
{
  Debye.SetTheta (Theta_V(Vval));
  return -dTheta_dV(Vval) * Debye.dF_dTheta(T);
}
  
double
DebyeFreeEnergy::P_FD (double Vval, double Tval)
{
  double eps = 1.0e-6;
  return ((F(Vval-eps,Tval)-F(Vval+eps,Tval))/(2.0*eps));
}  

double
DebyeFreeEnergy::dF_dT (double Vval, double Tval)
{
  Debye.SetTheta (Theta_V(Vval));
  return Debye.dF_dT (Tval);
}

double
DebyeFreeEnergy::dF_dT_FD (double Vval, double Tval)
{
  double eps = 1.0e-8;
  return ((F(Vval, Tval+eps) - F(Vval, Tval-eps))/(2.0*eps));
}
