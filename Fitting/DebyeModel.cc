#include "DebyeModel.h"

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

