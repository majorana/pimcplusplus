#ifndef MY_RAMAN_MODEL_H
#define MY_RAMAN_MODEL_H

#include "../Blitz.h"

class MyRamanModel
{
private:
  double b0, b1, b2, R0, R1, R2, nu0, nu1, nu2;
public:
  inline void SetParams (TinyVector<double,7> params)
  {
    b0    = params[0];  
    R0    = params[1];  
    R1    = params[2];  
    R2    = params[3];  
    nu0   = params[4];  
    nu1   = params[5];  
    nu2   = params[6];
  }

  inline void GetParams (TinyVector<double,7> &params)
  {
    params[0] = b0 ;  
    params[1] = R0 ;  
    params[2] = R1 ;  
    params[3] = R2 ;  
    params[4] = nu0;  
    params[5] = nu1;  
    params[6] = nu2;
  }

  
  inline TinyVector<double,7> 
  Grad_FD(TinyVector<double,2> PT)
  {
    TinyVector<double,7> save, grad, plus, minus;
    GetParams(save);
    for (int i=0; i<7; i++) {
      double eps = max(1.0e-6 * fabs(save[i]), 1.0e-10);
      plus = minus = save;
      plus[i]  = save[i] + eps;
      minus[i] = save[i] - eps;
      SetParams(plus);
      double nu_plus  = (*this)(PT);
      SetParams(minus);
      double nu_minus = (*this)(PT);
      grad[i] = (nu_plus - nu_minus) / (2.0*eps);
    }
    SetParams(save);
    return grad;
  }

  inline double operator() (TinyVector<double,2> PT)
  {
    double P = PT[0];
    double T = PT[1];

    double nu_0 = nu0 * pow(b0*P/(R0) + 1.0, 1.0/(b0));    
    double nu_1 = nu1 * pow(b0*P/(R1) + 1.0, 1.0/(b0));
    double nu_2 = nu2 * pow(b0*P/(R2) + 1.0, 1.0/(b0));

    return nu_0 + nu_1*exp(-nu_2/T);

    // double b  =  b0 +  b1*T + b2*T*T;
    // double R  =  R0 +  R1*exp(-R2/T);
    // double nubar = nu0 + (nu1+nu1_P*P)*exp(-(nu2+nu2_P)/T);
    // return nubar * pow(b*P/R + 1.0, 1.0/b);
  }

  inline TinyVector<double,7> 
  Grad(TinyVector<double,2> PT)
  {
    TinyVector<double,7> G, GFD;
    double P = PT[0];
    double T = PT[1];

    double nu_0 = nu0 * pow(b0*P/R0 + 1.0, 1.0/b0);    
    double nu_1 = nu1 * pow(b0*P/R1 + 1.0, 1.0/b0);
    double nu_2 = nu2 * pow(b0*P/R2 + 1.0, 1.0/b0);

    double dnu_0_db0 = nu_0 * (-1.0/(b0) * log(b0*P/R0 + 1.0) + (P/b0)/(b0*P+R0));
    double dnu_1_db0 = nu_1 * (-1.0/(b0) * log(b1*P/R1 + 1.0) + (P/b1)/(b1*P+R1));
    double dnu_2_db0 = nu_2 * (-1.0/(b0) * log(b2*P/R2 + 1.0) + (P/b2)/(b2*P+R2));

    double dnu_0_dR0     = -nu0*P/(R0) * pow(b0*P/R0 + 1.0, (1.0/b0-1.0));
    double dnu_1_dR1     = -nu1*P/(R1) * pow(b0*P/R1 + 1.0, (1.0/b0-1.0));
    double dnu_2_dR2     = -nu2*P/(R2) * pow(b0*P/R2 + 1.0, (1.0/b0-1.0));

    G[0] = dnu_0_db0 + dnu_1_db0 + dnu_2_db0;
    // G[1] = dnu_1_db1 * exp(-nu_2/T);
    // G[2] = dnu_2_db2 *(-nu_1/T)*exp(-nu_2/T);
    
    G[1] = dnu_0_dR0;
    G[2] = dnu_1_dR1 * exp(-nu_2/T);
    G[3] = dnu_2_dR2 *(-nu_1/T)*exp(-nu_2/T);

    G[4] = pow(b0*P/R0 + 1.0, 1.0/b0);
    G[5] = pow(b0*P/R1 + 1.0, 1.0/b0) * exp(-nu_2/T); 
    G[6] = pow(b0*P/R2 + 1.0, 1.0/b0) *(-nu_1/T)*exp(-nu_2/T);

    // fprintf (stderr, "\nb_n  = %10.5e %10.5e %10.5e\n", b0, b1, b2);
    // fprintf (stderr, "R_n  = %10.5e %10.5e %10.5e\n", R0, R1, R2);
    // fprintf (stderr, "nu_n = %10.5e %10.5e %10.5e\n", nu0, nu1, nu2);

    // cerr << "\nb  = " << b 
    // 	 << "   R  = " << R 
    // 	 << "   nu = " << nu << endl;
    
    // double dnu_db     = nu * (-1.0/(b*b) * log(b*P/R + 1.0) + (P/b)/(b*P+R));
    // double dnu_dR     = -nubar*P/(R*R) * pow(b*P/R + 1.0, (1.0/b-1.0));
    // double dnu_dnubar = pow(b*P/R + 1.0, 1.0/b);
    

    // fprintf (stderr, "b0     = %12.5e\n", b0);
    // fprintf (stderr, "b1     = %12.5e\n", b1);
    // fprintf (stderr, "b2     = %12.5e\n", b2);
    
    // fprintf (stderr, "R0     = %12.5e\n", R0);
    // fprintf (stderr, "R1     = %12.5e\n", R1);
    // fprintf (stderr, "R2     = %12.5e\n", R2);
    
    // fprintf (stderr, "nu0    = %12.5e\n", nu0);
    // fprintf (stderr, "nu1    = %12.5e\n", nu1);
    // fprintf (stderr, "nu2    = %12.5e\n\n", nu2);
    

     GFD = Grad_FD(PT);
    // fprintf (stderr, "\nGFD = %16.12e %16.12e %12.12e\n", GFD[0], GFD[1], GFD[2]);
    // fprintf (stderr, "G   = %16.12e %16.12e %12.12e\n",   G[0],   G[1],   G[2]);

    // fprintf (stderr, "GFD = %16.12e %16.12e %12.12e\n", GFD[3], GFD[4], GFD[5]);
    // fprintf (stderr, "G   = %16.12e %16.12e %12.12e\n",   G[3],   G[4],   G[5]);

    // fprintf (stderr, "GFD = %16.12e %16.12e %12.12e\n", GFD[6], GFD[7], GFD[8]);
    // fprintf (stderr, "G   = %16.12e %16.12e %12.12e\n",   G[6],   G[7],   G[8]);

    return GFD;
  }
};



#endif
