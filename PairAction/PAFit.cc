#include "PAFit.h"

void PAszFitClass::WriteBetaIndependentInfo (IOSectionClass &outSection)
{
  outSection.WriteVar ("Type", "szfit");
  outSection.WriteVar ("Order", Order);
  outSection.NewSection("qGrid");
  grid->Write(outSection);
  outSection.CloseSection();
}


void PAszFitClass::ReadParams(IOSectionClass &inSection)
{
  assert(inSection.OpenSection ("qGrid"));
  grid = ReadGrid (inSection);
  GridIsMine = true;
  assert(inSection.ReadVar ("Order", Order));
  UsePBC = inSection.ReadVar ("Box", Box);
}

void PAszFitClass::WriteFit (IOSectionClass &outSection, Rho &rho)
{
  Array<double,2> Coefs (grid->NumPoints, Order);
  Coefs = 0.0;
  outSection.WriteVar ("Coefs", Coefs);
}

bool PAszFitClass::Read (IOSectionClass &inSection,
			 double lowestBeta, int NumBetas)
{

}

double PAszFitClass::U(double r, double rp, double costheta, int level)
{
  return (0.0);
}


// void Rho::WriteFit (string fileName)
// {
//   FILE *fout = fopen (fileName.c_str(), "w");
//   FILE *fpout = fopen ("var.dat", "w");
//   assert (fout != NULL);
//   int lmax = U_ls.size()-1;
//   LinearGrid qgrid (1e-4, 7.5, 101);
  
//   // HACK
//   double myBeta = 0.125;
//   //double myBeta = Beta();

//   for (int qindex=0; qindex<qgrid.NumPoints; qindex++)
//     {
//       cerr << " q = " << qindex << endl;
//       double q = qgrid(qindex);
//       double smax = min (2.0*q, sqrt(-4.0*lambda*myBeta*log(1.0e-4)));
//       LinearGrid zgrid(-smax, smax, 101);
//       LinearGrid sgrid(0.0, smax, 101);
//       for (int zindex=0; zindex<zgrid.NumPoints; zindex++)
// 	{
// 	  double z = zgrid(zindex);
// 	  double r = q + 0.5*z;
// 	  double rp = q - 0.5*z;
// 	  if (r <= 0.0) r = 1.0e-5;
// 	  if (rp <= 0.0) rp = 1.0e-5;
// 	  double x = Transform.r2x(r);
// 	  double xp = Transform.r2x(rp);
// 	  Array<double,1> Ulvec(lmax+1), dUlvec(lmax+1);
// 	  // First, calculate the diagonal parts
// 	  U_lArray(r,r, Ulvec, dUlvec);
// 	  double Ur, Urp, dUr, dUrp;
// 	  UdU(r,r,1.0,Ulvec, dUlvec, Ur, dUr);
// 	  U_lArray(rp,rp, Ulvec, dUlvec);
// 	  UdU(r,r,1.0,Ulvec, dUlvec, Ur, dUr);
// // 	  for (int l=0; l<=lmax; l++)
// // 	    Ularray(l) = U_ls(l).U(x,x);
// // 	  double Udiag_r = U(r,r,1.0, Ularray);

// // 	  for (int l=0; l<=lmax; l++)
// // 	    Ularray(l) = U_ls(l).U(xp,xp);
// // 	  double Udiag_rp = U(rp,rp,1.0, Ularray);

// 	  // Now calculate off-diagonal elements
// // 	  for (int l=0; l<=lmax; l++)
// // 	    Ularray(l) = U_ls(l).U(x,xp);
// 	  U_lArray(r,rp, Ulvec, dUlvec);
// 	  //LinearGrid sgrid(fabs(r-rp)+1e-9, r+rp-1e-9, 101);
// 	  for (int sindex=0; sindex<sgrid.NumPoints; sindex++)
// 	    {
// 	      double s = sgrid(sindex);	      
// 	      double costheta;
// 	      if ((r<=0.0) || (rp<=0.0))
// 		costheta = 1.0;
// 	      else
// 		costheta = (r*r + rp*rp - s*s)/(2.0*r*rp);
// 	      //if (costheta > 1.0)  costheta = 1.0;
// 	      //if (costheta < -1.0) costheta = -1.0;
// 	      double theta = acos (costheta);
// 	      double Uval, dUval, Ufitval, dUfitval, Ufit2val, dUfit2val;
// 	      if (isnan(theta))
// 		{
// // 		  cerr << "NAN in theta.\n"
// // 		       << "costheta = " << costheta << endl
// // 		       << "r = " << r << " rp = " << rp << endl 
// // 		       << "q = " << q << " z = " << z << " s = " << s << endl;
// 		  Uval = NAN;	   dUval = NAN;
// 		  Ufitval = NAN;   dUfitval = NAN;
// 		  Ufit2val = NAN;  dUfit2val = NAN;
// 		}
// 	      else
// 		{
// 		  if ((r < grid->End) && (rp < grid->End)) {
// 		    UdU(r,rp,costheta, Ulvec, dUlvec, Uval, dUval);
// 		    //Uval = U(r,rp,costheta, Ularray);
// 		    //Uval -= 0.5*(Udiag_r + Udiag_rp);
// 		    Ufitval = Ufit(q, s, z);
// 		    //Ufitval -= 0.5*(Udiag_r + Udiag_rp);
// 		    Ufit2val = Ufit2(q,s,z);
// 		    dUfit2val = Ufit2.dU(q,s,z);
// 		    //Ufit2val -= 0.5*(Udiag_r + Udiag_rp);
// 		  }
// 		  else {
// 		    cerr << "r or rp outside grid.\n";
// 		    Uval = 0.0;
// 		  }
// 		}
// 	      fprintf (fout, "%1.14e %1.14e %1.14e %1.14e %1.14e %1.14e\n", 
// 		       Uval, Ufitval, Ufit2val, dUval, dUfitval, dUfit2val);
// 	      fprintf (fpout, "%1.6e %1.6e %1.6e %1.8e\n", q, z, s,
// 		       exp(-s*s/(2.0*beta)));
// 	    }
// 	}
//     }
//   fclose (fout);
//   fclose(fpout);
// }
