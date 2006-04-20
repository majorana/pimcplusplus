#include "PAFit.h"
#include "../IO/FileExpand.h"

void TestGradients (IOSectionClass &in)
{
  string PAfilename;
  IOSectionClass PAIO;
  assert (in.ReadVar("PAfilename", PAfilename));
  PAfilename = ExpandFileName(PAfilename);
  assert(PAIO.OpenFile(PAfilename));
  double tau;
  assert (in.ReadVar("tau", tau));
  PairActionFitClass &PA = *(ReadPAFit(PAIO, tau, 1));
  Vec3 Rion, r, rp;
  Array<double,1> tmp;
  assert(in.ReadVar("Rion", tmp));
  Rion[0]=tmp(0);   Rion[1]=tmp(1);   Rion[2]=tmp(2); 
  assert(in.ReadVar("r", tmp));
  r[0]=tmp(0);      r[1]=tmp(1);      r[2]=tmp(2); 
  assert(in.ReadVar("rp", tmp));
  rp[0]=tmp(0);      rp[1]=tmp(1);    rp[2]=tmp(2); 
  r = r - Rion;
  rp = rp - Rion;
  double rmag = sqrt(dot(r,r));
  double rpmag = sqrt(dot(rp,rp));
  double q = 0.5*(rmag + rpmag);
  double z = (rmag - rpmag);
  double s = sqrt (dot(r-rp,r-rp));
  cerr << "(q,z,s) = (" << q << ", " << z << ", " << s << ")\n";

  double d_dq, d_dz, d_ds, d_dqFD, d_dzFD, d_dsFD;
  PA.Derivs(q, z, s*s, 0, d_dq, d_dz, d_ds);
  cerr << "Analytic:\n";
  cerr << "  d_dq = " << d_dq << endl;
  cerr << "  d_dz = " << d_dz << endl;
  cerr << "  d_ds = " << d_ds << endl;

  d_dqFD = (PA.U(q+0.0001,z,s*s,0) - PA.U(q-0.0001,z,s*s,0))/0.0002;
  d_dzFD = (PA.U(q,z+0.0001,s*s,0) - PA.U(q,z-0.0001,s*s,0))/0.0002;
  d_dsFD = (PA.U(q,z,(s+0.0001)*(s+0.0001),0) - 
	    PA.U(q,z,(s-0.0001)*(s-0.0001),0))/0.0002;

  cerr << "Finite difference:\n";
  cerr << "  d_dq = " << d_dqFD << endl;
  cerr << "  d_dz = " << d_dzFD << endl;
  cerr << "  d_ds = " << d_dsFD << endl;

  PAIO.CloseFile();
}


main(int argc, char **argv)
{
  if (argc < 2) {
    cerr << "Usage:\n"
	 << "  TestGradients myfile.in\n";
    exit(1);
  }
  else {
    IOSectionClass in;
    assert (in.OpenFile(argv[1]));
    TestGradients(in);
  }
}
