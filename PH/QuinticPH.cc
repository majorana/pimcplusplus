#include "QuinticPH.h"

bool QuinticPH::IsPH()
{ 
  return true; 
}

double QuinticPH::A(double r)
{
  double val;
  if (r < CoreRadius) {
    double pa = pA(r);
    val = pa*pa + ABmin;
  }
  else
    val = 1.0;
  return val;
}


double QuinticPH::dAdr(double r)
{
  double val;
  if (r < CoreRadius) 
    val = 2.0*pA(r)*pA.Deriv(r);
  else
    val = 0.0;
  return val;
}


double QuinticPH::d2Adr2(double r)
{
  double val;
  if (r < CoreRadius) {
    double dpa = pA.Deriv(r);
    val = 2.0*(dpa*dpa + pA(r)*pA.Deriv2(r));
  }
  else
    val = 0.0;
  return val;
}


double QuinticPH::B(double r)
{
  double val;
  if (r < CoreRadius) {
    double pb = pB(r);
    val = pb*pb + ABmin;
  }
  else
    val = 0.0;
  return val;
}

double QuinticPH::V(double r)
{
  double val;
  if (r < CoreRadius && UseVcore)
    val = Vcore(r);
  else
    val = Vouter->V(r);
  return val;
}


double QuinticPH::dVdr(double r)
{
  double val;
  if (r < CoreRadius && UseVcore)
    val = Vcore.Deriv(r);
  else
    val = Vouter->dVdr(r);
  return val;
}


double QuinticPH::d2Vdr2(double r)
{
  double val;
  if (r < CoreRadius && UseVcore)
    val = Vcore.Deriv2(r);
  else
    val = Vouter->d2Vdr2(r);
  return val;
}


void QuinticPH::Write(IOSectionClass &out)
{
  out.NewSection("Agrid");  Agrid.Write(out);  out.CloseSection();
  out.NewSection("Bgrid");  Bgrid.Write(out);  out.CloseSection();
  out.NewSection("Vgrid");  Vgrid.Write(out);  out.CloseSection();

  out.WriteVar ("ABmin", ABmin);
  out.WriteVar ("pA", pA.Data());
  out.WriteVar ("pB", pB.Data());
  out.WriteVar ("Vcore", Vcore.Data());
}


/// Note:  Vouter must be set before calling "Read".
void QuinticPH::Read(IOSectionClass &in)
{
  assert (in.OpenSection ("Agrid")); Agrid.Read(in); in.CloseSection();
  assert (in.OpenSection ("Bgrid")); Bgrid.Read(in); in.CloseSection();
  assert (in.OpenSection ("Vgrid")); Vgrid.Read(in); in.CloseSection();
  CoreRadius = Agrid.End;

  assert(in.ReadVar("ABmin", ABmin));
  Array<double,1> temp;
  assert (in.ReadVar ("pA", temp));
  pA.Init (&Agrid, temp, 0.0, 0.0, 0.0, 0.0);
  assert (in.ReadVar ("pB", temp));
  pB.Init (&Bgrid, temp, 0.0, 0.0, 0.0, 0.0);
  double dVend  = Vouter->dVdr(CoreRadius);
  double d2Vend = Vouter->d2Vdr2(CoreRadius);
  assert (in.ReadVar ("Vcore", temp));
  Vcore.Init (&Bgrid, temp, NAN, dVend, NAN, d2Vend);

  CoreRadius = Agrid.End;
}

