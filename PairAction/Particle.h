#ifndef PARTICLE_H
#define PARTICLE_H

class ParticleClass
{
public:
  string Name;
  double lambda;
  double Charge;
  int Ndim;
  void Write (IOSectionClass &outSection)
  {
    outSection.WriteVar ("Name", Name);
    outSection.WriteVar ("lambda", lambda);
    outSection.WriteVar ("Charge", Charge);
    outSection.WriteVar ("Ndim", Ndim);
  }
  bool Read  (IOSectionClass &inSection)
  {
    bool success;
    success =  inSection.ReadVar ("Name", Name);
    success &= inSection.ReadVar ("lambda", lambda);
    success &= inSection.ReadVar ("Charge", Charge);
    success &= inSection.ReadVar ("Ndim", Ndim);
    return (success);
  }
};

#endif
