#ifndef NODAL_ACTION_H
#define NODAL_ACTION_H

#include "PathDataClass.h"

class NodalActionClass
{
protected:
  /// These hold the first and last particle to which this instance of
  /// the nodal action applies.  That is if we have 10 up electrons
  /// and 10 down electrons, then we would have one instance of the
  /// NodalActionClass for the ups and one for the down.  The first
  /// would have FirstPtcl=0 and LastPtcl=9, the second FirstPtcl=10,
  /// LastPtcl=19.
  int SpeciesNum;
  PathDataClass &PathData;
public:
  virtual double Action (int startSlice, int endSlice,
			 const Array<int,1> &changePtcls, int level) = 0;
  //virtual void Update (int slice, int ptcl, dVec newPos, dVec oldPos) = 0;
  NodalActionClass (PathDataClass &pathData, int speciesNum) : 
    PathData(pathData), SpeciesNum (speciesNum)
  {
    
  }
};

class FPNodalActionClass : public NodalActionClass
{
private:
  PathClass &Path;
  Array<double,2> DetMatrix, Cofactors;
  Array<dVec,1> GradVec;
  double Det();
  void GradientDet (int slice, double &det, Array<dVec,1> &gradient);
  void GradientDetFD (int slice, double &det, Array<dVec,1> &gradient);
  double NodalDist (int slice);
public:
  double Action (int startSlice, int endSlice,
		 const Array<int,1> &changedPtcls, int level);
  //void Update (int slice, int ptcl, dVec newPos, dVec oldPos);
  FPNodalActionClass(PathDataClass &pathData, int speciesNum) :
    NodalActionClass (pathData, speciesNum), Path(pathData.Path) 
  {
    int N = Path.Species(speciesNum).LastPtcl - 
      Path.Species(speciesNum).FirstPtcl+1;
    DetMatrix.resize(N,N);
    Cofactors.resize(N,N);
    GradVec.resize(N);
  }
};


#endif
