#ifndef SHORT_RANGE_CLASS_H
#define SHORT_RANGE_CLASS_H

class ShortRangeClass : public ActionBaseClass
{
protected:
  Array<PairActionClass*,1> PairActionVector;
  Array<PairActionClass*,2> PairMatrix;
  Array<bool,1> DoPtcl;
public:
  void Read (IOSectionClass &in);
  double Evaluate (int slice1, int slice2, 
		   const Array<int,1> &activeParticles, int level);
  ShortRangeClass (PathDataClass &PathData) : ActionBaseClass (pathData)
  {
    DoPtcl.resize(PathData.Path.NumParticles());
  }
};

#endif
