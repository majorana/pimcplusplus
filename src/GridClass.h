#ifndef ONGRIDCLASS_H
#define ONGRIDCLASS_H
#include "Common.h"

#include <list.h>


class PathClass;

class CellInfoClass
{
public:
  dVec left;
  dVec right;
  Array<int,1> MyLoc;
  list<CellInfoClass*> NeighborGrids;
  CellInfoClass()
  {
    MyLoc.resize(3);

  }
  Array<list<int>,1> Particles;
};

class CellMethodClass
{
public:
  int Xeffect;
  int Yeffect;
  int Zeffect;
  Array<TinyVector<int,3>,1> AffectedCells;
  PathClass &Path;
  double CutoffDistance;
  Array<int,1> NumGrid;
  Array<CellInfoClass,3> GridsArray;
  void Init(dVec box,Array<int,1> numGrid);
  bool GridsAffect(CellInfoClass &grid1,CellInfoClass &grid2);
  void BuildNeighborGrids();
  void PrintNeighborGrids();
  void BinParticles(int slice);
  void PrintParticles(int slice);
  void FindBox(dVec myPoint,int &x,int &y,int &z);
  void ReGrid(int slice,int ptcl);
  bool InBox(CellInfoClass &theGrid,dVec thePoint);
CellMethodClass(PathClass &path) : Path(path){
    CutoffDistance=8.0;
    //do nothing for now
}
  


  
};



#endif
