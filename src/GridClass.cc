#include "GridClass.h"
#include "PathClass.h"

void CellMethodClass::Init(dVec box,Array<int,1> numGrid)
{
  double actionDistance=8.0;
  cerr<<"Starting my initialization"<<endl;
  NumGrid.resize(3);
  NumGrid(0)=numGrid(0);
  NumGrid(1)=numGrid(1);
  NumGrid(2)=numGrid(2);
  cerr<<"Num grid 1 is"<<numGrid(0)<<endl;
  cerr<<"Num grid 2 is"<<numGrid(1)<<endl;
  cerr<<"Num grid 3 is"<<numGrid(2)<<endl;
  GridsArray.resize(numGrid(0),numGrid(1),numGrid(2));
  double xStart=-box[0]/2.0;
  double yStart=-box[1]/2.0;
  double zStart=-box[2]/2.0;
  double xSize=box[0]/numGrid(0);
  double ySize=box[1]/numGrid(1);
  double zSize=box[2]/numGrid(2);
  Xeffect=(int)(ceil(actionDistance/xSize));
  Yeffect=(int)(ceil(actionDistance/ySize));
  Zeffect=(int)(ceil(actionDistance/zSize));
  cerr<<"My effect is "<<Xeffect<<" "<<Yeffect<<" "<<Zeffect<<endl;
  for (int xCnt=0;xCnt<numGrid(0);xCnt++){
    for (int yCnt=0;yCnt<numGrid(1);yCnt++){
      for (int zCnt=0;zCnt<numGrid(2);zCnt++){
	//	cerr<<"We have "<<xCnt<<" "<<yCnt<<" "<<zCnt<<endl;
	GridsArray(xCnt,yCnt,zCnt).left[0]=xCnt*xSize+xStart;
	GridsArray(xCnt,yCnt,zCnt).right[0]=(xCnt+1)*xSize+xStart;
	GridsArray(xCnt,yCnt,zCnt).left[1]=yCnt*ySize+yStart;
	GridsArray(xCnt,yCnt,zCnt).right[1]=(yCnt+1)*ySize+yStart;
	GridsArray(xCnt,yCnt,zCnt).left[2]=zCnt*zSize+zStart;
	GridsArray(xCnt,yCnt,zCnt).right[2]=(zCnt+1)*zSize+zStart;
	GridsArray(xCnt,yCnt,zCnt).MyLoc(0)=xCnt;
	GridsArray(xCnt,yCnt,zCnt).MyLoc(1)=yCnt;
	GridsArray(xCnt,yCnt,zCnt).MyLoc(2)=zCnt;
	GridsArray(xCnt,yCnt,zCnt).Particles.resize(Path.NumTimeSlices());				      
	//	cerr<<"Done"<<endl;
      }
    }
  }
  cerr<<"Ending my initialization"<<endl;
}
bool CellMethodClass::GridsAffect(CellInfoClass &grid1,CellInfoClass &grid2)
{
  dVec diff1=grid1.left-grid2.right;
  dVec diff2=grid1.right-grid2.left;
  Path.PutInBox(diff1);
  Path.PutInBox(diff2);
  return (
	  (diff1[0]<CutoffDistance || 
	   diff1[1]<CutoffDistance ||
	   diff1[2]<CutoffDistance ||
	   diff2[0]<CutoffDistance ||
	   diff2[1]<CutoffDistance ||
	   diff2[2]<CutoffDistance)
	  );
}


void CellMethodClass::BuildNeighborGrids()
{
  cerr<<"begin"<<endl;
  cerr<<NumGrid(0)<<" "<<NumGrid(1)<<" "<<NumGrid(2)<<endl;
  for (int x=0;x<NumGrid(0);x++){
    for (int y=0;y<NumGrid(1);y++){
      for (int z=0;z<NumGrid(2);z++){
	for (int x2=0;x2<NumGrid(0);x2++){
	  for (int y2=0;y2<NumGrid(1);y2++){
	    for (int z2=0;z2<NumGrid(2);z2++){
	      //	      cerr<<"are we here"<<x<<" "<<y<<" "<<z<<" "<<endl;
	      //	      cerr<<x2<<" "<<y2<<" "<<z2<<endl;
	      if (GridsAffect(GridsArray(x,y,z),GridsArray(x2,y2,z2))){
		//cerr<<"Printing"<<x<<" "<<y<<" "<<z<<endl;
		GridsArray(x,y,z).NeighborGrids.push_back(&GridsArray(x2,y2,z2));
	      }
	    }
	  }
	}
      }
    }
  }
  cerr<<endl;
}

bool CellMethodClass::InBox(CellInfoClass &theGrid,dVec thePoint) 
{
  Path.PutInBox(thePoint);
  return (
  	  (theGrid.left[0]<=thePoint[0] && thePoint[0]<theGrid.right[0]) &&
  	  (theGrid.left[1]<=thePoint[1] && thePoint[1]<theGrid.right[1]) && 
  	  (theGrid.left[2]<=thePoint[2] && thePoint[2]<theGrid.right[2])
  	  );

}
///This needs to be rewritten to be a reasonable speed!!
void CellMethodClass::FindBox(dVec myPoint,int &x,int &y,int &z)
{
  x=0;
  y=0;
  z=0;
  //  cerr<<"Starting"<<endl;
  Path.PutInBox(myPoint);
  //  cerr<<myPoint<<endl;
  dVec box=Path.GetBox();
  //  cerr<<"My box is "<<box<<endl;
  double xStart=-box[0]/2.0;
  double yStart=-box[1]/2.0;
  double zStart=-box[2]/2.0;
  double xSize=box[0]/NumGrid(0);
  double ySize=box[1]/NumGrid(1);
  double zSize=box[2]/NumGrid(2);
  x=(int)floor((myPoint[0]-xStart-0.001)/xSize);
  y=(int)floor((myPoint[1]-yStart-0.001)/ySize);
  z=(int)floor((myPoint[2]-zStart-0.001)/zSize);
  //  cerr<<"Ending"<<endl;
}
  


void CellMethodClass::BinParticles(int slice)
{
  int x,y,z;
  for (int x=0;x<NumGrid(0);x++)
    for (int y=0;y<NumGrid(1);y++)
      for (int z=0;z<NumGrid(2);z++)
	GridsArray(x,y,z).Particles(slice).clear();
  //  cerr<<"My grid is"<<GridsArray.extent(0)<<" "<<GridsArray.extent(1)<<" "<<GridsArray.extent(2)<<endl;
  //  cerr<<"My grid has size "<<GridsArray(0,0,0).Particles.size()<<" "<<slice<<endl;
  for (int ptcl=0;ptcl<Path.NumParticles();ptcl++){
    //    cerr<<"about to find box"<<endl;
    FindBox(Path(slice,ptcl),x,y,z);
    //    cerr<<"I am particle "<<ptcl<<" and have decided to be in grid "<<x<<" "<<y<<" "<<z<<endl;
    GridsArray(x,y,z).Particles(slice).push_back(ptcl);
    //    cerr<<"And pushed back"<<endl;
  }	    
  //  cerr<<"Done"<<endl;
}

void CellMethodClass::PrintParticles(int slice)
{
  for (int x=0;x<NumGrid(0);x++){
    for (int y=0;y<NumGrid(1);y++){
      for (int z=0;z<NumGrid(2);z++){
	cerr<<"I am grid: "<<x<<" "<<y<<" "<<z<<": ";
	for (list<int>::iterator i=GridsArray(x,y,z).Particles(slice).begin();i!=GridsArray(x,y,z).Particles(slice).end();i++){
	  cerr<<*i<<" ";
	}
	cerr<<endl;
      }
    }
  }
}	



void CellMethodClass::PrintNeighborGrids()
{
  for (int x=0;x<NumGrid(0);x++){
    for (int y=0;y<NumGrid(1);y++){
      for (int z=0;z<NumGrid(2);z++){
	cerr<<"I am grid: "<<x<<" "<<y<<" "<<z;
	for (list<CellInfoClass*>::iterator i=GridsArray(x,y,z).NeighborGrids.begin();i!=GridsArray(x,y,z).NeighborGrids.end();i++){
	  cerr<<(*i)->MyLoc<<endl;
	}
      }
    }
  }
}	
								    
			     



