#include "WaterMoveRing.h"

dVec WaterRotateRing::Rotate(dVec coord,int u1,int u2,int u3,double theta)
{
  double x0 = coord[u1];
  double y0 = coord[u2];
  double costheta = cos(theta);
  double sintheta = sin(theta);
  dVec new_coord;
  new_coord[u1] = x0*costheta - y0*sintheta;
  new_coord[u2] = y0*costheta + x0*sintheta;
  new_coord[u3] = coord[u3];
  return new_coord;

}

dVec WaterTranslateRing::Translate(double epsilon)
{
	double x = (PathData.Path.Random.Local()-0.5)*2*epsilon;	
	double y = (PathData.Path.Random.Local()-0.5)*2*epsilon;	
	double z = (PathData.Path.Random.Local()-0.5)*2*epsilon;	
	dVec translate;
	translate[0] = x;
	translate[1] = y;
	translate[2] = z;
	return translate;
}

void WaterTranslateRing::Molecule2Atoms(int moleculeNum)
{

  coord_loc.resize(5);
  int speciesO = PathData.Path.SpeciesNum("O");
  int speciese = PathData.Path.SpeciesNum("e");
  int speciesp = PathData.Path.SpeciesNum("p");
  coord_loc(0)=PathData.Path.Species(speciesO).FirstPtcl+moleculeNum;
  coord_loc(1)=PathData.Path.Species(speciesp).FirstPtcl+moleculeNum;
  int numProtons=PathData.Path.Species(speciesp).LastPtcl-
    PathData.Path.Species(speciesp).FirstPtcl+1;
  coord_loc(2)=PathData.Path.Species(speciesp).FirstPtcl+moleculeNum+numProtons/2;
  coord_loc(3)=PathData.Path.Species(speciese).FirstPtcl+moleculeNum;
  int numElectrons=PathData.Path.Species(speciese).LastPtcl-
    PathData.Path.Species(speciese).FirstPtcl+1;
  coord_loc(4)=PathData.Path.Species(speciese).FirstPtcl+moleculeNum+numElectrons/2;
/*  cerr<<"The values of coord_loc are "<<coord_loc(0)<<" ";
  cerr<<coord_loc(1)<<" ";
  cerr<<coord_loc(2)<<" ";
  cerr<<coord_loc(3)<<" ";
  cerr<<coord_loc(4)<<" ";
  cerr<<endl;
  cerr<<"leave function"<<endl;
*/
}


void WaterRotateRing::Molecule2Atoms(int moleculeNum)
{
  coord_loc.resize(5);
  int speciesO = PathData.Path.SpeciesNum("O");
  int speciese = PathData.Path.SpeciesNum("e");
  int speciesp = PathData.Path.SpeciesNum("p");
  coord_loc(0)=PathData.Path.Species(speciesO).FirstPtcl+moleculeNum;
  coord_loc(1)=PathData.Path.Species(speciesp).FirstPtcl+moleculeNum;
  int numProtons=PathData.Path.Species(speciesp).LastPtcl-
    PathData.Path.Species(speciesp).FirstPtcl+1;
  coord_loc(2)=PathData.Path.Species(speciesp).FirstPtcl+moleculeNum+numProtons/2;
  coord_loc(3)=PathData.Path.Species(speciese).FirstPtcl+moleculeNum;
  int numElectrons=PathData.Path.Species(speciese).LastPtcl-
    PathData.Path.Species(speciese).FirstPtcl+1;
  coord_loc(4)=PathData.Path.Species(speciese).FirstPtcl+moleculeNum+numElectrons/2;
/*  cerr<<"The values water of coord_loc are "<<coord_loc(0)<<" ";
  cerr<<coord_loc(1)<<" ";
  cerr<<coord_loc(2)<<" ";
  cerr<<coord_loc(3)<<" ";
  cerr<<coord_loc(4)<<" ";
  cerr<<endl;
  cerr<<"leave function"<<endl;
*/

}

void WaterTranslateRing::MakeMove()
{
  double step = 0.2;
  int speciesO = PathData.Path.SpeciesNum("O");
  int speciesp = PathData.Path.SpeciesNum("p");
  int speciese = PathData.Path.SpeciesNum("e");
  int  numWater=PathData.Path.Species(speciesO).LastPtcl-
    PathData.Path.Species(speciesO).FirstPtcl+1;
  int choosemol = (int)floor(PathData.Path.Random.Local()*numWater);
  
// set up time slice parameters
  int numSlices = PathData.Path.TotalNumSlices;
  int endSlice = numSlices;
  int startSlice = 0;
/*  if(numSlices>1){
    int P_max = numSlices - 1;
    slice = (int)floor(P_max*PathData.Path.Random.Local()) + 1;
    startSlice = slice-1;
    endSlice = slice+1;
  }
*/


  Molecule2Atoms(choosemol);

  for (int i=0;i<coord_loc.size();i++){
    ActiveParticles(i)=coord_loc(i);
//    cerr<<"My coord_loc is "<<i<<" "<<coord_loc(i)<<endl;
  }
//  cerr << "OLD ";
  double oldAction = PathData.Actions.TIP5PWater.Action(startSlice,endSlice,ActiveParticles,0);
  oldAction += PathData.Actions.Kinetic.Action(startSlice,endSlice,ActiveParticles,0); 
  
  dVec move = Translate(step); 
  double move_mag = sqrt(move(0)*move(0) + move(1)*move(1) + move(2)*move(2));
  for(int i=0;i<coord_loc.size();i++){
//  So here loop over time slices
    for(int slice = 0; slice < numSlices; slice++){
      dVec old_coord=PathData.Path(slice,coord_loc(i));
      dVec new_coord = old_coord + move; 
//    cerr<<"New Coord is "<<new_coord<<endl;
      PathData.Path.SetPos(slice,coord_loc(i),new_coord);
    }
  }
 
//  cerr << "NEW " ; 
  double newAction = PathData.Actions.TIP5PWater.Action(startSlice,endSlice,ActiveParticles,0);
  newAction += PathData.Actions.Kinetic.Action(startSlice,endSlice,ActiveParticles,0); 
//  cerr<<"TRANSLATE:  The actions are "<<newAction<<" "<<oldAction<<endl;
  if (-(newAction-oldAction)>=log(PathData.Path.Random.Local())){
    PathData.AcceptMove(startSlice,endSlice,ActiveParticles);
    NumAccepted++;
    total_r_squared += move_mag*move_mag;
//    cerr << move << " has magnitude " << move_mag << " and mag squared " << move_mag*move_mag << endl;
//    cerr<<"TRANSLATE:  I've accepted " << NumAccepted << " " << NumMoves+1 <<endl;
  }
  else {
//    cerr<<"TRANSLATE:  I've rejected"<<endl;
    PathData.RejectMove(startSlice,endSlice,ActiveParticles);
  }
  NumMoves++;
  if (NumMoves%10000 == 0){
    cerr << NumMoves << " moves; current ring translate ratio is " << double(NumAccepted)/NumMoves << " with step size " << step << endl;
    double avg = sqrt(total_r_squared)/NumAccepted;
    double diff = total_r_squared/NumMoves;
    cerr << "TRANSLATE RING diffusion value is " << diff << " avg step is " << avg << endl;
  }
}

void WaterRotateRing::MakeMove()
{
  cerr << "ROTATE RING" << endl;
  double theta_step = 0.5;
  double dtheta = 2*M_PI*theta_step;
  int speciesO = PathData.Path.SpeciesNum("O");
  int  numWater=PathData.Path.Species(speciesO).LastPtcl-
    PathData.Path.Species(speciesO).FirstPtcl+1;
  int choosemol = (int)floor(PathData.Path.Random.Local()*numWater);

// choose a time slice to move
  int numSlices = PathData.Path.TotalNumSlices;
  int endSlice = numSlices;
  int startSlice = 0;
/*  if(numSlices>1){
    int P_max = numSlices - 1;
    slice = (int)floor(P_max*PathData.Path.Random.Local()) + 1;
    startSlice = slice-1;
    endSlice = slice+1;
  }
*/
  
  Molecule2Atoms(choosemol);

  for (int i=0;i<coord_loc.size();i++){
    ActiveParticles(i)=coord_loc(i);
  }

/*  these lines are just here for testing the rotation move ****
  int oxy = 0;
  int check = 3;
  double old_dist;
  dVec r;
  PathData.Path.DistDisp(slice,coord_loc(oxy),coord_loc(check),old_dist,r);
  cerr << "BEFORE " << check << " " << PathData.Path(slice,coord_loc(check)) << old_dist << " ";
*************************************************************/
//  cerr << "OLD " ;
  double oldAction = PathData.Actions.TIP5PWater.Action(startSlice,endSlice,ActiveParticles,0);
  oldAction += PathData.Actions.Kinetic.Action(startSlice,endSlice,ActiveParticles,0); 

  int x,y;
  int z = (int)floor(3*PathData.Path.Random.Local());
  double theta = (PathData.Path.Random.Local()-0.5)*dtheta;
  if (z == 0){
    x = 1;
    y = 2;
  }
  else if (z == 1){
    x = 2;
    y = 0;
  }
  else if (z == 2){
    x = 0;
    y = 1;
  }
//  cerr << "ROTATE " << theta/M_PI << " pi rad about " << z << endl;
  // error statement?
  for(int i=1;i<coord_loc.size();i++){
//  Loop over slices
    for(int slice=0;slice<endSlice;slice++){
      dVec old_coord=PathData.Path(slice,coord_loc(i));
      old_coord=old_coord-PathData.Path(slice,coord_loc(0));
      dVec new_coord = Rotate(old_coord,x,y,z,theta);
      new_coord=new_coord+PathData.Path(slice,coord_loc(0));
      PathData.Path.SetPos(slice,coord_loc(i),new_coord);
    }
  }
  
/*  these lines are just here for testing the rotation move ****
  double new_dist;
  dVec r2;
  PathData.Path.DistDisp(slice,coord_loc(oxy),coord_loc(check),new_dist,r2);
  cerr << "AFTER " << check << " " << PathData.Path(slice,coord_loc(check)) << old_dist << endl;
 *************************************************************/
//  cerr << "NEW " ;
  double newAction = PathData.Actions.TIP5PWater.Action(startSlice,endSlice,ActiveParticles,0);
  newAction += PathData.Actions.Kinetic.Action(startSlice,endSlice,ActiveParticles,0); 
//  cerr<<"ROTATE:  The actions are "<<newAction<<" "<<oldAction<<endl;
 
  if (-(newAction-oldAction)>=log(PathData.Path.Random.Local())){
    PathData.AcceptMove(startSlice,endSlice,ActiveParticles);
  //  cerr<<"ROTATE:  I've accepted " << NumAccepted << " " << NumMoves+1 <<endl;
    NumAccepted++;
    total_r_squared += 2*(1-cos(theta));
  }
  else {
    PathData.RejectMove(startSlice,endSlice,ActiveParticles);
//    cerr<<"ROTATE:  I've rejected"<<endl;
  }
  NumMoves++;
  if (NumMoves%10000 == 0){
    cerr << NumMoves << " moves; current ring rotate ratio is " << double(NumAccepted)/NumMoves << " with angle size " << dtheta << " (" << theta_step << ")" << endl;
    double avg = sqrt(total_r_squared)/NumAccepted;
    double diff = total_r_squared/NumMoves;
    cerr << "ROTATE RING diffusion value is " << diff << " avg angle is " << avg << endl;
  }
}
