#include "WaterMove.h"

//ArbitraryRotate
// rotates a coordinate, coord, a specified angle, phi, about a specified axis, axis.  axis and coord are specified as Cartesian coordinates with respect to a common origin; axis should be a unit vector
dVec WaterRareRotate::ArbitraryRotate(dVec axis,dVec coord, double phi){
  dVec coord_perp;
  dVec coord_aligned;
  Strip(axis,coord,coord_aligned,coord_perp);
  double perp_mag = PathData.Actions.TIP5PWater.Mag(coord_perp);
//cerr << "stripped vector " << coord << " off " << axis << " to give aligned " << coord_aligned << " and perp " << coord_perp << endl;
  coord =  PathData.Actions.TIP5PWater.Normalize(coord_perp);
  double cosphi = cos(phi);
  double sinphi = sin(phi);
  double gamma = coord(1)/coord(2) - axis(1)/axis(2);
  double eta = coord(0)/coord(2) - axis(0)/axis(2);
  double x = (cosphi/(gamma*coord(2)) - sinphi*axis(2)/coord(0))/(coord(1)/coord(0) + eta/gamma);
  double y = coord(1)*x/coord(0) + sinphi*axis(2)/coord(0);
  double z = -axis(0)*x/axis(2) - axis(1)*y/axis(2);
  dVec newcoord;
  newcoord(0) = x;
  newcoord(1) = y;
  newcoord(2) = z;
//cerr << "then new coord_perp is " << newcoord << endl;
  return (PathData.Actions.TIP5PWater.Scale(newcoord,perp_mag) + coord_aligned);
}

void WaterRareRotate::Strip(dVec R, dVec u,dVec &aligned, dVec &perp){
  //double phi = PathData.Actions.TIP5PWater.GetAngle(R,u);
  //aligned = PathData.Actions.TIP5PWater.Scale(R,PathData.Actions.TIP5PWater.Mag(u)*cos(phi));
  aligned = PathData.Actions.TIP5PWater.Scale(R,PathData.Actions.TIP5PWater.dotprod(R,u,1.0));
  perp = u - aligned;
}

int WaterRareRotate::FindSlot(double angularmag, int flag){
  double ratio;
  if (flag ==1)
    ratio = numSlots*angularmag/(sqrt(3.0)*2*M_PI);
  else if (flag ==2)
    ratio = numSlots*angularmag/(2*M_PI);
  int index = 0;
  while (index < ratio)
    index ++;
  return index; 
}

void WaterRareRotate::IntegrityCheck(int slice, Array<int,1> activeParticles){
  // constants and tolerances for comparison
  double Amargin = 0.0001;
  double Lmargin = 0.0001;
  double TetAngle = 1.9106;
  double Hlength = 1.0;
  double Elength = 0.8;
  bool pass = true;

  // get vectors WRT origin (oxygen) to determine molecule's geometry
  Array<dVec,1> V;
  V.resize(activeParticles.size());
  for (int i = 0; i < activeParticles.size(); i++){
    V(i) = PathData.Path(slice,activeParticles(i));
    if (i > 0){
      V(i) -= V(0);
    }
  }

  // test geometry
  // check bond lengths
  for (int i = 1; i < activeParticles.size(); i++){
    double length = PathData.Actions.TIP5PWater.Mag(V(i));
    double complength;
    if (i == 1 || i == 2)
      complength = Hlength;
    else if (i == 3 || i == 4)
      complength = Elength;
    if (abs(complength - length) > Lmargin){
      pass = false;
      cerr << activeParticles(0) << ": Integrity Check Failed." << endl;
      cerr << "  Bond length for ptcl " << activeParticles(i) << " " << length << " differs from " << complength << endl;
    }
  }

  // check angular separation 
  for (int i = 1; i < activeParticles.size(); i++){
    for (int j = i+1; j < activeParticles.size(); j++){
      double theta = PathData.Actions.TIP5PWater.GetAngle(V(i),V(j));
      if (abs(TetAngle - theta) > Amargin){
        pass = false;
        cerr << activeParticles(0) << ": Integrity Check Failed." << endl;
        cerr << "  Angle between  ptcls " << activeParticles(i) << "  and " << activeParticles(j) << " " << theta << " differs from " << TetAngle << endl;
      }
    }
  }
}

dVec WaterRotate::Rotate(dVec coord,int u1,int u2,int u3,double theta)
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

dVec WaterRareRotate::Rotate(dVec coord,int u1,int u2,int u3,double theta)
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

dVec WaterTranslate::Translate(double epsilon)
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

void 
WaterRotate::Read(IOSectionClass &moveInput)
{
  string typeCheck;
  assert(moveInput.ReadVar("type",typeCheck));
  assert(typeCheck=="WaterRotate");
  assert(moveInput.ReadVar("name",Name));
  assert(moveInput.ReadVar("step",Theta));  
}


/*
// MODIFIED VERSION
void WaterTranslate::Molecule2Atoms(int moleculeNum)
{

  coord_loc.resize(3);
  int speciesO = PathData.Path.SpeciesNum("O");
  int speciese = PathData.Path.SpeciesNum("e");
  int speciesp = PathData.Path.SpeciesNum("p");
  coord_loc(0)=PathData.Path.Species(speciesO).FirstPtcl+moleculeNum;
  coord_loc(1)=PathData.Path.Species(speciesp).FirstPtcl+moleculeNum;
  coord_loc(2)=PathData.Path.Species(speciese).FirstPtcl+moleculeNum;

  //cerr<<"The values of coord_loc are "<<coord_loc(0)<<" ";
  //cerr<<coord_loc(1)<<" ";
  //cerr<<coord_loc(2)<<" ";
  //cerr<<endl;
  //cerr<<"leave function"<<endl;

}



void WaterRotate::Molecule2Atoms(int moleculeNum)
{

  coord_loc.resize(3);
  int speciesO = PathData.Path.SpeciesNum("O");
  int speciese = PathData.Path.SpeciesNum("e");
  int speciesp = PathData.Path.SpeciesNum("p");
  coord_loc(0)=PathData.Path.Species(speciesO).FirstPtcl+moleculeNum;
  coord_loc(1)=PathData.Path.Species(speciesp).FirstPtcl+moleculeNum;
  coord_loc(2)=PathData.Path.Species(speciese).FirstPtcl+moleculeNum;

  //cerr<<"The values of coord_loc are "<<coord_loc(0)<<" ";
  //cerr<<coord_loc(1)<<" ";
  //cerr<<coord_loc(2)<<" ";
  //cerr<<endl;
  //cerr<<"leave function"<<endl;
}
*/

//ORIGINAL

void WaterTranslate::Molecule2Atoms(int moleculeNum)
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
  //cerr<<"The values of coord_loc are "<<coord_loc(0)<<" ";
  //cerr<<coord_loc(1)<<" ";
  //cerr<<coord_loc(2)<<" ";
  //cerr<<coord_loc(3)<<" ";
  //cerr<<coord_loc(4)<<" ";
  //cerr<<endl;
  //cerr<<"leave function"<<endl;
}


//ORIGINAL

void WaterRotate::Molecule2Atoms(int moleculeNum)
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
  //cerr<<"The values water of coord_loc are "<<coord_loc(0)<<" ";
  //cerr<<coord_loc(1)<<" ";
  //cerr<<coord_loc(2)<<" ";
  //cerr<<coord_loc(3)<<" ";
  //cerr<<coord_loc(4)<<" ";
  //cerr<<endl;
  //cerr<<"leave function"<<endl;
}

void WaterRareRotate::Molecule2Atoms(int moleculeNum)
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
}

void WaterTranslate::MakeMove()
{
  double step = 0.08;
  //double step = Step; // Using Step from input file
  int speciesO = PathData.Path.SpeciesNum("O");
  int speciesp = PathData.Path.SpeciesNum("p");
  int speciese = PathData.Path.SpeciesNum("e");
  int  numWater=PathData.Path.Species(speciesO).LastPtcl-
    PathData.Path.Species(speciesO).FirstPtcl+1;
  int choosemol = (int)floor(PathData.Path.Random.Local()*numWater);
  // We want to evaluate the kinetic action only between oxygens (i.e. COMs)
  Array<int,1> OActiveParticles;
  OActiveParticles.resize(1);
  
  // choose a time slice to move
  int numSlices = PathData.Path.TotalNumSlices;
//cerr << "numSlices is " << numSlices;
  int slice=0;
  int endSlice = 0;
  int startSlice = 0;
  if(numSlices>1){
    int P_max = numSlices - 1;
    slice = (int)floor(P_max*PathData.Path.Random.Local()) + 1;
    startSlice = slice-1;
    endSlice = slice+1;
//cerr << ".  We chose slice " << slice << endl;
  }

  Molecule2Atoms(choosemol);

  for (int i=0;i<coord_loc.size();i++){
    ActiveParticles(i)=coord_loc(i);
    //cerr<<"My coord_loc is "<<i<<" "<<coord_loc(i)<<endl;
  }
//cerr << "Assigned ActiveParticles " << ActiveParticles << endl;
  OActiveParticles(0) = ActiveParticles(0);

  double oldAction = 0.0;
// oldAction += PathData.Actions.TIP5PWater.Action(slice,slice,ActiveParticles,0);
 oldAction += PathData.Actions.ST2Water.Action(slice,slice,ActiveParticles,0);
 oldAction += PathData.Actions.Kinetic.Action(startSlice,endSlice,OActiveParticles,0); 
  dVec move = Translate(step); 
  double move_mag = sqrt(move(0)*move(0) + move(1)*move(1) + move(2)*move(2));
//cerr << "going to move " << move << " with magnitude " << move_mag << endl;
  for(int i=0;i<coord_loc.size();i++){
    dVec old_coord=PathData.Path(slice,coord_loc(i));
    dVec new_coord = old_coord + move; 
//cerr<< i << "; moved "<<old_coord << " to " << new_coord<<endl;
    PathData.Path.SetPos(slice,coord_loc(i),new_coord);
  }
 
  double newAction = 0.0;
// newAction += PathData.Actions.TIP5PWater.Action(slice,slice,ActiveParticles,0);
 newAction += PathData.Actions.ST2Water.Action(slice,slice,ActiveParticles,0);
 newAction += PathData.Actions.Kinetic.Action(startSlice,endSlice,OActiveParticles,0); 
//cerr << " and with kinetic the action is " << newAction << endl;
  //cerr<<"TRANSLATE:  The actions are "<<newAction<<" "<<oldAction<<endl;
  if (-(newAction-oldAction)>=log(PathData.Path.Random.Local())){
    PathData.AcceptMove(startSlice,endSlice,ActiveParticles);
    NumAccepted++;
    total_r_squared += move_mag*move_mag;
    //cerr << move << " has magnitude " << move_mag << " and mag squared " << move_mag*move_mag << endl;
    //cerr<<"TRANSLATE:  I've accepted " << NumAccepted << " " << TimesCalled+1 <<endl;
  }
  else {
   //cerr<<"TRANSLATE:  I've rejected"<<endl;
    PathData.RejectMove(startSlice,endSlice,ActiveParticles);
  }
  TimesCalled++;
  if (TimesCalled%10000 == 0){
    cerr << TimesCalled << " moves; current translate ratio is " << double(NumAccepted)/TimesCalled << " with step size " << step << endl;
    double avg = sqrt(total_r_squared)/NumAccepted;
    double diff = total_r_squared/TimesCalled;
    cerr << "TRANSLATE diffusion value is " << diff << " avg step is " << avg << endl;
  }
}

void WaterRotate::MakeMove()
{
  //double dtheta = 2*M_PI*Theta; // Using Theta from input file
  double dtheta = 2*M_PI*0.05;
  int speciesO = PathData.Path.SpeciesNum("O");
  int  numWater=PathData.Path.Species(speciesO).LastPtcl-
    PathData.Path.Species(speciesO).FirstPtcl+1;
  int choosemol = (int)floor(PathData.Path.Random.Local()*numWater);
  Array<int,1> HActiveParticles;
  HActiveParticles.resize(1);
  Array<int,1> OActiveParticles;
  OActiveParticles.resize(1);
  Array<int,1> p2ActiveParticles;
  p2ActiveParticles.resize(1);

// choose a time slice to move
  int numSlices = PathData.Path.TotalNumSlices;
  int slice=0;
  int endSlice = 0;
  int startSlice = 0;
  if(numSlices>1){
    int P_max = numSlices - 1;
    slice = (int)floor(P_max*PathData.Path.Random.Local()) + 1;
    startSlice = slice-1;
    endSlice = slice+1;
  }
  
  Molecule2Atoms(choosemol);

  for (int i=0;i<coord_loc.size();i++){
    ActiveParticles(i)=coord_loc(i);
  }
  OActiveParticles(0) = ActiveParticles(0);
  HActiveParticles(0) = ActiveParticles(1);
  //HActiveParticles(1) = ActiveParticles(2);
  p2ActiveParticles(0) = ActiveParticles(2);

/*  these lines are just here for testing the rotation move ****
  int oxy = 0;
  int check = 3;
  int check2 = 4;
  double old_dist;
  dVec r;
  PathData.Path.DistDisp(slice,coord_loc(oxy),coord_loc(check),old_dist,r);
cerr << "BEFORE " << check << " " << PathData.Path(slice,coord_loc(check)) << old_dist << " ";
dVec v1 = PathData.Path(slice,coord_loc(check)) - PathData.Path(slice,coord_loc(0));
dVec v2 = PathData.Path(slice,coord_loc(check2)) - PathData.Path(slice,coord_loc(0));
cerr << "test angle is " << PathData.Actions.TIP5PWater.GetAngle(v1,v2) << endl;
***********************************************************/
  double oldAction = 0.0;
//  oldAction += PathData.Actions.TIP5PWater.Action(slice,slice,ActiveParticles,0);
  oldAction += PathData.Actions.ST2Water.Action(slice,slice,ActiveParticles,0);
//  oldAction += PathData.Actions.Kinetic.Action(startSlice,endSlice,HActiveParticles,0); 
//  oldAction += PathData.Actions.TIP5PWater.RotationalKinetic(startSlice,endSlice,HActiveParticles,0);
//  oldAction += PathData.Actions.TIP5PWater.ProtonKineticAction(startSlice,endSlice,HActiveParticles,0);
//  oldAction += PathData.Actions.TIP5PWater.SecondProtonKineticAction(startSlice,endSlice,p2ActiveParticles,0);
//  oldAction += PathData.Actions.TIP5PWater.NewRotKinAction(startSlice,endSlice,HActiveParticles,0);
  oldAction += PathData.Actions.TIP5PWater.FixedAxisAction(startSlice,endSlice,HActiveParticles,0);
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
//cerr << "ROTATE " << theta/M_PI << " pi rad about " << z << endl;
  // error statement?
  for(int i=1;i<coord_loc.size();i++){
//cerr << "MOVING PTCL " << coord_loc(i) << " which is of species " <<  PathData.Path.ParticleSpeciesNum(coord_loc(i)) << endl;
    dVec old_coord=PathData.Path(slice,coord_loc(i));
    old_coord=old_coord-PathData.Path(slice,coord_loc(0));
    dVec new_coord = Rotate(old_coord,x,y,z,theta);
    new_coord=new_coord+PathData.Path(slice,coord_loc(0));
    PathData.Path.SetPos(slice,coord_loc(i),new_coord);
  }
  
/*  these lines are just here for testing the rotation move ****
  double new_dist;
  dVec r2;
  PathData.Path.DistDisp(slice,coord_loc(oxy),coord_loc(check),new_dist,r2);
  cerr << "AFTER " << check << " " << PathData.Path(slice,coord_loc(check)) << old_dist ;
v1 = PathData.Path(slice,coord_loc(check)) - PathData.Path(slice,coord_loc(0));
v2 = PathData.Path(slice,coord_loc(check2)) - PathData.Path(slice,coord_loc(0));
cerr << "test angle is " << PathData.Actions.TIP5PWater.GetAngle(v1,v2) << endl;
 *************************************************************/
  double newAction = 0.0;
//  newAction += PathData.Actions.TIP5PWater.Action(slice,slice,ActiveParticles,0);
  newAction += PathData.Actions.ST2Water.Action(slice,slice,ActiveParticles,0);
//  newAction += PathData.Actions.Kinetic.Action(startSlice,endSlice,HActiveParticles,0); 
//  newAction += PathData.Actions.TIP5PWater.RotationalKinetic(startSlice,endSlice,HActiveParticles,0);
//  newAction += PathData.Actions.TIP5PWater.ProtonKineticAction(startSlice,endSlice,HActiveParticles,0); //  newAction += PathData.Actions.TIP5PWater.SecondProtonKineticAction(startSlice,endSlice,p2ActiveParticles,0); //  newAction += PathData.Actions.TIP5PWater.NewRotKinAction(startSlice,endSlice,HActiveParticles,0);
  newAction += PathData.Actions.TIP5PWater.FixedAxisAction(startSlice,endSlice,HActiveParticles,0);

  //cerr<<"ROTATE:  The actions are "<<newAction<<" "<<oldAction<<endl;
 
  if (-(newAction-oldAction)>=log(PathData.Path.Random.Local())){
    PathData.AcceptMove(startSlice,endSlice,ActiveParticles);
    //cerr<<"ROTATE:  I've accepted " << NumAccepted << " " << TimesCalled+1 <<endl;
    NumAccepted++;
    total_r_squared += 2*(1-cos(theta));
  }
  else {
    PathData.RejectMove(startSlice,endSlice,ActiveParticles);
    //cerr<<"ROTATE:  I've rejected"<<endl;
  }
  TimesCalled++;
  newtime = clock();
  elapsed += (newtime - oldtime);
  oldtime = newtime;
  if (TimesCalled%10000 == 0){
    cerr << TimesCalled << " moves; current rotate ratio is " << double(NumAccepted)/TimesCalled << " with angle size " << dtheta << endl;
    double avg = sqrt(total_r_squared)/NumAccepted;
    double diff = total_r_squared/TimesCalled;
    cerr << "ROTATE diffusion value is " << diff << " avg angle is " << avg << endl;
  }
}


void WaterRareRotate::MakeMove()
{
  double dtheta = 2*M_PI;
  int speciesO = PathData.Path.SpeciesNum("O");
  int  numWater=PathData.Path.Species(speciesO).LastPtcl-
    PathData.Path.Species(speciesO).FirstPtcl+1;
  int choosemol = (int)floor(PathData.Path.Random.Local()*numWater);

// choose a time slice to move
  int numSlices = PathData.Path.TotalNumSlices;
  int slice=0;
  int endSlice = 0;
  int startSlice = 0;
  if(numSlices>1){
    int P_max = numSlices - 1;
    slice = (int)floor(P_max*PathData.Path.Random.Local()) + 1;
    startSlice = slice-1;
    endSlice = slice+1;
  }
  
  Molecule2Atoms(choosemol);

  for (int i=0;i<coord_loc.size();i++){
    ActiveParticles(i)=coord_loc(i);
  }

  double oldAction = 0.0;
  oldAction += PathData.Actions.ST2Water.Action(slice,slice,ActiveParticles,0);

  // generate three random angles
  double theta1 = (PathData.Path.Random.Local()-0.5)*dtheta;
  double theta2 = (PathData.Path.Random.Local()-0.5)*dtheta;
  double theta3 = (PathData.Path.Random.Local()-0.5)*dtheta;
  double angularmag = sqrt(theta1*theta1 + theta2*theta2 + theta3*theta3);

  // make three moves
  // Here we use ArbitraryRotate()
  // generate axes of rotation
  dVec O = PathData.Path(slice,coord_loc(0));
  dVec P1 = PathData.Path(slice,coord_loc(1)); 
  dVec P2 = PathData.Path(slice,coord_loc(2));
//cerr << "for ptcl " << coord_loc(0) << " we have protons " << coord_loc(1) << " and " << coord_loc(2) << endl;
//cerr << "coordinates are " << O << P1 << P2 << endl;
  P1 -= O;
  P2 -= O;
//cerr << "now coordinates are " << O << P1 << P2 << endl;
  dVec R1 = PathData.Actions.TIP5PWater.Normalize(PathData.Actions.TIP5PWater.GetBisector(P1,P2));
//cerr << "going to rotate about " << R1 << endl;
  dVec R2 = PathData.Actions.TIP5PWater.Normalize(PathData.Actions.TIP5PWater.CrossProd(P1,P2));
  dVec R3 = PathData.Actions.TIP5PWater.Normalize(PathData.Actions.TIP5PWater.CrossProd(R1,R2));
  // now do rotations
  for(int i=1;i<coord_loc.size();i++){
    dVec old_coord = PathData.Path(slice,coord_loc(i));
    old_coord = old_coord-PathData.Path(slice,coord_loc(0));
    //dVec new_coord = old_coord;
    dVec new_coord = ArbitraryRotate(R1,old_coord,theta1);
//cerr << "   " << i << ": now ptcl " << coord_loc(i) << " of mol " << PathData.Path.MolRef(coord_loc(i)) << " is " << new_coord << endl;
    new_coord = ArbitraryRotate(R2,new_coord,theta2);
    new_coord = ArbitraryRotate(R3,new_coord,theta3);
    new_coord = new_coord+PathData.Path(slice,coord_loc(0));
    PathData.Path.SetPos(slice,coord_loc(i),new_coord);
  }

// This is the original version which uses Rotate()
/*
  for(int i=1;i<coord_loc.size();i++){
    dVec old_coord=PathData.Path(slice,coord_loc(i));
    old_coord=old_coord-PathData.Path(slice,coord_loc(0));
    dVec new_coord = Rotate(old_coord,0,1,2,theta1);
    new_coord = Rotate(new_coord,1,2,0,theta2);
    new_coord = Rotate(new_coord,2,0,1,theta3);
    new_coord=new_coord+PathData.Path(slice,coord_loc(0));
    PathData.Path.SetPos(slice,coord_loc(i),new_coord);
  }*/
  
  // catalog accepted moves by magnitude
  double newAction = 0.0;
  newAction += PathData.Actions.ST2Water.Action(slice,slice,ActiveParticles,0);
 
  if (-(newAction-oldAction)>=log(PathData.Path.Random.Local())){
    PathData.AcceptMove(startSlice,endSlice,ActiveParticles);
    //cerr<<"ROTATE:  I've accepted " << NumAccepted << " " << TimesCalled+1 <<endl;
    NumAccepted++;
    total_r_squared += 2*(1-cos(angularmag));
    int slot = FindSlot(angularmag,1);
    AcceptLog(0,slot)++;
    AcceptLog(1,FindSlot(abs(theta1),2))++;
    AcceptLog(2,FindSlot(abs(theta2),2))++;
    AcceptLog(3,FindSlot(abs(theta3),2))++;
  }
  else {
    PathData.RejectMove(startSlice,endSlice,ActiveParticles);
    //cerr<<"ROTATE:  I've rejected"<<endl;
  }

  IntegrityCheck(slice,ActiveParticles);

  TimesCalled++;
  if (TimesCalled%10000 == 0){
    cerr << TimesCalled << " moves; current rare rotate ratio is " << double(NumAccepted)/TimesCalled << " with angle size " << dtheta << endl;
    double avg = sqrt(total_r_squared)/NumAccepted;
    double diff = total_r_squared/TimesCalled;
    cerr << "RARE ROTATE diffusion value is " << diff << " avg angle is " << avg << endl;
  }
  if (TimesCalled%200000 == 0){
    cerr << TimesCalled << ": Printing acceptance table of length " << AcceptLog.size() << endl;
    cerr << "printing histogram for angularmag" << endl;
    for (int i = 0; i < numSlots; i++)
      cerr << i*sqrt(3.0)*2*M_PI/numSlots << " " << AcceptLog(0,i) << endl;
    cerr << "printing histogram for theta1" << endl;
    for (int i = 0; i < numSlots; i++)
      cerr << i*2*M_PI/numSlots << " " << AcceptLog(1,i) << endl;
    cerr << "printing histogram for theta2" << endl;
    for (int i = 0; i < numSlots; i++)
      cerr << i*2*M_PI/numSlots << " " << AcceptLog(2,i) << endl;
    cerr << "printing histogram for theta3" << endl;
    for (int i = 0; i < numSlots; i++)
      cerr << i*2*M_PI/numSlots << " " << AcceptLog(3,i) << endl;
  }
}
