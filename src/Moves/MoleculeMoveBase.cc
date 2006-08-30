#include "MoleculeMoveBase.h"
#include "MoveUtils.h"

// default constructor loads MolMembers array of arrays
MolMoveClass::MolMoveClass(PathDataClass& myPathData, IO::IOSectionClass outSection):
  LocalStageClass (myPathData,outSection)
{
	cerr << "MolMoveClass constructor" << endl;
	cerr << "	numMol is " << numMol << endl;
	numAccepted = 0;
	numMoves = 0;

	mode = SINGLE;
 
	/* 
  cerr << "In MolMoveClass Constructor.  MolMembers loaded: " << endl;
  for(int i = 0; i < MolMembers.size() ; i++){
    cerr << i << ":  ";
    for(int j = 0; j < MolMembers(i).size(); j++) cerr << MolMembers(i)(j) << ", ";
    cerr << endl;
  }*/
}

void MolMoveClass::Read (IOSectionClass &in){
	/// Read in the molecule to move by string or int
	cerr << "MolMoveClass::Read" << endl;
	string setMolecule;
	int indexMolecule;
	if(in.ReadVar("Molecule",setMolecule)){
		cerr << "got Molecule " << setMolecule << endl;
		indexMolecule = 0;
		cerr << " index is " << indexMolecule << " which gets " << PathData.Path.MoleculeName[indexMolecule];
		cerr  << " out of total size " << PathData.Path.MoleculeName.size() << endl;
		while((PathData.Path.MoleculeName[indexMolecule] != setMolecule) && (indexMolecule<PathData.Path.MoleculeName.size())){
			cerr << "indexMolecule is " << indexMolecule << "... max is " << PathData.Path.MoleculeName.size() << endl;
			indexMolecule++;
		}
		molIndex = indexMolecule;
		numMol = PathData.Path.MoleculeNumber[molIndex];
		cerr << "Selected molecule " << PathData.Path.MoleculeName[molIndex] << endl;
	}
	else if (in.ReadVar("Molecule",indexMolecule)){
		molIndex = indexMolecule;
		numMol = PathData.Path.MoleculeNumber[molIndex];
		cerr << "Selected molecule " << PathData.Path.MoleculeName[molIndex] << endl;
	}
	else{
		molIndex = 0;
		numMol = PathData.Path.MoleculeNumber[molIndex];
		cerr << "Selected molecule " << PathData.Path.MoleculeName[molIndex] << " by default" << endl;
	}
	
/// Read in update mode: GLOBAL, SEQUENTIAL, SINGLE	
	string setMode;
	in.ReadVar("Mode", setMode);
	if(setMode == "GLOBAL"){
		mode = GLOBAL;
		MoveList.resize(numMol);
		for(int m=0; m<numMol; m++)
			MoveList(m) = m + PathData.Path.offset[molIndex];
		cerr << "Set for GLOBAL particle updates." << endl;
	}
	else if(setMode == "SEQUENTIAL"){
		mode = SEQUENTIAL;
		MoveList.resize(1);
		MoveList(0) = 0;
		cerr << "Set for SEQUENTIAL particle updates." << endl;
	}
	else{
		cerr << "Set for SINGLE particle updates." << endl;
		MoveList.resize(1);
	}
		
  Array<string,1> ActionList;
  assert (in.ReadVar ("Actions", ActionList));
	for(int a=0; a<ActionList.size(); a++){
		string setAction = ActionList(a);
	  if(setAction == "ST2Potential"){
  		Actions.push_back(&PathData.Actions.ST2Water);
			cerr << "  Added ST2 Water Potential" << endl;
		}else if(setAction == "TIP5PPotential"){
  		Actions.push_back(&PathData.Actions.TIP5PWater);
			cerr << "Added TIP5P Water Potential" << endl;
#ifdef USE_QMC
		}else if(setAction == "CEIMCAction"){
			//ActionBaseClass* newAction(PathData.Actions.CEIMCAction);
			PathData.Actions.CEIMCAction.Read(in);
  		Actions.push_back(&PathData.Actions.CEIMCAction);
			cerr << "Added CEIMC calculation of BO energy" << endl;
#endif
		}else if(setAction == "IonInteraction"){
  		Actions.push_back(&PathData.Actions.IonInteraction);
			cerr << "Added intermolecular ion-ion interaction" << endl;
		} else
    	cerr << "You specified " << setAction << ", which is not supported for this type of move" << endl;
	}

  MolMembers.resize(numMol);
  vector<int> catalog(numMol); // for storing info to load MolMembers
  // initialize catalog (all zeros)
  for(int m = 0; m < catalog.size(); m++)
    catalog[m] = 0;
	// get number of ptcls for each molecule m; store in catalog
  for(int p = 0; p <PathData.Path.NumParticles(); p++)
		catalog[PathData.Path.MolRef(p)]++;
	// resize MolMembers arrays appropritely; re-initialize catalog
  for(int m = 0; m < catalog.size(); m++){
    MolMembers(m).resize(catalog[m]);
		catalog[m] = 0;
  }
	// load ptcls into the MolMembers array; catalog indexes
  for(int p = 0; p <PathData.Path.NumParticles(); p++){
    int m = PathData.Path.MolRef(p);
		MolMembers(m)(catalog[m]) = p;
		catalog[m]++;
	}
}

dVec MolMoveClass::GetCOM(int slice, int mol){
  return PathData.Path(slice,MolMembers(mol)(0));
}

// Generates a translation vector and moves all
// particles in the molecule
dVec MolMoveClass::TranslateMol(int slice, Array<int,1>& activePtcls, double epsilon){
  dVec translate;
  for(int i = 0; i<3; i++) translate[i] = (PathData.Path.Random.Local()-0.5)*2*epsilon;	
//cerr << "  translate by " << translate << endl;
  for(int ptcl = 0; ptcl<activePtcls.size(); ptcl++){
    dVec newP = PathData.Path(slice,activePtcls(ptcl)) + translate;
    PathData.Path.SetPos(slice,activePtcls(ptcl),newP);
  }
	return translate;
}

// This will orchestrate the rotation of molecule mol
// about a specified axis, axis, by an angle theta
void MolMoveClass::RotateMol(int slice, Array<int,1>& activePtcls, dVec& axis, double theta){
	int mol = activePtcls(0);
  axis = Normalize(axis);
  // find COM vector
  dVec O = GetCOM(slice, mol);
	// nonsense to rotate the ptcl at the origin; "O", so loop starts from 1 not 0
  for(int ptcl = 1; ptcl<activePtcls.size(); ptcl++){
    dVec P = PathData.Path(slice, activePtcls(ptcl)) - O;
    dVec newP = ArbitraryRotate(axis, P, theta) + O;
    PathData.Path.SetPos(slice,activePtcls(ptcl),newP);
  }
}

// This will orchestrate and execute an arbitrary rotation
// of a molecule using ArbitraryRotate
// i.e. it generates its own axis of rotation
void MolMoveClass::RotateMol(int slice, Array<int,1>& activePtcls, double theta){
  // generate a unit vector
  dVec axis;
  for(int i = 0; i<3; i++) axis(i) = PathData.Path.Random.Local()-0.5;
  RotateMol(slice, activePtcls, axis, theta);
}

void MolMoveClass::Accept(){
	numAccepted++;
}

void MolMoveClass::Advance(){
	MoveList(0) = MoveList(0) + 1;
	if(MoveList(0) >= PathData.Path.numMol)
		MoveList(0) = 0;
}

// ArbitraryRotate takes a unit vector, axis, about
// which to rotate and coord, a coordinate WRT
// the origin of axis, and an angle phi by which
// to rotate.  Then, it returns the rotated coordinate

// Is this really the simplest algorithm to do a rotation
// about an arbitrary axis????
dVec ArbitraryRotate(dVec axis,dVec coord, double phi){
  dVec coord_perp;
  dVec coord_aligned;
  Strip(axis,coord,coord_aligned,coord_perp);
  double perp_mag = Mag(coord_perp);
  coord =  Normalize(coord_perp);
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
  return (Scale(newcoord,perp_mag) + coord_aligned);
}
