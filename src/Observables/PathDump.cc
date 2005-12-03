#include "PathDump.h"


void PathDumpClass::Accumulate()
{
  WriteBlock();
}

void PathDumpClass::Read(IOSectionClass &in)
{
  ObservableClass::Read(in);
  assert (in.ReadVar ("dumpFreq", DumpFreq));
  DumpNodes = in.OpenSection("NodeDump");
  if (DumpNodes) {
    assert(in.ReadVar("Ptcl", NodePtcl));
    assert(in.ReadVar("Slice", NodeSlice));
    int nx, ny, nz;
    assert(in.ReadVar("Nx", nx));
    assert(in.ReadVar("Ny", ny));
    assert(in.ReadVar("Nz", nz));
    dVec box = PathData.Path.GetBox();
    Xgrid.Init(-0.5*box[0], 0.5*box[0], nx);
    Ygrid.Init(-0.5*box[1], 0.5*box[1], ny);
    Zgrid.Init(-0.5*box[2], 0.5*box[2], nz);
    if (PathData.Path.Communicator.MyProc()==0) {
      IOSection.WriteVar("NodeSlice", NodeSlice);
      IOSection.WriteVar("NodePtcl", NodePtcl);
      IOSection.NewSection("Xgrid");
      Xgrid.Write(IOSection);
      IOSection.CloseSection();
      IOSection.NewSection("Ygrid");
      Ygrid.Write(IOSection);
      IOSection.CloseSection();
      IOSection.NewSection("Zgrid");
      Zgrid.Write(IOSection);
      IOSection.CloseSection();
    }
    in.CloseSection();
  }
}

void PathDumpClass::WriteBlock()
{
  PathClass &Path = PathData.Path;
  int start, end, numProcs, myProc;
  numProcs = Path.Communicator.NumProcs();
  myProc   = Path.Communicator.MyProc();

  int refSave = PathData.Path.GetRefSlice();
  PathData.MoveRefSlice(0);

  if (DumpNodes) 
    NodeDump();


  if (PathData.Path.OpenPaths){
    Array<double,1> tailLoc(NDIM);
    OpenLinkVar.Write((int)PathData.Path.OpenLink);
    for (int dim=0;dim<NDIM;dim++){
      tailLoc(dim)=PathData.Path(PathData.Path.OpenLink,
				 PathData.Path.NumParticles())[dim];
    }
    TailLocVar.Write(tailLoc);
    RefLinkVar.Write(PathData.Path.RefSlice);
    OpenLinkPtclVar.Write(PathData.Path.OpenPtcl);
  }

  Array<double,3> pathArray;
  int totalSlices = Path.TotalNumSlices;
  Path.SliceRange(numProcs-1, start,end);
  int maxShift = end-start;
  int slicesLeft = totalSlices;
  int offset = 0;
  int numPtcls = PathData.NumParticles();

  if (myProc == 0)
    pathArray.resize(numPtcls,totalSlices,NDIM);
  while (slicesLeft > maxShift) {
    int relNodeSlice = NodeSlice - offset;
    if ((relNodeSlice>=0) && (relNodeSlice<Path.NumTimeSlices())
	&& (myProc == 0))
      NodeDump ();
    // First copy
    PathData.MoveJoin(maxShift);
    if (myProc == 0)
      for (int i=0; i<maxShift; i++)
	for (int ptcl=0; ptcl < numPtcls; ptcl++) 
	  for (int dim=0; dim<NDIM; dim++) 
	    pathArray(ptcl, i+offset, dim) = Path(i,ptcl)[dim];
    
    // Now shift
    PathData.ShiftData(-maxShift);
    PathData.Join = 0;
    slicesLeft -= maxShift;
    offset += maxShift;
  }
  // Move join out of the way
  PathData.MoveJoin (Path.NumTimeSlices()-1);
  // Copy last part
  if (myProc == 0)
    for (int i=0; i<slicesLeft; i++)
      for (int ptcl=0; ptcl < numPtcls; ptcl++)
	for (int dim=0; dim<NDIM; dim++) 
	  pathArray(ptcl, i+offset, dim) = Path(i,ptcl)[dim];
  
  // Reset path to original position
  PathData.MoveRefSlice (refSave);

  Array<int,1> permVec(numPtcls);
  Path.TotalPermutation(permVec);

  // Now write
  PathVar.Write (pathArray);
  PermVar.Write (permVec);

  if (FirstTime && (myProc == 0)) {
    WriteInfo();
    Array<string,1> speciesNames(numPtcls);
    for (int spIndex=0;spIndex<PathData.NumSpecies();spIndex++){
      SpeciesClass &species = PathData.Path.Species(spIndex);
      for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++)
	speciesNames(ptcl)=species.Name;
    }
    IOSection.WriteVar("SpeciesNames", speciesNames);
    IOSection.WriteVar("Type","Path");
  }
  PathVar.Flush();
  FirstTime = false;
}


void
PathDumpClass::NodeDump()
{
  int speciesNum = PathData.Path.ParticleSpeciesNum(NodePtcl);
  if (PathData.Actions.NodalActions(speciesNum) != NULL) {
    NodeType nodeType = PathData.Actions.NodalActions(speciesNum)->Type();
    if (nodeType == GROUND_STATE)
      GroundStateNodeDump();
    else if (nodeType == GROUND_STATE_FP)
      FixedPhaseNodeDump();
    else if (nodeType == FREE_PARTICLE)
      FreeParticleNodeDump();
  }
}

void
PathDumpClass::FixedPhaseNodeDump()
{
  if (PathData.Path.UseCorrelatedSampling()) 
    PathData.Path.SetIonConfig(0);

  /// HACK HACK HACK
  /// Find ptcl closest to ion we moved.
  double closestDist = 1.0e100;
  int closestPtcl = 16;
  for (int ptcl=16; ptcl < 32; ptcl++) {
    dVec disp = PathData.Path(NodeSlice,ptcl) - PathData.Path(NodeSlice,0);
    PathData.Path.PutInBox(disp);
    double dist = sqrt(dot(disp,disp));
    if (dist < closestDist) {
      closestDist = dist;
      closestPtcl = ptcl;
    }
  }
  NodePtcl = closestPtcl;



  int nx = Xgrid.NumPoints;
  int ny = Ygrid.NumPoints;
  int nz = Zgrid.NumPoints;

  PathClass &Path = PathData.Path;
  Array<double,3> nodeDump(nx,ny,nz);
  int speciesNum = PathData.Path.ParticleSpeciesNum(NodePtcl);
  FixedPhaseActionClass &FP = 
    *((FixedPhaseActionClass*)PathData.Actions.NodalActions(speciesNum));

  /// This must be fixed!
  /// HACK HACK
  int relNodeSlice = NodeSlice; // + Path.GetRefSlice();
  if (relNodeSlice > Path.TotalNumSlices)
    relNodeSlice -= Path.TotalNumSlices;

  dVec savePos = Path(relNodeSlice,NodePtcl);
  dVec newPos;
  for (int ix=0; ix<nx; ix++) {
    newPos[0] = Xgrid(ix);
    for (int iy=0; iy<ny; iy++) {
      newPos[1] = Ygrid(iy);
      for (int iz=0; iz<nz; iz++) {
	newPos[2] = Zgrid(iz);
	Path(relNodeSlice,NodePtcl) = newPos;
	nodeDump(ix,iy,iz) = log10(FP.CalcGrad2(relNodeSlice));
      }
    }
  }
  Path(relNodeSlice,NodePtcl) = savePos;
  NodeVarA.Write(nodeDump);
  if (PathData.Path.UseCorrelatedSampling()) {
    PathData.Path.SetIonConfig(1);
    for (int ix=0; ix<nx; ix++) {
      newPos[0] = Xgrid(ix);
      for (int iy=0; iy<ny; iy++) {
	newPos[1] = Ygrid(iy);
	for (int iz=0; iz<nz; iz++) {
	  newPos[2] = Zgrid(iz);
	  Path(relNodeSlice,NodePtcl) = newPos;
	  nodeDump(ix,iy,iz) = log10(FP.CalcGrad2(relNodeSlice));
	}
      }
    }
    Path(relNodeSlice,NodePtcl) = savePos;
    NodeVarB.Write(nodeDump);
    PathData.Path.SetIonConfig(0);
  }
}


void
PathDumpClass::GroundStateNodeDump()
{
  int nx = Xgrid.NumPoints;
  int ny = Ygrid.NumPoints;
  int nz = Zgrid.NumPoints;

  PathClass &Path = PathData.Path;
  Array<double,3> nodeDump(nx,ny,nz);
  int speciesNum = Path.ParticleSpeciesNum(NodePtcl);
  GroundStateNodalActionClass &GS = 
    *((GroundStateNodalActionClass*)PathData.Actions.NodalActions(speciesNum));

  int relNodeSlice = NodeSlice + Path.GetRefSlice();
  if (relNodeSlice > PathData.Path.TotalNumSlices)
    relNodeSlice -= Path.TotalNumSlices;

  dVec savePos = Path(relNodeSlice,NodePtcl);
  dVec newPos;
  if (PathData.Path.UseCorrelatedSampling()) 
    PathData.Path.SetIonConfig(0);
  for (int ix=0; ix<nx; ix++) {
    newPos[0] = Xgrid(ix);
    for (int iy=0; iy<ny; iy++) {
      newPos[1] = Ygrid(iy);
      for (int iz=0; iz<nz; iz++) {
	newPos[2] = Zgrid(iz);
	Path(relNodeSlice,NodePtcl) = newPos;
	nodeDump(ix,iy,iz) = GS.Det(relNodeSlice);
      }
    }
  }
  Path(relNodeSlice,NodePtcl) = savePos;
  NodeVarA.Write(nodeDump);
  if (PathData.Path.UseCorrelatedSampling()) {
    PathData.Path.SetIonConfig(1);
    for (int ix=0; ix<nx; ix++) {
      newPos[0] = Xgrid(ix);
      for (int iy=0; iy<ny; iy++) {
	newPos[1] = Ygrid(iy);
	for (int iz=0; iz<nz; iz++) {
	  newPos[2] = Zgrid(iz);
	  Path(relNodeSlice,NodePtcl) = newPos;
	  nodeDump(ix,iy,iz) = GS.Det(relNodeSlice);
	}
      }
    }
    Path(relNodeSlice,NodePtcl) = savePos;
    NodeVarB.Write(nodeDump);
    PathData.Path.SetIonConfig(0);
  }
}


void
PathDumpClass::FreeParticleNodeDump()
{
  int nx = Xgrid.NumPoints;
  int ny = Ygrid.NumPoints;
  int nz = Zgrid.NumPoints;

  PathClass &Path = PathData.Path;
  Array<double,3> nodeDump(nx,ny,nz);
  int speciesNum = PathData.Path.ParticleSpeciesNum(NodePtcl);
  FreeNodalActionClass &GS = 
    *((FreeNodalActionClass*)PathData.Actions.NodalActions(speciesNum));

  int relNodeSlice = NodeSlice + Path.GetRefSlice();
  if (relNodeSlice > Path.TotalNumSlices)
    relNodeSlice -= Path.TotalNumSlices;

  if (PathData.Path.UseCorrelatedSampling())
    PathData.Path.SetIonConfig(0);
  dVec savePos = Path(relNodeSlice,NodePtcl);
  dVec newPos;
  for (int ix=0; ix<nx; ix++) {
    newPos[0] = Xgrid(ix);
    for (int iy=0; iy<ny; iy++) {
      newPos[1] = Ygrid(iy);
      for (int iz=0; iz<nz; iz++) {
	newPos[2] = Zgrid(iz);
	Path(relNodeSlice,NodePtcl) = newPos;
	nodeDump(ix,iy,iz) = GS.Det(relNodeSlice);
      }
    }
  }
  Path(relNodeSlice,NodePtcl) = savePos;
  NodeVarA.Write(nodeDump);
}
