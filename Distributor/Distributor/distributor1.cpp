// ============================================================================
// THIS FILE IS COPYRIGHTED - SEE Kernel/COPYRIGHT.txt
// ============================================================================
// File   : Distributor/distributor1.cpp
// SVN    : $Rev$ $Date$
// Cksum  : 
// ============================================================================

//============================================================================
// General routines
// ===========================================================================
# include "distributor.h"

# include "Pcm/Base/PcmTask.h"
# include "Pcm/Base/PcmDefine.h"

using namespace E_FUNC;

// ------------------------------------------------------------------------
SplDistributor::SplDistributor() :
  _meanPointsPerProc(-1)
{
}
// ------------------------------------------------------------------------
const list<DesBlock*> 
SplDistributor::getSplittedBlocks()  const
{
  return getBlocksToSplit();
}
// ------------------------------------------------------------------------
const list<DesMesh*>  
SplDistributor::getSplittedMeshes()  const
{
  return getMeshesToSplit();
}
// ------------------------------------------------------------------------
const list<DesMask*>  
SplDistributor::getSplittedMasks()   const
{
  return getMasksToSplit();
}
// ------------------------------------------------------------------------
const list<DesBoundary*> 
SplDistributor::getSplittedBnds() const
{
  return getBndsToSplit();
}
// ------------------------------------------------------------------------
void  
SplDistributor::printSplitInfo()
{
  if ( _meanPointsPerProc >= 0 )
  {
    cerr << "Distributions of meshes amongst computing nodes : \n";
    list<DesMesh*> splMsh = _meshesToSplit;
    splMsh.sort();
    splMsh.unique();
    FldArrayI nbPtsPerProcs(PcmTask::getTask()->numberOfProc());
    nbPtsPerProcs.setAllValuesAt(0);
    for ( list<DesMesh*>::iterator itMesh = _meshesToConserve.begin();
          itMesh != _meshesToConserve.end(); itMesh++ )
    {
      E_Int node = (*itMesh)->getI(KEY_NODE);
      nbPtsPerProcs[node] = nbPtsPerProcs[node] +
        (*itMesh)->getI(KEY_IM)*
        (*itMesh)->getI(KEY_JM)*(*itMesh)->getI(KEY_KM);
    }

    for ( list<DesMesh*>::iterator it = splMsh.begin();
          it != splMsh.end(); it++ )
    {
      E_Int im = (*it)->getI(KEY_IM);
      E_Int jm = (*it)->getI(KEY_JM);
      E_Int km = (*it)->getI(KEY_KM);
      E_Int proc = (*it)->getI(KEY_NODE);
      cerr << "Mesh " << (*it)->getName() << " : \n";
      cerr << "\tMesh dimension (" 
           << im << ", "
           << jm << ", "
           << km << " ) on node "
           << proc << endl;
      nbPtsPerProcs[proc] += im*jm*km;
    }
    cerr << "Mean points per computing node : "
         << _meanPointsPerProc << endl;
    for ( E_Int p = 0; p < PcmTask::getTask()->numberOfProc(); p++ )
      cerr << "\tNumber of points on computing node " << p << " -> " << nbPtsPerProcs[p] << endl;
  }
  else
  {
    cerr << "No splitting yet !\n";
  }
}
// ------------------------------------------------------------------------
E_Int
SplDistributor::getIdOfOriginalBlock( E_Int idSpltBlock )
{
  return idSpltBlock;
}
// ------------------------------------------------------------------------
SplImpl*
SplDistributor::clone() const
{
  // Fonction provisoire :
  return new SplDistributor;
}
