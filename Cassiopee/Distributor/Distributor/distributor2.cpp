// ============================================================================
// THIS FILE IS COPYRIGHTED - SEE Kernel/COPYRIGHT.txt
// ============================================================================
// File   : Distributor/distributor2.cpp
// SVN    : $Rev$ $Date$
// Cksum  : 
// ============================================================================

// ===========================================================================
// Part specific to split
// ===========================================================================

# include "distributor.h"
# include "Pcm/Base/PcmTask.h"
# include "Pcm/Base/PcmDefine.h"

using namespace E_FUNC;

// ---------------------------------------------------------------------------
void
SplDistributor::split()
{
  if (PcmTask::getTask()->myNode() == 0) printf("Distributing...\n");

  //==========================================================================
  // Pre-treatment:
  // Create: locMeshes and glbMeshes
  // Create: sortedLocMeshes and sortedGlbMeshes
  // Unique par nom des listes _meshesToConserve
  //==========================================================================
  list<DesMesh*>& locMeshes = _meshesToConserve;
  list<DesMesh*> sortedLocMeshes;
  for (list<DesMesh*>::iterator itMesh = locMeshes.begin();
       itMesh != locMeshes.end(); itMesh++)
  {
    E_Boolean found = E_False;
    for ( list<DesMesh*>::iterator itMesh2 = sortedLocMeshes.begin();
          itMesh2 != sortedLocMeshes.end(); itMesh2++)
      if ( (*itMesh)->getNameK() == (*itMesh2)->getNameK() )
        found = E_True;

    if (found == E_False)
      sortedLocMeshes.push_back(*itMesh);
  }
  locMeshes = sortedLocMeshes;
  
  list<DesMesh*>& glbMeshes = _meshesToSplit;
  list<DesMesh*> sortedGlbMeshes;
  for ( list<DesMesh*>::iterator itMesh = glbMeshes.begin();
        itMesh != glbMeshes.end(); itMesh++ )
  {
    E_Boolean found = E_False;
    for ( list<DesMesh*>::iterator itMesh2 = sortedGlbMeshes.begin();
          itMesh2 != sortedGlbMeshes.end(); itMesh2++)
      if ( (*itMesh)->getNameK() == (*itMesh2)->getNameK() )
        found = E_True;

    if (found == E_False)
      sortedGlbMeshes.push_back(*itMesh);
  }
  glbMeshes = sortedGlbMeshes;
  
  //----------------------------
  // Data for genetic algorithm
  //----------------------------
  // Size of population
  E_Int sizeOfPopulation = E_max(E_Int(30), 2*_meshesToSplit.size());
  sizeOfPopulation = E_min( sizeOfPopulation, 200);
  // Number max of generations
  //const E_Int nitMax = E_max(E_Int(50), 2*_meshesToSplit.size());
  const E_Int nitMax = 50;
  if (PcmTask::getTask()->myNode() == 0)
    printf("Population size = %d, %d.\n", sizeOfPopulation, nitMax);
  // latence : incompressible time to establish a communication
  const E_Float latence = 1.e-2;
  // Communication speed
  const E_Float comSpeed = 1.e-2;
  
  // Local data
  srand48(2003);
  E_Int nbProc = PcmTask::getTask()->numberOfProc();
  E_Int me     = PcmTask::getTask()->myNode();
  
  // Nb de noeuds des maillages imposes par processeurs
  // OUT : nbNodePerProc
  FldArrayI nbNodePerProc(nbProc);
  nbNodePerProc.setAllValuesAt(0);
  E_Int nbLocNodes = 0;
  for ( list<DesMesh*>::iterator itMesh = locMeshes.begin();
        itMesh != locMeshes.end(); itMesh++ )
  {
    if ( (*itMesh)->getI(KEY_NODE) == me )
    {
      nbLocNodes += (*itMesh)->getI(KEY_IM)*
        (*itMesh)->getI(KEY_JM)*(*itMesh)->getI(KEY_KM);
    }
  }
  E_PCM::pgall( nbLocNodes, nbNodePerProc );
  
  // Calcul du nombre total de noeuds des blocs imposes :
  nbLocNodes = 0;
  for ( E_Int i = 0; i < nbProc; i++ )
    nbLocNodes += nbNodePerProc[i];
  
  for ( list<DesMesh*>::iterator itMesh = glbMeshes.begin();
        itMesh != glbMeshes.end(); itMesh++ )
  {
    nbLocNodes += (*itMesh)->getI(KEY_IM)*
      (*itMesh)->getI(KEY_JM)*(*itMesh)->getI(KEY_KM);
  }
  
  // On en deduit le nombre de noeuds moyen devant etre contenus dans
  // chaque processeur :
  E_Float meanPtsPerProc = E_Float(nbLocNodes)/E_Float(nbProc);
  _meanPointsPerProc = meanPtsPerProc;

  //----------------------------------------------------------
  // Algorithme genetique pour essayer d'equilibrer au mieux :
  //----------------------------------------------------------
  E_Int   nbItWithFail = 0;
  E_Float bestEval = 1.E+12;
  const unsigned long  MAXRAND = 1UL<<31;
  FldArrayI nicePeople( glbMeshes.size(), sizeOfPopulation );
  FldArrayF evalPeople( sizeOfPopulation );
  E_Float* evalp = evalPeople.begin();
  FldArrayI survivors( sizeOfPopulation );
  E_Int* survivorsp = survivors.begin();
  
  // Initialization :
  // Pour le premier people ( of nice ! ),
  // on essaie de mettre un truc pas trop con.
  FldArrayI nbPtsPerProcs(nbNodePerProc);
  E_Int* nbPtsPerProcsp = nbPtsPerProcs.begin();
  E_Int ind = 0;
  E_Int indBest = 0;
  for ( list<DesMesh*>::iterator itMesh = glbMeshes.begin();
        itMesh != glbMeshes.end(); itMesh++ )
  {
    E_Float bestEval = meanPtsPerProc;
    E_Int nbPtsForMesh = (*itMesh)->getI(KEY_IM)*
      (*itMesh)->getI(KEY_JM)*(*itMesh)->getI(KEY_KM);
    for ( E_Int p = 0; p < nbProc; p++ )
    {
      E_Float eval = E_abs(meanPtsPerProc-nbPtsForMesh-nbPtsPerProcsp[p]);
      if ( eval < bestEval )
      {
	bestEval = eval;
	indBest  = p; 
      }
    }
    nbPtsPerProcsp[indBest] += nbPtsForMesh;
    nicePeople(ind,1) = indBest;
    ind++;
  }
  
  // Pour le deuxieme, c'est une bete humaine:
  ind = 0;
  for (unsigned int i = 0; i < glbMeshes.size(); i++)
    nicePeople(i,2) = ind%nbProc;
  
  // Pour les autres, c'est completement au hasard :
  for (E_Int j = 3; j <= sizeOfPopulation; j++)
    for (size_t i = 0; i < glbMeshes.size(); i++)
    {
      E_Int targetProc = E_Int((nbProc*E_Float(lrand48()))/(E_Float(MAXRAND)));
      nicePeople(i,j) = targetProc;
    }
  
  // First evaluation for nice people :
  for (E_Int j = 1; j <= sizeOfPopulation; j++)
  {
    nbPtsPerProcs = nbNodePerProc;
    nbPtsPerProcsp = nbPtsPerProcs.begin();
    E_Int ind = 0;
    for ( list<DesMesh*>::iterator itMesh = glbMeshes.begin();
          itMesh != glbMeshes.end(); itMesh++ )
    {
      nbPtsPerProcsp[nicePeople(ind,j)] += (*itMesh)->getI(KEY_IM)*
        (*itMesh)->getI(KEY_JM)*(*itMesh)->getI(KEY_KM);
      ind++;
    }
    evalp[j-1] = 0.;
    for ( E_Int p = 0; p < nbProc; p++ )
      evalp[j-1] += (nbPtsPerProcsp[p]-meanPtsPerProc)*
        (nbPtsPerProcsp[p]-meanPtsPerProc);
    evalp[j-1] = sqrt(evalp[j-1]);
    ind = 0;
    for (list<SplApi::ConnectInfo*>::iterator itrConnect = _connect.begin();
         itrConnect != _connect.end(); itrConnect++)
    {
      SplApi::ConnectInfo* c = *itrConnect;
      while (c != NULL)
      {
        // le voisin est-il sur le meme processeur?
        if ( nicePeople(c->neighbour-1,j) == nicePeople(ind,j) )
        {
          evalp[j-1] += latence;
          evalp[j-1] += comSpeed*c->size;
        }
        c = c->next;
      }
      ind++;
    }
    // Check for empty processors
    for (E_Int i = 0; i < nbProc; i++)
    {
      if (nbPtsPerProcsp[i] == 0) evalp[j-1] += 1.e6;
    }
  }
  bestEval = 1.E+12;
  E_Int jBest = 1;
  
  while (nbItWithFail < nitMax)
  {
    // Phase of selection :
    // 1st : Search best evaluation in the population :
    jBest = 1;
    for ( E_Int j = 1; j <= sizeOfPopulation; j++ )
      jBest = ( evalp[j-1] < evalp[jBest-1] ? j : jBest );
    
    // On tire au hasard pour voir quelle partie de la population survit :
    // La probabilite est de : evalBest/evalCurr. Le personne la mieux
    // adaptee a sa survie assuree :
    E_Float oldBest = bestEval;
    bestEval = evalp[jBest-1];
    
    if ( oldBest <= bestEval )
    {
      // Dans ce cas, on a degenere la population
      nbItWithFail++;
    }
    else
    {
      // Ah, on a une meilleurs population :
      nbItWithFail = 0;
    }
    E_Int nbSurvivors = 0;
    if ( bestEval == 0. ) break;
    for ( E_Int j = 1; j <= sizeOfPopulation; j++ )
    {
      if ( drand48() > bestEval/evalp[j-1] )
      { 
        evalp[j-1] = -1;
      }
      else
      {
        survivorsp[nbSurvivors++] = j;
      }
    }
    
    // Les participants vires sont remplaces par des clones
    // des survivants :
    for ( E_Int j = 1; j <= sizeOfPopulation; j++ )
      if ( evalp[j-1] == -1 )
      {
	E_Int indClone  = survivorsp[E_Int((1.*nbSurvivors*lrand48())/
                                           (E_Float(MAXRAND)))];
	for ( E_Int i = 0; i < glbMeshes.size(); i++ )
	{
	  nicePeople(i,j) = nicePeople(i,indClone);
	}
      }
    FldArrayI uglyPeople(nicePeople);
    // Croisement de la population de nice people
    for (E_Int j = 1; j <= sizeOfPopulation; j++)
    {
      if ( j != jBest )
      {
        // tirage au hasard de papa et de maman
        E_Int indPapa  =
        survivorsp[E_Int((1.*nbSurvivors*lrand48())/(MAXRAND+1.0))];
        E_Int indMaman =
        survivorsp[E_Int((1.*nbSurvivors*lrand48())/(MAXRAND+1.0))];
	for ( E_Int i = 0; i < glbMeshes.size(); i++ )
	  if ( drand48() >= 0.5 )
	  {
	    nicePeople(i,j) = uglyPeople(i,indPapa);
	  }
	  else
	  {
	    nicePeople(i,j) = uglyPeople(i,indMaman);
	  }
      }
      else
      {
	for ( E_Int i = 0; i < glbMeshes.size(); i++ )
	{
	  nicePeople(i,j) = uglyPeople(i,j);
	}
      }
    }
    // Malheureusement, la piscine etait radioactive. Papa et Maman
    // sont irradies, et vont peut-etre muter. A remarquer que le plus
    // adapte ne pourra pas muter, ce qui est souhaitable...
    for ( E_Int j = 1; j <= sizeOfPopulation; j++ )
    {
      if ( j != jBest )
      {
	for ( E_Int i = 0; i < glbMeshes.size(); i++ )
	{
	  // 10% chance de muter :
	  if ( drand48() < 0.1 )
	  {
	    E_Int targetProc = E_Int((1.*nbProc*lrand48())/(MAXRAND+1.0));
	    nicePeople(i,j) = targetProc;
	  }
	}
      }
    }
    
    // On reevalue l'adaptation des survivants :
    for ( E_Int j = 1; j <= sizeOfPopulation; j++ )
    {
      nbPtsPerProcs = nbNodePerProc;
      nbPtsPerProcsp = nbPtsPerProcs.begin();
      E_Int ind = 0;
      for ( list<DesMesh*>::iterator itMesh = glbMeshes.begin();
            itMesh != glbMeshes.end(); itMesh++ )
      {
        nbPtsPerProcsp[nicePeople(ind,j)] += (*itMesh)->getI(KEY_IM)*
          (*itMesh)->getI(KEY_JM)*(*itMesh)->getI(KEY_KM);
        ind++;
      }
      evalp[j-1] = 0.;
      for ( E_Int p = 0; p < nbProc; p++ )
	evalp[j-1] += (nbPtsPerProcsp[p]-meanPtsPerProc)*
          (nbPtsPerProcsp[p]-meanPtsPerProc);
      evalp[j-1] = sqrt(evalp[j-1]);
      ind = 0;
      for (list<SplApi::ConnectInfo*>::iterator itrConnect = _connect.begin();
           itrConnect != _connect.end(); itrConnect++)
      {
        SplApi::ConnectInfo* c = *itrConnect;
        while (c != NULL)
        {
          // le voisin est-il sur le meme processeur?
          if ( nicePeople(c->neighbour-1,j) != nicePeople(ind,j) )
          {
            evalp[j-1] += latence;
            evalp[j-1] += comSpeed*c->size;
          }
          c = c->next;
        }
        ind++;
      }

      // Check for empty processors
      for (E_Int i = 0; i < nbProc; i++)
      {
        if (nbPtsPerProcsp[i] == 0) evalp[j-1] += 1.e6;
      }
    }
   
  } // loop
  
  
  ind = 0;
  for ( list<DesMesh*>::iterator itMesh = glbMeshes.begin();
        itMesh != glbMeshes.end(); itMesh++ )
  {
    (*itMesh)->setI(KEY_NODE, nicePeople(ind,jBest));
    ind++;
  }
  
  // Calcul du nombre de blocs par processeur
  FldArrayI NbBlocksPerProc(nbProc);
  NbBlocksPerProc.setAllValuesAt(0);
  for ( list<DesMesh*>::iterator itMesh = glbMeshes.begin();
        itMesh != glbMeshes.end(); itMesh++ )
  {
    NbBlocksPerProc[(*itMesh)->getI(KEY_NODE)]++;
  }
  for ( list<DesMesh*>::iterator itMesh = locMeshes.begin();
        itMesh != locMeshes.end(); itMesh++ )
  {
    NbBlocksPerProc[(*itMesh)->getI(KEY_NODE)]++;
  }

  E_Int number = 0;
  for ( list<DesBlock*>::iterator itBlock = _blocksToSplit.begin();
        itBlock != _blocksToSplit.end(); itBlock++ )
  {
    if ((*itBlock)->getNumber() > number)
      number = (*itBlock)->getNumber();
  }

  // Check for empty processors
  for (E_Int i = 0; i < nbProc; i++)
  {
    if (NbBlocksPerProc[i] == 0)
    {
      if (PcmTask::getTask()->myNode() == 0)
        cerr << "Warning: the processor "<<i<<" is empty."<<endl;
    }
  }

  // Calcul des volumes de communication
  nbPtsPerProcs = nbNodePerProc;
  nbPtsPerProcsp = nbPtsPerProcs.begin();
  ind = 0;
  for ( list<DesMesh*>::iterator itMesh = glbMeshes.begin();
        itMesh != glbMeshes.end(); itMesh++ )
  {
    nbPtsPerProcsp[nicePeople(ind,jBest)] += (*itMesh)->getI(KEY_IM)*
      (*itMesh)->getI(KEY_JM)*(*itMesh)->getI(KEY_KM);
    ind ++;
  }
  
  E_Float charge = 0.;
  for ( E_Int p = 0; p < nbProc; p++ )
    charge += (nbPtsPerProcsp[p]-meanPtsPerProc)*
      (nbPtsPerProcsp[p]-meanPtsPerProc);
  charge = sqrt(charge);
  ind = 0;
  E_Float comNumber = 0;
  E_Float comVolume = 0.;
  for (list<SplApi::ConnectInfo*>::iterator itrConnect = _connect.begin();
       itrConnect != _connect.end(); itrConnect++)
  {
    SplApi::ConnectInfo* c = *itrConnect;
    while (c != NULL)
    {
      // le voisin est-il sur le meme processeur?
      if ( nicePeople(c->neighbour-1,jBest) != nicePeople(ind,jBest) )
      {
        comNumber += latence;
        comVolume += comSpeed*c->size;
      }
      c = c->next;
    }
    ind++;
  }

  // Ecriture sortie
  if (PcmTask::getTask()->myNode() == 0)
    printf("Best config: charge = %f, latence = %f, comVol = %f, F = %f.\n",
           charge, comNumber/latence, comVolume/comSpeed, evalp[jBest-1]);
}

