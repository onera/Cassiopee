/*    
    Copyright 2013-2025 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/

# include "Algorithm.h"
# include <stdlib.h>

using namespace std;
using namespace K_FLD;
using namespace K_CONST;
//=============================================================================
// Compute iblank field on grids
//=============================================================================
void computeIBlank(vector<StructBlock*>& structBlocks)
{
  E_Int nzone = structBlocks.size();
  printf("Eliminating overlapped regions..."); fflush(stdout);
  FldArrayIS dejaVu(nzone, nzone);
  dejaVu.setAllValuesAtNull();
  
  // Construit la liste des blocs recouvrants le bloc i
  for (E_Int i = 0; i < nzone; i++)
    structBlocks[i]->selectTypeOfBlks(structBlocks);

  // Construit la liste des points "interpolable" du bloc i
  // ainsi que les cellules d'interpolation
  for (E_Int i = 0; i < nzone; i++)
    structBlocks[i]->compInterpolationInformation(structBlocks);
 
  // Verifie certaines cellules pour certains cas particuliers
  //for (E_Int i = 0; i < nzone; i++)
  //  structBlocks[i]->checkValidityOfInterpolationCell(i, structBlocks);

  for (E_Int i = 0; i < nzone ; i++)
    structBlocks[i]->cleanListsOfBlks(structBlocks);
  
  for (E_Int i = 0; i < nzone; i++)
  {
    for ( E_Int j = 0; j < nzone; j++)
    {
      if ( i != j && dejaVu(i,j+1) == 0 )
      {
        E_Int Nij = structBlocks[i]->getNbOfInterpolablePtsForBlk(j);
        E_Int Nji = structBlocks[j]->getNbOfInterpolablePtsForBlk(i);
        
        // sufficient to say if it is not its own block 
        if ( Nij != 0 && Nji != 0 )
        {
          printf(".");
          
          if ( Nij <= Nji ) // blanking pts of the coarser mesh i
          {
            structBlocks[i]->compIBlank(structBlocks, j);
            structBlocks[i]->addOverlappingBlk(structBlocks[j]);
            structBlocks[j]->addOverlappingBlk(structBlocks[i]);
          }
          else if ( Nij > Nji ) // blanking pts of the coarser mesh j
          {
            structBlocks[j]->compIBlank( structBlocks, i);
            structBlocks[i]->addOverlappingBlk(structBlocks[j]);
            structBlocks[j]->addOverlappingBlk(structBlocks[i]);
          }
         
          // Tag that blanking is already done for the pair ij
          dejaVu(i,j+1) = 1;
          dejaVu(j,i+1) = 1;
        }
      }
    }
  }
  printf("done.\n");

  /*-------------------------------------*/
  /* Temporary save : overlap eliminated */
  /*-------------------------------------*/
  printf("Saving blanked file blanked.tp..."); fflush(stdout);
  E_Bool add0 = false;
  for (E_Int i = 0; i < nzone; i++)
  {
    printf("."); fflush(stdout);
    structBlocks[i]->write("blanked.tp", add0);
    add0 = true;
  }
}

//=============================================================================
void computeStringing(vector<StructBlock*>& structBlocks,
                      vector<CString*>& strings)
{
  E_Int nzone = structBlocks.size();

  printf("Identifying stringed points..."); fflush(stdout);
  for (E_Int i = 0; i < nzone; i++)
  {
    structBlocks[i]->compIBlankCC();
    structBlocks[i]->updateIBlankArray();
  }
  
  for (E_Int i = 0; i < nzone; i++)
    structBlocks[i]->identifyGapBndPts( i, structBlocks);
  printf("done.\n");
  
  /*-----------*/
  /* Stringing */
  /*-----------*/
  printf("Ordering stringed points..."); fflush(stdout);
  for (E_Int i = 0; i < nzone; i++)
  {
    printf("."); fflush(stdout);
    structBlocks[i]->stringing(structBlocks[i]->getStrings());
  }
  printf("done.\n");
  
  E_Bool add = false;
  for (E_Int i = 0; i < nzone; i++)
  {
    list<FldArrayI*>& st = structBlocks[i]->getStrings();
    for (list<FldArrayI*>::iterator itr = st.begin();
         itr != st.end(); itr++)
    {
      structBlocks[i]->writeLine(**itr, "strings.tp", add);
      add = true;
    }
  }
  printf("File strings.tp written.\n");

  /*-------------------------------*/
  /* Creation des strings globales */
  /*-------------------------------*/
  printf("Creating global strings..."); fflush(stdout); 
  for (E_Int i = 0; i < nzone; i++)
  {
    list<FldArrayI*>& st = structBlocks[i]->getStrings();    
    for (list<FldArrayI*>::iterator itr = st.begin();
         itr != st.end(); itr++)
    {
      CString* cs = new CString(**itr, structBlocks[i]);
      strings.push_back(cs);
      delete *itr;
    }
    st.clear();
  }
  printf("done.\n");
}

//=============================================================================
void computeMatchingSegments(vector<CString*>& strings,
                             vector<SegmentPair*>& segPairs)
{
  SegmentPair* onePair;
  E_Int ns = strings.size();
  
  printf("Searching for matching segment..."); fflush(stdout);
  for (E_Int i = 0 ; i < ns; i++)
    strings[i]->selectStrings(strings, i);
 
  for ( E_Int i = 0; i < ns; i++)
  {
    printf("."); fflush(stdout);
    vector<CString*> strOut1;
    vector<CString*> strOut2;
    strings[i]->searchForMatchingStrings( strings, strOut1, strOut2);
    
    if (strOut1.size() != strOut2.size())
    {
      printf("Error: zip: segments of strings are not set by pairs.\n"); 
      exit(0);
    }
    
    // Build segment pairs
    E_Int strOut1Size = strOut1.size();
    for (E_Int j = 0; j < strOut1Size; j++)
    {
      onePair = new SegmentPair( strOut1[j], strOut2[j]);
      segPairs.push_back(onePair);
    }

    for(E_Int v = 0; v < strOut1Size; v++)
    {
      delete strOut1[v];
      delete strOut2[v];
    }
  }
  
  printf("done.\n");

  E_Int nSeg = segPairs.size();
  E_Bool add = false;
  for (E_Int i = 0; i < nSeg; i++)
  {
    FldArrayI& s1 = segPairs[i]->getIndArray1();
    FldArrayI& s2 = segPairs[i]->getIndArray2();
  
    segPairs[i]->getBlock1()->writeLine(s1, "segments.tp", add);
    add = true;
    segPairs[i]->getBlock2()->writeLine(s2, "segments.tp", add);
  }
}

//=============================================================================
// Zip zones between pairs of segments.
//In: segPairs
//Out : field : champ (coord+solution) sur les triangles
//      triconnect : connectivite triangles  
//=============================================================================
void zipInbetweenPairs(vector<SegmentPair*>& segPairs,
                       vector<FldArrayF*>& field,
                       vector<FldArrayI*>& triConnect,
                       FldArrayB& isZipped)
{
  printf("Zipping..."); fflush(stdout);
  vector<FldArrayF*> listOfEdges;
  E_Int nSeg = segPairs.size();
  E_Int nfieldTot = 0;
  if ( nSeg > 0 )
    nfieldTot = segPairs[0]->getBlock1()->getNbOfFields();
  else 
    return;

  for (E_Int i = 0; i < nSeg; i++)
  {
    FldArrayF* tf = new FldArrayF();
    FldArrayI* ti = new FldArrayI();
    field.push_back(tf);
    triConnect.push_back(ti);
  }
  
  for (E_Int i = 0; i < nSeg; i++)
  {
    isZipped[i] =
      segPairs[i]->computeZipper( nfieldTot, *field[i], *triConnect[i] );
    printf("."); fflush(stdout);
  }
  printf("done.\n");
}

//=============================================================================
void closePockets(vector<CString*>& strings, 
                  vector<SegmentPair*>& segPairs,
                  vector<FldArrayF*>& field,
                  vector<FldArrayI*>& triConnect)
{
  if (strings.size() == 0 )
    return;

  E_Int nSeg = segPairs.size();
  E_Float matchTol = strings[0]->getBlock()->getMatchTol();
    
  printf("Closing pockets..."); fflush(stdout);

  // 0-Creation
  vector<SingleSegment*> singleSegments;
  vector<FldArrayF*> vectOfLastEdges;

  // 1-Update dejaVu to begin search for remaining segments
  for (E_Int i = 0; i < nSeg; i++)
  {
    segPairs[i]->updateFlagForStrings(segPairs);  
    segPairs[i]->identifyLastEdges(vectOfLastEdges);
  }

  // 2-Create list of remaining segments for non treated points
  E_Int stringsSize = strings.size();
  for (E_Int i = 0; i < stringsSize; i++)
  {
    list<FldArrayI*> remainingSegs;
    strings[i]->compRemainingSegments(remainingSegs);
    StructBlock* blk = strings[i]->getBlock();
    if ( remainingSegs.size() != 0 )
    {
      //build list of remaining seg
      for (list<FldArrayI*>::iterator itr = remainingSegs.begin();
           itr != remainingSegs.end(); itr++)
      {
        SingleSegment* ss = new SingleSegment(**itr, blk);
        singleSegments.push_back(ss);
        delete *itr;
      }
    }
  }
  
  // 3-Create list of remaining segments for extremities of segment pairs
  FldArrayIS dejaVu;
  dejaVu.malloc(vectOfLastEdges.size());
  dejaVu.setAllValuesAtNull();
  E_Int vectOfLastEdgesSize = vectOfLastEdges.size();
  for (E_Int v1 = 0; v1 < vectOfLastEdgesSize; v1++)
  {
    if (dejaVu[v1] == 0)
    {
      FldArrayF& field1 = *(vectOfLastEdges[v1]);     
      for (E_Int v2 = 0; v2 < vectOfLastEdgesSize; v2++)
      {
        if (dejaVu[v2] == 0 && v2 != v1)
        {      
          FldArrayF& field2 = *(vectOfLastEdges[v2]);     
          E_Bool isId = testIfEdgesAreMatching(matchTol, field1, field2);
          if (isId == true)
          {
            dejaVu[v1] = 1;
            dejaVu[v2] = 1;
            goto next;
          }
        }
      }
      SingleSegment* ss = new SingleSegment(field1);
      singleSegments.push_back(ss);
    }
    next:;
  }
  dejaVu.malloc(0);
  
  for (E_Int v = 0; v < vectOfLastEdgesSize; v++)
    delete vectOfLastEdges[v];
  vectOfLastEdges.clear();

  // 4-write the last edges in strings2.tp
  writeSingleSegments(singleSegments);

  vector<Pocket*> pockets;
  buildPocket(matchTol, singleSegments, pockets);
  
  // 5-triangulate pockets
  E_Int pocketsSize = pockets.size();
  for (E_Int i = 0; i < pocketsSize; i++)
  {
    FldArrayF* tf = new FldArrayF();
    FldArrayI* ti = new FldArrayI();
    field.push_back(tf);
    triConnect.push_back(ti);
  }
  
  FldArrayB closed(pocketsSize);
  for (E_Int p = 0; p < pocketsSize; p++)
    closed[p] = pockets[p]->closePocket( *field[p+nSeg], 
                                         *triConnect[p+nSeg]);

  printf("done.\n");
  
  // write pockets
  E_Bool add = false;
  for (E_Int p = 0; p < pocketsSize; p++)
  {
    pockets[p]->writeLine((char*)"pockets.tp", add);
    delete pockets[p];
    add = true;
  }
  pockets.clear();
  printf("File pockets.tp written.\n"); 
}

//=============================================================================
void mergeAllZonesInOne(vector<StructBlock*>& structBlocks,
                        vector<FldArrayF*>& field,
                        vector<FldArrayI*>& idgT, 
                        FldArrayI& idgG, FldArrayF& fieldG,
                        FldArrayI& FirstPt)
{
  E_Int nzone = structBlocks.size();
  FldArrayI NbBndPt(nzone);
  FldArrayI FirstBndPt(nzone);
  E_Int NbBndPtTot ;
  FldArrayI NbPt(nzone);
  E_Int NbPtTot;
  E_Int size = field.size();
  assert( nzone > 0); 
  E_Int nfieldTot = structBlocks[0]->getNbOfFields();
  E_Float matchTol = structBlocks[0]->getMatchTol();

  printf("Computing zone boundaries, degenerations and connectivity...");
  fflush(stdout);
  for (E_Int i = 0 ; i < nzone ; i++)
  {
    structBlocks[i]->compIBndIDgUnsConnectEN();
    printf("."); fflush(stdout);
  }
  printf("done.\n");
  
  /* Calcul de la taille du tableau global des noeuds
     et de l'indice du premier noeud de chaque zone dans ce tableau */
  NbPtTot = 0;
  for (E_Int i = 0; i < nzone; i++)
  {
    FirstPt[i] = NbPtTot;
    NbPt[i] = structBlocks[i]->getNbOfMeshPts();
    NbPtTot = NbPtTot + NbPt[i];
  }
  
  /*-----------------------------*/
  /* TABLEAU DES NOEUDS A MERGER */
  /*-----------------------------*/
  /* Calcul de la taille du tableau des noeuds a merger
     et de l'indice du premier noeud de chaque zone dans ce tableau */
  NbBndPtTot = 0;
  for (E_Int i = 0; i < nzone; i++)
  {
    FirstBndPt[i] = NbBndPtTot;
    NbBndPt[i] = structBlocks[i]->getIBndSize();
    NbBndPtTot = NbBndPtTot + NbBndPt[i];
  }
  /* Numerotation globale des tableaux _ibnd[i]*/
  for (E_Int i = 0; i < nzone; i++)
    structBlocks[i]->incrIBnd(FirstPt[i]);
  
  /* Creation du tableau global des points a merger */
  FldArrayI ibndG(NbBndPtTot);
  for (E_Int i = 0; i < nzone; i++)
    ibndG.copyFrom(FirstBndPt[i], structBlocks[i]->getIBnd());

  /* TABLEAU DES DEGENERESCENCES
     Numerotation globale des tableaux _idg[i]*/
  for (E_Int i = 0; i < nzone; i++)
    structBlocks[i]->incrIDg(FirstPt[i]);
  
  /* Creation du tableau de degenerescence global*/
  idgG.malloc(NbPtTot);
  for (E_Int i = 0; i < nzone; i++)
    idgG.copyFrom(FirstPt[i], structBlocks[i]->getIDg());

  /* Creation du tableau de field (coord + CdfField) global : fieldG */
  fieldG.malloc(NbPtTot, nfieldTot);
  for (E_Int i = 0; i < nzone; i++)
    fieldG.copyFrom(FirstPt[i], structBlocks[i]->getGlobalField());

  /*--------------------------*/
  /* MERGE DES POINTS DE IBND */
  /*--------------------------*/
  printf("Merging zones..."); fflush(stdout);
  E_Float xref, yref, zref;
  E_Float dx, dy, dz;
  E_Float D;
  E_Int ji, test;
  E_Float* fieldGx = fieldG.begin(1);
  E_Float* fieldGy = fieldG.begin(2);
  E_Float* fieldGz = fieldG.begin(3);
  
  for (E_Int i = 1; i < NbBndPtTot; i++)
  {
    E_Int indi = ibndG[i];
    xref = fieldGx[indi];
    yref = fieldGy[indi];
    zref = fieldGz[indi];

    ji = 0;
    test = 0;
    while (ji < i && test == 0)
    {
      E_Int indji = ibndG[ji];
      dx = fieldGx[indji] - xref;
      dy = fieldGy[indji] - yref;
      dz = fieldGz[indji] - zref;
      D = sqrt( dx*dx + dy*dy + dz*dz );
     
      if ( D < matchTol)
      {
        idgG[indi] = idgG[indji];
        test = 1;
        if ((i%10) == 0) { printf("."); fflush(stdout); }
      }
      ji++;
    }
  }
  printf("done.\n");
  
  /*---------------------------------------*/
  /* Merge des points des zones de raccord */
  /*---------------------------------------*/
  printf("Merging zipped points..."); fflush(stdout);
  // Remember : the kth elements are segments (of strings) or  pockets
 
  
  for (E_Int k = 0; k < size; k++)
  {
    FldArrayI* td = new FldArrayI();  
    idgT.push_back(td); 
    E_Int s = field[k]->getSize();
    (*idgT[k]).malloc(s);

    E_Float* xt = field[k]->begin(1);
    E_Float* yt = field[k]->begin(2);
    E_Float* zt = field[k]->begin(3);
    
    for (E_Int i = 0; i < s; i++)
    {
      xref = xt[i];
      yref = yt[i];
      zref = zt[i];
      ji = 0;
      test = 0;
     
      while ( ji < NbBndPtTot && test == 0 )
      {
        E_Int indji = ibndG[ji];
        dx = fieldGx[indji] - xref;
        dy = fieldGy[indji] - yref;
        dz = fieldGz[indji] - zref;
        D = sqrt( dx*dx + dy*dy + dz*dz );
 
        if (D < matchTol)
        {
          (*idgT[k])[i] = idgG[indji];
          test = 1;
          if ((i%100) == 0) { printf("."); fflush(stdout); }
        }
        ji++;
      }
    }
    delete field[k];
  }
  field.clear();
  printf("done.\n");  
}

//=============================================================================
void computeConnectivity(vector<StructBlock*>& structBlocks,
                         vector<FldArrayI*>& triConnect,
                         vector<FldArrayI*>& idgT,
                         FldArrayI& idgG, FldArrayI& unsConnectENG,
                         FldArrayI& FirstPt)
{
  E_Int nzone = structBlocks.size();
  E_Int size = triConnect.size();
  // Construction de la connectivite globale a partir des StructBlocks
  printf("Building connectivity..."); fflush(stdout);
  FldArrayI NbElmt(nzone);
  FldArrayI FirstElmt(nzone);
  E_Int NbElmtTot = 0;

  for (E_Int i = 0; i < nzone; i++)
  {
    FirstElmt[i] = NbElmtTot;
    NbElmt[i] = structBlocks[i]->getUnsConnectENSize();
    NbElmtTot = NbElmtTot + NbElmt[i];
  }
  
  unsConnectENG.malloc(NbElmtTot, 3);
  for (E_Int i = 0; i < nzone; i++)
  {
    structBlocks[i]->incrUnsConnectEN(FirstPt[i]);
    unsConnectENG.copyFrom(FirstElmt[i], structBlocks[i]->getUnsConnectEN());
  }

  // Prise en compte des degenerescences globales
  // ATTENTION : incrementation des indices pour Tecplot NON effectuee
  for (E_Int j = 1; j <= 3; j++)
    for (E_Int i = 0; i < NbElmtTot; i++)
      unsConnectENG(i,j) = idgG[unsConnectENG(i,j)];
  
  // Construction de la connectivite globale : ajout zones de raccord
  // ATTENTION : incrementation des indices pour Tecplot NON effectuee
  FldArrayI NbElmtTri(size);
  FldArrayI FirstElmtTri(size);
  for (E_Int k = 0; k < size; k++)
  {
    FirstElmtTri[k] = NbElmtTot;
    NbElmtTri[k] = (*triConnect[k]).getSize();
    NbElmtTot = NbElmtTot + NbElmtTri[k];
  }
  unsConnectENG.reAllocMat(NbElmtTot, 3);
  E_Int* cn1 = unsConnectENG.begin(1);
  E_Int* cn2 = unsConnectENG.begin(2);
  E_Int* cn3 = unsConnectENG.begin(3);

  for (E_Int k = 0; k < size; k++)
  {
    E_Int firstTri = FirstElmtTri[k];
    E_Int* tricn1 = triConnect[k]->begin(1);
    E_Int* tricn2 = triConnect[k]->begin(2);
    E_Int* tricn3 = triConnect[k]->begin(3);
    
    for (E_Int i = 0; i < NbElmtTri[k]; i++)
    {
      E_Int noTri = firstTri + i;
      E_Int noTri1 = tricn1[i]-1;
      E_Int noTri2 = tricn2[i]-1;
      E_Int noTri3 = tricn3[i]-1;
      
      cn1[noTri] = (*idgT[k])[noTri1];
      cn2[noTri] = (*idgT[k])[noTri2];
      cn3[noTri] = (*idgT[k])[noTri3];
    }
    delete triConnect[k];
    delete idgT[k];
  }
  
  triConnect.clear();
  idgT.clear();
  printf("done.\n");
}

//=============================================================================
void deletingUnusedNodes(FldArrayI& unsConnectENG, 
                         FldArrayF& fieldG)
{
  E_Int NbPtTot = fieldG.getSize();
  E_Int nfieldTot = fieldG.getNfld();
  E_Int NbElmtTot = unsConnectENG.getSize();
  printf("Compacting nodes numbers..."); fflush(stdout);
  // Numerotation des noeuds utilises
  // + Incrementation des noeuds pour Tecplot OK
  FldArrayI compactNodes(NbPtTot, 2);
  compactNodes.setAllValuesAtNull();
  E_Int* compactN1 = compactNodes.begin(1);
  E_Int* compactN2 = compactNodes.begin(2);

  for (E_Int j = 1; j <= 3; j++)
  {
    E_Int* cnj = unsConnectENG.begin(j);
    for (E_Int i = 0; i < NbElmtTot; i++)
      compactN1[cnj[i]] = 1;
  }  
  compactN2[0] = compactN1[0];

  for (E_Int i = 1; i < NbPtTot; i++)
    compactN2[i] = compactN2[i-1] + compactN1[i];

  E_Int jj = 0;
  // Suppression des champs des noeuds inutilises
  for (E_Int k = 1; k <= nfieldTot; k++)
  {
    jj = 0;
    E_Float* fieldGk = fieldG.begin(k);
    
    for (E_Int i = 0; i < NbPtTot; i++ )
    {
      if (compactN1[i] == 1)
      {
        fieldGk[jj] = fieldGk[i];
        jj++;
      }
    }
  }
  NbPtTot = jj ;
  fieldG.reAllocMat(NbPtTot, nfieldTot);
  // MAJ de la table de connectivite
  for (E_Int i = 0; i < NbElmtTot; i++)
    for (E_Int j = 1; j <= 3; j++)
      unsConnectENG(i, j) = compactN2[unsConnectENG(i,j)];
  printf("done.\n");
}
