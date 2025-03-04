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

#include "Connect/connect.h"

using namespace K_FLD;
using namespace std;
//=============================================================================
/* orderBAR2Struct: ordonne une BAR selon un i-array structure. La BAR
   ne doit pas etre ramifiee et ne doit pas contenir de doublons (i.e.
   cleanConnectivity doit etre faite avant).
   Attention: fout doit etre allouee en externe */
//=============================================================================
void K_CONNECT::orderBAR2Struct(E_Int posx, E_Int posy, E_Int posz, 
                                FldArrayF& f, FldArrayI& cn, 
                                FldArrayF& fout)
{
  E_Int nelts = cn.getSize(); // connectivite BAR recuperee avec getConnect
  E_Int nfld = f.getNfld();
  FldArrayIS dejaVu(nelts); dejaVu.setAllValuesAtNull();
  E_Int n=0;
  E_Int indp=0, indn=0;
  // si pas de boucle: recherche du pt de depart
  short found1=0, found2=0;
  // recherche d'un pt n appartenant pas a deux elements
  for (E_Int et1 = 0; et1 < nelts; et1++)
  {
    found1 = 0; found2 = 0;
    E_Int ind11 = cn(et1,1); E_Int ind21 = cn(et1,2);
    for (E_Int et2 = 0; et2 < nelts; et2++)
    {
      if (et1 != et2) 
      {
        E_Int ind12 = cn(et2,1); E_Int ind22 = cn(et2,2);
        if (ind11 == ind12 || ind11 == ind22) found1 = 1;   
        if (ind21 == ind12 || ind21 == ind22) found2 = 1;   
      }
    }
    
    if (found1 == 0){indp = ind11-1; indn = ind21-1; dejaVu[et1]=1; break;}
    else if (found2 == 0) {indp = ind21-1; indn = ind11-1; dejaVu[et1]=1; break;}
  }
  if (found1 == 1 && found2 == 1) {indp = cn(0,1)-1;indn = cn(0,2)-1; dejaVu[0] = 1;}// boucle
  for (E_Int eq = 1; eq <= nfld; eq++)
  {fout(n,eq) = f(indp,eq); fout(n+1,eq) = f(indn,eq);}
  n = n+2;
  indp = indn;
  E_Int npts2 = fout.getSize();
  short found = 1;

  while (found == 1 && n < npts2) 
  {
    found = 0;
    for (E_Int et = 0; et < nelts; et++)
    {
      if (dejaVu[et] == 0) 
      {
        if (cn(et,1) == indp+1)  
        {
          indn = cn(et,2)-1; indp = indn; found = 1; dejaVu[et]=1; 
          for (E_Int eq = 1; eq <= nfld; eq++) fout(n,eq) = f(indn,eq); 
          n++;
          goto next;
        }
        else if (cn(et,2) == indp+1)
        {
          indn = cn(et,1)-1; indp = indn; found = 1; dejaVu[et]=1; 
          for (E_Int eq = 1; eq <= nfld; eq++) fout(n,eq) = f(indn,eq); 
          n++;
          goto next;
        }
      }
    }
    next:;
  }
  fout.reAllocMat(n, nfld);
}

//=============================================================================
/* Ordonnancement des elements de type bar
   IN: cBN: connectivite elt->noeuds non ordonnee
   IN: field: champ defini aux noeuds  non ordonnee
   OUT: fieldout: field ordonne
   OUT: cBNout: connectivite elt->noeud ordonne
   les tableaux de sortie doivent etre alloues en externe
*/
//=============================================================================
void K_CONNECT::orderBAR(FldArrayF& field, FldArrayI& cBN,
                         FldArrayF& fieldout, FldArrayI& cBNout)
{
  E_Int nnodes = field.getSize();

  //creation de la connectivite noeud -> element
  vector< vector<E_Int> > cNB(nnodes);
  connectEV2VE(cBN, cNB);

  E_Int nbars = cBN.getSize();
  FldArrayIS dejaVu(nbars); dejaVu.setAllValuesAtNull();
  E_Int nfld = field.getNfld();

  //initialisation
  E_Int ind0 = 0; 
  E_Int elt0 = (cNB[ind0])[0];
  for (E_Int eq = 1 ; eq <= nfld; eq++)  
    fieldout(0,eq) = field(0,eq);
  cBNout(0,1) = 1;

  E_Int n1, n2;

  E_Int ib = 0; // compteur bar
  E_Int in = 1; // compteur noeuds

  while ( dejaVu[elt0] == 0 ) 
  {
    //recherche du 2eme noeud de elt0 
    n1 = cBN(elt0,1)-1; n2 = cBN(elt0,2)-1;
    if (ind0 == n1) ind0 = n2;
    else if (ind0 == n2) ind0 = n1;

    //mise a jour des tableaux
    cBNout(ib,2) = ind0+1;
    for (E_Int eq = 1; eq <= nfld; eq++)
      fieldout(in,eq) = field(ind0,eq);
    ib++; in++;
    dejaVu[elt0] = 1;

    if (ib >= nbars) return;
    
    // recherche du suivant
    vector<E_Int>& elts = cNB[ind0];
    elt0 = -1;
    E_Int eltsSize = elts.size();
    for (E_Int ie = 0; ie < eltsSize; ie++)
    {
      E_Int eltt = elts[ie];
      if (dejaVu[eltt] == 0) 
      {
        elt0 = eltt;
        cBNout(ib,1) = ind0;
        break;
      }
    }
    
    if (elt0 == -1) // pas d element suivant trouve
    {
      // il y a une ramification 
      //recherche d une ramification
      for (E_Int ni = 0; ni < nnodes; ni++)
      {
        vector<E_Int>& eltsni = cNB[ni];
        E_Int prod = 1;
        E_Int sum = 0;
        E_Int eltn = -1;
        E_Int eltsniSize = eltsni.size();
        for (E_Int ei = 0; ei < eltsniSize; ei++)
        {
          E_Int elti = eltsni[ei]; 
          E_Int dj = dejaVu[elti];
          prod = prod * dj;
          sum = sum + dj;
          if (eltn == -1 && dj == 0) // stockage du premier non deja traite
            eltn = elti;
        }
        
        if (prod == 0 && sum != 0) //ramification
        {
          ind0 = ni;
          elt0 = eltn;
          cBNout(ib,1) = ind0;
          for (E_Int eq = 1; eq <= nfld; eq++)
            fieldout(in,eq) = field(ind0,eq);
          in++;
          dejaVu[elt0] = 1;
        }
      }
    } // fin de recherche de ramification
  }
}