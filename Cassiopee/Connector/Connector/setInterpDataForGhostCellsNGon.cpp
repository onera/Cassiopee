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
# include <unordered_map>
# include "connector.h"
using namespace K_FLD;

//=============================================================================
/* Returns the indices of ghost points (cells) and the indices of 
   donor points according to the 1-to-1 grid connectivity for NGon grids */
//=============================================================================
PyObject* K_CONNECTOR::setInterpDataForGCNGon(PyObject* self, PyObject* args)
{
  PyObject* FL; // face list of match
  PyObject* FLDonor; // face list donor
  PyObject* rindElt; // rind of elements (couche 0, 1,2)
  PyObject* rindEltDonor; // rind of elements donor (couche 0, 1,2)
  
  PyObject* PE; // Parent-element
  PyObject* PEDonor; // PE donor

  PyObject* array; // array2 of zone
  PyObject* arrayDonor; // array2 of zone donor

  if (!PYPARSETUPLE_(args, OOOO_ OOOO_,
                    &FL, &FLDonor,
                    &rindElt, &rindEltDonor,
                    &array, &arrayDonor, 
                    &PE, &PEDonor))
  {
    return NULL;
  }
  
  FldArrayI* FLI;
  E_Int res = K_NUMPY::getFromNumpyArray(FL, FLI);
  if (res == 0) 
  {    
    PyErr_SetString(PyExc_TypeError, 
                    "setInterpDataForGhostCells: 1st arg is not a valid numpy array.");
    return NULL;
  }
  E_Int* FLp = FLI->begin();
  FldArrayI* FLDonorI;
  res = K_NUMPY::getFromNumpyArray(FLDonor, FLDonorI);
  if (res == 0)
  {    
    RELEASESHAREDN(FL, FLI);
    PyErr_SetString(PyExc_TypeError, 
                    "setInterpDataForGhostCells: 2nd arg is not a valid numpy array.");
    return NULL;
  }
  E_Int* FLDonorp = FLDonorI->begin();

  FldArrayI* Intext;
  res = K_NUMPY::getFromNumpyArray(rindElt, Intext);
  if (res == 0)
  {    
    RELEASESHAREDN(FL, FLI);
    RELEASESHAREDN(FLDonor, FLDonorI);
    PyErr_SetString(PyExc_TypeError, 
                    "setInterpDataForGhostCells: 3nd arg is not a valid numpy array.");
    return NULL;
  }
  E_Int* iptIntExt = Intext->begin();

  FldArrayI* IntextD;
  res = K_NUMPY::getFromNumpyArray(rindEltDonor, IntextD);
  if (res == 0)
  {    
    RELEASESHAREDN(FL, FLI);
    RELEASESHAREDN(FLDonor, FLDonorI);
    RELEASESHAREDN(rindElt, Intext);
    PyErr_SetString(PyExc_TypeError, 
                    "setInterpDataForGhostCells: 4nd arg is not a valid numpy array.");
    return NULL;
  }
  E_Int* iptIntExtD= IntextD->begin();


  K_FLD::FldArrayF* f; K_FLD::FldArrayI* c;
  char* varString; char* eltType; E_Int ni,nj,nk;
  E_Int res1 = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, c, eltType);
  if (res1 == 0 || res1 == 1 || strcmp(eltType, "NGON") != 0)
  {    
    RELEASESHAREDN(FL, FLI);
    RELEASESHAREDN(FLDonor, FLDonorI);
    RELEASESHAREDN(rindElt, Intext);
    RELEASESHAREDN(rindEltDonor, IntextD);
    RELEASESHAREDB(res1, array, f, c);
    PyErr_SetString(PyExc_TypeError, 
                    "setInterpDataForGhostCells: 5th arg must be an NGON array.");
    return NULL;
  }
  E_Int* nface = c->getNFace();
  E_Int* indPH = c->getIndPH();

  K_FLD::FldArrayF* fd; K_FLD::FldArrayI* cd;
  char* varStringd; char* eltTyped; E_Int nid,njd,nkd;
  E_Int res2 = K_ARRAY::getFromArray3(arrayDonor, varStringd, fd, nid, njd, nkd, cd, eltTyped);
  if (res2 == 0 || res2 == 1 || strcmp(eltTyped, "NGON") != 0)
  {    
    RELEASESHAREDN(FL, FLI);
    RELEASESHAREDN(FLDonor, FLDonorI);
    RELEASESHAREDN(rindElt, Intext);
    RELEASESHAREDN(rindEltDonor, IntextD);
    RELEASESHAREDB(res1, array, f, c);
    RELEASESHAREDB(res2, arrayDonor, fd, cd);
    PyErr_SetString(PyExc_TypeError, 
                    "setInterpDataForGhostCells: 6th arg must be an NGON array.");
    return NULL;
  }
  E_Int* nfaceDonor = cd->getNFace();
  E_Int* indPHDonor = cd->getIndPH();

  FldArrayI* PEI;
  res = K_NUMPY::getFromNumpyArray(PE, PEI);
  if (res == 0) 
  {    
    RELEASESHAREDN(FL, FLI);
    RELEASESHAREDN(FLDonor, FLDonorI);
    RELEASESHAREDN(rindElt, Intext);
    RELEASESHAREDN(rindEltDonor, IntextD);
    RELEASESHAREDB(res1, array, f, c);
    RELEASESHAREDB(res2, arrayDonor, fd, cd);
    PyErr_SetString(PyExc_TypeError, 
                    "setInterpDataForGhostCells: 7th arg is not a valid numpy array.");
    return NULL;
  }
  E_Int* PEp1 = PEI->begin(1);
  E_Int* PEp2 = PEI->begin(2);
  FldArrayI* PEDonorI;
  res = K_NUMPY::getFromNumpyArray(PEDonor, PEDonorI);
  if (res == 0) 
  {    
    RELEASESHAREDN(FL, FLI);
    RELEASESHAREDN(FLDonor, FLDonorI);
    RELEASESHAREDN(rindElt, Intext);
    RELEASESHAREDN(rindEltDonor, IntextD);
    RELEASESHAREDN(PE, PEI);
    RELEASESHAREDB(res1, array, f, c);
    RELEASESHAREDB(res2, arrayDonor, fd, cd);
    PyErr_SetString(PyExc_TypeError, 
                    "setInterpDataForGhostCells: 8th arg is not a valid numpy array.");
    return NULL;
  }
  E_Int* PEDonorp1 = PEDonorI->begin(1);
  E_Int* PEDonorp2 = PEDonorI->begin(2);


  /*
  intext[0]  #element couche zero
  intext[1]  #element couche 1 de type raccord
  intext[2]  #element couche 2 de type raccord
  intext[3]  #element couche 1 et 2 de type BC
  */

  E_Int cell_layer0 = iptIntExt[ 0];
  E_Int cell_layer0D= iptIntExtD[0];

                         // El0             El1 rac        EL1 et 2  BC
  E_Int cell_match2_deb = iptIntExt[ 0] + iptIntExt[ 1] + iptIntExt[3];    //EL2 rac
  E_Int cell_match2_fin = iptIntExt[ 0] + iptIntExt[ 1] + iptIntExt[3] + iptIntExt[ 2];

  //E_Int cell_match2_debD= iptIntExtD[ 0] + iptIntExtD[ 1] + iptIntExtD[3];    //EL2 rac
  //E_Int cell_match2_finD= iptIntExtD[ 0] + iptIntExtD[ 1] + iptIntExtD[3] + iptIntExtD[ 2];

  printf("Elt0 = " SF_D_ ",  Elt1= " SF_D_ ", Elt2 = " SF_D_ ", Elt3= " SF_D_ ", Elt4= " SF_D_ " \n", iptIntExt[0], iptIntExt[1], iptIntExt[2], iptIntExt[3], iptIntExt[4] );
  //printf("match :  deb= %d,  fin= %d \n", cell_match2_deb, cell_match2_fin);
  //printf("matchD:  deb= %d,  fin= %d \n", cell_match2_debD, cell_match2_finD);
  // map
  std::unordered_map<E_Int, E_Int> map;
  E_Int ind, indopp, e1, e2, eopp1, eopp2, e1g, e1oppg;
  E_Int indf, ne1, ne2, e2g, indfd, ne1d, ne2d, e2oppg;
  for (E_Int i = 0; i < FLI->getSize(); i++)
  {
    ind    = FLp[i]-1;        //FL starts from 1 
    indopp = FLDonorp[i]-1;

    e1    = PEp1[ind];
    e2    = PEp2[ind];
    eopp1 = PEDonorp1[indopp];
    eopp2 = PEDonorp2[indopp];
    
    //printf("face = %d, eltG = %d, eltD = %d \n", ind, e1,e2);
    //printf("faceD= %d, eltGD= %d, eltDD= %d \n", indopp, eopp1,eopp2);

    e1g = -1; e1oppg = -1; //indice element ghost 1ere couche
    if (e1 > 0 && e1 > cell_layer0   ) e1g = e1;
    if (e2 > 0 && e2 > cell_layer0   ) e1g = e2;
    if (eopp1 > 0 && eopp1 <= cell_layer0D) e1oppg = eopp1;
    if (eopp2 > 0 && eopp2 <= cell_layer0D) e1oppg = eopp2;
    if (e1g > 0){ map[e1g] = e1oppg; }//both start from 1
    //if (e1g > 0){ map[e1g] = e1oppg; printf("e1g= %d, e1opp= %d \n", e1g, e1oppg); }//both start from 1
    //if (e1g > -1){ map[e1g] = e1oppg; printf("e1g= %d, e1opp= %d \n", e1g, e1oppg); }//both start from 1

    E_Int* pt  = &(nface[indPH[e1g-1]]);
    E_Int* ptd = &(nfaceDonor[indPHDonor[e1oppg-1]]);

    E_Int nf = pt[0];
    for (E_Int j = 0; j < nf; j++) 
    {
      indf = pt[j+1]-1;
      ne1  = PEp1[indf];
      ne2  = PEp2[indf];
      if (ne1 == e1g) e2g = ne2;
      else            e2g = ne1;

      indfd = ptd[j+1]-1;
      ne1d  = PEDonorp1[indfd];
      ne2d  = PEDonorp2[indfd];
      if (ne1d == e1oppg) e2oppg = ne2d;
      else                e2oppg = ne1d;

      if (e2g >= cell_match2_deb && e2g <= cell_match2_fin && e2oppg <= cell_layer0D){ map[e2g] = e2oppg;}
      //if (e2g >= cell_match2_deb && e2g <= cell_match2_fin && e2oppg <= cell_layer0D){ map[e2g] = e2oppg;  printf("e2g= %d, e2opp= %d \n", e2g, e2oppg);}
 

      //if (e2g > cell_layer0 && e2g < cell_match && e2oppg < cell_layer0D){ map[e2g] = e2oppg;  printf("e2g= %d, e2opp= %d \n", e2g, e2oppg);}
      //if ( e2oppg <= cell_layer0D){ map[e2g] = e2oppg;  printf("e2g= %d, e2opp= %d \n", e2g, e2oppg);}
    }

  }
    //printf("FIN RAC \n");

  // Libere la memoire
  RELEASESHAREDN(FL, FLI);
  RELEASESHAREDN(FLDonor, FLDonorI);
  RELEASESHAREDN(rindElt, Intext);
  RELEASESHAREDN(rindEltDonor, IntextD);
  RELEASESHAREDN(PE, PEI);
  RELEASESHAREDN(PEDonor, PEDonorI);
  RELEASESHAREDB(res1, array, f, c);
  RELEASESHAREDB(res2, arrayDonor, fd, cd);

  // Remet la map a plat
  E_Int size = map.size();

  PyObject* PL = K_NUMPY::buildNumpyArray(size,1,1,1);
  PyObject* PLD = K_NUMPY::buildNumpyArray(size,1,1,1);
  PyObject* itype = K_NUMPY::buildNumpyArray(size,1,1,1);
  PyObject* coeff = K_NUMPY::buildNumpyArray(size,1,0,1);
  E_Int* PLp = K_NUMPY::getNumpyPtrI(PL);
  E_Int* PLDp = K_NUMPY::getNumpyPtrI(PLD);
  E_Int* itypep = K_NUMPY::getNumpyPtrI(itype);
  E_Float* coeffp = K_NUMPY::getNumpyPtrF(coeff);


  E_Int cpt = 0;
  for (std::pair<E_Int,E_Int> elt : map)
  {
    E_Int k = elt.first;
    E_Int ind = elt.second;
    PLp[cpt] = k-1;
    PLDp[cpt] = ind-1;
    cpt++;
  }

  for (E_Int i = 0; i < size; i++) 
  {
    itypep[i] = 1;
    coeffp[i] = 1.;
  }

  PyObject* tpl = Py_BuildValue("[OOOO]", PL, PLD, itype, coeff);
  Py_DECREF(PL); Py_DECREF(PLD); Py_DECREF(itype); Py_DECREF(coeff);
  return tpl;
}
