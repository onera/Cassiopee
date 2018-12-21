/*    
    Copyright 2013-2019 Onera.

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

// convertit un maillage prismatique (penta) en maillage tetraedrique
// voir article Dompierre et al 1999
// How to subdivide Pyramids, Prisms, and Hexahedra into Tetrahedra

# include "converter.h"
# include "kcore.h"
# include <string.h>
# include <stdio.h>

using namespace K_FLD;
using namespace std;
using namespace K_FUNC;

//=============================================================================
/* Conversion du maillage prismatique en maillage tetraedrique. Chaque prisme
   est decompose en 3 tetraedres. Pour chaque prisme, on determine l'indice
   minimum parmi ts les sommets. On determine ensuite les 2 diagonales passant
   par ce point et les deux facettes quad adjacentes. Sur la 3eme facette,
   on determine aussi le pt d'indice global minimum, ce qui permet de 
   determiner la derniere diagonale. De proche en proche, on construit 
   ainsi directement tous les tetraedres. */
//=============================================================================
PyObject* K_CONVERTER::convertPenta2Tetra(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;
  
  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, ni, nj, nk, cn, 
                                    eltType, true);

  // Test non structure ?
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "convertPenta2Tetra: invalid array.");
    return NULL;
  }
  if (res == 1)
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError,
                    "convertPenta2Tetra: input array must be unstructured.");
    return NULL;
  }

  // Test penta type ?
  if (strcmp(eltType, "PENTA") != 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "convertPenta2Tetra: unstructured array must be PENTA.");
    return NULL;
  }
  
  // Chaque prisme se decompose en 3 tetraedres
  E_Int neltsp = cn->getSize();
  E_Int nelts = 3*neltsp;
  FldArrayI& cn0 = *cn;
  E_Int eltt = 4; //TETRA 
  PyObject* tpl = K_ARRAY::buildArray(f->getNfld(), varString, f->getSize(), nelts, eltt, NULL);
  E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF fn(f->getSize(), f->getNfld(), fnp, true); fn = *f;
  E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
  FldArrayI ct(nelts, 4, cnnp, true);

  E_Int* ct1 = ct.begin(1);
  E_Int* ct2 = ct.begin(2);
  E_Int* ct3 = ct.begin(3);
  E_Int* ct4 = ct.begin(4);

#pragma omp parallel default(shared)
  {
  E_Int cnt;
  E_Int indir[6];
  E_Int diag;
  E_Int i1,i2,i3,i4,i5,i6;
  E_Int *cn01, *cn02, *cn03, *cn04, *cn05, *cn06;

#pragma omp for
  for (E_Int elt = 0; elt < neltsp; elt++)
  {
    /* determination du prisme tourne: imin -> pt 1
       determination du second point min sur la facette quad opposee: diag
       retour du tableau d'indirection I1I2I3I4I5I6 */
    diag = 0;
    buildSortedPrism(elt, *cn, diag, indir);
   
    i1 = indir[0]; cn01 = cn0.begin(i1);
    i2 = indir[1]; cn02 = cn0.begin(i2);
    i3 = indir[2]; cn03 = cn0.begin(i3);
    i4 = indir[3]; cn04 = cn0.begin(i4);
    i5 = indir[4]; cn05 = cn0.begin(i5);
    i6 = indir[5]; cn06 = cn0.begin(i6);
    if (diag == -1) //config 2-6
    {
      // build tetras: I1I2I3I6,I1I2I6I5,I1I5I6I4 
      // t1: I1I2I3I6
      cnt = 3*elt;
      ct1[cnt] = cn01[elt];
      ct2[cnt] = cn02[elt];
      ct3[cnt] = cn03[elt];
      ct4[cnt] = cn06[elt];

      // t2: I1I2I6I5
      cnt = 3*elt+1;
      ct1[cnt] = cn01[elt];
      ct2[cnt] = cn02[elt];
      ct3[cnt] = cn06[elt];
      ct4[cnt] = cn05[elt];

      // t3: I1I5I6I4
      cnt = 3*elt+2;
      ct1[cnt] = cn01[elt];
      ct2[cnt] = cn05[elt];
      ct3[cnt] = cn06[elt];
      ct4[cnt] = cn04[elt];
    }
    else if (diag == 1)//config 3-5
    {
      // build tetras: I1I2I3I5, I1I5I3I6, I1I5I6I4
      // t1: I1I2I3I5
      cnt = 3*elt;
      ct1[cnt] = cn01[elt];
      ct2[cnt] = cn02[elt];
      ct3[cnt] = cn03[elt];
      ct4[cnt] = cn05[elt];

      // t2: I1I5I3I6
      cnt = 3*elt+1;
      ct1[cnt] = cn01[elt];
      ct2[cnt] = cn05[elt];
      ct3[cnt] = cn03[elt];
      ct4[cnt] = cn06[elt];

      // t3: I1I5I6I4
      cnt = 3*elt+2;
      ct1[cnt] = cn01[elt];
      ct2[cnt] = cn05[elt];
      ct3[cnt] = cn06[elt];
      ct4[cnt] = cn04[elt];
    }
    //else //erreur de codage
    //{
    //  printf("Error: convertPenta2Tetra: bad value for diag.\n");
    //  RELEASESHAREDU(array, f, cn);
    //  return NULL;
    //}
  }
  }
  
  RELEASESHAREDU(array, f, cn);
  return tpl;
}

//=============================================================================
/* Determination du prisme ayant le sommet de plus petit indice en bas a gauche
   diag vaut -1 si deuxieme min est I2 ou I6
   diag vaut 1 si deuxieme min est I3 ou I5
 */
//=============================================================================
void K_CONVERTER::buildSortedPrism(E_Int elt, FldArrayI& cn, E_Int& diag,
                                   E_Int* indir)
{
  // Determination de indmin
  E_Int indmin = cn(elt,1);
  E_Int imin = 1;
  E_Int ind;
  for (E_Int i = 2; i <= 6; i++)
  {
    ind = cn(elt,i);
    if (ind < indmin)
    {
      indmin = ind;
      imin = i;
    }
  }
  switch (imin)
  {
    case 1 :
      indir[0] = 1;
      indir[1] = 2;
      indir[2] = 3;
      indir[3] = 4;
      indir[4] = 5;
      indir[5] = 6;
      break;

    case 2 : 
      indir[0] = 2;
      indir[1] = 3;
      indir[2] = 1;
      indir[3] = 5;
      indir[4] = 6;
      indir[5] = 4;
      break;

    case 3 : 
      indir[0] = 3;
      indir[1] = 1;
      indir[2] = 2;
      indir[3] = 6;
      indir[4] = 4;
      indir[5] = 5;
      break;

    case 4 : 
      indir[0] = 4;
      indir[1] = 6;
      indir[2] = 5;
      indir[3] = 1;
      indir[4] = 3;
      indir[5] = 2;
      break;

    case 5 : 
      indir[0] = 5;
      indir[1] = 4;
      indir[2] = 6;
      indir[3] = 2;
      indir[4] = 1;
      indir[5] = 3;
      break;

    case 6 : 
      indir[0] = 6;
      indir[1] = 5;
      indir[2] = 4;
      indir[3] = 3;
      indir[4] = 2;
      indir[5] = 1;
      break;

    default ://erreur de codage
      printf("Error: code error in function buildSortedPrism.\n");
      diag = 0;
      return;
  }
  //determination de l indice min sur la 3eme facette quad
  // soit I2, I6, I3, I5 
  E_Int indI2 = cn(elt,indir[1]);
  E_Int indI3 = cn(elt,indir[2]);
  E_Int indI5 = cn(elt,indir[4]);
  E_Int indI6 = cn(elt,indir[5]);
  
  indmin = indI2;
  diag = -1;

  if (indI6 < indmin)
  {
    indmin = indI6;
    diag = -1;
  } 
  if (indI3 < indmin) 
  {
    indmin = indI3;
    diag = 1;
  }
  if (indI5 < indmin)
  {
    indmin = indI5;
    diag = 1;
  }
}
