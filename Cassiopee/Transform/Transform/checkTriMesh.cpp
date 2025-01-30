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
// ct: connectivite TRI
# include "transform.h" 
# include "kcore.h"

using namespace K_FLD;
using namespace std;
using namespace K_FUNC;

//=============================================================================
/* Verifie un maillage TRI */
//=============================================================================
PyObject* K_TRANSFORM::checkTriMesh(PyObject* self, PyObject* args)
{
  PyObject* o; E_Int mode;
  if (!PYPARSETUPLE_(args, O_ I_, &o, &mode)) return NULL;
   
  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(o, varString,
                                    f, ni, nj, nk, cn, eltType, true);
  // Test non structure ?
  if (res != 2)
  {
    if (res == 1) RELEASESHAREDS(o,f);
    PyErr_SetString(PyExc_TypeError,
                    "checkTriMesh: input array must be unstructured.");
    return NULL;
  }

  // Test TRI
  if (strcmp(eltType, "TRI") != 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "checkTriMesh: unstructured array must be TRI.");
    RELEASESHAREDU(o, f, cn); return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "checkTriMesh: coord must be present in array.");
    RELEASESHAREDU(o, f, cn); return NULL;
  }
  posx++; posy++; posz++;

  E_Float* x = f->begin(posx);
  E_Float* y = f->begin(posy);
  E_Float* z = f->begin(posz);
  FldArrayI ct = *cn;

  E_Int nec, ninv;
  checkTriMesh(ct, f->getSize(), x,y,z, nec, ninv);

  RELEASESHAREDU(o, f, cn);
  return Py_BuildValue("(l,l)", nec, ninv);
}

//================================================================================
// check mesh
// verifie s'il y a des mailles ecrasees ou inversees dans un maillage TRI
//=================================================================================
void K_TRANSFORM::checkTriMesh(FldArrayI& ct, E_Int np, 
                               E_Float* x, E_Float* y, E_Float* z,
                               E_Int& ne, E_Int& ni)
{
  E_Int ntr = ct.getSize();
  E_Int* ct1 = ct.begin(1);
  E_Int* ct2 = ct.begin(2);
  E_Int* ct3 = ct.begin(3);
  E_Int ind1, ind2, ind3, ind4;
  // Connectivite elts-elts voisins
  vector< vector<E_Int> > cEEN(ntr);
  K_CONNECT::connectEV2EENbrs("TRI", np, ct, cEEN);
  ct1 = ct.begin(1); ct2 = ct.begin(2); ct3 = ct.begin(3);
  E_Float ndir1, ndir2;
  E_Float ptA[3], ptB[3], ptC[3], dir1[3];
  E_Float ptD[3], dir2[3];
  E_Float inverse1, rad1, rad2, ndirl;
  E_Int indA, indB, indC, indD, ind5, ind6, swap;

  E_Int maillesEcrasees = 0;
  E_Int maillesInversees = 0;  

  for (E_Int i = 0; i < ntr; i++)
  {
    ind1 = ct1[i]-1; ind2 = ct2[i]-1; ind3 = ct3[i]-1;
    vector<E_Int>& voisins = cEEN[i];
    swap = -1;

    /* verification du triangle initial */
    ptA[0] = x[ind1]; ptA[1] = y[ind1]; ptA[2] = z[ind1];
    ptB[0] = x[ind2]; ptB[1] = y[ind2]; ptB[2] = z[ind2];
    ptC[0] = x[ind3]; ptC[1] = y[ind3]; ptC[2] = z[ind3];
    dir1[0] = (ptB[1]-ptA[1])*(ptC[2]-ptA[2])-(ptB[2]-ptA[2])*(ptC[1]-ptA[1]);
    dir1[1] = (ptB[2]-ptA[2])*(ptC[0]-ptA[0])-(ptB[0]-ptA[0])*(ptC[2]-ptA[2]);
    dir1[2] = (ptB[0]-ptA[0])*(ptC[1]-ptA[1])-(ptB[1]-ptA[1])*(ptC[0]-ptA[0]);
    ndirl = sqrt(dir1[0]*dir1[0]+dir1[1]*dir1[1]+dir1[2]*dir1[2]);
    if (ndirl < 1.e-11) 
    {
      printf("check: " SF_D_ ": " SF_F_ " maille ecrase.\n", i, ndirl); 
      maillesEcrasees += 1;
    }

    for (size_t v = 0; v < voisins.size(); v++)
    {
      E_Int ie = voisins[v];
      ind4 = ct1[ie]-1; ind5 = ct2[ie]-1; ind6 = ct3[ie]-1;
      if      (ind1 == ind4 && ind2 == ind6) {indA = ind3; indB = ind1; indC = ind2; indD = ind5;}
      else if (ind1 == ind4 && ind2 == ind5) {indA = ind3; indB = ind1; indC = ind2; indD = ind6;}
      else if (ind1 == ind5 && ind2 == ind4) {indA = ind3; indB = ind1; indC = ind2; indD = ind6;}
      else if (ind1 == ind5 && ind2 == ind6) {indA = ind3; indB = ind1; indC = ind2; indD = ind4;}
      else if (ind1 == ind6 && ind2 == ind4) {indA = ind3; indB = ind1; indC = ind2; indD = ind5;}
      else if (ind1 == ind6 && ind2 == ind5) {indA = ind3; indB = ind1; indC = ind2; indD = ind4;}

      else if (ind1 == ind4 && ind3 == ind6) {indA = ind2; indB = ind3; indC = ind1; indD = ind5;}
      else if (ind1 == ind4 && ind3 == ind5) {indA = ind2; indB = ind3; indC = ind1; indD = ind6;}
      else if (ind1 == ind5 && ind3 == ind4) {indA = ind2; indB = ind3; indC = ind1; indD = ind6;}
      else if (ind1 == ind5 && ind3 == ind6) {indA = ind2; indB = ind3; indC = ind1; indD = ind4;}
      else if (ind1 == ind6 && ind3 == ind4) {indA = ind2; indB = ind3; indC = ind1; indD = ind5;}
      else if (ind1 == ind6 && ind3 == ind5) {indA = ind2; indB = ind3; indC = ind1; indD = ind4;}

      else if (ind2 == ind4 && ind3 == ind6) {indA = ind1; indB = ind2; indC = ind3; indD = ind5;}
      else if (ind2 == ind4 && ind3 == ind5) {indA = ind1; indB = ind2; indC = ind3; indD = ind6;}
      else if (ind2 == ind5 && ind3 == ind4) {indA = ind1; indB = ind2; indC = ind3; indD = ind6;}
      else if (ind2 == ind5 && ind3 == ind6) {indA = ind1; indB = ind2; indC = ind3; indD = ind4;}
      else if (ind2 == ind6 && ind3 == ind4) {indA = ind1; indB = ind2; indC = ind3; indD = ind5;}
      else if (ind2 == ind6 && ind3 == ind5) {indA = ind1; indB = ind2; indC = ind3; indD = ind4;}
      else 
      {
        indA = 0; indB = 0; indC = 0; indD = 0;
        printf("what?? problem\n");
      }

      ptA[0] = x[indA]; ptA[1] = y[indA]; ptA[2] = z[indA];
      ptB[0] = x[indB]; ptB[1] = y[indB]; ptB[2] = z[indB];
      ptC[0] = x[indC]; ptC[1] = y[indC]; ptC[2] = z[indC];
      ptD[0] = x[indD]; ptD[1] = y[indD]; ptD[2] = z[indD];

      // AB ^ AC  
      dir1[0] = (ptB[1]-ptA[1])*(ptC[2]-ptA[2])-(ptB[2]-ptA[2])*(ptC[1]-ptA[1]);
      dir1[1] = (ptB[2]-ptA[2])*(ptC[0]-ptA[0])-(ptB[0]-ptA[0])*(ptC[2]-ptA[2]);
      dir1[2] = (ptB[0]-ptA[0])*(ptC[1]-ptA[1])-(ptB[1]-ptA[1])*(ptC[0]-ptA[0]);
      ndir1 = sqrt(dir1[0]*dir1[0]+dir1[1]*dir1[1]+dir1[2]*dir1[2]);
      rad1 = K_COMPGEOM::circumCircleRadius(ptA[0], ptA[1], ptA[2],
                                            ptB[0], ptB[1], ptB[2],
                                            ptC[0], ptC[1], ptC[2]);

      // DC ^ DB
      dir2[0] = (ptC[1]-ptD[1])*(ptB[2]-ptD[2])-(ptC[2]-ptD[2])*(ptB[1]-ptD[1]);
      dir2[1] = (ptC[2]-ptD[2])*(ptB[0]-ptD[0])-(ptC[0]-ptD[0])*(ptB[2]-ptD[2]);
      dir2[2] = (ptC[0]-ptD[0])*(ptB[1]-ptD[1])-(ptC[1]-ptD[1])*(ptB[0]-ptD[0]);
      ndir2 = sqrt(dir2[0]*dir2[0]+dir2[1]*dir2[1]+dir2[2]*dir2[2]);
      inverse1 = dir1[0]*dir2[0]+dir1[1]*dir2[1]+dir1[2]*dir2[2];
      if (ndir1 > 1.e-12 && ndir2 > 1.e-12) inverse1 = inverse1/(ndir1*ndir2);
      rad2 = K_COMPGEOM::circumCircleRadius(ptB[0], ptB[1], ptB[2],
                                            ptC[0], ptC[1], ptC[2],
                                            ptD[0], ptD[1], ptD[2]);

      if (inverse1 < -0.9)
      {
        printf("check: " SF_D_ ": " SF_F_ " maille inversee.\n", i, inverse1); 
        maillesInversees += 1;
      }
    }

  }

  printf("Check: Mailles inversees=" SF_D_ " - mailles ecrasees=" SF_D_ "\n", maillesInversees, maillesEcrasees);
  ne = maillesEcrasees;
  ni = maillesInversees;
}
