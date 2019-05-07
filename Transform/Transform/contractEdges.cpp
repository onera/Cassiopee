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
// ct: connectivite TRI
# include "transform.h" 
# include "kcore.h"

using namespace K_FLD;
using namespace std;
using namespace K_FUNC;

//=============================================================================
/* Contract les edges d'un maillage TRI suivant differents criteres
  Les points sont modifies, le nbre d'elements est modifie */
//=============================================================================
PyObject* K_TRANSFORM::contractEdges(PyObject* self, PyObject* args)
{
  PyObject* o; E_Int mode;
  if (!PYPARSETUPLEI(args, "Ol", "Oi", &o, &mode)) return NULL;
   
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
                    "contractEdges: input array must be unstructured.");
    return NULL;
  }

  // Test TRI
  if (strcmp(eltType, "TRI") != 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "contractEdges: unstructured array must be TRI.");
    RELEASESHAREDU(o, f, cn); return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "contract: coord must be present in array.");
    RELEASESHAREDU(o, f, cn); return NULL;
  }
  posx++; posy++; posz++;

  E_Float* x = f->begin(posx);
  E_Float* y = f->begin(posy);
  E_Float* z = f->begin(posz);
  FldArrayI ct = *cn;

  contractEdges(ct, f->getSize(), x,y,z, mode);

  PyObject* tpl = K_ARRAY::buildArray(*f, varString, ct, 2, NULL);
  RELEASESHAREDU(o, f, cn);
  return tpl;
}

//=======================================================================
// IN: ct: connectivite TRI
// IN: np: nbre de noeuds dans le maillage TRI
// IN: x,y,z: ptrs sur les coords du maillage
// IN: mode: 1 contract les mailles ecrasees
//           2 contract les mailles de gros radius
//           3 contracte les mailles si inverse
//=======================================================================
void K_TRANSFORM::contractEdges(FldArrayI& ct, E_Int np, 
                                E_Float* x, E_Float* y, E_Float* z, 
                                E_Int mode)
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

  // swap edges dans les triangles
  E_Float ndir1, ndir2;
  E_Float ptA[3], ptB[3], ptC[3], dir1[3];
  E_Float ptD[3], dir2[3];
  E_Float inverse1, rad1, rad2, ndirl;
  E_Int indA, indB, indC, indD, ind5, ind6, swap;
  E_Float xx, yy, zz, e1, e2, e3, rad;
  E_Int neigh, mat, v;
  //E_Int maillesEcrasees = 0;
  //E_Int maillesInversees = 0;

  short* decim = new short [ntr];
  for (E_Int i = 0; i < ntr; i++) decim[i] = 1;
  short* fixed = new short [np];
  for (E_Int i = 0; i < np; i++) fixed[i] = 0;

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
    if (ndirl < 1.e-11) printf("contractEdges: %d: %f init ecrase.\n", i, ndirl); 

    if (mode == 1) // decimate si degenere (ecrase)
    {
      if (ndirl < 1.e-11 && decim[i] == 1)
      {
        // Contracte le plus petit edge
        e1 = sqrt((ptB[0]-ptA[0])*(ptB[0]-ptA[0])+
                  (ptB[1]-ptA[1])*(ptB[1]-ptA[1])+
                  (ptB[2]-ptA[2])*(ptB[2]-ptA[2]));
        e2 = sqrt((ptC[0]-ptB[0])*(ptC[0]-ptB[0])+
                  (ptC[1]-ptB[1])*(ptC[1]-ptB[1])+
                  (ptC[2]-ptB[2])*(ptC[2]-ptB[2]));
        e3 = sqrt((ptA[0]-ptC[0])*(ptA[0]-ptC[0])+
                  (ptA[1]-ptC[1])*(ptA[1]-ptC[1])+
                  (ptA[2]-ptC[2])*(ptA[2]-ptC[2]));

        neigh = -1;
        if (e1 <= e2 && e1 <= e3 && fixed[ind1] == 0 && fixed[ind2] == 0) // contracte e1
        {
          xx = 0.5*(x[ind1]+x[ind2]);
          yy = 0.5*(y[ind1]+y[ind2]);
          zz = 0.5*(z[ind1]+z[ind2]);
          x[ind1] = xx; x[ind2] = xx;
          y[ind1] = yy; y[ind2] = yy;
          z[ind1] = zz; z[ind2] = zz;
          fixed[ind1] = 1; fixed[ind2] = 1;
          for (size_t k = 0; k < voisins.size(); k++)
          {
            v = voisins[k]; mat = 0;
            if (ct1[v] == ind1+1 || ct1[v] == ind2+1) mat++;
            if (ct2[v] == ind1+1 || ct2[v] == ind2+1) mat++;
            if (ct3[v] == ind1+1 || ct3[v] == ind2+1) mat++;
            if (mat == 2) { neigh = v; break; }
          }
        }
        else if (e2 <= e1 && e2 <= e3 && fixed[ind2] == 0 && fixed[ind3] == 0)
        {
          xx = 0.5*(x[ind2]+x[ind3]);
          yy = 0.5*(y[ind2]+y[ind3]);
          zz = 0.5*(z[ind2]+z[ind3]);
          x[ind2] = xx; x[ind3] = xx;
          y[ind2] = yy; y[ind3] = yy;
          z[ind2] = zz; z[ind3] = zz;
          fixed[ind2] = 1; fixed[ind3] = 1;
          for (size_t k = 0; k < voisins.size(); k++)
          {
            v = voisins[k]; mat = 0;
            if (ct1[v] == ind2+1 || ct1[v] == ind3+1) mat++;
            if (ct2[v] == ind2+1 || ct2[v] == ind3+1) mat++;
            if (ct3[v] == ind2+1 || ct3[v] == ind3+1) mat++;
            if (mat == 2) { neigh = v; break; }
          }
        }
        else if (fixed[ind1] == 0 && fixed[ind3] == 0)
        {
          xx = 0.5*(x[ind1]+x[ind3]);
          yy = 0.5*(y[ind1]+y[ind3]);
          zz = 0.5*(z[ind1]+z[ind3]);
          x[ind1] = xx; x[ind3] = xx;
          y[ind1] = yy; y[ind3] = yy;
          z[ind1] = zz; z[ind3] = zz;
          fixed[ind1] = 1; fixed[ind3] = 1;
          for (size_t k = 0; k < voisins.size(); k++)
          {
            v = voisins[k]; mat = 0;
            if (ct1[v] == ind1+1 || ct1[v] == ind3+1) mat++;
            if (ct2[v] == ind1+1 || ct2[v] == ind3+1) mat++;
            if (ct3[v] == ind1+1 || ct3[v] == ind3+1) mat++;
            if (mat == 2) { neigh = v; break; }
          }
        }
        if (neigh > 0 && decim[neigh] == 1) { decim[i] = 0; decim[neigh] = 0; }
      }
    }
    else if (mode == 3) // decimate si radius trop grand
    {
      rad = K_COMPGEOM::circumCircleRadius(ptA, ptB, ptC);
      e1 = sqrt((ptB[0]-ptA[0])*(ptB[0]-ptA[0])+
                (ptB[1]-ptA[1])*(ptB[1]-ptA[1])+
                (ptB[2]-ptA[2])*(ptB[2]-ptA[2]));
      e2 = sqrt((ptC[0]-ptB[0])*(ptC[0]-ptB[0])+
                (ptC[1]-ptB[1])*(ptC[1]-ptB[1])+
                (ptC[2]-ptB[2])*(ptC[2]-ptB[2]));
      e3 = sqrt((ptA[0]-ptC[0])*(ptA[0]-ptC[0])+
                (ptA[1]-ptC[1])*(ptA[1]-ptC[1])+
                (ptA[2]-ptC[2])*(ptA[2]-ptC[2]));
      // max radius

    }
    else if (mode == 2) // decimate edge si inverse avec voisin
    {
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
        //printf("ind1 %d %d %d\n", ind1, ind2, ind3);
        //printf("ind4 %d %d %d\n", ind4, ind5, ind6);
        //printf("indA %d %d %d %d\n", indA, indB, indC, indD);

        ptA[0] = x[indA]; ptA[1] = y[indA]; ptA[2] = z[indA];
        ptB[0] = x[indB]; ptB[1] = y[indB]; ptB[2] = z[indB];
        ptC[0] = x[indC]; ptC[1] = y[indC]; ptC[2] = z[indC];
        ptD[0] = x[indD]; ptD[1] = y[indD]; ptD[2] = z[indD];

        // AB ^ AC  
        dir1[0] = (ptB[1]-ptA[1])*(ptC[2]-ptA[2])-(ptB[2]-ptA[2])*(ptC[1]-ptA[1]);
        dir1[1] = (ptB[2]-ptA[2])*(ptC[0]-ptA[0])-(ptB[0]-ptA[0])*(ptC[2]-ptA[2]);
        dir1[2] = (ptB[0]-ptA[0])*(ptC[1]-ptA[1])-(ptB[1]-ptA[1])*(ptC[0]-ptA[0]);
        ndir1 = sqrt(dir1[0]*dir1[0]+dir1[1]*dir1[1]+dir1[2]*dir1[2]);
        rad1 = K_COMPGEOM::circumCircleRadius(ptA, ptB, ptC);

        // DC ^ DB
        dir2[0] = (ptC[1]-ptD[1])*(ptB[2]-ptD[2])-(ptC[2]-ptD[2])*(ptB[1]-ptD[1]);
        dir2[1] = (ptC[2]-ptD[2])*(ptB[0]-ptD[0])-(ptC[0]-ptD[0])*(ptB[2]-ptD[2]);
        dir2[2] = (ptC[0]-ptD[0])*(ptB[1]-ptD[1])-(ptC[1]-ptD[1])*(ptB[0]-ptD[0]);
        ndir2 = sqrt(dir2[0]*dir2[0]+dir2[1]*dir2[1]+dir2[2]*dir2[2]);
        inverse1 = dir1[0]*dir2[0]+dir1[1]*dir2[1]+dir1[2]*dir2[2];
        if (ndir1 > 1.e-12 && ndir2 > 1.e-12) inverse1 = inverse1/(ndir1*ndir2);
        rad2 = K_COMPGEOM::circumCircleRadius(ptB, ptC, ptD);

        if (inverse1 < -0.9) printf("contractEdges: %d: %f inverse.\n", i, inverse1); 

        if (inverse1 < -0.9 && decim[i] == 1 && decim[ie] == 1 && fixed[indB] == 0 && fixed[indC] == 0) 
        {
          decim[i] = 0; decim[ie] = 0;
          fixed[indB] = 1; fixed[indC] = 1;
          xx = 0.5*(x[indB]+x[indC]);
          yy = 0.5*(y[indB]+y[indC]);
          zz = 0.5*(z[indB]+z[indC]);
          x[indB] = xx; x[indC] = xx;
          y[indB] = yy; y[indC] = yy;
          z[indB] = zz; z[indC] = zz;
          break;
        }
      }
    }
  }

  // Supprime les mailles decim
  E_Int ne = 0;
  for (E_Int i = 0; i < ntr; i++) ne += decim[i];

  printf("contractEdges: final number of elements=%d\n", ne);
  FldArrayI ctn(ne, 3);
  E_Int* ctn1 = ctn.begin(1);
  E_Int* ctn2 = ctn.begin(2);
  E_Int* ctn3 = ctn.begin(3);
  E_Int j = 0;
  for (E_Int i = 0; i < ntr; i++)
  {
    if (decim[i] == 1)
    {
      ctn1[j] = ct1[i]; ctn2[j] = ct2[i]; ctn3[j] = ct3[i]; j++; 
    }
  }

  ct = ctn;

  //printf("Mailles inversees=%d - mailles ecrasees=%d\n", maillesInversees, maillesEcrasees);
  delete [] decim;
  delete [] fixed;
}


