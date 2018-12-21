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
/* Flip les edges d'un maillage TRI */
//=============================================================================
PyObject* K_TRANSFORM::flipEdges(PyObject* self, PyObject* args)
{
  PyObject* o; E_Int mode; E_Int nit;
  if (!PYPARSETUPLEI(args, "Oll", "Oii", &o, &mode, &nit)) return NULL;
   
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
                    "flipEdges: input array must be unstructured.");
    return NULL;
  }

  // Test TRI
  if (strcmp(eltType, "TRI") != 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "flipEdges: unstructured array must be TRI.");
    RELEASESHAREDU(o, f, cn); return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "flipEdges: coord must be present in array.");
    RELEASESHAREDU(o, f, cn); return NULL;
  }
  posx++; posy++; posz++;

  E_Int posi = K_ARRAY::isNamePresent("indic", varString); posi++;
  E_Float* indic = NULL;
  if (posi > 0) indic = f->begin(posi);

  E_Float* x = f->begin(posx);
  E_Float* y = f->begin(posy);
  E_Float* z = f->begin(posz);
  FldArrayI ct = *cn;

  for (E_Int i = 0; i < nit; i++)
  {
    printf("iteration %d=mode=%d================\n",i,mode);
    flipEdges(ct, f->getSize(), x,y,z, indic, mode);
  }

  PyObject* tpl = K_ARRAY::buildArray(*f, varString, ct, 2, NULL);
  RELEASESHAREDU(o, f, cn);
  return tpl;
}

//=======================================================================
// IN: ct: connectivite TRI
// IN: np: nbre de noeuds dans le maillage TRI
// IN: x,y,z: ptrs sur les coords du maillage
// IN: mode: 1 flip mailles ecrasees avec check inverse
//           2 flip mailles ecrasees sans check inverse
//           3 flip toutes mailles avec optimisation du radius (qualite)
//=======================================================================
void K_TRANSFORM::flipEdges(FldArrayI& ct, E_Int np, 
                            E_Float* x, E_Float* y, E_Float* z, 
                            E_Float* indic, E_Int mode)
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
  E_Float ndir1, ndir2, ndir3, ndir4;
  E_Float ptA[3], ptB[3], ptC[3], dir1[3];
  E_Float ptD[3], dir2[3], dir3[3], dir4[3];
  E_Float inverse1, inverse2, rad1, rad2, rad3, rad4, ndirl;
  E_Int indA, indB, indC, indD, ind5, ind6, swap, ie, iv1, iv2, iv, pos1, pos2;
  E_Int tA, tB, tC, tD;
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
    if (ndirl < 1.e-11) printf("flipEdges: %d: %f init ecrase.\n", i, ndirl); 

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

      // BD ^ BA
      dir3[0] = (ptD[1]-ptB[1])*(ptA[2]-ptB[2])-(ptD[2]-ptB[2])*(ptA[1]-ptB[1]);
      dir3[1] = (ptD[2]-ptB[2])*(ptA[0]-ptB[0])-(ptD[0]-ptB[0])*(ptA[2]-ptB[2]);
      dir3[2] = (ptD[0]-ptB[0])*(ptA[1]-ptB[1])-(ptD[1]-ptB[1])*(ptA[0]-ptB[0]);
      ndir3 = sqrt(dir3[0]*dir3[0]+dir3[1]*dir3[1]+dir3[2]*dir3[2]);
      rad3 = K_COMPGEOM::circumCircleRadius(ptA, ptB, ptD);

      // CA ^ CD
      dir4[0] = (ptA[1]-ptC[1])*(ptD[2]-ptC[2])-(ptA[2]-ptC[2])*(ptD[1]-ptC[1]);
      dir4[1] = (ptA[2]-ptC[2])*(ptD[0]-ptC[0])-(ptA[0]-ptC[0])*(ptD[2]-ptC[2]);
      dir4[2] = (ptA[0]-ptC[0])*(ptD[1]-ptC[1])-(ptA[1]-ptC[1])*(ptD[0]-ptC[0]);
      ndir4 = sqrt(dir4[0]*dir4[0]+dir4[1]*dir4[1]+dir4[2]*dir4[2]);
      inverse2 = dir3[0]*dir4[0]+dir3[1]*dir4[1]+dir3[2]*dir4[2];
      if (ndir3 > 1.e-12 && ndir4 > 1.e-12) inverse2 = inverse2/(ndir3*ndir4);
      rad4 = K_COMPGEOM::circumCircleRadius(ptA, ptC, ptD);

      if (indic != NULL)
      { 
        tA = floor(indic[indA]+0.5); tB = floor(indic[indB]+0.5); 
        tC = floor(indic[indC]+0.5); tD = floor(indic[indD]+0.5);
        //printf("%d %d %d %d\n", tA, tB, tC, tD);
      }
      else {tA = 0; tB = 0; tC = 0; tD = 0;}
      if (tB == 1 && tC == 1 && ndir1 > 1.e-12) // no swap
      {
        goto next;
      }

      if (ndir1 < 1.e-11) printf("flipEdges: %d: %f ecrase. inv1=%f,inv2=%f,ndirs=%f %f %f\n", i, ndir1,inverse1,inverse2,ndir2,ndir3,ndir4); 

      if (mode == 2)
      {
        if (ndir1 < 1.e-12 && ndir3 > 1.e-12 && ndir4 > 1.e-12)
        {
          printf("flipEdges: %d: swap maille ecrasee1\n", i);
          printf("flipEdges: %d: improving %f %f -> %f %f\n", i, ndir1, ndir2, ndir3, ndir4);
          swap = ie; goto next;
        }
      }

      // swappable?
      if (mode == 1 && inverse2 > -0.9)
      {
        if (inverse1 < -0.9)
        {
          // corrige l'inverse
          printf("flipEdges: %d: swap inverse\n", i);
          printf("flipEdges: %d: improving %f %f -> %f %f\n", i, ndir1, ndir2, ndir3, ndir4);
          swap = ie; goto next;
        }
        else if (ndir1 < 1.e-12 && ndir3 > 1.e-12 && ndir4 > 1.e-12)
        {
          // corrige mailles ecrasees 
          printf("flipEdges: %d: swap maille ecrasee1\n", i);
          printf("flipEdges: %d: improving %f %f -> %f %f\n", i, ndir1, ndir2, ndir3, ndir4);
          swap = ie; goto next;
        }
      }

      if (mode == 3 && inverse2 > -0.9)
      {
        /*
        if (ndir1 < ndir3 && ndir1 < ndir4 && ndir2 < ndir3 && ndir2 < ndir4) // swap
        {
          // Ameliore les deux surfaces
          printf("flipEdges: swap sur\n");
          swap = ie; goto next;
        }
        else if (K_FUNC::E_abs(ndir3-ndir4) < K_FUNC::E_abs(ndir2-ndir1)) // swap
        {
          // Ameliore l'ecart entre les 2 surfaces
          printf("flipEdges: swap opt surfaces\n");
          swap = ie; goto next;
        }*/
        if (rad3 < rad1 && rad3 < rad2 && rad4 < rad1 && rad4 < rad2)
        {
          // optimisation suivant equilaterite
          printf("flipEdges: %d: swap opt equi\n",i );
          printf("flipEdges: %d: improving rad %f %f -> %f %f\n", i, rad1, rad2, rad3, rad4);
          swap = ie; goto next;
        }
      }
    }
    
    next:;
    // DBX
    //swap = -1;
    //if (swap != -1) 
    //{ printf("I can swap %d %d\n", i, swap); printf("improving %f %f -> %f %f\n", ndir1, ndir2, ndir3, ndir4);}
    //if (i != 0 && i != 1 && i != 2 && i != 3) swap = -1;
    //if (i != 3) swap = -1;
    
    if (swap < 0)
    {
      if (inverse1 < -0.9) maillesInversees++;
      if (ndir1 < 1.e-12) maillesEcrasees++;
    }
    else
    {
      if (inverse2 < -0.9) maillesInversees++;
      if (ndir3 < 1.e-12) maillesEcrasees++; 
    }

    // update cEEN
    if (swap != -1)
    {
      ie = swap;
      // voisins de i
      vector< E_Int >& voisins = cEEN[i];
      iv1 = -1;
      for (size_t v = 0; v < voisins.size(); v++)
      {
        iv = voisins[v];
        ind1 = ct1[iv]-1; ind2 = ct2[iv]-1; ind3 = ct3[iv]-1;
        //printf("vois i %d %d %d\n", ind1,ind2,ind3);
        if      (ind1 == indA && ind2 == indC) { iv1 = iv; pos1 = v; break; }
        else if (ind1 == indA && ind3 == indC) { iv1 = iv; pos1 = v; break; }
        else if (ind2 == indA && ind3 == indC) { iv1 = iv; pos1 = v; break; }
        else if (ind2 == indA && ind1 == indC) { iv1 = iv; pos1 = v; break; }
        else if (ind3 == indA && ind1 == indC) { iv1 = iv; pos1 = v; break; }
        else if (ind3 == indA && ind2 == indC) { iv1 = iv; pos1 = v; break; }
      }
      // voisins de ie
      vector< E_Int >& voisins2 = cEEN[ie];
      iv2 = -1;
      for (size_t v = 0; v < voisins2.size(); v++)
      {
        iv = voisins2[v];
        ind1 = ct1[iv]-1; ind2 = ct2[iv]-1; ind3 = ct3[iv]-1;
        //printf("vois ie %d %d %d\n", ind1,ind2,ind3);
      
        if      (ind1 == indB && ind2 == indD) { iv2 = iv; pos2 = v; break; }
        else if (ind1 == indB && ind3 == indD) { iv2 = iv; pos2 = v; break; }
        else if (ind2 == indB && ind3 == indD) { iv2 = iv; pos2 = v; break; }
        else if (ind2 == indB && ind1 == indD) { iv2 = iv; pos2 = v; break; }
        else if (ind3 == indB && ind1 == indD) { iv2 = iv; pos2 = v; break; }
        else if (ind3 == indB && ind2 == indD) { iv2 = iv; pos2 = v; break; }
      }
      // swap
      if (iv1 != -1 && iv2 != -1)
      {
        //printf("ech: %d %d - %d %d %d %d\n",i,ie,iv1,iv2,pos1,pos2);
        cEEN[i][pos1] = iv2;
        cEEN[ie][pos2] = iv1;
      }
      else if (iv1 != -1)
      {
        //printf("erase: %d %d - %d %d %d\n",i,ie,iv1,pos1,cEEN[i].size());
        cEEN[i].erase(cEEN[i].begin()+pos1);
        cEEN[ie].push_back(iv1);
      }
      else if (iv2 != -1)
      {
        //printf("erase: %d %d - %d %d %d\n",i,ie,iv2,pos2,cEEN[ie].size());
        cEEN[i].push_back(iv2);
        cEEN[ie].erase(cEEN[ie].begin()+pos2);
      }
   
      if (iv1 != -1)
      {
        vector< E_Int >& voisins = cEEN[iv1];
        for (size_t v = 0; v < voisins.size(); v++)
        {
           if (voisins[v] == i) { cEEN[iv1][v] = ie; break; }
        }
      }
      if (iv2 != -1)
      {
        vector< E_Int >& voisins = cEEN[iv2];
        for (size_t v = 0; v < voisins.size(); v++)
        {
           if (voisins[v] == ie) { cEEN[iv2][v] = i; break; }
        }
      }

      ct1[i] = indA+1; ct2[i] = indB+1; ct3[i] = indD+1;
      ct1[ie] = indA+1; ct2[ie] = indD+1; ct3[ie] = indC+1;
    }
  }

  printf("Mailles inversees=%d - mailles ecrasees=%d\n", maillesInversees, maillesEcrasees);

}
