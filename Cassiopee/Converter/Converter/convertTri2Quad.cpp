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

// convertit un maillage triangulaire en maillage quad

# include "converter.h" 
# include "kcore.h"
# include <string.h>
# include <stdio.h>

using namespace K_FLD;
using namespace std;
using namespace K_FUNC;

E_Float getMaxAngle(E_Float* A, E_Float* B, E_Float* C, E_Float* D);

//=============================================================================
/* Conversion du maillage triangulaire en maillage quad.
   On rend un maillage quad et un maillage tri. Le maillage tri contient
   les elements du maillage de depart qui n'ont pas pu etre fusionnes. */
//=============================================================================
PyObject* K_CONVERTER::convertTri2Quad(PyObject* self, PyObject* args)
{
  PyObject* pArray;
  E_Float angle;
  if (!PYPARSETUPLE_(args, O_ R_, &pArray, &angle)) return NULL;
 
  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(pArray, varString,
                                     f, ni, nj, nk, cn, eltType);

  // Test non structure ?
  if (res != 2)
  {
    if (res == 1) RELEASESHAREDS(pArray,f);
    PyErr_SetString(PyExc_TypeError,
                    "convertTri2Quad: input array must be unstructured.");
    return NULL;
  }

  // Test TRI
  if (strcmp(eltType, "TRI") != 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "convertTri2Quad: unstructured array must be TRI.");
    RELEASESHAREDU(pArray, f, cn); return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "convertTri2Quad: coord must be present in array.");
    RELEASESHAREDU(pArray, f, cn); return NULL;
  }
  posx++; posy++; posz++;

  // Connectivite
  E_Int nv = f->getSize();
  E_Int ne = cn->getSize();
  vector< vector<E_Int> > cEEN(ne);
  K_CONNECT::connectEV2EENbrs(eltType, nv, *cn, cEEN);

  // Tableaux de tags
  FldArrayIS dejaVu(ne); dejaVu.setAllValuesAtNull();
  short* dejaVup = dejaVu.begin();
  FldArrayIS merged(ne); merged.setAllValuesAtNull();
  short* mergedp = merged.begin();

  // QUAD de sortie
  E_Int nq = 0;
  FldArrayI cq(ne, 4);
  E_Int* cq1 = cq.begin(1);
  E_Int* cq2 = cq.begin(2);
  E_Int* cq3 = cq.begin(3);
  E_Int* cq4 = cq.begin(4);

  E_Int elt, i, ind, eltv, eltvF, size;
  E_Float angleB=0, angleBF, angleA;
  E_Float ptA1[3], ptB1[3], ptC1[3];
  E_Float ptA2[3], ptB2[3], ptC2[3];
  E_Float* x = f->begin(posx);
  E_Float* y = f->begin(posy);
  E_Float* z = f->begin(posz);
  FldArrayI& cnp = *cn;
  E_Int* cnp1 = cnp.begin(1);
  E_Int* cnp2 = cnp.begin(2);
  E_Int* cnp3 = cnp.begin(3);

  debut: ;

  // Recherche d'un element de frontiere : elt
  for (i = 0; i < ne; i++)
  { if (dejaVup[i] == 0 && cEEN[i].size() < 3) break; }

  if (i < ne) elt = i; // elt frontiere trouve
  else
  { // recherche du premier elt non deja vu
    for (i = 0; i < ne; i++) 
    { if (dejaVup[i] == 0) break; }
    elt = i;
    if (i == ne) goto fin;
  }
  dejaVup[elt] = 1;

  // Cherche le meilleur voisin avec qui fusionner
  eltvF = -1; angleA = 0.; angleBF = 360.;
  size = cEEN[elt].size();
  for (i = 0; i < size; i++)
  {
    eltv = cEEN[elt][i];
    
    if (dejaVup[eltv] == 0)
    {
      ind = cnp1[elt]-1;
      ptA1[0] = x[ind]; ptA1[1] = y[ind]; ptA1[2] = z[ind];
      ind = cnp2[elt]-1;
      ptB1[0] = x[ind]; ptB1[1] = y[ind]; ptB1[2] = z[ind];
      ind = cnp3[elt]-1;
      ptC1[0] = x[ind]; ptC1[1] = y[ind]; ptC1[2] = z[ind];
      ind = cnp1[eltv]-1;
      ptA2[0] = x[ind]; ptA2[1] = y[ind]; ptA2[2] = z[ind];
      ind = cnp2[eltv]-1;
      ptB2[0] = x[ind]; ptB2[1] = y[ind]; ptB2[2] = z[ind];
      ind = cnp3[eltv]-1;
      ptC2[0] = x[ind]; ptC2[1] = y[ind]; ptC2[2] = z[ind];
      
      angleA = 
        K_COMPGEOM::getAlphaAngleBetweenTriangles(ptA1, ptB1, ptC1,
                                                  ptA2, ptB2, ptC2);

      if (E_abs(angleA - 180.) < angle)
      {
        if (cnp1[eltv] != cnp1[elt] &&
            cnp1[eltv] != cnp2[elt] && cnp1[eltv] != cnp3[elt] ) 
        {
          if (cnp2[eltv] == cnp3[elt])
          {
            angleB = getMaxAngle(ptA1, ptB1, ptA2, ptC1); // OK
          }
          else if (cnp2[eltv] == cnp1[elt])
          {
            angleB = getMaxAngle(ptA1, ptB1, ptC1, ptA2);
          }
          else
          {
            angleB = getMaxAngle(ptA1, ptA2, ptB1, ptC1);
          }
        }
        else if (cnp2[eltv] != cnp1[elt] && 
                 cnp2[eltv] != cnp2[elt] && cnp2[eltv] != cnp3[elt] ) 
        {
          if (cnp3[eltv] == cnp3[elt])
          {
            angleB = getMaxAngle(ptA1, ptB1, ptB2, ptC1);
          }
          else if (cnp3[eltv] == cnp1[elt])
          {
            angleB = getMaxAngle(ptA1, ptB1, ptC1, ptB2);
          }
          else
          {
            angleB = getMaxAngle(ptA1, ptB2, ptB1, ptC1);
          }
        }
        else if (cnp3[eltv] != cnp1[elt] && 
                 cnp3[eltv] != cnp2[elt] && cnp3[eltv] != cnp3[elt] )
        {
          if (cnp1[eltv] == cnp3[elt])
          { 
            angleB = getMaxAngle(ptA1, ptB1, ptC2, ptC1);
          }
          else if (cnp1[eltv] == cnp1[elt])
          {
            angleB = getMaxAngle(ptA1, ptB1, ptC1, ptC2);
          }
          else
          {
            angleB = getMaxAngle(ptA1, ptC2, ptB1, ptC1);
          }
        }
        if (angleB < angleBF) { eltvF = eltv; angleBF = angleB; }
      }
    }
  }
  
  // Tres mauvais quad. On ne le fusionne pas.
  //printf("Je regarde l'element %d\n", elt);
  //printf("J'ai trouve %d : %f\n", eltvF, angleBF);
  if (angleBF > 80.) eltvF = -1;
  //eltvF = -1; // DEBUG
  eltv = eltvF;

  // Fusionne avec eltv
  if (eltv != -1)
  { 
    mergedp[elt] = 1; mergedp[eltv] = 1; dejaVup[eltv] = 1;
        
    if (cnp1[eltv] != cnp1[elt] && 
        cnp1[eltv] != cnp2[elt] && cnp1[eltv] != cnp3[elt] ) 
    {
      if (cnp2[eltv] == cnp3[elt])
      {
        cq1[nq] = cnp1[elt];
        cq2[nq] = cnp2[elt];
        cq3[nq] = cnp1[eltv];
        cq4[nq] = cnp3[elt];
      }
      else if (cnp2[eltv] == cnp1[elt])
      {
        cq1[nq] = cnp1[elt];
        cq2[nq] = cnp2[elt];
        cq3[nq] = cnp3[elt];
        cq4[nq] = cnp1[eltv];
      }
      else
      {
        cq1[nq] = cnp1[elt];
        cq2[nq] = cnp1[eltv];
        cq3[nq] = cnp2[elt];
        cq4[nq] = cnp3[elt];
      }
      nq++;
    }
    else if (cnp2[eltv] != cnp1[elt] && 
             cnp2[eltv] != cnp2[elt] && cnp2[eltv] != cnp3[elt] ) 
    {
      if (cnp3[eltv] == cnp3[elt])
      {
        cq1[nq] = cnp1[elt];
        cq2[nq] = cnp2[elt];
        cq3[nq] = cnp2[eltv];
        cq4[nq] = cnp3[elt];
      }
      else if (cnp3[eltv] == cnp1[elt])
      {
        cq1[nq] = cnp1[elt];
        cq2[nq] = cnp2[elt];
        cq3[nq] = cnp3[elt];
        cq4[nq] = cnp2[eltv];
      }
      else
      {
        cq1[nq] = cnp1[elt];
        cq2[nq] = cnp2[eltv];
        cq3[nq] = cnp2[elt];
        cq4[nq] = cnp3[elt];
      }
      nq++;
    }
    else if (cnp3[eltv] != cnp1[elt] && 
             cnp3[eltv] != cnp2[elt] && cnp3[eltv] != cnp3[elt] )
    {
      if (cnp1[eltv] == cnp3[elt])
      {
        cq1[nq] = cnp1[elt];
        cq2[nq] = cnp2[elt];
        cq3[nq] = cnp3[eltv];
        cq4[nq] = cnp3[elt];
      }
      else if (cnp1[eltv] == cnp1[elt])
      {
        cq1[nq] = cnp1[elt];
        cq2[nq] = cnp2[elt];
        cq3[nq] = cnp3[elt];
        cq4[nq] = cnp3[eltv];
      }
      else
      {
        cq1[nq] = cnp1[elt];
        cq2[nq] = cnp3[eltv];
        cq3[nq] = cnp2[elt];
        cq4[nq] = cnp3[elt];
      }
      nq++;
    }
  }

  goto debut;

  fin: ;
  // Retourne une liste faite de 2 arrays (un QUAD et un TRI)
  // Array QUAD
  cq.reAllocMat(nq, 4);
  FldArrayF fq(*f);
  K_CONNECT::cleanConnectivity(posx, posy, posz, 1.e-12, "QUAD", fq, cq);
  
  // Array TRI
  E_Int nt = 0;
  FldArrayI ct(ne, 3);
  E_Int* ct1 = ct.begin(1);
  E_Int* ct2 = ct.begin(2);
  E_Int* ct3 = ct.begin(3);
  FldArrayF ft(*f);
  for (i = 0; i < ne; i++)
  {
    if (mergedp[i] == 0) 
    { ct1[nt] = cnp1[i]; ct2[nt] = cnp2[i]; ct3[nt] = cnp3[i]; nt++; }
  }
  ct.reAllocMat(nt, 3);
  K_CONNECT::cleanConnectivity(posx, posy, posz, 1.e-12, "TRI", ft, ct);
  
  PyObject* tpl = PyList_New(0);
  PyObject* o;
  o = K_ARRAY::buildArray(fq, varString, cq, 3, NULL);
  PyList_Append(tpl, o); Py_DECREF(o);
  o = K_ARRAY::buildArray(ft, varString, ct, 2, NULL);
  PyList_Append(tpl, o); Py_DECREF(o);
  
  return tpl;
}

//=============================================================================
/* Retourne l'ecart max a 90 degres dans un Quad */
//=============================================================================
E_Float getMaxAngle(E_Float* ptA, E_Float* ptB, E_Float* ptC, E_Float* ptD)
{
  E_Int ret = K_COMPGEOM::checkQuadConvexity(ptA, ptB, ptC, ptD);
  if (ret != 0) return 370.; // very bad
  E_Float dirVect[3];
  // produit vectoriel AB x AC pour calculer dirVect
  E_Float dx1 = ptB[0]-ptA[0];
  E_Float dy1 = ptB[1]-ptA[1];
  E_Float dz1 = ptB[2]-ptA[2];
  E_Float dx2 = ptC[0]-ptA[0];
  E_Float dy2 = ptC[1]-ptA[1];
  E_Float dz2 = ptC[2]-ptA[2];
  E_Float p1 = dy1*dz2-dz1*dy2;
  E_Float p2 = dz1*dx2-dz2*dx1;
  E_Float p3 = dx1*dy2-dy1*dx2;
  E_Float n = sqrt(p1*p1+p2*p2+p3*p3);
  if (n > 1.e-12) 
  {
    dirVect[0] = p1/n; dirVect[1] = p2/n; dirVect[2] = p3/n; 
  }
  else
  {
    dirVect[0] = 0.; dirVect[1] = 0.; dirVect[2] = 1.;
  }
  E_Float angle1 = K_COMPGEOM::getAlphaAngleBetweenBars(ptA, ptD,
                                                        ptA, ptB, dirVect);
  E_Float angle2 = K_COMPGEOM::getAlphaAngleBetweenBars(ptB, ptA,
                                                        ptB, ptC, dirVect);
  E_Float angle3 = K_COMPGEOM::getAlphaAngleBetweenBars(ptC, ptB,
                                                        ptC, ptD, dirVect);
  E_Float angle4 = K_COMPGEOM::getAlphaAngleBetweenBars(ptD, ptC,
                                                        ptD, ptA, dirVect);
  //printf("A %f %f %f\n",ptA[0], ptA[1], ptA[2]);
  //printf("B %f %f %f\n",ptB[0], ptB[1], ptB[2]);
  //printf("C %f %f %f\n",ptC[0], ptC[1], ptC[2]);
  //printf("D %f %f %f\n",ptD[0], ptD[1], ptD[2]);
  //printf("%f %f %f %f\n", angle1, angle2, angle3, angle4);

  if (angle1 > 180.) angle1 = 360.-angle1;
  if (angle2 > 180.) angle2 = 360.-angle2;
  if (angle3 > 180.) angle3 = 360.-angle3;
  if (angle4 > 180.) angle4 = 360.-angle4;
  
  E_Float angle = 0.;
  angle = K_FUNC::E_max(angle, E_abs(angle1-90.));
  angle = K_FUNC::E_max(angle, E_abs(angle2-90.));
  angle = K_FUNC::E_max(angle, E_abs(angle3-90.));
  angle = K_FUNC::E_max(angle, E_abs(angle4-90.));

  return angle;
}
