/*    
    Copyright 2013-2026 ONERA.

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

# include "generator.h"
# include <math.h>       // for the maths

using namespace K_FLD;
using namespace std;

//===========================================================================
/* Calcul si les Bounding box de 2 arrays s'intersectent */
//===========================================================================
PyObject* K_GENERATOR::bboxIntersection(PyObject* self, PyObject* args)
{
  PyObject* array1; PyObject* array2;
  E_Float tol;
  if (!PYPARSETUPLE_(args, OO_ R_,
                    &array1, &array2, &tol)) return NULL;

  // Check array1
  E_Int ni1, nj1, nk1;
  FldArrayF* f1; FldArrayI* cn1;
  char* varString1; char* eltType1;
  E_Int res1 = K_ARRAY::getFromArray3(array1, varString1, f1, ni1, nj1, nk1, cn1, eltType1);
    
  if (res1 != 1 && res1 != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "bboxIntersection: array1 is invalid.");
    return NULL;
  }
  // Check array2
  E_Int ni2, nj2, nk2;
  FldArrayF* f2;
  FldArrayI* cn2;
  char* varString2;
  char* eltType2;
  E_Int res2 = K_ARRAY::getFromArray3(array2, varString2, f2, ni2, nj2, nk2, cn2, eltType2);
  
  if (res2 != 1 && res2 != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "bboxIntersection: array2 is invalid.");
    return NULL;
    RELEASESHAREDB(res1, array1, f1, cn1); 
  }

  E_Int isIntersect = 0;
  // determination de x,y et z ds les arrays
  E_Int posx1 = K_ARRAY::isCoordinateXPresent(varString1);
  E_Int posy1 = K_ARRAY::isCoordinateYPresent(varString1);
  E_Int posz1 = K_ARRAY::isCoordinateZPresent(varString1);
  E_Int posx2 = K_ARRAY::isCoordinateXPresent(varString2);
  E_Int posy2 = K_ARRAY::isCoordinateYPresent(varString2);
  E_Int posz2 = K_ARRAY::isCoordinateZPresent(varString2);
  if (posx1 == -1 || posy1 == -1 || posz1 == -1)
  {
    RELEASESHAREDB(res1, array1, f1, cn1); 
    RELEASESHAREDB(res2, array2, f2, cn2); 
    PyErr_SetString(PyExc_TypeError,
                    "bboxIntersection: can't find coordinates in array1.");
    return NULL;
  }
  if (posx2 == -1 || posy2 == -1 || posz2 == -1)
  {
    RELEASESHAREDB(res1, array1, f1, cn1); 
    RELEASESHAREDB(res2, array2, f2, cn2); 
    PyErr_SetString(PyExc_TypeError,
                    "bboxIntersection: can't find coordinates in array2.");
    return NULL;
  }
  posx1++; posy1++; posz1++; 
  posx2++; posy2++; posz2++;

  // bbox de a1
  E_Int npts1 = f1->getSize();
  E_Float xmin1, ymin1, zmin1, xmax1, ymax1, zmax1;
  K_COMPGEOM::boundingBoxUnstruct(npts1, f1->begin(posx1), f1->begin(posy1), f1->begin(posz1), 
                                  xmin1, ymin1, zmin1, xmax1, ymax1, zmax1);

  // bbox de a2
  E_Int npts2 = f2->getSize();
  E_Float xmin2, ymin2, zmin2, xmax2, ymax2, zmax2;
  K_COMPGEOM::boundingBoxUnstruct(npts2, f2->begin(posx2), f2->begin(posy2), f2->begin(posz2), 
                                  xmin2, ymin2, zmin2, xmax2, ymax2, zmax2);
  
  RELEASESHAREDB(res1, array1, f1, cn1); 
  RELEASESHAREDB(res2, array2, f2, cn2); 
  if ((xmax1 > xmin2-tol && xmin1 < xmax2+tol) &&
      (ymax1 > ymin2-tol && ymin1 < ymax2+tol) &&
      (zmax1 > zmin2-tol && zmin1 < zmax2+tol)) 
    isIntersect = 1;

  return Py_BuildValue(I_, isIntersect);
}

//===========================================================================
/* Determine if an Axis-Aligned Bounding Box (zone1) intersects an
   Oriented Bounding Box (zone2)
   (zone, in place) */
//===========================================================================
PyObject* K_GENERATOR::_bboxIntersectionZ(PyObject* self, PyObject* args)
{
  char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;
  PyObject* zone1; PyObject* zone2; E_Float tol;
  
  if (!PYPARSETUPLE_(args, OO_ R_ SSS_,
                    &zone1, &zone2, &tol, &GridCoordinates, &FlowSolutionNodes,
                    &FlowSolutionCenters)) return NULL;

  // Checks coordinates of zone 1
  vector<PyArrayObject*> hook1;
  E_Int im1, jm1, km1, cnSize1, cnNfld1;
  char* varString1; char* eltType1;
  vector<E_Float*> fields1; vector<E_Int> locs1;
  vector<E_Int*> cn1;
  K_PYTREE::getFromZone(zone1, 1, -1, varString1, fields1, locs1, im1, jm1, km1, 
                        cn1, cnSize1, cnNfld1, eltType1, hook1,
                        GridCoordinates, FlowSolutionNodes, FlowSolutionCenters);
  E_Int posx1 = K_ARRAY::isCoordinateXPresent(varString1);
  E_Int posy1 = K_ARRAY::isCoordinateYPresent(varString1);
  E_Int posz1 = K_ARRAY::isCoordinateZPresent(varString1);
  delete [] varString1;
  if (posx1 == -1 || posy1 == -1 || posz1 == -1)
  {
    RELEASESHAREDZ(hook1, (char*)NULL, (char*)NULL);
    PyErr_SetString(PyExc_TypeError,
                    "bboxIntersection: cannot find coordinates in zone1.");
    return NULL;
  }

  // Checks coordinates of zone 2
  vector<PyArrayObject*> hook2;
  E_Int im2, jm2, km2, cnSize2, cnNfld2;
  char* varString2; char* eltType2;
  vector<E_Float*> fields2; vector<E_Int> locs2;
  vector<E_Int*> cn2;
  K_PYTREE::getFromZone(zone2, 1, -1, varString2, fields2, locs2, im2, jm2, km2, 
                        cn2, cnSize2, cnNfld2, eltType2, hook2,
                        GridCoordinates, FlowSolutionNodes, FlowSolutionCenters);
  E_Int posx2 = K_ARRAY::isCoordinateXPresent(varString2);
  E_Int posy2 = K_ARRAY::isCoordinateYPresent(varString2);
  E_Int posz2 = K_ARRAY::isCoordinateZPresent(varString2);
  delete [] varString2;
  if (posx2 == -1 || posy2 == -1 || posz2 == -1)
  {
    RELEASESHAREDZ(hook1, (char*)NULL, (char*)NULL);
    RELEASESHAREDZ(hook2, (char*)NULL, (char*)NULL);
    PyErr_SetString(PyExc_TypeError,
                    "bboxIntersection: cannot find coordinates in zone2.");
    return NULL;
  }

  E_Float* xt1 = fields1[posx1]; 
  E_Float* yt1 = fields1[posy1]; 
  E_Float* zt1 = fields1[posz1]; 
  E_Int nt1 = im1 * jm1 * km1;
  
  E_Float* xt2 = fields2[posx2]; 
  E_Float* yt2 = fields2[posy2]; 
  E_Float* zt2 = fields2[posz2];
  E_Int nt2 = im2 * jm2 * km2;
  
  // AABB of zone 1:
  E_Float xmin1 = xt1[0];
  E_Float xmax1 = xt1[0];
  E_Float ymin1 = yt1[0];
  E_Float ymax1 = yt1[0];
  E_Float zmin1 = zt1[0];
  E_Float zmax1 = zt1[0];
  
  for (E_Int i=0; i < nt1; i++)
  {
     if (xt1[i] >= xmax1)  {xmax1 = xt1[i];}
  }
  for (E_Int i=0; i < nt1; i++)
  {
     if (xt1[i] <= xmin1)  {xmin1 = xt1[i];}
  }
  for (E_Int i=0; i < nt1; i++)
  {
     if (yt1[i] >= ymax1)  {ymax1 = yt1[i];}
  }

  for (E_Int i=0; i < nt1; i++)
  {
     if (yt1[i] <= ymin1) {ymin1 = yt1[i];}
  }

  for (E_Int i=0; i < nt1; i++)
  {
     if (zt1[i] >= zmax1) {zmax1 = zt1[i];}
  }
  for (E_Int i=0; i < nt1; i++)
  {
     if (zt1[i] <= zmin1) {zmin1 = zt1[i];}
  }  
  
  // AABB of zone 2:
  E_Float xmin2 = xt2[0];
  E_Float xmax2 = xt2[0];
  E_Float ymin2 = yt2[0];
  E_Float ymax2 = yt2[0];
  E_Float zmin2 = zt2[0];
  E_Float zmax2 = zt2[0];
  
  for (E_Int i=0; i < nt2; i++)
  {
     if (xt2[i] >= xmax2)  {xmax2 = xt2[i];}
  }
  for (E_Int i=0; i < nt2; i++)
  {
     if (xt2[i] <= xmin2)  {xmin2 = xt2[i];}
  }
  for (E_Int i=0; i < nt2; i++)
  {
     if (yt2[i] >= ymax2)  {ymax2 = yt2[i];}
  }

  for (E_Int i=0; i < nt2; i++)
  {
     if (yt2[i] <= ymin2) {ymin2 = yt2[i];}
  }

  for (E_Int i=0; i < nt2; i++)
  {
     if (zt2[i] >= zmax2) {zmax2 = zt2[i];}
  }
  for (E_Int i=0; i < nt2; i++)
  {
     if (zt2[i] <= zmin2) {zmin2 = zt2[i];}
  }

  RELEASESHAREDZ(hook1, (char*)NULL, (char*)NULL);
  RELEASESHAREDZ(hook2, (char*)NULL, (char*)NULL);
  
  E_Int isIntersect = 0;
  if ((xmax1 > xmin2-tol && xmin1 < xmax2+tol) &&
      (ymax1 > ymin2-tol && ymin1 < ymax2+tol) &&
      (zmax1 > zmin2-tol && zmin1 < zmax2+tol)) 
     {isIntersect = 1;}

  return Py_BuildValue(I_, isIntersect);
}
