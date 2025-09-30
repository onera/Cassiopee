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

#include "generator.h"

using namespace K_FLD;
using namespace std;
using namespace K_CONST;

//=============================================================================
/* Generates a set of Cartesian grids from an octree HEXA mesh */
//=============================================================================
PyObject* K_GENERATOR::octree2Struct(PyObject* self, PyObject* args)
{
  PyObject* octree; PyObject* listOfVmin;
  if (!PYPARSETUPLE_(args, OO_, &octree, &listOfVmin)) return NULL;
  if (PyList_Check(listOfVmin) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "octree2Struct: 2nd argument must be a list.");
    return NULL;
  }
  E_Int nvmin = PyList_Size(listOfVmin);
  if (nvmin == 0)
  { 
    PyErr_SetString(PyExc_TypeError, 
                    "octree2Struct: vmin list is empty.");
    return NULL;
  }
  // Check array
  E_Int ni, nj, nk;
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(octree, varString, f, ni, nj, nk, 
                                     cn, eltType);
  E_Int dim = 3;
  if (res != 2) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "octree2Struct: array must be unstructured.");
    if (res == 1) RELEASESHAREDS(octree, f);
    return NULL;   
  }
  if (strcmp(eltType, "HEXA") == 0) dim = 3;
  else if (strcmp(eltType, "QUAD") == 0) dim = 2;
  else
  {
    PyErr_SetString(PyExc_TypeError, 
                    "octree2Struct: array must be HEXA or QUAD.");
    RELEASESHAREDU(octree,f,cn);
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "octree2Struct: coordinates not found in array.");
    RELEASESHAREDU(octree,f,cn);
    return NULL;
  }
  posx++; posy++; posz++;

  // determination du niveau le plus fin
  E_Float dhmin = 1.e10;
  E_Int nelts = cn->getSize();
  E_Float* xo = f->begin(posx);
  E_Float* yo = f->begin(posy);
  E_Float* zo = f->begin(posz);
  E_Int* cn1 = cn->begin(1);
  E_Int* cn2 = cn->begin(2);

  FldArrayI levels(nelts); E_Int* lp = levels.begin();
  E_Int levelMax = 0;

#pragma omp parallel default(shared)
  {
    E_Float dhminLoc = 1.e10;
    E_Float dhloc;
#pragma omp for
    for (E_Int et = 0; et < nelts; et++)
    {
      E_Int v1 = cn1[et]-1; E_Int v2 = cn2[et]-1;  
      dhloc = K_FUNC::E_abs(xo[v2]-xo[v1]);
      dhminLoc = K_FUNC::E_min(dhminLoc, dhloc);
    }
#pragma omp critical
    {
      if (dhminLoc < dhmin) { dhmin = dhminLoc; }
    }
  }
  
  //printf("dhmin = %f\n", dhmin);
#pragma omp parallel default(shared)
  { 
    E_Float dhloc;
    E_Int levelMaxLoc = 0;
#pragma omp for
    for (E_Int et = 0; et < nelts; et++)
    {
      E_Int v1 = cn1[et]-1; E_Int v2 = cn2[et]-1;  
      dhloc = K_FUNC::E_abs(xo[v2]-xo[v1]); 
      E_Int lloc = E_Int(dhloc/dhmin+1.e-6);
      E_Int found = 1; E_Int levelloc = 0;
      while (found < lloc){ found *= 2; levelloc++;}
      lp[et] = levelloc;
      levelMaxLoc = K_FUNC::E_max(levelloc, levelMaxLoc);
    }
#pragma omp critical
    {
      if (levelMaxLoc > levelMax) { levelMax = levelMaxLoc; }
    }
  }
  vector<E_Int> vmint;
  E_Int vminlok;
  if (nvmin <= levelMax)
  {
    for (E_Int v = 0; v < nvmin; v++)
    {
      PyObject* tpl0 = PyList_GetItem(listOfVmin,v);
      vminlok = PyLong_AsLong(tpl0);
      vmint.push_back(vminlok);
    }
    vminlok = vmint[nvmin-1];
    for (E_Int v = nvmin; v <= levelMax; v++)
      vmint.push_back(vminlok);
  }
  else // nvmin > levelMax
  {
    for (E_Int v = 0; v <= levelMax; v++)
    {
      PyObject* tpl0 = PyList_GetItem(listOfVmin,v);
      vminlok = PyLong_AsLong(tpl0);
      vmint.push_back(vminlok);
    }
  }
  //printf("levelMax %d - %d\n",levelMax,vmint.size());
  //for (E_Int i = 0; i < vmint.size(); i++) printf("%d\n",vmint[i]);
  //for (E_Int e = 0; e < nelts; e++) printf(" %d ",lp[e]);

  // decoupage de chaque QUAD ou HEXA en vmin points par direction
  // attention: QUAD en (x,y): z = 0.
  E_Int api = f->getApi();
  E_Int nfld = f->getNfld();
  PyObject* l = PyList_New(nelts);

  if (dim == 2) // octree QUAD 
  {
    E_Float** xtt = new E_Float* [nelts];
    for (E_Int et = 0; et < nelts; et++) // can not thread because of python
    {
      E_Int lloc = lp[et]; 
      E_Int vminloc = vmint[lloc];
      PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, vminloc, vminloc, 1, api);
      FldArrayF* coords;
      K_ARRAY::getFromArray3(tpl, coords);
      xtt[et] = coords->begin();
      PyList_SET_ITEM(l, et, tpl);
      RELEASESHAREDS(tpl, coords);
    }

#pragma omp parallel for default(shared)
    for (E_Int et = 0; et < nelts; et++)
    {
      E_Int lloc = lp[et];
      E_Int vminloc = vmint[lloc];
      E_Int npts = vminloc*vminloc;
      E_Float* xt = xtt[et];
      E_Float* yt = xt+npts;
      E_Float* zt = yt+npts;
      E_Int ind1 = cn1[et]-1; E_Int ind2 = cn2[et]-1;
      E_Float xmin = xo[ind1]; E_Float xmax = xo[ind2];
      E_Float ymin = yo[ind1]; 
      E_Float vmini = 1./(vminloc-1);
      E_Float dh = (xmax-xmin)*vmini;
      E_Int ind;

      for (E_Int j = 0; j < vminloc; j++)
        for (E_Int i = 0; i < vminloc; i++)
        {
          ind = i + j*vminloc; 
          xt[ind] = xmin + i*dh; yt[ind] = ymin + j*dh; zt[ind] = 0.;
        }
    }
    delete [] xtt;
  }
  else // octree HEXA ->3D structure
  {
    E_Float** xtt = new E_Float* [nelts];

    for (E_Int et = 0; et < nelts; et++)
    {
      E_Int lloc = lp[et]; 
      E_Int vminloc = vmint[lloc];
      PyObject* tpl = K_ARRAY::buildArray3(3, "x,y,z", vminloc, vminloc, vminloc, api);
      FldArrayF* coords;
      K_ARRAY::getFromArray3(tpl, coords);
      xtt[et] = coords->begin();
      PyList_SET_ITEM(l, et, tpl);
      RELEASESHAREDS(tpl, coords);
    }

#pragma omp parallel for default(shared)
    for (E_Int et = 0; et < nelts; et++)
    {
      E_Int lloc = lp[et]; 
      E_Int vminloc = vmint[lloc];
      E_Int vminloc2 = vminloc*vminloc;
      E_Int npts = vminloc2*vminloc;
      E_Float* xt = xtt[et];
      E_Float* yt = xt+npts;
      E_Float* zt = yt+npts;
      E_Int ind1 = cn1[et]-1; E_Int ind2 = cn2[et]-1;
      E_Float xmin = xo[ind1]; E_Float xmax = xo[ind2]; 
      E_Float ymin = yo[ind1]; E_Float zmin = zo[ind1];
      E_Float vmini = 1./(vminloc-1);
      E_Float dh = (xmax-xmin)*vmini;
      E_Int ind;

      for (E_Int k = 0; k < vminloc; k++)
        for (E_Int j = 0; j < vminloc; j++)
          for (E_Int i = 0; i < vminloc; i++)
          {
            ind = i + j*vminloc + k*vminloc2; 
            xt[ind] = xmin + i*dh; yt[ind] = ymin + j*dh; zt[ind] = zmin + k*dh;
          }
    }
    delete [] xtt;
  }
  RELEASESHAREDU(octree, f, cn);
  return l;
}
