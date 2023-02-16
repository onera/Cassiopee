/*    
    Copyright 2013-2023 Onera.

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
# include "connector.h"

using namespace std;
using namespace K_FUNC;
using namespace K_FLD;

//=============================================================================
/*  Nouvel algo changeWall (in place) */
//=============================================================================
PyObject* K_CONNECTOR::changeWall2(PyObject* self, PyObject* args)
{
  PyObject* array; //domaine a adapter avec cellN
  PyObject* walls1; // liste walls de la zone arrayCenters correspondant au mismatch
  PyObject* walls2; // liste walls en mismatch avec normales
  E_Float hCL; // hauteur de la couche limite (limite le recollement)
  PyObject* dirList; // direction du walls1 (1: i, 2: j; 3: k)

  if (!PYPARSETUPLE(args, "OOOdO", "OOOdO", "OOOfO", "OOOfO",
                    &array, &walls1, &walls2, &hCL, &dirList))
     return NULL;

  // Check array - zone a modifier
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray2(array, varString, f, nil, njl, nkl, 
                                     cn, eltType);
  
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "changeWall2: invalid array.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
   
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "changeWall2: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  E_Int posCellN = K_ARRAY::isCellNatureField1Present(varString);
  if (posCellN == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "changeWall2: can't find cellN in array.");
    return NULL;
  }
  posCellN++;

  // Check walls 1 (liste TRI)
  E_Int nw1 = PyList_Size(walls1);
  vector<E_Int> dirs(nw1);
  vector<FldArrayF*> fw1(nw1);
  vector<FldArrayI*> cw1(nw1);
  vector<E_Int> posx1(nw1);
  vector<E_Int> posy1(nw1);
  vector<E_Int> posz1(nw1);
  for (E_Int i = 0; i < nw1; i++)
  {
    PyObject* o = PyList_GetItem(walls1, i);
    E_Int ni2,nj2,nk2;
    char* eltType2; char* varString2;
    E_Int res = K_ARRAY::getFromArray2(o, varString2, fw1[i], ni2, nj2, nk2, 
                                       cw1[i], eltType2);
    PyObject* v = PyList_GetItem(dirList, i);
    dirs[i] = PyLong_AsLong(v);
    posx1[i] = K_ARRAY::isCoordinateXPresent(varString2)+1;
    posy1[i] = K_ARRAY::isCoordinateYPresent(varString2)+1;
    posz1[i] = K_ARRAY::isCoordinateZPresent(varString2)+1;
  }

  // Check walls 2 (liste TRI)
  E_Int nw2 = PyList_Size(walls2);
  vector<FldArrayF*> fw2(nw2);
  vector<FldArrayI*> cw2(nw2);
  vector<E_Int> posx2(nw2);
  vector<E_Int> posy2(nw2);
  vector<E_Int> posz2(nw2);
  for (E_Int i = 0; i < nw2; i++)
  {
    PyObject* o = PyList_GetItem(walls2, i);
    E_Int ni2,nj2,nk2;
    char* eltType2; char* varString2;
    E_Int res = K_ARRAY::getFromArray2(o, varString2, fw2[i], ni2, nj2, nk2, 
                                       cw2[i], eltType2);
    posx2[i] = K_ARRAY::isCoordinateXPresent(varString2)+1;
    posy2[i] = K_ARRAY::isCoordinateYPresent(varString2)+1;
    posz2[i] = K_ARRAY::isCoordinateZPresent(varString2)+1;                                   
  }

  // Modification des points de array
  E_Int ind, noet, ind1;
  E_Float* fx = f->begin(posx);
  E_Float* fy = f->begin(posy);
  E_Float* fz = f->begin(posz);
  E_Float* fcellN = f->begin(posCellN);
  E_Float x,y,z,xp,yp,zp;
  E_Float dirx,diry,dirz;
  E_Float h1,h2;

  for (E_Int k = 0; k < nkl; k++)
  for (E_Int j = 0; j < njl; j++)
  for (E_Int i = 0; i < nil; i++)
  {
    ind = i + j*nil + k*nil*njl;
    
    if (fcellN[ind] == 2)
    {
      x = fx[ind]; y = fy[ind]; z = fz[ind];

      xp = K_CONST::E_MAX_FLOAT; yp = K_CONST::E_MAX_FLOAT; zp = K_CONST::E_MAX_FLOAT;

      // project on all walls1
      for (E_Int iw = 0; iw < nw1; iw++)
      {
        if (dirs[iw] == 1) { ind1 = ind1+1; }
        else if (dirs[iw] == 2) { ind1 = ind1+nil; }
        else { ind1 = ind1+nil*njl; }
        dirx = x - (*fw1[iw])(ind1,posx1[iw]);
        diry = y - (*fw1[iw])(ind1,posy1[iw]);
        dirz = z - (*fw1[iw])(ind1,posz1[iw]);
        
        noet = K_COMPGEOM::projectDir(x, y, z, dirx, diry, dirz, 
                                      fw1[iw]->begin(posx1[iw]), 
                                      fw1[iw]->begin(posy1[iw]),
                                      fw1[iw]->begin(posz1[iw]),
                                      *(cw1[iw]), xp, yp, zp);
        h1 = (x-xp)*(x-xp)+(y-yp)*(y-yp)+(z-zp)*(z-zp);
      }
      // project on all walls2
      for (E_Int iw = 0; iw < nw2; iw++)
      {
        if (dirs[iw] == 1) { ind1 = ind1+1; }
        else if (dirs[iw] == 2) { ind1 = ind1+nil; }
        else { ind1 = ind1+nil*njl; }
        dirx = x - (*fw1[iw])(ind1,posx1[iw]);
        diry = y - (*fw1[iw])(ind1,posy1[iw]);
        dirz = z - (*fw1[iw])(ind1,posz1[iw]);
        
        noet = K_COMPGEOM::projectDir(x, y, z, dirx, diry, dirz, 
                                      fw1[iw]->begin(posx1[iw]), 
                                      fw1[iw]->begin(posy1[iw]),
                                      fw1[iw]->begin(posz1[iw]),
                                      *(cw1[iw]), xp, yp, zp);
        h1 = (x-xp)*(x-xp)+(y-yp)*(y-yp)+(z-zp)*(z-zp);
      }

    }
  }
  

  RELEASESHAREDB(res, array, f, cn);
  for (E_Int i = 0; i < nw1; i++) RELEASESHAREDU(PyList_GetItem(walls1, i), fw1[i], cw1[i]);
  for (E_Int i = 0; i < nw2; i++) RELEASESHAREDU(PyList_GetItem(walls2, i), fw2[i], cw2[i]);
  
  Py_INCREF(Py_None);
  return Py_None;

}