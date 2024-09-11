/*    
    Copyright 2013-2024 Onera.

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

using namespace K_FLD;
using namespace std;

//=============================================================================
// Selectionne les elements interieurs a des courbes donnees
// IN: TRI array
// IN: liste de courbes BAR-arrays
//=============================================================================
PyObject* K_GENERATOR::selectInsideElts(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* curvesList;
  if ( !PyArg_ParseTuple(args, "OO", &array, &curvesList ) )
  {
    return NULL;
  }

  // Check array
  E_Int im, jm, km;
  FldArrayF* f;
  char* varString;
  char* eltType;
  FldArrayI* cn;
  E_Int res = 
    K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType);

  if (res != 2)
  {
    if ( res == 1) {delete f;}
    PyErr_SetString(PyExc_TypeError,
                    "selectInsideElts: array must be a TRI-array.");
    return NULL;
  }
  if (strcmp(eltType, "TRI") != 0)
  {
    delete f; delete cn;
    PyErr_SetString(PyExc_TypeError,
                    "selectInsideElts: array must be a TRI-array.");
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if ( posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(
      PyExc_TypeError,
      "selectInsideElts: array must contain variables (x,y,z).");
    delete f; delete cn;
    return NULL;
  }
  posx++; posy++; posz++;

  // Check every array in curvesList
  if (PyList_Check(curvesList) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "selectInsideElts: second argument must be a list.");
    return NULL;
  }

  // Extract infos from arrays
  vector<E_Int> resl;
  vector<char*> structVarString;
  vector<char*> unstrVarString;
  vector<FldArrayF*> structF;
  vector<FldArrayF*> unstrF;
  vector<E_Int> nit;
  vector<E_Int> njt; 
  vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltTypet;
  vector<PyObject*> objst, objut;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = false;
  E_Boolean skipDiffVars = true;

  E_Int isOk = K_ARRAY::getFromArrays(
    curvesList, resl, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltTypet, objst, objut, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured);
  
  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "selectInsideElts: invalid list of arrays.");
    for (size_t i = 0; i < structF.size(); i++) delete structF[i];
    for (size_t i = 0; i < unstrF.size(); i++) delete unstrF[i];
    for (size_t i = 0; i < cnt.size(); i++) delete cnt[i];
    for (size_t i = 0; i < structVarString.size(); i++) delete[] structVarString[i];
    for (size_t i = 0; i < unstrVarString.size(); i++) delete[] unstrVarString[i];
    for (size_t i = 0; i < eltTypet.size(); i++) delete[] eltTypet[i];
    for (size_t i = 0; i < objst.size(); i++) Py_DECREF(objst[i]);
    for (size_t i = 0; i < objut.size(); i++) Py_DECREF(objut[i]);
    delete f; delete cn;
    return NULL;
  }
  if (structF.size() > 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "selectInsideElts: arrays must be unstructured.");
    for (size_t i = 0; i < structF.size(); i++) delete structF[i];
    for (size_t i = 0; i < unstrF.size(); i++) delete unstrF[i];
    for (size_t i = 0; i < cnt.size(); i++) delete cnt[i];
    for (size_t i = 0; i < structVarString.size(); i++) delete[] structVarString[i];
    for (size_t i = 0; i < unstrVarString.size(); i++) delete[] unstrVarString[i];
    for (size_t i = 0; i < eltTypet.size(); i++) delete[] eltTypet[i];
    for (size_t i = 0; i < objst.size(); i++) Py_DECREF(objst[i]);
    for (size_t i = 0; i < objut.size(); i++) Py_DECREF(objut[i]);
    delete f; delete cn;
    return NULL;
  }
  for (size_t i = 0; i < eltTypet.size(); i++)
  {
    if (strcmp(eltTypet[i], "BAR") != 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "selectInsideElts: arrays must be BAR-arrays.");
      for (size_t i = 0; i < structF.size(); i++) delete structF[i];
      for (size_t i = 0; i < unstrF.size(); i++) delete unstrF[i];
      for (size_t i = 0; i < cnt.size(); i++) delete cnt[i];
      for (size_t i = 0; i < structVarString.size(); i++) delete[] structVarString[i];
      for (size_t i = 0; i < unstrVarString.size(); i++) delete[] unstrVarString[i];
      for (size_t i = 0; i < eltTypet.size(); i++) delete[] eltTypet[i];
      for (size_t i = 0; i < objst.size(); i++) Py_DECREF(objst[i]);
      for (size_t i = 0; i < objut.size(); i++) Py_DECREF(objut[i]);
      delete f; delete cn;
      return NULL;
    }
  }

  // Construction des courbes
  E_Int ncurves = unstrF.size();
  vector<FldArrayF*> curves;
  for (E_Int i = 0; i < ncurves; i++)
  {
    FldArrayF& tp = *unstrF[i];
    FldArrayF* t = new FldArrayF(tp.getSize(), 3);
    E_Int posxt = K_ARRAY::isCoordinateXPresent(unstrVarString[i])+1;
    E_Int posyt = K_ARRAY::isCoordinateYPresent(unstrVarString[i])+1;
    E_Int poszt = K_ARRAY::isCoordinateZPresent(unstrVarString[i])+1;
    t->setOneField(tp, posxt, 1);
    t->setOneField(tp, posyt, 2);
    t->setOneField(tp, poszt, 3);
    curves.push_back(t);
  }
  
  // liberation de la memoire
  for (E_Int i = 0; i < ncurves; i++) delete unstrF[i];

  // Calcul des centres des elements
  FldArrayF& fp = *f;
  FldArrayI& cnp = *cn;
  FldArrayF centers(cnp.getSize(), 3);
  E_Int sizeElt = cnp.getNfld();
  E_Float iElt = 1./sizeElt;
  E_Float xc, yc, zc;
  E_Int ind;
  for (E_Int i = 0; i < cnp.getSize(); i++)
  {
    xc = 0.; yc = 0.; zc = 0.;
    for (E_Int j = 1; j <= sizeElt; j++)
    {
      ind = cnp(i, j)-1;
      xc = xc + fp(ind,posx);
      yc = yc + fp(ind,posy);
      zc = zc + fp(ind,posz);
    }
    xc = xc*iElt; yc = yc*iElt; zc = zc*iElt;
    centers(i, 1) = xc; centers(i, 2) = yc; centers(i, 3) = zc;
  }

  // Test des centres des elements pour chaque courbe
  FldArrayI* connect = new FldArrayI(cn->getSize(), cn->getNfld());
  FldArrayI& connectp = *connect;
  
  E_Int test, ne;
  E_Float p[3];
  
  ne = 0;
  for (E_Int ind = 0; ind < cnp.getSize(); ind++)
  {
    p[0] = centers(ind,1); p[1] = centers(ind,2); p[2] = centers(ind,3);
    test = 0;
    for (E_Int i = 0; i < ncurves; i++)
    {
      test = test + K_COMPGEOM::pointInPolygon2D(curves[i]->begin(1), curves[i]->begin(2), curves[i]->begin(3),
                                                 *cnt[i],
                                                 p,
                                                 NULL, NULL);
    }
    if (test - (test/2)*2. != 0) // dedans
    {
      for (E_Int i = 1; i <= sizeElt; i++)
        connectp(ne,i) = cnp(ind,i);
      ne++;
    }
  }
  connectp.reAllocMat(ne, sizeElt);

  // sortie
  PyObject* tpl;
  tpl = K_ARRAY::buildArray(*f, varString, connectp, -1, eltType);
  
  for (size_t i = 0; i < cnt.size(); i++) delete cnt[i];
  //for (size_t i = 0; i < objut.size(); i++) Py_DECREF(objut[i]);
  for (E_Int i = 0; i < ncurves; i++) delete curves[i];
  delete f; delete cn; delete connect;
  
  return tpl;
}
