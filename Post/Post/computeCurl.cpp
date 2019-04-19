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
# include <string.h>
# include "post.h"

using namespace K_FLD;
using namespace std;

extern "C"
{
  void k6conv2center1_(const E_Int& ni, const E_Int& nj, const E_Int& nk, 
                       const E_Int& nfld, E_Float* fieldnode, 
                       E_Float* fieldcenter);

  void k6compstructcurlt_(
    const E_Int& ni, const E_Int& nj, const E_Int& nk, 
    const E_Int& ncells, const E_Int& nint, 
    E_Float* xt, E_Float* yt, E_Float* zt,
    E_Float* ux, E_Float* uy, E_Float* uz, 
    E_Float* rotx, E_Float* roty, E_Float* rotz,
    E_Float* surf, E_Float* snorm, E_Float* centerInt, 
    E_Float* vol, E_Float* uintx, E_Float* uinty, E_Float* uintz);

  void k6compstructcurl2dt_(
    const E_Int& ni, const E_Int& nj, const E_Int& ncells, 
    const E_Float* xt, const E_Float* yt, const E_Float* zt,  
    const E_Float* ux, const E_Float* uy, const E_Float* uz,
    E_Float* rotx, E_Float* roty, E_Float* rotz);

  void k6compunstrcurl_(E_Int& dim,E_Int& npts,E_Int& nelts,E_Int& nedges, 
                        E_Int& nnodes,E_Int* cn,E_Float* xt,E_Float* yt,
                        E_Float* zt,E_Float* snx,E_Float* sny,E_Float* snz,
                        E_Float* surf, E_Float* vol, 
                        E_Float* uintx, E_Float* uinty, E_Float* uintz, 
                        E_Float* ux, E_Float* uy, E_Float* uz, 
                        E_Float* xint, E_Float* yint, E_Float* zint, 
                        E_Float* rotx, E_Float* roty, E_Float* rotz);
  
  void k6compunstrcurl2d_(E_Int& npts, E_Int& nelts, E_Int&nedges, 
                          E_Int& nnodes, E_Int* cn, E_Float* xt,
                          E_Float* yt, E_Float* zt, 
                          E_Float* snx, E_Float* sny, E_Float* snz, 
                          E_Float* surf, E_Float* vol,
                          E_Float* ux, E_Float* uy, E_Float* uz, 
                          E_Float* rotx, E_Float* roty, E_Float* rotz);
}

//=============================================================================
/* Calcul du rotationnel d'un champ defini par un vecteur (u,v,w) en noeuds. 
   Le rotationnel est fourni aux centres des cellules */
//=============================================================================
PyObject* K_POST::computeCurl(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* vars0;
  if (!PyArg_ParseTuple(args, "OO", &array, &vars0)) return NULL;
  
  //extraction des variables constituant le vecteur dont le rot est calcule
  vector<char*> vars;
  if (PyList_Check(vars0) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeCurl: a list of 3 variables for curl computation must be defined.");
    return NULL; 
  }
  if (PyList_Size(vars0) != 3)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeCurl: 3 variables must be defined to extract the curl.");
    return NULL;
  }
  for (int i = 0; i < PyList_Size(vars0); i++)
  {
    PyObject* tpl0 = PyList_GetItem(vars0, i);
    if (PyString_Check(tpl0))
    {
      char* str = PyString_AsString(tpl0);
      vars.push_back(str);
    }
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(tpl0))
    {
      char* str = PyBytes_AsString(PyUnicode_AsUTF8String(tpl0));
      vars.push_back(str);  
    }
#endif
    else  
    {
      PyErr_SetString(PyExc_TypeError,
                      "computeCurl: varname must be a string.");
      return NULL;
    }
  }
  // Check array
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk;// number of points of array
  E_Int res = K_ARRAY::getFromArray(array, varString, f, ni, nj, nk, cn, 
                                    eltType, true);
  if (res != 1 && res != 2) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeCurl: invalid array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }

  E_Int posu = K_ARRAY::isNamePresent(vars[0], varString);
  E_Int posv = K_ARRAY::isNamePresent(vars[1], varString);
  E_Int posw = K_ARRAY::isNamePresent(vars[2], varString);
  if (posu == -1 || posv == -1 || posw == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeCurl: one variable was not found in array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  posu++; posv++; posw++;
  
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeCurl: coordinates not found in array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  
  posx++; posy++; posz++;
  E_Int nfld = f->getNfld();
  if (nfld < 6)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeCurl: no field to compute.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }

  // calcul du rotationnel 
  PyObject* tpl;
  char* varStringOut = new char[15];
  strcpy(varStringOut, "rotx,roty,rotz"); 
  if (res == 1)
  {
    E_Int ni1 = K_FUNC::E_max(1,ni-1);
    E_Int nj1 = K_FUNC::E_max(1,nj-1);
    E_Int nk1 = K_FUNC::E_max(1,nk-1);
    E_Int ncells = ni1*nj1*nk1;
    tpl = K_ARRAY::buildArray(3, varStringOut, ni1, nj1, nk1);   
    E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF fp(ncells, 3, fnp, true); 
    computeCurlStruct(ni, nj, nk, 
                      f->begin(posx), f->begin(posy), f->begin(posz),
                      f->begin(posu), f->begin(posv), f->begin(posw),
                      fp.begin(1), fp.begin(2), fp.begin(3));
  }
  else // non structure
  {
    if (strcmp(eltType, "TRI")   != 0 &&
        strcmp(eltType, "QUAD")  != 0 &&
        strcmp(eltType, "TETRA") != 0 &&
        strcmp( eltType, "HEXA") != 0 &&
        strcmp(eltType, "PENTA") != 0) 
    {
      PyErr_SetString(PyExc_TypeError,
                      "computeNormCurl: not a valid element type.");
      RELEASESHAREDU(array,f, cn); return NULL;
    }
    
    E_Int npts = f->getSize();
    tpl = K_ARRAY::buildArray(3, varStringOut, npts, cn->getSize(), -1, eltType, true, cn->getSize()*cn->getNfld());
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    K_KCORE::memcpy__(cnnp, cn->begin(), cn->getSize()*cn->getNfld());
    E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
    E_Int nelts = cn->getSize();    
    FldArrayF fp(nelts, 3, fnp, true);
  
    // calcul du rotationnel aux centres des elements
    computeCurlNS(eltType, npts, *cn, 
                  f->begin(posx), f->begin(posy), f->begin(posz),
                  f->begin(posu), f->begin(posv), f->begin(posw),
                  fp.begin(1), fp.begin(2), fp.begin(3));              
  } 
  RELEASESHAREDB(res, array, f, cn);
  delete [] varStringOut;
  return tpl;
}
//==============================================================================
E_Int K_POST::computeCurlStruct(E_Int ni, E_Int nj, E_Int nk, 
                                E_Float* xt, E_Float* yt, E_Float* zt,
                                E_Float* ux, E_Float* uy, E_Float* uz,
                                E_Float* rotx, E_Float* roty, E_Float* rotz)
{
  if ((ni == 1 && nj == 1) || (ni == 1 && nk == 1) || (nj == 1 && nk == 1))
    return -1;
     
  E_Int dim = 3;
  if (ni == 1)
  {
    ni = nj; nj = nk; nk = 2; dim = 2;
  }
  else if (nj == 1)
  {
    nj = nk; nk = 2; dim = 2;
  }
  else if (nk == 1)
  {
    nk = 2; dim = 2;
  }
  // Calcul du rotationnel aux centres
  E_Int ni1 = K_FUNC::E_max(1,ni-1);
  E_Int nj1 = K_FUNC::E_max(1,nj-1);
  E_Int nk1 = K_FUNC::E_max(1,nk-1);
  E_Int ncells = ni1*nj1*nk1;
  E_Int nint = ni*nj1*nk1 + ni1*nj*nk1 + ni1*nj1*nk;
  
  FldArrayF surf(nint, 3);
  FldArrayF snorm(nint);
  FldArrayF centerInt(nint,3);
  FldArrayF vol(ncells);
  FldArrayF uint(nint,3);
  if (dim == 2) 
    k6compstructcurl2dt_(ni, nj, ncells, 
                         xt, yt, zt, ux, uy, uz,
                         rotx, roty, rotz);  
  else 
    k6compstructcurlt_(ni, nj, nk, ncells, nint, 
                       xt, yt, zt, ux, uy, uz,
                       rotx, roty, rotz,
                       surf.begin(), snorm.begin(),
                       centerInt.begin(), vol.begin(), 
                       uint.begin(1), uint.begin(2), uint.begin(3)); 
  return 1;
}
//=============================================================================
E_Int K_POST::computeCurlNS(char* eltType, E_Int npts, FldArrayI& cn, 
                            E_Float* xt, E_Float* yt, E_Float* zt,
                            E_Float* ux, E_Float* uy, E_Float* uz,
                            E_Float* rotx, E_Float* roty, E_Float* rotz)
{
  E_Int nelts = cn.getSize();
  E_Int nnodes; //nb de noeuds par elts
  E_Int nedges; //nb de facettes par elts
  E_Int dim = 3; 
  if (strcmp(eltType, "TRI") == 0) 
  {
    nnodes = 3; nedges = 3; dim = 2;
  }
  else if (strcmp(eltType, "QUAD") == 0) 
  {
    nnodes = 4; nedges = 4; dim = 2;
  }
  else if (strcmp(eltType, "TETRA") == 0)
  {
    nedges = 4; nnodes = 4;
  }
  else if (strcmp( eltType, "HEXA") == 0) 
  {
    nedges = 6; nnodes = 8;
  }
  else if (strcmp(eltType, "PENTA") == 0) 
  {
    nedges = 5; nnodes = 6;
  }

  if (dim == 2)
  {
    FldArrayF snx(nelts,1);
    FldArrayF sny(nelts,1);
    FldArrayF snz(nelts,1);
    FldArrayF surf(nelts,1);
    FldArrayF vol(nelts);
    k6compunstrcurl2d_(npts, nelts, nedges, nnodes, cn.begin(),
                       xt, yt, zt, 
                       snx.begin(), sny.begin(), snz.begin(), 
                       surf.begin(), vol.begin(),
                       ux, uy, uz, rotx, roty, rotz);
  }
  else //dim = 3
  {
    FldArrayF snx(nelts, nedges);
    FldArrayF sny(nelts, nedges);
    FldArrayF snz(nelts, nedges);
    FldArrayF surf(nelts, nedges);
    FldArrayF vol(nelts);
    FldArrayF uintx(nelts, nedges);
    FldArrayF uinty(nelts, nedges);
    FldArrayF uintz(nelts, nedges);
    FldArrayF xint(nelts, nedges);
    FldArrayF yint(nelts, nedges);
    FldArrayF zint(nelts, nedges);
    k6compunstrcurl_(dim, npts, nelts, 
                     nedges, nnodes, cn.begin(),  
                     xt, yt, zt, 
                     snx.begin(), sny.begin(), snz.begin(), 
                     surf.begin(), vol.begin(),
                     uintx.begin(), uinty.begin(), uintz.begin(), 
                     ux, uy, uz, xint.begin(), yint.begin(), zint.begin(), 
                     rotx, roty, rotz);
  }
  return 1;
}
