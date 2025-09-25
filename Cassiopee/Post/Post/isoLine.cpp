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

// Build an isoline on a surface triangular mesh

# include "post.h"
using namespace std;
using namespace K_FLD;

//=============================================================================
/* Construit une isoline sur un maillage triangulaire surfacique */
//=============================================================================
PyObject* K_POST::isoLine(PyObject* self, PyObject* args)
{
  // surf: maillage surfacique (x,y,z+sol)
  // field: nom du field dont on cherche l'iso
  // value: valeur de l'iso
  PyObject* surf;
  char* field; E_Float value;
  if (!PYPARSETUPLE_(args, O_ S_ R_, &surf, &field, &value))
  {
    return NULL;
  }

  /*-----------------------------------------------*/
  /* Extraction des donnees du maillage surfacique */ 
  /*-----------------------------------------------*/
  char* varString0; char* eltType0;
  FldArrayF* f;
  E_Int nil, njl, nkl;
  FldArrayI* cn = NULL;
  E_Int res = 
    K_ARRAY::getFromArray3(surf, varString0, f, nil, njl, nkl, 
                           cn, eltType0);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "isoLine: input array is invalid.");
    return NULL;
  }
  if (res != 2 || (strcmp(eltType0, "TRI") != 0 && 
        strcmp(eltType0, "BAR") != 0))
  {
    RELEASESHAREDB(res, surf, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "isoLine: input array must be TRI or BAR.");
    return NULL;
  }
  
  // Check size of array
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString0);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString0);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString0);

  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "isoLine: coordinates not found in array.");
    RELEASESHAREDU(surf, f, cn);
    return NULL;
  }
  posx++; posy++; posz++;

  // position de la variable iso
  E_Int posf = K_ARRAY::isNamePresent(field, varString0);
  if (posf == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "isoLine: variable doesn't exist in array.");
    RELEASESHAREDU(surf, f, cn);
    return NULL;
  }
  posf++;

  // Position du cellN eventuellement
  E_Int poscellN = K_ARRAY::isCellNatureField2Present(varString0)+1;

  FldArrayF fiso; FldArrayI ciso;
  if (strcmp(eltType0, "TRI") == 0)
  {
    doIsoLine(*f, *cn, posf, value, poscellN, fiso, ciso);
    RELEASESHAREDU(surf, f, cn);
    E_Float tolc = 1.e-12; 
    K_CONNECT::cleanConnectivity(posx, posy, posz, tolc, 
                                 "BAR", fiso, ciso);

  }
  else
  {
    doIsoNode(*f, *cn, posf, value, poscellN, fiso, ciso);
    RELEASESHAREDU(surf, f, cn);
  }
  
  if (fiso.getSize() == 0 || ciso.getSize() == 0)
  {
    PyErr_SetString(PyExc_ValueError,
                    "isoLine: isoline is empty.");
    return NULL;
  }

  PyObject* t;
  E_Int api = 1;  // f->getApi();
  if (strcmp(eltType0, "TRI") == 0)
    t = K_ARRAY::buildArray3(fiso, varString0, ciso, "BAR", api);
  else
    t = K_ARRAY::buildArray3(fiso, varString0, ciso, "NODE", api);
  return t;
}

//=============================================================================
/* IN: f, cn: surface
   IN: posf: position de la variable iso
   IN: value: valeur de l'iso
   IN: poscellN: position de cellN (eventuellement, sinon=0)
   OUT: fiso, ciso: l'iso BAR.
*/
//=============================================================================
void K_POST::doIsoLine(FldArrayF& f, FldArrayI& cn, E_Int posf, E_Float value,
                       E_Int poscellN,
                       FldArrayF& fiso, FldArrayI& ciso)
{
  E_Int nelts = cn.getSize();
  E_Int nfld = f.getNfld();
  E_Float eps = 1.e-10;
  E_Float valuep = value + eps;
  E_Float valuem = value - eps;
  E_Float sign1, sign2, sign3, f1, f2, f3;
  E_Float amb1, amb2, amb3;
  E_Int ind1, ind2, ind3;

  fiso.malloc(2*nelts, nfld);
  ciso.malloc(nelts, 2);
  E_Int* ciso1 = ciso.begin(1); 
  E_Int* ciso2 = ciso.begin(2);
  E_Int nbars = 0; // nombre de bars dans l'iso
  E_Int npts = 0; // nombre de pts dans l'iso
  E_Float* fp = f.begin(posf);
  E_Int* cn1 = cn.begin(1);
  E_Int* cn2 = cn.begin(2);
  E_Int* cn3 = cn.begin(3);

  for (E_Int i = 0; i < nelts; i++)
  {
    ind1 = cn1[i]-1; ind2 = cn2[i]-1; ind3 = cn3[i]-1;
    f1 = fp[ind1]; f2 = fp[ind2]; f3 = fp[ind3];
    if (f1 > valuep && f2 > valuep && f3 > valuep) goto next;
    if (f1 < valuem && f2 < valuem && f3 < valuem) goto next;
    
    sign1 = K_FUNC::E_sign(f1-value);
    sign2 = K_FUNC::E_sign(f2-value);
    sign3 = K_FUNC::E_sign(f3-value);

    amb1 = K_FUNC::E_abs(f1-value);
    amb2 = K_FUNC::E_abs(f2-value);
    amb3 = K_FUNC::E_abs(f3-value);

    if (amb1 < eps && amb2 < eps && amb3 < eps) goto next;
    if (amb1 < eps) sign1 = 0.;
    if (amb2 < eps) sign2 = 0.;
    if (amb3 < eps) sign3 = 0.;

    // P12 et P13
    if (sign1 != sign2 && sign1 != sign3)
    {
      vertexInterp(nfld, value, f, poscellN,
                   f1, f2, ind1, ind2,
                   fiso, npts);
      vertexInterp(nfld, value, f, poscellN,
                   f1, f3, ind1, ind3,
                   fiso, npts);
      ciso1[nbars] = npts-1;
      ciso2[nbars] = npts;
      nbars++;
    }
    else if (sign1 != sign2 && sign2 != sign3)
    {
      vertexInterp(nfld, value, f, poscellN,
                   f1, f2, ind1, ind2,
                   fiso, npts);
      vertexInterp(nfld, value, f, poscellN,
                   f2, f3, ind2, ind3,
                   fiso, npts);
      ciso1[nbars] = npts-1;
      ciso2[nbars] = npts;
      nbars++;
    }
    else if (sign1 != sign3 && sign2 != sign3)
    {
      vertexInterp(nfld, value, f, poscellN,
                   f1, f3, ind1, ind3,
                   fiso, npts);
      vertexInterp(nfld, value, f, poscellN,
                   f2, f3, ind2, ind3,
                   fiso, npts);
      ciso1[nbars] = npts-1;
      ciso2[nbars] = npts;
      nbars++;
    }
    
    next: ;
  }

  fiso.reAllocMat(npts, nfld);
  ciso.reAllocMat(nbars, 2);
}

//=============================================================================
/* IN: f, cn: BAR
   IN: posf: position de la variable iso
   IN: value: valeur de l'iso
   IN: poscellN: position de cellN (eventuellement, sinon=0)
   OUT: fiso, ciso: l'iso NODE.
*/
//=============================================================================
void K_POST::doIsoNode(FldArrayF& f, FldArrayI& cn, E_Int posf, E_Float value,
                       E_Int poscellN,
                       FldArrayF& fiso, FldArrayI& ciso)
{
  E_Int nelts = cn.getSize();
  E_Int nfld = f.getNfld();
  E_Float eps = 1.e-10;
  E_Float valuep = value + eps;
  E_Float valuem = value - eps;
  E_Float sign1, sign2, f1, f2;
  E_Float amb1, amb2;
  E_Int ind1, ind2;

  fiso.malloc(nelts, nfld);
  ciso.malloc(nelts, 1);
  E_Int* ciso1 = ciso.begin(1); 
  E_Int nnodes = 0; // nombre de nodes dans l'iso
  E_Int npts = 0; // nombre de pts dans l'iso
  E_Float* fp = f.begin(posf);
  E_Int* cn1 = cn.begin(1);
  E_Int* cn2 = cn.begin(2);

  for (E_Int i = 0; i < nelts; i++)
  {
    ind1 = cn1[i]-1; ind2 = cn2[i]-1;
    f1 = fp[ind1]; f2 = fp[ind2];
    if (f1 > valuep && f2 > valuep) goto next;
    if (f1 < valuem && f2 < valuem) goto next;
    
    sign1 = K_FUNC::E_sign(f1-value);
    sign2 = K_FUNC::E_sign(f2-value);

    amb1 = K_FUNC::E_abs(f1-value);
    amb2 = K_FUNC::E_abs(f2-value);
    
    if (amb1 < eps && amb2 < eps) goto next;
    if (amb1 < eps) sign1 = 0.;
    if (amb2 < eps) sign2 = 0.;

    if (sign1 != sign2)
    {
      vertexInterp(nfld, value, f, poscellN,
                   f1, f2, ind1, ind2,
                   fiso, npts);
      ciso1[nnodes] = npts;
      nnodes++;
    }    
    next: ;
  }

  fiso.reAllocMat(npts, nfld);
  ciso.reAllocMat(nnodes, 1);
}
