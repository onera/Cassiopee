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

// grow an array of one layer

# include "generator.h"
using namespace K_FLD;

//=============================================================================
// Grow an array of one layer
//=============================================================================
PyObject* K_GENERATOR::growMesh(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* vect;

  if (!PYPARSETUPLE_(args, OO_, &array, &vect)) return NULL;

  // Check arrays
  E_Int ni, nj, nk, niv, njv, nkv;
  FldArrayF* f; FldArrayF* fv;
  char* varString; char* eltType;
  FldArrayI* cn;
  char* varStringv; char* eltTypev;
  FldArrayI* cnv;
  
  // Extraction des infos sur le maillage
  E_Int res =  K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, 
                                      cn, eltType);

  E_Int api = f->getApi();

  // Check data
  if (res == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "grow: array is invalid.");
    return NULL;
  }
  // Extraction des infos sur le vecteur
  E_Int resv = K_ARRAY::getFromArray3(vect, varStringv, fv, niv, njv, nkv, 
                                      cnv, eltTypev);
  if (resv == -1)
  {
    RELEASESHAREDB(res, array, f,cn);
    PyErr_SetString(PyExc_TypeError,
                    "grow: vector array is invalid.");
    return NULL;
  }
  if (res == 2)
  {
    if (strcmp(eltType, "TRI") != 0 && strcmp(eltType, "QUAD") != 0 
        && strcmp(eltType, "BAR") != 0)
    {
      RELEASESHAREDU(array, f, cn); RELEASESHAREDB(resv, vect, fv, cnv);
      PyErr_SetString(PyExc_TypeError,
                      "grow: array must be a surface array.");
      return NULL;
    }
  }
  if (res == 1)
  {
    if (ni != 1 && nj != 1 && nk != 1)
    {
      RELEASESHAREDS(array, f);  RELEASESHAREDB(resv, vect, fv, cnv);

      PyErr_SetString(PyExc_TypeError,
                      "grow: array must be a surface array.");
      return NULL;
    }
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f,cn); RELEASESHAREDB(resv, vect, fv, cnv);
    PyErr_SetString(PyExc_TypeError,
                    "grow: coordinates must be present in array.");
    return NULL;
  }
  posx++; posy++; posz++;
  
  if (f->getSize() != fv->getSize())
  {
    RELEASESHAREDB(res, array, f,cn); RELEASESHAREDB(resv, vect, fv, cnv);
    PyErr_SetString(PyExc_TypeError,
                    "grow: array and vector must have the same size.");
    return NULL;
  }
  if (fv->getNfld() != 3)
  {
    RELEASESHAREDB(res, array, f,cn); RELEASESHAREDB(resv, vect, fv, cnv);
    PyErr_SetString(PyExc_TypeError,
                    "grow: vector must have 3 variables.");
    return NULL;
  } 

  if (res == 1)
  {
    E_Int nic, njc, nkc;
    E_Int nfld = f->getNfld();
    E_Int dir = 1;
    if (ni == 1) dir = 1;
    else if (nj == 1) dir = 2;
    else if (nk == 1) dir = 3;
    
    E_Int step, ind, indp;
    switch (dir)
    {
      case 1:
        nic = 2; njc = nj; nkc = nk;
        step = 1;
        break;

      case 2:
        nic = ni; njc = 2; nkc = nk;
        step = ni;
        break;

      case 3:
        nic = ni; njc = nj; nkc = 2;
        step = ni*nj; 
        break;
    }
    PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, nic, njc, nkc, api);
    FldArrayF* coord;
    K_ARRAY::getFromArray3(tpl, coord);

    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      E_Float* coordp0 = coord->begin(nv);
      E_Float* fp = f->begin(nv);
      
      for (E_Int k = 0; k < nk; k++)
        for (E_Int j = 0; j < nj; j++)
          for (E_Int i = 0; i < ni; i++)
          {
            ind = i + j*ni + k*ni*nj;
            indp = i + j*nic + k*nic*njc;
            coordp0[indp] = fp[ind];
            coordp0[indp+step] = fp[ind];
          }
    }

    E_Float* coordx = coord->begin(posx);
    E_Float* coordy = coord->begin(posy);
    E_Float* coordz = coord->begin(posz);
    E_Float* fx = f->begin(posx);
    E_Float* fy = f->begin(posy);
    E_Float* fz = f->begin(posz);
    E_Float* fvx = fv->begin(1);
    E_Float* fvy = fv->begin(2);
    E_Float* fvz = fv->begin(3);
  
    for (E_Int k = 0; k < nk; k++)
      for (E_Int j = 0; j < nj; j++)
        for (E_Int i = 0; i < ni; i++)
        {
          ind = i + j*ni + k*ni*nj;
          indp = i + j*nic + k*nic*njc;        
          coordx[indp+step] = fx[ind] + fvx[ind];
          coordy[indp+step] = fy[ind] + fvy[ind];
          coordz[indp+step] = fz[ind] + fvz[ind];
        }
 
    RELEASESHAREDS(tpl, coord);
    RELEASESHAREDS(array, f); RELEASESHAREDB(resv, vect, fv, cnv);
    return tpl;
  }
  else
  {
    E_Int nt = f->getSize(); E_Int nfld = f->getNfld();
    E_Int nelts = cn->getSize(); E_Int nvert = cn->getNfld()*2;
    E_Int npts = 2*nt;// nbr of pts in the new zone

    char* eltType2 = new char[K_ARRAY::VARSTRINGLENGTH];
    eltType2[0] = '\0';
    if (nvert == 6) strcpy(eltType2, "PENTA");
    else if (nvert == 8) strcpy(eltType2, "HEXA");
    else strcpy(eltType2, "QUAD");

    PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts, nelts, eltType2, 0, api);
    FldArrayF* coord; FldArrayI* cn2;
    K_ARRAY::getFromArray3(tpl, coord, cn2);
    K_FLD::FldArrayI& cm2 = *(cn2->getConnect(0));

    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      E_Float* coordp0 = coord->begin(nv);
      E_Float* fp = f->begin(nv);
      for (E_Int ind = 0; ind < nt; ind++)
      {
        coordp0[ind] = fp[ind];
        coordp0[ind+nt] = fp[ind];
      }
    }
    E_Float* coordx = coord->begin(posx);
    E_Float* coordy = coord->begin(posy);
    E_Float* coordz = coord->begin(posz);
    E_Float* fx = f->begin(posx);
    E_Float* fy = f->begin(posy);
    E_Float* fz = f->begin(posz);
    E_Float* fvx = fv->begin(1);
    E_Float* fvy = fv->begin(2);
    E_Float* fvz = fv->begin(3);
    for (E_Int ind = 0; ind < nt; ind++)
    {
      coordx[ind+nt] = fx[ind] + fvx[ind];
      coordy[ind+nt] = fy[ind] + fvy[ind];
      coordz[ind+nt] = fz[ind] + fvz[ind];
    }

    // Connectivite
    if (nvert == 6) // PENTA
    {
      E_Int* cn10 = cn->begin(1);
      E_Int* cn20 = cn->begin(2);
      E_Int* cn30 = cn->begin(3);

      for (E_Int ind = 0; ind < nelts; ind++)
      {
        cm2(ind, 1) = cn10[ind];
        cm2(ind, 2) = cn20[ind];
        cm2(ind, 3) = cn30[ind];
        cm2(ind, 4) = cn10[ind] + nt;
        cm2(ind, 5) = cn20[ind] + nt;
        cm2(ind, 6) = cn30[ind] + nt;
      }
    }
    else if (nvert == 8) // HEXA
    {
      E_Int* cn10 = cn->begin(1);
      E_Int* cn20 = cn->begin(2);
      E_Int* cn30 = cn->begin(3);
      E_Int* cn40 = cn->begin(4);

      for (E_Int ind = 0; ind < nelts; ind++)
      {
        cm2(ind, 1) = cn10[ind];
        cm2(ind, 2) = cn20[ind];
        cm2(ind, 3) = cn30[ind];
        cm2(ind, 4) = cn40[ind];
        cm2(ind, 5) = cn10[ind] + nt;
        cm2(ind, 6) = cn20[ind] + nt;
        cm2(ind, 7) = cn30[ind] + nt;
        cm2(ind, 8) = cn40[ind] + nt;
      }
    }
    else if (nvert == 4) // QUAD
    {
      E_Int* cn10 = cn->begin(1);
      E_Int* cn20 = cn->begin(2);
      for (E_Int ind = 0; ind < nelts; ind++)
      {
        cm2(ind, 1) = cn10[ind];
        cm2(ind, 2) = cn20[ind];
        cm2(ind, 3) = cn10[ind] + nt;
        cm2(ind, 4) = cn20[ind] + nt;
      }
    }

    delete[] eltType2;
    RELEASESHAREDU(tpl, coord, cn2);
    RELEASESHAREDB(res, array, f, cn); RELEASESHAREDB(resv, vect, fv, cnv);
    return tpl;
  }
}
