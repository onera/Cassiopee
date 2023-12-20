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

// Initialization by vortices

# include "initiator.h"

using namespace K_FLD;
using namespace std;

extern "C"
{
  void k6lamb_(const E_Float& x0, const E_Float& y0,
               const E_Float& Gamma, const E_Float& MInf,
               const E_Int& npts,
               const E_Float* xc, const E_Float* yc, const E_Float* zc,
               E_Float* ro, E_Float* rou, E_Float* rov, E_Float* row, E_Float* roe);

  void k6scully_(const E_Float& x0, const E_Float& y0,
                 const E_Float& Gamma, const E_Float& a, const E_Float& MInf,
                 const E_Int& npts,
                 const E_Float* xc, const E_Float* yc, const E_Float* zc,
                 E_Float* field);

  void k6visbal_(const E_Float& x0, const E_Float& y0,
                 const E_Float& Gamma, const E_Float& MInf,
                 const E_Int& npts,
                 const E_Float* xc, const E_Float* yc, const E_Float* zc,
                 E_Float* field);

  void k6yee_(const E_Float& x0, const E_Float& y0,
              const E_Float& Gamma, const E_Float& Minf,
              const E_Int& npts,
              const E_Float* xc, const E_Float* yc, const E_Float* zc,
              E_Float* field);

  void k6scully2_(const E_Float& x0, const E_Float& y0,
                  const E_Float& Gamma, const E_Float& a, const E_Float& MInf,
                  const E_Int& npts,
                  const E_Float* xc, const E_Float* yc, const E_Float* zc,
                  E_Float* field);

  void k6wissocq_(const E_Float& x0, const E_Float& y0,
                  const E_Float& Gamma, const E_Float& MInf,
                  const E_Int& npts,
                  const E_Float* xc, const E_Float* yc, const E_Float* zc,
                  E_Float* field);
}
// ============================================================================
/* Init by a lamb vortex */
// ============================================================================
PyObject* K_INITIATOR::initLamb(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float x0, y0, Gam, MInf;
  
  if (!PYPARSETUPLE_(args, O_ TRR_ RR_, &array, &x0, &y0, &Gam, &MInf))
    return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);

  if (res != 1 && res !=2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "initLamb: invalid type of array.");
    return NULL;
  }
  char varStringOut[K_ARRAY::VARSTRINGLENGTH];

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_ValueError,
                    "initLamb: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  E_Int defVars = 0; //says if output array is built using default variables or not
  E_Int neqTot = f->getNfld();
  E_Int npts = f->getSize();
  E_Int posro = K_ARRAY::isDensityPresent(varString);
  E_Int posrou = K_ARRAY::isMomentumXPresent(varString);
  E_Int posrov = K_ARRAY::isMomentumYPresent(varString);
  E_Int posrow = K_ARRAY::isMomentumZPresent(varString);
  E_Int posroe = K_ARRAY::isEnergyStagnationDensityPresent(varString);
  if (posro == -1 || posrou == -1 || posrov == -1 || posrow == -1 || posroe == -1)
  {
    defVars = 1;
    neqTot = 8;
    strcpy(varStringOut, "x,y,z,ro,rou,rov,row,roE");
    posro = 4; posrou = 5; posrov = 6; posrow = 7; posroe = 8;
  }
  else
  {
    posro++; posrou++; posrov++; posrow++; posroe++;
    strcpy(varStringOut, varString);
  }

  PyObject* tpl;
  if (res == 1)
    tpl = K_ARRAY::buildArray3(neqTot, varStringOut, im, jm, km, f->getApi());
  else
    tpl = K_ARRAY::buildArray3(neqTot, varStringOut, *f, *cn, eltType, 0, -1, true);

  FldArrayF* f2;  
  K_ARRAY::getFromArray3(tpl, f2);
  
  // copy previous data
  E_Float* fp; E_Float* f2p;
  if (defVars == 1)
  {
    fp = f->begin(posx); f2p = f2->begin(1);
    for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i]; 
    fp = f->begin(posy); f2p = f2->begin(2);
    for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i]; 
    fp = f->begin(posz); f2p = f2->begin(3);
    for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
  }
  else
  {
    for (E_Int n = 1; n <= neqTot; n++)
    {
      fp = f->begin(n); f2p = f2->begin(n);
      for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
    }
  }
  // init with lamb
  k6lamb_(x0, y0, Gam, MInf, npts,
          f->begin(posx), f->begin(posy), f->begin(posz),
          f2->begin(posro), f2->begin(posrou), f2->begin(posrov), 
          f2->begin(posrow), f2->begin(posroe));

  RELEASESHAREDB(res, array, f, cn);
  RELEASESHAREDS(tpl, f2);
  return tpl;
}

// ===========================================================================
/* Init by a Visbal vortex */
// ============================================================================
PyObject* K_INITIATOR::initVisbal(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float x0, y0, Gam, MInf;
  
  if (!PYPARSETUPLE_(args, O_ TRR_ RR_, &array, &x0, &y0, &Gam, &MInf))
    return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  FldArrayF cfdField, fieldTot;
  E_Int res =
    K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType);
  char varStringOut[K_ARRAY::VARSTRINGLENGTH];

  if (res == 1 || res == 2)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      delete f;
      PyErr_SetString(PyExc_ValueError,
                      "initVisbal: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;

    E_Int defVars = 0;//says if output array is built using default variables or not
    //E_Int nCfdFld;
    E_Int neqTot = f->getNfld();

    E_Int posro = K_ARRAY::isDensityPresent(varString);
    E_Int posrou = K_ARRAY::isMomentumXPresent(varString);
    E_Int posrov = K_ARRAY::isMomentumYPresent(varString);
    E_Int posrow = K_ARRAY::isMomentumZPresent(varString);
    E_Int posroe = K_ARRAY::isEnergyStagnationDensityPresent(varString);
    vector<E_Int> posTot;

    if (posro == -1 || posrou == -1 || posrov == -1 || posrow == -1 || posroe == -1)
    {
      defVars = 1;
      //nCfdFld = 5;
      neqTot = 8;
      strcpy(varStringOut, "x,y,z,ro,rou,rov,row,roE");
    }
    else
    {
      posro++; posrou++; posrov++; posrow++; posroe++;
      posTot.push_back(posro);
      posTot.push_back(posrou);
      posTot.push_back(posrov);
      posTot.push_back(posrow);
      posTot.push_back(posroe);
      //nCfdFld = neqTot-3;
      strcpy(varStringOut, varString);
    }
    E_Int npts = f->getSize();

    FldArrayF cfdField(npts, 5);
    cfdField.setAllValuesAtNull();

    k6visbal_(x0, y0, Gam, MInf, npts,
              f->begin(posx), f->begin(posy), f->begin(posz),
              cfdField.begin());

    fieldTot.malloc(npts, neqTot);
    E_Float* xt = f->begin(posx);
    E_Float* yt = f->begin(posy);
    E_Float* zt = f->begin(posz);
    E_Float* xt2 = fieldTot.begin(1);
    E_Float* yt2 = fieldTot.begin(2);
    E_Float* zt2 = fieldTot.begin(3);
    E_Int posTotsize = posTot.size();

    switch (defVars)
    {
      case 0 :
        for (E_Int i = 0; i < npts; i++)
        {
          for (E_Int eq = 1; eq <= neqTot; eq++)
            fieldTot(i,eq) = (*f)(i,eq);

          for (E_Int c = 0; c < posTotsize; c++)
          {
            E_Int pos0 = posTot[c];
            fieldTot(i,pos0) = cfdField(i,c+1);
          }
        }
        break;

      default:
        for (E_Int i = 0; i < npts; i++)
        {
          xt2[i] = xt[i]; yt2[i] = yt[i]; zt2[i] = zt[i];

          for (E_Int eq = 1; eq <= 5; eq++)
            fieldTot(i,eq+3) = cfdField(i,eq);
        }
        break;
    }
    // build array
    delete f;
    FldArrayF* an = new FldArrayF(fieldTot);
    PyObject* tpl;
    if (res == 1)
    {
      tpl = K_ARRAY::buildArray(*an, varStringOut, im, jm, km);
      delete an;
    }
    else
    {
      tpl = K_ARRAY::buildArray(*an, varStringOut, *cn, -1, eltType);
      delete an; delete cn;
    }
    return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "initVisbal : invalid type of array.");
    return NULL;
  }
}

// ============================================================================
/* Init by a scully vortex */
// ============================================================================
PyObject* K_INITIATOR::initScully(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float x0, y0, Gam, coreRadius, MInf;
  E_Int model=0;
  
  if (!PYPARSETUPLE_(args, O_ TRR_ RRR_ I_,
                     &array, &x0, &y0, &Gam, &coreRadius, &MInf, &model))
    return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  FldArrayF cfdField, fieldTot;

  E_Int res =
    K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType);
  char varStringOut[K_ARRAY::VARSTRINGLENGTH];

  if (res == 1 || res == 2)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      delete f;
      PyErr_SetString(PyExc_ValueError,
                      "initScully: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;

    E_Int defVars = 0;//says if output array is built using default variables or not
    //E_Int nCfdFld;
    E_Int neqTot = f->getNfld();

    E_Int posro = K_ARRAY::isDensityPresent(varString);
    E_Int posrou = K_ARRAY::isMomentumXPresent(varString);
    E_Int posrov = K_ARRAY::isMomentumYPresent(varString);
    E_Int posrow = K_ARRAY::isMomentumZPresent(varString);
    E_Int posroe = K_ARRAY::isEnergyStagnationDensityPresent(varString);
    vector<E_Int> posTot;

    if (posro == -1 || posrou == -1 || posrov == -1 || posrow == -1 || posroe == -1)
    {
      defVars = 1;
      //nCfdFld = 5;
      neqTot = 8;
      strcpy(varStringOut, "x,y,z,ro,rou,rov,row,roE");
    }
    else
    {
      posro++; posrou++; posrov++; posrow++; posroe++;
      posTot.push_back(posro);
      posTot.push_back(posrou);
      posTot.push_back(posrov);
      posTot.push_back(posrow);
      posTot.push_back(posroe);
      //nCfdFld = neqTot-3;
      strcpy(varStringOut, varString);
    }
    E_Int npts = f->getSize();

    FldArrayF cfdField(npts, 5);
    cfdField.setAllValuesAtNull();

    switch (model)
    {
      case 0:
        k6scully_(x0, y0, Gam, coreRadius, MInf, npts,
                  f->begin(posx), f->begin(posy), f->begin(posz),
                  cfdField.begin());
      default :
        k6scully2_(x0, y0, Gam, coreRadius, MInf, npts,
                   f->begin(posx), f->begin(posy), f->begin(posz),
                   cfdField.begin());
    }

    fieldTot.malloc(npts, neqTot);
    E_Float* xt = f->begin(posx);
    E_Float* yt = f->begin(posy);
    E_Float* zt = f->begin(posz);
    E_Float* xt2 = fieldTot.begin(1);
    E_Float* yt2 = fieldTot.begin(2);
    E_Float* zt2 = fieldTot.begin(3);
    E_Int posTotsize = posTot.size();

    switch (defVars)
    {
      case 0 :
        for (E_Int i = 0; i < npts; i++)
        {
          for (E_Int eq = 1; eq <= neqTot; eq++)
            fieldTot(i,eq) = (*f)(i,eq);

          for (E_Int c = 0; c < posTotsize; c++)
          {
            E_Int pos0 = posTot[c];
            fieldTot(i,pos0) = cfdField(i,c+1);
          }
        }
        break;

      default:
        for (E_Int i = 0; i < npts; i++)
        {
          xt2[i] = xt[i]; yt2[i] = yt[i]; zt2[i] = zt[i];
          for (E_Int eq = 1; eq <= 5; eq++)
            fieldTot(i,eq+3) = cfdField(i,eq);
        }
        break;
    }
    // build array
    delete f;
    FldArrayF* an = new FldArrayF(fieldTot);
    PyObject* tpl;
    if (res == 1)
    {
      tpl = K_ARRAY::buildArray(*an, varStringOut, im, jm, km);
      delete an;
    }
    else
    {
      tpl = K_ARRAY::buildArray(*an, varStringOut, *cn, -1, eltType);
      delete an; delete cn;
    }
    return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "initScully: unknown type of array.");
    return NULL;
  }
}

// ===========================================================================
/* Init by a Yee vortex */
// ============================================================================
PyObject*
K_INITIATOR::initYee(PyObject* self,
                     PyObject* args)
{
  PyObject* array;
  E_Float x0, y0, Gam, Minf;
  
  if (!PYPARSETUPLE_(args, O_ TRR_ RR_, &array, &x0, &y0, &Gam, &Minf))
    return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f;
  char* varString; char* eltType;
  FldArrayI* cn;
  FldArrayF cfdField, fieldTot;
  E_Int res =
    K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType);
  char varStringOut[K_ARRAY::VARSTRINGLENGTH];

  if (res == 1 || res == 2)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      delete f;
      PyErr_SetString(PyExc_ValueError,
                      "initYee: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;

    E_Int defVars = 0;//says if output array is built using default variables or not
    //E_Int nCfdFld;
    E_Int neqTot = f->getNfld();

    E_Int posro = K_ARRAY::isDensityPresent(varString);
    E_Int posrou = K_ARRAY::isMomentumXPresent(varString);
    E_Int posrov = K_ARRAY::isMomentumYPresent(varString);
    E_Int posrow = K_ARRAY::isMomentumZPresent(varString);
    E_Int posroe = K_ARRAY::isEnergyStagnationDensityPresent(varString);
    vector<E_Int> posTot;

    if (posro == -1 || posrou == -1 || posrov == -1 || posrow == -1 || posroe == -1)
    {
      defVars = 1;
      //nCfdFld = 5;
      neqTot = 8;
      strcpy(varStringOut, "x,y,z,ro,rou,rov,row,roE");
    }
    else
    {
      posro++; posrou++; posrov++; posrow++; posroe++;
      posTot.push_back(posro);
      posTot.push_back(posrou);
      posTot.push_back(posrov);
      posTot.push_back(posrow);
      posTot.push_back(posroe);
      //nCfdFld = neqTot-3;
      strcpy(varStringOut, varString);
    }
    E_Int npts = f->getSize();

    FldArrayF cfdField(npts, 5);
    cfdField.setAllValuesAtNull();

    k6yee_(x0, y0, Gam, Minf, npts,
              f->begin(posx), f->begin(posy), f->begin(posz),
              cfdField.begin());

    fieldTot.malloc(npts, neqTot);
    E_Float* xt = f->begin(posx);
    E_Float* yt = f->begin(posy);
    E_Float* zt = f->begin(posz);
    E_Float* xt2 = fieldTot.begin(1);
    E_Float* yt2 = fieldTot.begin(2);
    E_Float* zt2 = fieldTot.begin(3);
    E_Int posTotsize = posTot.size();

    switch (defVars)
    {
      case 0 :
        for (E_Int i = 0; i < npts; i++)
        {
          for (E_Int eq = 1; eq <= neqTot; eq++)
            fieldTot(i,eq) = (*f)(i,eq);

          for (E_Int c = 0; c < posTotsize; c++)
          {
            E_Int pos0 = posTot[c];
            fieldTot(i,pos0) = cfdField(i,c+1);
          }
        }
        break;

      default:
        for (E_Int i = 0; i < npts; i++)
        {
          xt2[i] = xt[i]; yt2[i] = yt[i]; zt2[i] = zt[i];

          for (E_Int eq = 1; eq <= 5; eq++)
            fieldTot(i,eq+3) = cfdField(i,eq);
        }
        break;
    }
    // build array
    delete f;
    FldArrayF* an = new FldArrayF(fieldTot);
    PyObject* tpl;
    if (res == 1)
    {
      tpl = K_ARRAY::buildArray(*an, varStringOut, im, jm, km);
      delete an;
    }
    else
    {
      tpl = K_ARRAY::buildArray(*an, varStringOut, *cn, -1, eltType);
      delete an; delete cn;
    }
    return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "initYee : invalid type of array.");
    return NULL;
  }
}

// ============================================================================
/* Init by Wissocq's vortex */
// ============================================================================
PyObject* K_INITIATOR::initWissocq(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float x0, y0, Gam, MInf;
  
  if (!PYPARSETUPLE_(args, O_ TRR_ RR_, &array, &x0, &y0, &Gam, &MInf))
    return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  FldArrayF cfdField, fieldTot;
  E_Int res =
    K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType);

  char varStringOut[K_ARRAY::VARSTRINGLENGTH];

  if (res == 1 || res == 2)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      delete f;
      PyErr_SetString(PyExc_ValueError,
                      "initWissocq: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;

    E_Int defVars = 0; //says if output array is built using default variables or not
    //E_Int nCfdFld;
    E_Int neqTot = f->getNfld();
    E_Int posro = K_ARRAY::isDensityPresent(varString);
    E_Int posrou = K_ARRAY::isMomentumXPresent(varString);
    E_Int posrov = K_ARRAY::isMomentumYPresent(varString);
    E_Int posrow = K_ARRAY::isMomentumZPresent(varString);
    E_Int posroe = K_ARRAY::isEnergyStagnationDensityPresent(varString);
    vector<E_Int> posTot;
    if (posro == -1 || posrou == -1 || posrov == -1 || posrow == -1 || posroe == -1)
    {
      defVars = 1;
      //nCfdFld = 5;
      neqTot = 8;
      strcpy(varStringOut, "x,y,z,ro,rou,rov,row,roE");
    }
    else
    {
      posro++; posrou++; posrov++; posrow++; posroe++;
      posTot.push_back(posro);
      posTot.push_back(posrou);
      posTot.push_back(posrov);
      posTot.push_back(posrow);
      posTot.push_back(posroe);
      //nCfdFld = neqTot-3;
      strcpy(varStringOut, varString);
    }
    E_Int npts = f->getSize();

    FldArrayF cfdField(npts, 5);
    cfdField.setAllValuesAtNull();

    k6wissocq_(x0, y0, Gam, MInf, npts,
               f->begin(posx), f->begin(posy), f->begin(posz),
               cfdField.begin());

    fieldTot.malloc(npts, neqTot);
    E_Float* xt = f->begin(posx);
    E_Float* yt = f->begin(posy);
    E_Float* zt = f->begin(posz);
    E_Float* xt2 = fieldTot.begin(1);
    E_Float* yt2 = fieldTot.begin(2);
    E_Float* zt2 = fieldTot.begin(3);
    E_Int posTotsize = posTot.size();
    switch (defVars)
    {
      case 0 :
        for (E_Int i = 0; i < npts; i++)
        {
          for (E_Int eq = 1; eq <= neqTot; eq++)
            fieldTot(i,eq) = (*f)(i,eq);

          for (E_Int c = 0; c < posTotsize; c++)
          {
            E_Int pos0 = posTot[c];
            fieldTot(i,pos0) = cfdField(i,c+1);
          }
        }
        break;

      default:
        for (E_Int i = 0; i < npts; i++)
        {
          xt2[i] = xt[i]; yt2[i] = yt[i]; zt2[i] = zt[i];
          for (E_Int eq = 1; eq <= 5; eq++)
            fieldTot(i,eq+3) = cfdField(i,eq);
        }
        break;
    }

    // build array
    delete f;
    FldArrayF* an = new FldArrayF(fieldTot);
    PyObject* tpl;
    if (res == 1)
    {
      tpl = K_ARRAY::buildArray(*an, varStringOut, im, jm, km);
      delete an;
    }
    else
    {
      tpl = K_ARRAY::buildArray(*an, varStringOut, *cn, -1, eltType);
      delete an; delete cn;
    }
    return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "initWissocq: invalid type of array.");
    return NULL;
  }
}
