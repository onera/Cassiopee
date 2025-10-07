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

// Overlay two fields

# include "initiator.h"

using namespace K_FLD;
using namespace std;

// ============================================================================
/* Field overlay */
// ============================================================================
PyObject* K_INITIATOR::overlayField(PyObject* self, PyObject* args)
{
  PyObject* array1;
  PyObject* array2;
  E_Float MInf;

  if (!PYPARSETUPLE_(args, OO_ R_, &array1, &array2, &MInf)) return NULL;

  // Check array
  E_Int im1, jm1, km1, im2, jm2, km2;
  FldArrayF* f1; FldArrayF* f2;
  FldArrayI* cn1; FldArrayI* cn2;
  char *varString1; char *varString2;
  char* eltType1; char* eltType2;

  E_Int res1 = K_ARRAY::getFromArray3(array1, varString1, f1, im1, jm1, km1,
                                      cn1, eltType1);
  E_Int res2 = K_ARRAY::getFromArray3(array2, varString2, f2, im2, jm2, km2,
                                      cn2, eltType2);

  E_Int res0;
  E_Int posr1, posru1, posrv1, posrw1, posre1;
  E_Int posr2, posru2, posrv2, posrw2, posre2;
  vector<E_Int> pos1; vector<E_Int> pos2;

  FldArrayF cfdField;

  if ((res1 == 1 && res2 == 1) || (res1 == 2 && res2 == 2))
  {
    // check if there are at least 5 variables in cfd field
    if (f1->getNfld() < 8 || f2->getNfld() < 8)
    {
      RELEASESHAREDB(res1, array1, f1, cn1);
      RELEASESHAREDB(res2, array2, f2, cn2);
      PyErr_SetString(
        PyExc_ValueError,
        "overlayField: number of variables is too small.");
      return NULL;
    }
    char varString[K_ARRAY::VARSTRINGLENGTH];
    res0 = K_ARRAY::getPosition(varString1, varString2, pos1,  pos2, varString);
    strcpy(varString, "x,y,z,ro,rou,rov,row,roE");
    if (res0 == -1)
    {
      RELEASESHAREDB(res1, array1, f1, cn1);
      RELEASESHAREDB(res2, array2, f2, cn2);
      PyErr_SetString(PyExc_ValueError,
                      "overlayField: no common variables.");
      return NULL;
    }
    else if (res0 == 0)
    {
      printf("Warning: overlayField: some variables are different...\n");
      printf(" Only common variables are kept.\n");
    }
    // Search for location of conservative variables in f1
    posr1 = K_ARRAY::isDensityPresent(varString1);
    posru1 = K_ARRAY::isMomentumXPresent(varString1);
    posrv1 = K_ARRAY::isMomentumYPresent(varString1);
    posrw1 = K_ARRAY::isMomentumZPresent(varString1);
    posre1 = K_ARRAY::isEnergyStagnationDensityPresent(varString1);

    if (posr1 == -1 || posru1 == -1 || posrv1 == -1 || posrw1 == -1 || posre1 == -1)
    {
      RELEASESHAREDB(res1, array1, f1, cn1);
      RELEASESHAREDB(res2, array2, f2, cn2);
      PyErr_SetString(PyExc_ValueError,
                      "overlayField: conservative variables not found in array1");
      return NULL;
    }

    // Search for location of conservative variables in f2
    posr2 = K_ARRAY::isDensityPresent(varString2);
    posru2 = K_ARRAY::isMomentumXPresent(varString2);
    posrv2 = K_ARRAY::isMomentumYPresent(varString2);
    posrw2 = K_ARRAY::isMomentumZPresent(varString2);
    posre2 = K_ARRAY::isEnergyStagnationDensityPresent(varString2);

    if (posr2 == -1 || posru2 == -1 || posrv2 == -1 || posrw2 == -1 || posre2 == -1)
    {
      RELEASESHAREDB(res1, array1, f1, cn1);
      RELEASESHAREDB(res2, array2, f2, cn2);
      PyErr_SetString(PyExc_ValueError,
                      "overlayField: conservative variables not found in array2.");
      return NULL;
    }

    posr1++; posru1++; posrv1++; posrw1++; posre1++;
    posr2++; posru2++; posrv2++; posrw2++; posre2++;

    if (f1->getSize() != f2->getSize())
    {
      RELEASESHAREDB(res1, array1, f1, cn1);
      RELEASESHAREDB(res2, array2, f2, cn2);
      PyErr_SetString(PyExc_ValueError,
                      "overlayField: the two meshes are different.");
      return NULL;
    }
    E_Int api = f1->getApi();

    E_Float ro1, ro2, u1, u2, v1, v2, w1, w2, roe1, roe2, p1, p2;
    E_Float ro, p, u, v, w;
    E_Float roinf, pinf, uinf;
    E_Int ncells1 = f1->getSize();

    // Adimensionnement infini
    roinf = 1.;
    if (MInf != 0)
    {
      uinf = 1.;
      pinf = 1./(1.4*MInf*MInf);
    }
    else
    {
      uinf = 0.;
      pinf = 1./1.4;
    }
    cfdField.malloc(ncells1, 8);

    E_Float* fr1 = f1->begin(posr1); E_Float* fr2 = f2->begin(posr2);
    E_Float* fru1 = f1->begin(posru1); E_Float* fru2 = f2->begin(posru2);
    E_Float* frv1 = f1->begin(posrv1); E_Float* frv2 = f2->begin(posrv2);
    E_Float* frw1 = f1->begin(posrw1); E_Float* frw2 = f2->begin(posrw2);
    E_Float* fre1 = f1->begin(posre1); E_Float* fre2 = f2->begin(posre2);
    E_Float invr1, invr2;

    E_Float* xt1 = f1->begin(pos1[0]);
    E_Float* yt1 = f1->begin(pos1[1]);
    E_Float* zt1 = f1->begin(pos1[2]);
    E_Float invgam1 = 1/0.4;

    E_Float* xt = cfdField.begin(1);
    E_Float* yt = cfdField.begin(2);
    E_Float* zt = cfdField.begin(3);
    E_Float* rot = cfdField.begin(4); E_Float* rout = cfdField.begin(5);
    E_Float* rovt = cfdField.begin(6); E_Float* rowt = cfdField.begin(7);
    E_Float* roet = cfdField.begin(8);

    for (E_Int i = 0; i < ncells1; i++)
    {
      ro1 = fr1[i];
      invr1 = 1./ro1;
      u1 = fru1[i] * invr1;
      v1 = frv1[i] * invr1;
      w1 = frw1[i] * invr1;
      roe1 = fre1[i]-0.5*ro1*(u1*u1+v1*v1+w1*w1);
      p1 = 0.4*roe1;

      ro2= fr2[i];
      invr2 = 1./ro2;
      u2 = fru2[i] * invr2;
      v2 = frv2[i] * invr2;
      w2 = frw2[i] * invr2;
      roe2 = fre2[i]-0.5*ro2*(u2*u2+v2*v2+w2*w2);
      p2 = 0.4*roe2;

      ro = ro1+ro2-roinf;
      p = p1+p2-pinf;
      u = u1+u2-uinf; // on suppose 0 degre d'incidence
      v = v1+v2;
      w = w1+w2;

      xt[i] = xt1[i];
      yt[i] = yt1[i];
      zt[i] = zt1[i];
      rot[i] = ro;
      rout[i] = ro*u;
      rovt[i] = ro*v;
      rowt[i] = ro*w;
      roet[i] = p*invgam1  + 0.5 * ro * ( u*u + v*v + w*w );
    }

    // build array
    PyObject* tpl;
    if (res1 == 1)
    {
      tpl = K_ARRAY::buildArray3(cfdField, varString, im2, jm2, km2, api);
    }
    else
    {
      tpl = K_ARRAY::buildArray3(cfdField, varString, *cn2, eltType2, api);
    }
    RELEASESHAREDB(res1, array1, f1, cn1);
    RELEASESHAREDB(res2, array2, f2, cn2);
    return tpl;
  }
  else if (res1 == 2 || res2 == 2)
  {
    RELEASESHAREDB(res1, array1, f1, cn1);
    RELEASESHAREDB(res2, array2, f2, cn2);
    PyErr_SetString(PyExc_TypeError,
                    "overlayField: not used for unstructured arrays.");
    return NULL;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "overlayField: one or both arrays are invalid.");
    return NULL;
  }
}
