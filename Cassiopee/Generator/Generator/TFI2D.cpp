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

// TFI generator - 2D struct

# include "generator.h"

using namespace std;
using namespace K_FLD;

//===========================================================================
/* TFI 2D : structure 2D */
//===========================================================================
PyObject* K_GENERATOR::TFI2D(PyObject* arrays)
{
  // Extract infos from arrays
  vector<E_Int> res;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> fields; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt; vector<char*> eltType;
  vector<PyObject*> objs, obju;
  E_Bool skipNoCoord = true;
  E_Bool skipStructured = false;
  E_Bool skipUnstructured = true; 
  E_Bool skipDiffVars = true;

  E_Int isOk = K_ARRAY::getFromArrays(
    arrays, res, structVarString, unstrVarString,
    fields, unstrF, nit, njt, nkt, cnt, eltType, objs, obju, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);

  E_Int nzones = fields.size();
  E_Int nfld = 0;
  if (nzones != 0) nfld = fields[0]->getNfld();
  if (isOk == -1 || nfld < 3)
  {
    PyErr_SetString(PyExc_TypeError,
                    "TFI: invalid list of arrays.");
    for (E_Int nos = 0; nos < nzones; nos++)
      RELEASESHAREDS(objs[nos], fields[nos]);
    return NULL;
  }
  
  // verification que les arrays sont bien 1D
  for (E_Int v = 0; v < nzones; v++)
  {
    E_Int ni = nit[v]; E_Int nj = njt[v]; E_Int nk = nkt[v];
    if (ni < 2 || nj != 1 || nk != 1)
    {       
      for (E_Int nos = 0; nos < nzones; nos++)
        RELEASESHAREDS(objs[nos], fields[nos]);
      PyErr_SetString(PyExc_TypeError,
                      "TFI: one array is not valid: must be i-varying only.");
      return NULL;
    }
  }
  // coordonnees deja verifiees dans getFromArrays
  char* varString = structVarString[0];
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString); 
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString); 
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  posx++; posy++; posz++;

  // reordonner les arrays le cas echeant
  FldArrayIS newOrder(4);
  E_Int ok = reorderTFI2D(posx, posy, posz, nit, fields, newOrder);
  if (ok == 0) 
  {
    for (E_Int nos = 0; nos < nzones; nos++)
      RELEASESHAREDS(objs[nos], fields[nos]);
    PyErr_SetString(PyExc_TypeError,
                    "TFI: input arrays must be C0.");
    return NULL;
  }
  E_Int imin = newOrder[0]; E_Int imax = newOrder[2];
  E_Int jmin = newOrder[3]; E_Int jmax = newOrder[1];
  E_Int nj = fields[imin]->getSize();
  E_Int ni = fields[jmin]->getSize();
  E_Int nk = 1;

  E_Int api = fields[0]->getApi();
  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, ni, nj, nk, api);
  FldArrayF* coord;
  K_ARRAY::getFromArray3(tpl, coord);

  //short isok = TFIstruct2D(ni, nj, nfld, imin, imax, jmin, jmax, fields, coord);
  short isok = TFIstruct2D2(ni, nj, nfld, posx, posy, posz, imin, imax, jmin, jmax, fields, *coord);

  if (isok < 1) //echec
  {
    for (E_Int nos = 0; nos < nzones; nos++)
      RELEASESHAREDS(objs[nos], fields[nos]);
    if ( isok == -1)
      PyErr_SetString(PyExc_TypeError,
                      "TFI: imin and imax borders are not of same size ni.");
    else if ( isok == -2 )
            PyErr_SetString(PyExc_TypeError,
                      "TFI: jmin and jmax borders are not of same size nj.");
    else       
      PyErr_SetString(PyExc_TypeError,
                      "TFI: input arrays are not valid.");
    return NULL;
  }

  RELEASESHAREDS(tpl, coord);
  for (E_Int nos = 0; nos < nzones; nos++)
    RELEASESHAREDS(objs[nos], fields[nos]);
  return tpl;
}
 
//===========================================================================
/* TFI 2D structure.  */
//===========================================================================
short K_GENERATOR::TFIstruct2D(E_Int ni, E_Int nj, E_Int nfld, 
                               E_Int imin, E_Int imax, E_Int jmin, E_Int jmax, 
                               vector<FldArrayF*>& fields,
                               FldArrayF& coords)
{ 
  E_Int ni1 = ni-1;
  E_Int nj1 = nj-1;
  E_Int i, j, ind;
  E_Float ip, jp, ip1, jp1;
  E_Float invni1 = 1./ni1;
  E_Float invnj1 = 1./nj1;
  for (E_Int eq = 1; eq <= nfld; eq++)
  {
    E_Float* ximin = fields[imin]->begin(eq);
    E_Float* ximax = fields[imax]->begin(eq);
    E_Float* xjmin = fields[jmin]->begin(eq);
    E_Float* xjmax = fields[jmax]->begin(eq);
    E_Float* xt = coords.begin(eq);
    for (j = 0; j < nj; j++)
      for (i = 0; i < ni; i++)
      {
        ind = i+ j*ni;
        ip = i*invni1; ip1 = 1.-ip;
        jp = j*invnj1; jp1 = 1.-jp;
       
        xt[ind] = jp * xjmax[i] + jp1 * xjmin[i] + 
          ip * (ximax[j] - jp*ximax[nj1] - jp1*ximax[0] ) + 
          ip1* (ximin[j] - jp*ximin[nj1] - jp1*ximin[0] );
      }
  }
  return 1;
}

//===========================================================================
/* TFI 2D structure avec ponderation  */
//===========================================================================
short K_GENERATOR::TFIstruct2D2(E_Int ni, E_Int nj, E_Int nfld,
                                E_Int posx, E_Int posy, E_Int posz,
                                E_Int imin, E_Int imax, E_Int jmin, E_Int jmax,
                                vector<FldArrayF*>& fields,
                                FldArrayF& coords)
{
  E_Int ni1 = ni-1;
  E_Int nj1 = nj-1;
  E_Int i, j, ind;

  if (fields[imin]->getSize() != fields[imax]->getSize()) return -1;
  if (fields[jmin]->getSize() != fields[jmax]->getSize()) return -2;

  E_Float* pondximin = fields[imin]->begin(posx);
  E_Float* pondximax = fields[imax]->begin(posx);
  E_Float* pondxjmin = fields[jmin]->begin(posx);
  E_Float* pondxjmax = fields[jmax]->begin(posx);
  E_Float* pondyimin = fields[imin]->begin(posy);
  E_Float* pondyimax = fields[imax]->begin(posy);
  E_Float* pondyjmin = fields[jmin]->begin(posy);
  E_Float* pondyjmax = fields[jmax]->begin(posy);
  E_Float* pondzimin = fields[imin]->begin(posz);
  E_Float* pondzimax = fields[imax]->begin(posz);
  E_Float* pondzjmin = fields[jmin]->begin(posz);
  E_Float* pondzjmax = fields[jmax]->begin(posz);

  E_Float pondximindist = pow( pow((pondximin[nj1] - pondximin[0]), 2.0) + pow((pondyimin[nj1] - pondyimin[0]), 2.0) + pow((pondzimin[nj1] - pondzimin[0]) , 2.0) , 0.5);
  E_Float pondximaxdist = pow( pow((pondximax[nj1] - pondximax[0]), 2.0) + pow((pondyimax[nj1] - pondyimax[0]), 2.0) + pow((pondzimax[nj1] - pondzimax[0]) , 2.0) , 0.5);
  E_Float pondxjmindist = pow( pow((pondxjmin[ni1] - pondxjmin[0]), 2.0) + pow((pondyjmin[ni1] - pondyjmin[0]), 2.0) + pow((pondzjmin[ni1] - pondzjmin[0]) , 2.0) , 0.5);
  E_Float pondxjmaxdist = pow( pow((pondxjmax[ni1] - pondxjmax[0]), 2.0) + pow((pondyjmax[ni1] - pondyjmax[0]), 2.0) + pow((pondzjmax[ni1] - pondzjmax[0]) , 2.0) , 0.5);

  E_Float invpondximindist = 1.0 / pondximindist;
  E_Float invpondximaxdist = 1.0 / pondximaxdist;
  E_Float invpondxjmindist = 1.0 / pondxjmindist;
  E_Float invpondxjmaxdist = 1.0 / pondxjmaxdist;

  for (E_Int eq = 1; eq <= nfld; eq++)
  {
    E_Float* ximin = fields[imin]->begin(eq);
    E_Float* ximax = fields[imax]->begin(eq);
    E_Float* xjmin = fields[jmin]->begin(eq);
    E_Float* xjmax = fields[jmax]->begin(eq);
    E_Float* xt = coords.begin(eq);
    for (j = 0; j < nj; j++)
      for (i = 0; i < ni; i++)
      {
        ind = i+ j*ni;
        E_Float pondximindistpartial = pow( pow((pondximin[j] - pondximin[0]), 2.0) + pow((pondyimin[j] - pondyimin[0]), 2.0) + pow((pondzimin[j] - pondzimin[0]), 2.0) , 0.5);
        E_Float pondximaxdistpartial = pow( pow((pondximax[j] - pondximax[0]), 2.0) + pow((pondyimax[j] - pondyimax[0]), 2.0) + pow((pondzimax[j] - pondzimax[0]), 2.0) , 0.5);
        E_Float pondxjmindistpartial = pow( pow((pondxjmin[i] - pondxjmin[0]), 2.0) + pow((pondyjmin[i] - pondyjmin[0]), 2.0) + pow((pondzjmin[i] - pondzjmin[0]), 2.0) , 0.5);
        E_Float pondxjmaxdistpartial = pow( pow((pondxjmax[i] - pondxjmax[0]), 2.0) + pow((pondyjmax[i] - pondyjmax[0]), 2.0) + pow((pondzjmax[i] - pondzjmax[0]), 2.0) , 0.5);

        E_Float pondximindistratio = pondximindistpartial * invpondximindist;
        E_Float pondximaxdistratio = pondximaxdistpartial * invpondximaxdist;
        E_Float pondxjmindistratio = pondxjmindistpartial * invpondxjmindist;
        E_Float pondxjmaxdistratio = pondxjmaxdistpartial * invpondxjmaxdist;

        E_Float pondxidistratio = (pondximindistratio + pondximaxdistratio)*0.5;
        E_Float pondxjdistratio = (pondxjmindistratio + pondxjmaxdistratio)*0.5;

        E_Float pondxidistratio1 = 1.0 - pondxidistratio;
        E_Float pondxjdistratio1 = 1.0 - pondxjdistratio;

        xt[ind] = pondxidistratio1 * xjmin[i] + pondxidistratio * xjmax[i]
                + pondxjdistratio1 * ximin[j] + pondxjdistratio * ximax[j]
                - pondxjdistratio1* pondxidistratio1 * ximin[0]
                - pondxjdistratio1* pondxidistratio* ximin[nj1]
                - pondxjdistratio * pondxidistratio1 * ximax[0]
                - pondxjdistratio * pondxidistratio * ximax[nj1];
      }
  }
  return 1;
}

//=========================================================================
/* Reordonne les arrays d'entree - cas 2D */
//=========================================================================
E_Int K_GENERATOR::reorderTFI2D(
  E_Int posx, E_Int posy, E_Int posz,
  vector<E_Int>& nit, vector<FldArrayF*>& fields, 
  FldArrayIS& newOrder)
{
  E_Float eps = 1.e-6;
  E_Int jm = 1; E_Int km = 1;
  FldArrayIS dejaVu(4); dejaVu.setAllValuesAtNull();

  E_Int ok = 0; E_Int i2 = -1; E_Int j2 = -1; E_Int k2 = -1;
  // 1- depart : f1 -> detection raccord en imax1
  E_Int vp = 0; E_Int vpp = 0;
  E_Int c = 1;
  while (c < 5)
  {
    i2 = -1; j2 = -1; k2 = -1;
    for (E_Int vn = 0; vn < 4; vn++)
    {
      if (dejaVu[vn] == 0 && vp != vn && vpp != vn) 
      {
        i2 = -1; j2 = -1; k2 = -1;       
        ok = K_CONNECT::getCoincidentVertexIndices(
          nit[vp], jm, km, nit[vp], jm, km, nit[vn], jm, km,
          posx, posy, posz, posx, posy, posz,
          *fields[vp], *fields[vn], i2, j2, k2, eps);
        if (ok > 0)
        { 
          newOrder[c-1] = vp;
          dejaVu[vn] = 1; vpp = vp; vp = vn;
          goto next;
        }

        if (c != 1) 
        {
          ok = K_CONNECT::getCoincidentVertexIndices(
            1, jm, km, nit[vp], jm, km, nit[vn], jm, km,
            posx, posy, posz, posx, posy, posz,
            *fields[vp], *fields[vn], i2, j2, k2, eps);
          if (ok > 0)
          { 
            newOrder[c-1] = vp;
            dejaVu[vn] = 1; vpp = vp; vp = vn; 
            goto next;
          }
        }
      }
    }
    printf("Error: TFI: some arrays are not continuous.\n");
    return 0;

    next:;
    if ( ( c == 1 && i2 != 1 ) || 
         ( c == 2 && i2 != nit[vp]) || 
         ( c == 3 && i2 != nit[vp]))
    {
      K_CONNECT::reorderStructField(nit[vp], jm, km, *fields[vp],-1,2,3);
    }
    c++;
  }
  return 1;
}
