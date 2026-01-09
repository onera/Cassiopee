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
# include <stdio.h>
# include "transform.h"

using namespace std;
using namespace K_FLD;
using namespace K_FUNC;
using namespace K_CONST;

// ============================================================================
/* Join two arrays */
// ============================================================================
PyObject* K_TRANSFORM::join(PyObject* self, PyObject* args)
{
  PyObject* array1; PyObject* array2;
  E_Float tol;
  if (!PYPARSETUPLE_(args, OO_ R_,
                    &array1, &array2, &tol))
  {
      return NULL;
  }

  // Check array
  E_Int im1, jm1, km1, im2, jm2, km2;
  FldArrayF* f1; FldArrayF* f2;
  FldArrayI* cn1; FldArrayI* cn2;
  char *varString1; char *varString2;
  char* eltType1; char* eltType2;
  E_Int res1 = K_ARRAY::getFromArray3(
    array1, varString1, f1, im1, jm1, km1, cn1, eltType1);
  E_Int res2 = K_ARRAY::getFromArray3(
    array2, varString2, f2, im2, jm2, km2, cn2, eltType2);

  vector<E_Int> pos1; vector<E_Int> pos2;
  E_Int im, jm, km;
  E_Int res = 0;

  if ((res1 == 1 && res2 == 2) || (res1 == 2 && res2 == 1))
  {
    RELEASESHAREDS(array1,f1); RELEASESHAREDU(array2,f2,cn2);
    PyErr_SetString(PyExc_TypeError,
                    "join: can not be used with one structured and one unstructured array.");
    return NULL;
  }
  else if ( res1 < 1 || res1 > 2 || res2 < 1 || res2 > 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "join: one array is invalid.");
    return NULL;
  }
  else
  {
    char* varString = new char [strlen(varString1) + strlen(varString2) + 4];
    E_Int res0 = K_ARRAY::getPosition(varString1, varString2, pos1, pos2, varString);
    if (res0 == -1)
    {
      RELEASESHAREDB(res1, array1, f1, cn1); RELEASESHAREDB(res2, array2, f2, cn2);
      PyErr_SetString(PyExc_ValueError,
                      "join: one array is empty.");
      delete [] varString; return NULL;
    }
    if (pos1.size() != (size_t)f1->getNfld() || pos2.size() != (size_t)f2->getNfld())
      printf("Warning: join: some variables are different. Only variables %s are kept.\n", varString);

    E_Int posx1 = K_ARRAY::isCoordinateXPresent(varString1);
    E_Int posy1 = K_ARRAY::isCoordinateYPresent(varString1);
    E_Int posz1 = K_ARRAY::isCoordinateZPresent(varString1);
    E_Int posx2 = K_ARRAY::isCoordinateXPresent(varString2);
    E_Int posy2 = K_ARRAY::isCoordinateYPresent(varString2);
    E_Int posz2 = K_ARRAY::isCoordinateZPresent(varString2);
    if (posx1 == -1 || posy1 == -1 || posz1 == -1 ||
        posx2 == -1 || posy2 == -1 || posz2 == -1)
    {
      RELEASESHAREDB(res1, array1, f1, cn1); RELEASESHAREDB(res2, array2, f2, cn2);
      PyErr_SetString(PyExc_ValueError,
                      "join: coordinates not found in arrays.");
      delete [] varString;
      return NULL;
    }
    posx1++; posy1++; posz1++;
    posx2++; posy2++; posz2++;

    PyObject* tpl = NULL;
    E_Int api = f1->getApi();
    if (res1 == 1 && res2 == 1) //structured
    {
      FldArrayF* an = new FldArrayF();
      FldArrayF& field = *an;
      res = joinStructured(*f1, im1, jm1, km1, posx1, posy1, posz1,
                           *f2, im2, jm2, km2, posx2, posy2, posz2,
                           pos1, pos2, field, im, jm, km, tol);

      if (res == 0)
        PyErr_SetString(PyExc_TypeError,
                        "join: failed to join structured connectivities");
      else tpl = K_ARRAY::buildArray3(*an, varString, im, jm, km, api);

      delete an; delete [] varString;
      RELEASESHAREDS(array1,f1); RELEASESHAREDS(array2,f2);
      return tpl;
    }
    else if (res1 == 2 && res2 == 2)
    {
      if (K_STRING::cmp(eltType1, "NGON") == 0)
      {
        if (K_STRING::cmp(eltType1, eltType2) != 0)
        {
          PyErr_SetString(PyExc_TypeError,
                          "join: can only join NGON array with another "
                          "NGON array");
        }
        else
        {
          tpl = joinNGON(*f1, *cn1, *f2, *cn2, posx1, posy1, posz1,
                        varString, tol);
        }
      }
      else
      {
        tpl = joinUnstructured(*f1, *cn1, *f2, *cn2, posx1, posy1, posz1,
                               eltType1, eltType2, varString1, tol);
      }

      RELEASESHAREDU(array1, f1, cn1); RELEASESHAREDU(array2, f2, cn2);
      if (tpl == NULL)
        PyErr_SetString(PyExc_TypeError,
                        "join: failed to join unstructured or NGON "
                        "connectivities");
      delete [] varString;
      return tpl;
    }
  }
  return NULL;
}
//=============================================================================
/* field est alloue ici */
//=============================================================================
E_Int
K_TRANSFORM::joinStructured(FldArrayF f1, E_Int im1, E_Int jm1, E_Int km1,
                            E_Int posx1, E_Int posy1, E_Int posz1,
                            FldArrayF f2, E_Int im2, E_Int jm2, E_Int km2,
                            E_Int posx2, E_Int posy2, E_Int posz2,
                            vector<E_Int>& pos1, vector<E_Int>& pos2,
                            FldArrayF& field,
                            E_Int& im, E_Int& jm, E_Int& km, E_Float tol)
{
  if (im1 != 1 && jm1 != 1 && km1 != 1)
  {
    if (im2 != 1 && jm2 != 1 && km2 != 1) // 3D
    {
      return joinstructured3d(f1, im1, jm1, km1, posx1, posy1, posz1,
                              f2, im2, jm2, km2, posx2, posy2, posz2,
                              pos1, pos2, field, im, jm, km, tol);
    }
    else
    {
      printf("Warning: join: arrays must be both 3D.\n");
      return 0;
    }
  }
  else if ((im1 == 1 && jm1 == 1) ||
           (jm1 == 1 && km1 == 1) ||
           (im1 == 1 && km1 == 1))
  {
    if ((im2 == 1 && jm2 == 1) ||
        (jm2 == 1 && km2 == 1) ||
        (im2 == 1 && km2 == 1))
    {
      return joinstructured1d(f1, im1, jm1, km1, posx1, posy1, posz1,
                              f2, im2, jm2, km2, posx2, posy2, posz2,
                              pos1, pos2, field, im, jm, km, tol);
    }
    else
    {
      printf("Warning: join: arrays must be both 1D.\n");
      return 0;
    }
  }
  else if (im1 == 1 || jm1 == 1 || km1 == 1)
  {
    if (im2 == 1 || jm2 == 1 || km2 == 1)
    {
      return joinstructured2d(f1, im1, jm1, km1, posx1, posy1, posz1,
                              f2, im2, jm2, km2, posx2, posy2, posz2,
                              pos1, pos2, field, im, jm, km, tol);
    }
    else
    {
      printf("Warning: join: arrays must be both 2D.\n");
      return 0;
    }
  }
  else
  {
    printf("Warning: join: invalid dimensions.\n");
    return 0;
  }
}
//=============================================================================
/* Join 3d */
//=============================================================================
E_Int
K_TRANSFORM::joinstructured3d(FldArrayF& f1, E_Int im1, E_Int jm1, E_Int km1,
                              E_Int posx1, E_Int posy1, E_Int posz1,
                              FldArrayF& f2, E_Int im2, E_Int jm2, E_Int km2,
                              E_Int posx2, E_Int posy2, E_Int posz2,
                              vector<E_Int>& pos1, vector<E_Int>& pos2,
                              FldArrayF& field,
                              E_Int& im, E_Int& jm, E_Int& km, E_Float tol)
{
  E_Int nof1, nof2;
  E_Int ok = K_CONNECT::detectMatchInterface(im1, jm1, km1,
                                             posx1, posy1, posz1,
                                             im2, jm2, km2,
                                             posx2, posy2, posz2,
                                             f1, f2, nof1, nof2, tol);
  if (ok == 0) return 0;
  switch (nof1)
  {
    case 1:
      K_CONNECT::reorderStructField(im1, jm1, km1, f1, -1, -2, 3);
      break;
    case 2:
      break;
    case 3:
      K_CONNECT::reorderStructField(im1, jm1, km1, f1, 2, -1, 3);
      break;
    case 4:
      K_CONNECT::reorderStructField(im1, jm1, km1, f1, -2, 1, 3);
      break;
    case 5:
      K_CONNECT::reorderStructField(im1, jm1, km1, f1, -2, 3, -1);
      break;
    case 6:
      K_CONNECT::reorderStructField(im1, jm1, km1, f1, 2, 3, 1);
      break;
    default:
      return 0;
  }
  switch (nof2)
  {
    case 1:
      break;
    case 2:
      K_CONNECT::reorderStructField(im2, jm2, km2, f2, -1, -2, 3);
      break;
    case 3:
      K_CONNECT::reorderStructField(im2, jm2, km2, f2, -2, 1, 3);
      break;
    case 4:
      K_CONNECT::reorderStructField(im2, jm2, km2, f2, 2, -1, 3);
      break;
    case 5:
      K_CONNECT::reorderStructField(im2, jm2, km2, f2, 2, 3, 1);
      break;
    case 6:
      K_CONNECT::reorderStructField(im2, jm2, km2, f2, -2, 3, -1);
      break;
    default:
      return 0;
  }

  // determination des pts coincidents avec A1(im1,1,1) et B1(im1,jm1,1)
  E_Int iA2, jA2, kA2, iB2, jB2, kB2;
  ok = K_CONNECT::getCoincidentVertexIndices(
    im1, 1, 1, im1, jm1, km1, im2, jm2, km2,
    posx1, posy1, posz1, posx2, posy2, posz2,
    f1, f2, iA2, jA2, kA2, tol);
  if (ok == 0) return 0;
  ok = K_CONNECT:: getCoincidentVertexIndices(
    im1, jm1, 1, im1, jm1, km1, im2, jm2, km2,
    posx1, posy1, posz1, posx2, posy2, posz2,
    f1, f2, iB2, jB2, kB2, tol);
  if (ok == 0) return 0;
  if (ok > 1)// ambigu : test un pt a cote
  {
    //Point A1(im1,2,2) a tester avec (1,2,2),(1,im2-1,2),(1,2,km2-1),(1,jm2-1,km2-1)
    ok = nextCornerMatchingIndices(im1, 2, 2, im1, jm1, km1, im2, jm2, km2,
                                   posx1, posy1, posz1, posx2, posy2, posz2,
                                   f1, f2, iA2, jA2, kA2, tol);
    if (ok == 0) return 0;
    //Point B1(im1,jm1-1,2) a tester avec les memes points
    ok = nextCornerMatchingIndices(im1, jm1-1, 2, im1, jm1, km1, im2, jm2, km2,
                                   posx1, posy1, posz1, posx2, posy2, posz2,
                                   f1, f2, iB2, jB2, kB2, tol);
    if (ok == 0) return 0;
    if (jA2 == 2) jA2 = 1;
    else jA2 = jm2;
    if (jB2 == 2) jB2 = 1;
    else jB2 = jm2;
    if (kA2 == 2) kA2 = 1;
    else kA2 = km2;
    if (kB2 == 2) kB2 = 1;
    else kB2 = km2;
  }
  if (iA2 != 1 || iB2 != 1) return 0;

  if (jA2 == 1 && kA2 == 1)
  {
    if (jB2 == jm2 && kB2 == 1)
    {
      ;
    }
    else if (jB2 == 1 && kB2 == km2)
      K_CONNECT::reorderStructField(im2, jm2, km2, f2, 1, 3, 2);
    else return 0;
  }
  else if (jA2 == jm2 && kA2 == 1)
  {
    if (jB2 == jm2 && kB2 == km2)
      K_CONNECT::reorderStructField(im2, jm2, km2, f2, 1, -3, 2);

    else if (jB2 == 1 && kB2 == 1)
      K_CONNECT::reorderStructField(im2, jm2, km2, f2, 1, -2, 3);

    else return 0;
  }
  else if (jA2 == jm2 && kA2 == km2)
  {
    if (jB2 == jm2 && kB2 == 1)
      K_CONNECT::reorderStructField(im2, jm2, km2, f2, 1, -3, -2);
    else if (jB2 == 1 && kB2 == km2)
      K_CONNECT::reorderStructField(im2, jm2, km2, f2, 1, -2, -3);
    else return 0;
  }
  else if (jA2 == 1 && kA2 == km2)
  {
    if (jB2 == 1 && kB2 == 1)
      K_CONNECT::reorderStructField(im2, jm2, km2, f2, 1, 3, -2);
    else if (jB2 == jm2 && kB2 == km2)
      K_CONNECT::reorderStructField(im2, jm2, km2, f2, 1, 2, -3);
    else return 0;
  }
  else return 0;

  // assemblage des arrays
  im = im1+im2-1;
  jm = jm2;
  km = km2;
  E_Int nfld = pos1.size();
  field.malloc(im*jm*km, nfld);

  E_Int im1jm1 = im1*jm1;
  E_Int im2jm2 = im2*jm2;

  #pragma omp parallel
  {
    E_Int ind, ind1, ind2;
    for (E_Int eq = 0 ; eq < nfld; eq++)
    {
      E_Int eq1 = pos1[eq];
      E_Int eq2 = pos2[eq];
      E_Float* floc1 = f1.begin(eq1);
      E_Float* floc2 = f2.begin(eq2);
      E_Float* fcnt = field.begin(eq+1);
      for (E_Int k = 0; k < km; k++)
        for (E_Int j = 0; j < jm; j++)
        {
          #pragma omp for
          for (E_Int i = 0; i < im1-1; i++)
          {
            ind1 = i + j * im1 + k * im1jm1;
            ind = i + j * (im1+im2-1) + k*(im1+im2-1)*jm;
            fcnt[ind] = floc1[ind1];
          }
          #pragma omp for
          for (E_Int i = 0; i < im2; i++)
          {
            ind2 = i + j * im2 + k * im2jm2;
            ind = i + im1-1 + j * (im1+im2-1) + k*(im1+im2-1)*jm;
            fcnt[ind] = floc2[ind2];
          }
        }
      }
    }

  // Remet l'array final avec la numerotation de f1
  /*
  switch (nof1)
  {
    case 1:
      K_CONNECT::reorderStructField(im, jm, km, field, -1, -2, 3);
      break;
    case 2:
      break;
    case 3:
      K_CONNECT::reorderStructField(im, jm, km, field, -2, 1, 3);
      break;
    case 4:
      K_CONNECT::reorderStructField(im, jm, km, field, 2, -1, 3);
      break;
    case 5:
      K_CONNECT::reorderStructField(im, jm, km, field, -3, -1, 2);
      break;
    case 6:
      K_CONNECT::reorderStructField(im, jm, km, field, 3, 1, 2);
      break;
    default:
      return 0;
  }
  */

  return 1;
}
//=============================================================================
/* Join 2d */
//=============================================================================
E_Int
K_TRANSFORM::joinstructured2d(FldArrayF& f1, E_Int im1, E_Int jm1, E_Int km1,
                              E_Int posx1, E_Int posy1, E_Int posz1,
                              FldArrayF& f2, E_Int im2, E_Int jm2, E_Int km2,
                              E_Int posx2, E_Int posy2, E_Int posz2,
                              std::vector<E_Int>& pos1,
                              std::vector<E_Int>& pos2,
                              FldArrayF& field,
                              E_Int& im, E_Int& jm, E_Int& km, E_Float tol)
{
  E_Int nof1, nof2;
  if (im1 == 1) K_CONNECT::reorderStructField(im1, jm1, km1, f1, 3, 1, 2);
  if (jm1 == 1) K_CONNECT::reorderStructField(im1, jm1, km1, f1, 1, 3, 2);
  if (im2 == 1) K_CONNECT::reorderStructField(im2, jm2, km2, f2, 3, 1, 2);
  if (jm2 == 1) K_CONNECT::reorderStructField(im2, jm2, km2, f2, 1, 3, 2);
  E_Int isok = K_CONNECT::detectMatchInterface(im1, jm1, km1,
                                               posx1, posy1, posz1,
                                               im2, jm2, km2,
                                               posx2, posy2, posz2,
                                               f1, f2, nof1, nof2, tol);
  if (isok == 0) return 0;
  //printf("2D, %d: %d %d\n", isok, nof1, nof2);

  switch (nof1)
  {
    case 1:
      K_CONNECT::reorderStructField(im1, jm1, km1, f1, -1, -2, 3);
      break;
    case 2:
      break;
    case 3:
      K_CONNECT::reorderStructField(im1, jm1, km1, f1, 2, -1, 3);
      break;
    case 4:
      K_CONNECT::reorderStructField(im1, jm1, km1, f1, 2, 1, 3); // -2,1,3
      break;
    default:
      return 0;
  }
  switch (nof2)
  {
    case 1:
      break;
    case 2:
      K_CONNECT::reorderStructField(im2, jm2, km2, f2, -1, -2, 3);
      break;
    case 3:
      K_CONNECT::reorderStructField(im2, jm2, km2, f2, 2, 1, 3); // -2,1,3
      break;
    case 4:
      K_CONNECT::reorderStructField(im2, jm2, km2, f2, 2, -1, 3);
      break;
    default:
      return 0;
  }
  // bloc en O ?
  E_Int nfld = pos1.size();
  // imax de blk1 est connecte a imin de blk2
  E_Int ind11 = im1-1; E_Int ind12 = im1-1 + (jm1-1)*im1;
  E_Int ind21 = 0; E_Int ind22 = (jm2-1)*im2;
  E_Float dx1 = f1(ind11, posx1)-f1(ind12, posx1);
  E_Float dy1 = f1(ind11, posy1)-f1(ind12, posy1);
  E_Float dz1 = f1(ind11, posz1)-f1(ind12, posz1);
  E_Float dx2 = f2(ind21, posx2)-f2(ind22, posx2);
  E_Float dy2 = f2(ind21, posy2)-f2(ind22, posy2);
  E_Float dz2 = f2(ind21, posz2)-f2(ind22, posz2);
  if (K_FUNC::E_abs(dx1) < tol && K_FUNC::E_abs(dy1) < tol &&
      K_FUNC::E_abs(dz1) < tol &&
      K_FUNC::E_abs(dx2) < tol && K_FUNC::E_abs(dy2) < tol &&
      K_FUNC::E_abs(dz2) < tol)
  {
    // reconnecter les pts du bloc 2 dans l'ordre par rapport au bloc 1
    FldArrayF tmp(f2.getSize(), f2.getNfld()); tmp.setAllValuesAtNull();
    E_Int ind2n;
    for (E_Int j1 = 1; j1 <= jm1; j1++)
    {
      ind11 = im1-1 + (j1-1)*im1;
      for (E_Int j2 = 1; j2 <= jm2; j2++)
      {
        ind21 = (j2-1)*im2;
        ind2n = (j1-1)*im2;
        dx1 = f1(ind11, posx1)-f2(ind21, posx1);
        dy1 = f1(ind11, posy1)-f2(ind21, posy1);
        dz1 = f1(ind11, posz1)-f2(ind21, posz1);
        if (K_FUNC::E_abs(dx1) < tol && K_FUNC::E_abs(dy1) < tol &&
            K_FUNC::E_abs(dz1) < tol )
        {
          for (E_Int eq = 1; eq <= nfld; eq++)
          {
            tmp(ind2n, eq) = f2(ind21, eq);
            for (E_Int i2 = 1; i2 < im2; i2++)
              tmp(ind2n+i2, eq) = f2(ind21+i2, eq);
          }
          break;
        }
      }
    }
    f2 = tmp;
  }
  else
  {
    // test match A1(im1,1) avec  (1,1) ou (1,jm2) ?
    E_Int indA1 = im1-1;
    E_Float x1 = f1(indA1, posx1);
    E_Float y1 = f1(indA1, posy1);
    E_Float z1 = f1(indA1, posz1);

    E_Float x2 = f2(0, posx2);  // (1,1) de mesh2
    E_Float y2 = f2(0, posy2);
    E_Float z2 = f2(0, posz2);

    E_Float dx = x2-x1;
    E_Float dy = y2-y1;
    E_Float dz = z2-z1;

    if (K_FUNC::fEqualZero(dx, tol) == false ||
        K_FUNC::fEqualZero(dy, tol) == false ||
        K_FUNC::fEqualZero(dz, tol) == false)
      K_CONNECT::reorderStructField(im2, jm2, km2, f2, 1,-2, 3);
  }
  /*-----------------------*/
  /* assemblage des arrays */
  /*-----------------------*/
  E_Int ind, ind1, ind2;

  im = im1+im2-1;
  jm = jm1;
  km = 1;
  field.malloc(im*jm, nfld);
  for (E_Int eq = 0 ; eq < nfld; eq++)
  {
    E_Int eq1 = pos1[eq];
    E_Int eq2 = pos2[eq];
    E_Float* floc1 = f1.begin(eq1);
    E_Float* floc2 = f2.begin(eq2);
    E_Float* fcnt = field.begin(eq+1);

    ind = 0;
    for (E_Int j = 0; j < jm; j++)
    {
      for (E_Int i = 0; i < im1-1; i++)
      {
        ind1 = i + j * im1;
        fcnt[ind] = floc1[ind1];
        ind++;
      }
      for (E_Int i = 0; i < im2; i++)
      {
        ind2 = i + j * im2;
        fcnt[ind] = floc2[ind2];
        ind++;
      }
    }
  }
  return 1;
}
//=============================================================================
/* Join 1d */
//=============================================================================
E_Int
K_TRANSFORM::joinstructured1d(FldArrayF& f1, E_Int im1, E_Int jm1, E_Int km1,
                              E_Int posx1, E_Int posy1, E_Int posz1,
                              FldArrayF& f2, E_Int im2, E_Int jm2, E_Int km2,
                              E_Int posx2, E_Int posy2, E_Int posz2,
                              std::vector<E_Int>& pos1,
                              std::vector<E_Int>& pos2,
                              FldArrayF& field,
                              E_Int& im, E_Int& jm, E_Int& km, E_Float tol)
{
  if (jm1 != 1) K_CONNECT::reorderStructField(im1, jm1, km1, f1, 3, 1, 2);
  if (km1 != 1) K_CONNECT::reorderStructField(im1, jm1, km1, f1, 2, 3, 1);
  if (jm2 != 1) K_CONNECT::reorderStructField(im2, jm2, km2, f2, 3, 1, 2);
  if (km2 != 1) K_CONNECT::reorderStructField(im2, jm2, km2, f2, 2, 3, 1);

  E_Int nof1, nof2;
  E_Int isok = K_CONNECT::detectMatchInterface(
    im1, jm1, km1, posx1, posy1, posz1,
    im2, jm2, km2, posx2, posy2, posz2,
    f1, f2, nof1, nof2, tol);
  if (isok == 0) return 0;

  switch (nof1)
  {
    case 1:
      K_CONNECT::reorderStructField(im1, jm1, km1, f1, -1, 2, 3);
      break;
    case 2: // deja ordonne
      break;
    default:
      return 0;
  }
  switch (nof2)
  {
    case 1: //deja ordonne
      break;
    case 2:
      K_CONNECT::reorderStructField(im2, jm2, km2, f2, -1, 2, 3);
      break;
    default:
      return 0;
  }

  // assemblage des arrays
  E_Int nfld = pos1.size();
  field.malloc(im1+im2-1, nfld);
  E_Int inc = im1-1;

  E_Int pos1Size = pos1.size();
  for (E_Int eq = 0; eq < pos1Size;  eq++)
  {
    E_Int eq1 = pos1[eq];
    E_Int eq2 = pos2[eq];
    E_Float* floc1 = f1.begin(eq1);
    E_Float* floc2 = f2.begin(eq2);
    E_Float* fcnt = field.begin(eq+1);
    for (E_Int i = 0; i < im1; i++) fcnt[i] = floc1[i];
    for (E_Int i = 1; i < im2; i++) fcnt[i+inc] = floc2[i];
  }

  // taille du maillage
  im = field.getSize(); jm = 1; km = 1;
  return 1;
}

//=============================================================================
/*
   Join topologique : somme des vertex et de la connectivite
   + cleanConnectivity
   Gere tous types d'elements (BE + BE -> ME, BE + ME -> ME, ME + ME -> ME)
   Il faut que les champs soient ranges dans le meme ordre.
*/
//=============================================================================
PyObject* K_TRANSFORM::joinUnstructured(FldArrayF& f1, FldArrayI& cn1,
                                        FldArrayF& f2, FldArrayI& cn2,
                                        E_Int posx, E_Int posy, E_Int posz,
                                        char* eltType1, char* eltType2,
                                        char* varString, E_Float tol)
{
  // Acces universel sur BE/ME
  E_Int nc1 = cn1.getNConnect(), nc2 = cn2.getNConnect();
  E_Int nc = nc1;
  E_Int nfld = f1.getNfld();
  E_Int npts1 = f1.getSize(), npts2 = f2.getSize();
  E_Int npts = npts1 + npts2;

  // Acces universel aux eltTypes
  vector<char*> eltTypes, eltTypes2;
  K_ARRAY::extractVars(eltType2, eltTypes2);
  // Concatenate elttypes, discard duplicates
  char eltType[K_ARRAY::VARSTRINGLENGTH]; eltType[0] = '\0';
  strcpy(eltType, eltType1);
  for (E_Int ic = 0; ic < nc2; ic++)
  {
    char* eltTypConn2 = eltTypes2[ic];
    if (strstr(eltType, eltTypConn2) == NULL)
    {
      strcat(eltType, ",");
      strcat(eltType, eltTypConn2);
      nc += 1;
    }
  }

  // ME
  E_Int api = f1.getApi();

  K_ARRAY::extractVars(eltType, eltTypes);
  vector<E_Int> nelts(nc, 0);
  // Table d'indirection pour la seconde connectivite ME
  vector<E_Int> indir2(nc2, -1);

  // Calcule le nombre d'elements par type de connectivite et
  // rempli la table d'indirection
  for (E_Int ic = 0; ic < nc; ic++)
  {
    char* eltTypConn = eltTypes[ic];
    if (ic < nc1)
    {
      FldArrayI& cm1 = *(cn1.getConnect(ic));
      nelts[ic] += cm1.getSize();
    }

    for (E_Int ic2 = 0; ic2 < nc2; ic2++)
    {
      char* eltTypConn2 = eltTypes2[ic2];
      if (K_STRING::cmp(eltTypConn, eltTypConn2) == 0)
      {
        FldArrayI& cm2 = *(cn2.getConnect(ic2));
        nelts[ic] += cm2.getSize();
        indir2[ic2] = ic;
      }
    }
  }

  // Fusion des connectivites
  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts, nelts,
                                       eltType, false, api);
  FldArrayF* f; FldArrayI* cn;
  K_ARRAY::getFromArray3(tpl, f, cn);

  #pragma omp parallel
  {
    // Copie des champs aux noeuds
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* f1n = f1.begin(n);
      E_Float* f2n = f2.begin(n);
      E_Float* fn = f->begin(n);
      #pragma omp for
      for (E_Int i = 0; i < npts1; i++) fn[i] = f1n[i];
      #pragma omp for
      for (E_Int i = 0; i < npts2; i++) fn[i+npts1] = f2n[i];
    }

    // Boucle sur toutes les connectivites du premier ME
    for (E_Int ic = 0; ic < nc1; ic++)
    {
      FldArrayI& cm1 = *(cn1.getConnect(ic));
      FldArrayI& cm = *(cn->getConnect(ic));
      #pragma omp for
      for (E_Int i = 0; i < cm1.getSize(); i++)
        for (E_Int j = 1; j <= cm1.getNfld(); j++)
          cm(i,j) = cm1(i,j);
    }

    // Boucle sur toutes les connectivites du second ME
    for (E_Int ic = 0; ic < nc2; ic++)
    {
      FldArrayI& cm2 = *(cn2.getConnect(ic));
      FldArrayI& cm = *(cn->getConnect(indir2[ic]));
      E_Int nelts2 = cm2.getSize();
      E_Int offsetElts = cm.getSize() - nelts2;
      #pragma omp for
      for (E_Int i = 0; i < nelts2; i++)
        for (E_Int j = 1; j <= cm2.getNfld(); j++)
          // Add offsets
          cm(i+offsetElts,j) = cm2(i,j) + npts1;
    }
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
  for (size_t ic = 0; ic < eltTypes2.size(); ic++) delete [] eltTypes2[ic];
  // Clean connectivity
  if (posx > 0 && posy > 0 && posz > 0)
  {
    PyObject* tpl2 = K_CONNECT::V_cleanConnectivity(varString, *f, *cn, eltType, tol);
    RELEASESHAREDU(tpl, f, cn); Py_DECREF(tpl);
    return tpl2;
  }
  else
  {
    RELEASESHAREDU(tpl, f, cn);
    return tpl;
  }
}
//=============================================================================
/* Join topologique: somme des vertex et de la connectivite
   + cleanConnectivity
   Il faut que les champs soient ranges dans le meme ordre */
//=============================================================================
PyObject* K_TRANSFORM::joinNGON(FldArrayF& f1, FldArrayI& cn1,
                                FldArrayF& f2, FldArrayI& cn2,
                                E_Int posx, E_Int posy, E_Int posz,
                                char* varString, E_Float tol)
{
  E_Int nfaces1 = cn1.getNFaces(), sizeFN1 = cn1.getSizeNGon();
  E_Int nfaces2 = cn2.getNFaces(), sizeFN2 = cn2.getSizeNGon();
  E_Int nelts1 = cn1.getNElts(), sizeEF1 = cn1.getSizeNFace();
  E_Int nelts2 = cn2.getNElts(), sizeEF2 = cn2.getSizeNFace();

  E_Int nfld = f1.getNfld();
  E_Int npts1 = f1.getSize();
  E_Int npts2 = f2.getSize();

  // Fusion des connectivites
  E_Int nfaces = nfaces1 + nfaces2;
  E_Int sizeFN = sizeFN1 + sizeFN2;
  E_Int nelts = nelts1 + nelts2;
  E_Int sizeEF = sizeEF1 + sizeEF2;
  E_Int npts = npts1 + npts2;

  E_Int api = f1.getApi();
  E_Int ngonType = cn1.getNGonType();
  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts, nelts,
                                       nfaces, "NGON", sizeFN, sizeEF,
                                       ngonType, false, api);
  FldArrayF* f; FldArrayI* cn;
  K_ARRAY::getFromArray3(tpl, f, cn);

  // Acces non universel sur les ptrs
  E_Int *ngon1 = cn1.getNGon(), *ngon2 = cn2.getNGon(), *ngon = cn->getNGon();
  E_Int *nface1 = cn1.getNFace(), *nface2 = cn2.getNFace(), *nface = cn->getNFace();
  E_Int *indPG1 = NULL, *indPG2 = NULL, *indPG = NULL;
  E_Int *indPH1 = NULL, *indPH2 = NULL, *indPH = NULL;
  if (ngonType == 2 || ngonType == 3)
  {
    indPG1 = cn1.getIndPG(); indPG2 = cn2.getIndPG(); indPG = cn->getIndPG();
    indPH1 = cn1.getIndPH(); indPH2 = cn2.getIndPH(); indPH = cn->getIndPH();
  }

  #pragma omp parallel
  {
    // Copie des champs aux noeuds
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* f1n = f1.begin(n);
      E_Float* f2n = f2.begin(n);
      E_Float* fn = f->begin(n);
      #pragma omp for nowait
      for (E_Int i = 0; i < npts1; i++) fn[i] = f1n[i];
      #pragma omp for nowait
      for (E_Int i = 0; i < npts2; i++) fn[i+npts1] = f2n[i];
    }

    // Copie des connectivites (add offset to all elements of the second
    // connectivity and correct outside the parallel block)
    #pragma omp for nowait
    for (E_Int i = 0; i < sizeFN1; i++) ngon[i] = ngon1[i];
    #pragma omp for nowait
    for (E_Int i = 0; i < sizeFN2; i++) ngon[i+sizeFN1] = ngon2[i] + npts1;

    #pragma omp for nowait
    for (E_Int i = 0; i < sizeEF1; i++) nface[i] = nface1[i];
    #pragma omp for nowait
    for (E_Int i = 0; i < sizeEF2; i++) nface[i+sizeEF1] = nface2[i] + nfaces1;

    if (ngonType == 2 || ngonType == 3)
    {
        #pragma omp for nowait
        for (E_Int i = 0; i < nfaces1; i++) indPG[i] = indPG1[i];
        #pragma omp for nowait
        for (E_Int i = 0; i < nfaces2; i++) indPG[i+nfaces1] = indPG2[i] + nfaces1;

        #pragma omp for nowait
        for (E_Int i = 0; i < nelts1; i++) indPH[i] = indPH1[i];
        #pragma omp for
        for (E_Int i = 0; i < nelts2; i++) indPH[i+nelts1] = indPH2[i] + nelts1;
      }
  }

  // Correction for number of vertices per face and number of faces per element
  // for the second connectivity
  if (api != 3)
  {
    E_Int ind = 0;
    for (E_Int i = 0; i < nfaces2; i++)
    {
      ngon[sizeFN1+ind] = ngon2[ind];
      ind += ngon2[ind]+1;
    }
    ind = 0;
    for (E_Int i = 0; i < nelts2; i++)
    {
      nface[sizeEF1+ind] = nface2[ind];
      ind += nface2[ind]+1;
    }
  }

  if (posx > 0 && posy > 0 && posz > 0)
  {
    PyObject* tpl2 = K_CONNECT::V_cleanConnectivity(varString, *f, *cn, "NGON", tol);
    RELEASESHAREDU(tpl, f, cn); Py_DECREF(tpl);
    return tpl2;
  }
  else
  {
    RELEASESHAREDU(tpl, f, cn);
    return tpl;
  }
}
//=============================================================================
/*  Test si un pt(i,j,k) coincident avec un des coins definis par
    (1,2,2), (1,im2-1,2),(1,2,km2-1),(1,jm2-1,km2-1) de f2
    attention: les indices demarrent a 1 en i j et k  */
//=============================================================================
E_Int K_TRANSFORM::nextCornerMatchingIndices(
  E_Int i1, E_Int j1, E_Int k1,
  E_Int im1, E_Int jm1, E_Int km1,
  E_Int im2, E_Int jm2, E_Int km2,
  E_Int posx1, E_Int posy1, E_Int posz1,
  E_Int posx2, E_Int posy2, E_Int posz2,
  FldArrayF& f1, FldArrayF& f2,
  E_Int& i2, E_Int& j2, E_Int& k2,
  E_Float eps)
{
  E_Int im1jm1 = im1*jm1;
  E_Int ind1 = (i1-1) + (j1-1) * im1 + (k1-1) * im1jm1;
  E_Float* xt1 = f1.begin(posx1);
  E_Float* xt2 = f2.begin(posx2);
  E_Float* yt1 = f1.begin(posy1);
  E_Float* yt2 = f2.begin(posy2);
  E_Float* zt1 = f1.begin(posz1);
  E_Float* zt2 = f2.begin(posz2);
  E_Float dx, dy, dz;
  E_Float xt11 = xt1[ind1];
  E_Float yt11 = yt1[ind1];
  E_Float zt11 = zt1[ind1];
  E_Int im2jm2 = im2*jm2;
  E_Float eps2 = 3*eps*eps;
  E_Int count = 0; E_Int isav=-1, jsav=-1, ksav=-1;
  // test (1,2,2)
  i2 = 1; j2 = 2; k2 = 2;
  E_Int ind2 = i2-1 + (j2-1)*im2+(k2-1)*im2jm2;
  dx = xt2[ind2] - xt11;
  dy = yt2[ind2] - yt11;
  dz = zt2[ind2] - zt11;
  if (dx*dx+dy*dy+dz*dz < eps2) {isav=i2;jsav=j2; ksav=k2;count++;}

  // test (1,jm2-1,2)
  i2 = 1; j2 = jm2-1; k2 = 2;
  ind2 = i2-1 + (j2-1)*im2+(k2-1)*im2jm2;
  dx = xt2[ind2] - xt11;
  dy = yt2[ind2] - yt11;
  dz = zt2[ind2] - zt11;
  if (dx*dx+dy*dy+dz*dz < eps2) {isav=i2;jsav=j2; ksav=k2;count++;}

  // test (1,2,km2-1)
  i2 = 1; j2 = 2; k2 = km2-1;
  ind2 = i2-1 + (j2-1)*im2+(k2-1)*im2jm2;
  ind2 = (km2-1)*im2jm2;
  dx = xt2[ind2] - xt11;
  dy = yt2[ind2] - yt11;
  dz = zt2[ind2] - zt11;
  if (dx*dx+dy*dy+dz*dz < eps2) {isav=i2;jsav=j2; ksav=k2;count++;}

  // test (1,jm2-1,km2-1)
  i2 = 1; j2 = jm2-1; k2 = km2-1;
  ind2 = i2-1 + (j2-1)*im2+(k2-1)*im2jm2;
  dx = xt2[ind2] - xt11;
  dy = yt2[ind2] - yt11;
  dz = zt2[ind2] - zt11;
  if (dx*dx+dy*dy+dz*dz < eps2) {isav=i2;jsav=j2; ksav=k2;count++;}
  i2 = isav; j2 = jsav; k2 = ksav;
  return count;
}
