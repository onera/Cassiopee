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

# include "transform.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* splitCurvatureRadius
   Decoupe un i-array aux pts de courbure importante de la
   courbe. Retourne la liste des morceaux. */
//=============================================================================
PyObject* K_TRANSFORM::splitCurvatureRadius(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float Rs; // courbure seuil
  if (!PYPARSETUPLE_(args, O_ R_, &array, &Rs))
  {
    return NULL;
  }

  E_Float cvs = 1./K_FUNC::E_max(Rs, 1.e-10);
  E_Float ds = 0.;
  // extraction de l'array 1D
  E_Int im, jm, km;
  FldArrayF* f;
  FldArrayI* cn;
  char* varString;
  char* et;
  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, et);
  if ( res != 1 && res != 2 )
  {
    PyErr_SetString(PyExc_TypeError,
                    "splitCurvatureRadius: input array is not valid.");
    return NULL;
  }
  if ( res == 2 )
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "splitCurvatureRadius: array must be an i-array.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if ( posx == -1 || posy == -1 || posz == -1 )
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError,
                    "splitCurvatureRadius: coordinates not found in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  if ( im < 2 || jm != 1 || km != 1 )
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError,
                    "splitCurvatureRadius: structured array must be an i-array.");
    return NULL;
  }

  E_Int api = f->getApi();
  E_Int npts = f->getSize();
  if ( npts < 6 )
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError,
                    "splitCurvatureRadius: not enough points in i-array.");
    return NULL;
  }

  // Insertion dans une liste des morceaux de courbe splittee
  PyObject* l = PyList_New(0);

  vector<FldArrayF*> fsplit;
  splitSplineStruct(ds, cvs, posx, posy, posz, f, fsplit);

  PyObject* tpl;
  E_Int fsplitSize = fsplit.size();
  for (E_Int v = 0; v < fsplitSize; v++)
  {
    FldArrayF& f0 = *fsplit[v];
    E_Int ni = f0.getSize();
    tpl = K_ARRAY::buildArray3(f0, varString, ni, 1, 1, api);
    delete &f0;
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
  return l;
}

//===========================================================================
/* Split d'une courbe structuree 1D  */
//===========================================================================
void K_TRANSFORM::splitSplineStruct(E_Float dmax, E_Float cmax,
                                    E_Int posx, E_Int posy, E_Int posz, FldArrayF* f,
                                    vector<FldArrayF*>& fsplit)
{
  vector<E_Int> splitPts;
  E_Int npts = f->getSize();

  E_Float* xt = f->begin(posx);
  E_Float* yt = f->begin(posy);
  E_Float* zt = f->begin(posz);

  // Calcul de la courbure aux pts de la courbe
  FldArrayF curv(npts);
  K_COMPGEOM::compCurvature(npts, xt, yt, zt, curv);

  // Determination des indices des points de split
  // on splitte la ou la courbure est forte
  E_Float eps = 1.e-6;
  E_Int N = npts-1;
  for (E_Int i = 1; i < npts-1; i++)
  {
    E_Float cv = K_FUNC::E_abs(curv[i]);
    E_Int im = i-1;
    E_Int ip = i+1;
    if ( i == 0 ) im = 0;
    if ( i == N ) ip = N;

    if ( cv > K_FUNC::E_abs(curv[im])+eps
         && cv > K_FUNC::E_abs(curv[ip])+eps ) //max local
    {
      if ( cv > cmax )
        splitPts.push_back(i);
    }
  }
  splitPts.push_back(N);

  // split de la spline aux pts splitPts
  E_Int iprev = 0;
  E_Int splitSize = splitPts.size();

  for (E_Int i = 0; i < splitSize; i++)
  {
    E_Int ii = splitPts[i];
    E_Int size = ii + 1 - iprev;
    E_Int nfld = f->getNfld();
    FldArrayF* f0 = new FldArrayF(size, nfld);
    E_Int cc = 0;
    for (E_Int ind = iprev; ind <= ii; ind++)
    {
      for (E_Int fld=1;fld<=nfld;fld++)
      {
	E_Float* f0fld = f0->begin(fld);
	E_Float* ffld  = f->begin(fld);
	f0fld[cc] = ffld[ind];
      }
      cc++;
    }
    fsplit.push_back(f0);
    iprev = ii;
  }

  delete f;
}
// //===========================================================================
// /* Split d'une courbe structuree 1D  */
// //===========================================================================
// void splitSplineStruct(E_Float dmax, E_Float cmax,
//                        E_Int posx, E_Int posy, E_Int posz, FldArrayF* f,
//                        vector<FldArrayF*>& fsplit)
// {
//   vector<E_Int> splitPts;
//   E_Int npts = f->getSize();

//   E_Float* xt = f->begin(posx);
//   E_Float* yt = f->begin(posy);
//   E_Float* zt = f->begin(posz);
//   FldArrayF curv(npts);

//   K_COMPGEOM::compCurvature(npts, xt, yt, zt, curv);

// //   FldArrayF cvt(npts);
// //   FldArrayF dpt(npts);
// //   E_Float inv6 = 1./6;
// //   E_Float inv3 = 1./3;
// //   E_Float inv12 = 1./12;
// //   E_Float inv36 = inv12*inv3;
// //   E_Float inv9 = inv3*inv3;
// //   for (E_Int i = 2; i <= N+2; i++)
// //   {
// //     E_Int im2 = node[i-2];
// //     E_Int im1 = node[i-1];
// //     E_Int ip2 = node[i+2];
// //     E_Int ip1 = node[i+1];
// //     E_Int ii = i-2;

// //     //calcul du deplacement
// //     E_Float dtx =
// //       inv36*(xt[im2]+xt[ip2])+2*inv9*(xt[im1]+xt[ip1])-0.5*xt[ii];
// //     E_Float dty =
// //       inv36*(yt[im2]+yt[ip2])+2*inv9*(yt[im1]+yt[ip1])-0.5*yt[ii];

// //     //calcul de la courbure avec deplacement
// //     E_Float bp1 =
// //       inv12*(xt[im2]+xt[ip2])+inv6*(xt[im1]+xt[ip1])-0.5*xt[ii];
// //     E_Float bp2 =
// //       inv12*(yt[im2]+yt[ip2])+inv6*(yt[im1]+yt[ip1])-0.5*yt[ii];
// //     E_Float cp1 =
// //       inv3*(xt[ip1]-xt[im1]) + inv12*(xt[ip2]-xt[im2]);
// //     E_Float cp2 =
// //       inv3*(yt[ip1]-yt[im1]) + inv12*(yt[ip2]-yt[im2]);

// //     E_Float coef = cp1*cp1+cp2*cp2;
// //     E_Float coef32 = sqrt(coef*coef*coef);
// //     cvt[ii] = K_FUNC::E_abs(2 * (cp1*bp2-cp2*bp1)/coef32);

// //     // calcul de la courbure locale
// //     E_Float dxpm = xt[ip1]-xt[im1];
// //     E_Float dypm = yt[ip1]-yt[im1];

// //     dpt[ii] = K_FUNC::E_abs(dtx+dty);
// //     curv[ii] =
// //       4*( dxpm*(yt[ip1]-2*yt[ii]+yt[im1]) - dypm*(xt[ip1]-2*xt[ii]+xt[im1]) ) /
// //       (pow( dxpm*dxpm + dypm*dypm, 3./2. ));
// //   }


//   // Determination des indices des points de split
//   // on splitte la ou la courbure est forte
//   E_Float eps = 1.e-6;
//   for (E_Int i = 0; i < npts; i++)
//   {
//     E_Float cv = K_FUNC::E_abs(curv[i]);
//     //E_Float cv = cvt[i];
// //     E_Float dt = dpt[i];
//     E_Int im = i-1;
//     E_Int ip = i+1;
//     if ( i == 0 ) im = 0;
//     if ( i == N ) ip = N;

//     //if ( cv > cvt[im]+eps && cv > cvt[ip]+eps )//max local
//     if ( cv > K_FUNC::E_abs(curv[im])+eps
//          && cv > K_FUNC::E_abs(curv[ip])+eps )//max local
//     {
//       if ( cv > cmax )//&& dt > dmax)
//         splitPts.push_back(i);
//     }
//   }
//   splitPts.push_back(N);

// //   // On split aussi ensuite suivant la courbure cumulee
// //   vector<E_Int> splitPtsF;
// //   E_Int splitSize = splitPts.size();
// //   E_Int istart = 0;
// //   for (E_Int i = 0; i < splitSize; i++)
// //   {
// //     E_Int ii = splitPts[i];
// //     // Calcul de l'integrale
// //     E_Float I = 0.;
// //     for (E_Int j = istart; j < ii; j++)
// //       I = I + curv[j];
// //     E_Int N = K_FUNC::E_abs(I) / (2000.*cmax) + 1;
// //     E_Float ctmax = K_FUNC::E_abs(I / (1.*N));
// //     printf("I : %f, N : %d, ctmax : %f\n", I, N, ctmax);
// //     E_Float ct = 0.;
// //     for (E_Int j = istart; j < ii; j++)
// //     {
// //       ct = ct + curv[j];
// //       printf("j %d %f %f\n", j, ct, curv[j]);
// //       if ( (ct > ctmax || ct < -ctmax) && j < ii-5)
// //       {
// //         printf("-- split -- %d\n", j);
// //         splitPtsF.push_back(j); ct = 0.;
// //       }
// //     }
// //     splitPtsF.push_back(ii);
// //     istart = ii;
// //   }

// //   splitPts = splitPtsF;
//   E_Int splitSize = splitPts.size();

//   // split de la spline aux pts splitPts
//   E_Int iprev = 0;
//   for (E_Int i = 0; i < splitSize; i++)
//   {
//     E_Int ii = splitPts[i];
//     E_Int size = ii + 1 - iprev;
//     FldArrayF* f0 = new FldArrayF(size, 3);
//     FldArrayF& f0p = *f0;
//     E_Int cc = 0;
//     for (E_Int ind = iprev; ind <= ii; ind++)
//     {
//       f0p(cc,1) = xt[ind];
//       f0p(cc,2) = yt[ind];
//       f0p(cc,3) = zt[ind];
//       cc++;
//     }
//     fsplit.push_back(f0);
//     iprev = ii;
//   }

//   delete f;
// }
