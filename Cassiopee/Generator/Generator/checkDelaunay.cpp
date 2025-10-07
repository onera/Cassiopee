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
# include "generator.h"

using namespace std;
using namespace K_FLD;
using namespace K_CONST;

// ============================================================================
/* Check if the Delaunay triangulation is inside the contours, defined 
   by a set of BARS. If not, a point is inserted in the contours where 
   the triangulation is out of the contours. */
// ============================================================================
PyObject* K_GENERATOR::checkDelaunay(PyObject* self, PyObject* args)
{
  E_Float eps = 1.e-12;
  PyObject* arrayc;
  PyObject* arrayd;  
  if (!PYPARSETUPLE_(args, OO_, &arrayc, &arrayd)) return NULL;
  
  // check array defining the contours
  E_Int imc, jmc, kmc;
  FldArrayF* fc; FldArrayI* cnc;
  char* varStringc; char* eltTypec;
  E_Int resc = K_ARRAY::getFromArray3(arrayc, varStringc, fc, imc, jmc, kmc,
                                      cnc, eltTypec);
  if (resc == 2)
  {
    if (strcmp(eltTypec, "BAR") != 0) 
    {
      PyErr_SetString(PyExc_ValueError,
                      "checkDelaunay: array defining contours must be BAR type.");
      RELEASESHAREDU(arrayc, fc, cnc); return NULL;
    }
  }
  else
  {
    if (resc == 1) RELEASESHAREDS(arrayc, fc);
    PyErr_SetString(PyExc_ValueError,
                    "checkDelaunay: array defining contours must be a BAR-array.");
    return NULL;
  }
  
  // check triangulation array 
  E_Int imd, jmd, kmd;
  FldArrayF* fd; FldArrayI* cnd;
  char* varStringd; char* eltTyped;
  E_Int resd =  K_ARRAY::getFromArray3(arrayd, varStringd, fd, imd, jmd, kmd,
                                       cnd, eltTyped);
  if (resd == 2)
  {
    if ( strcmp(eltTyped, "TRI") != 0 ) 
    {
      RELEASESHAREDU(arrayc, fc, cnc);
      RELEASESHAREDU(arrayd, fd, cnd);
      PyErr_SetString(PyExc_ValueError,
                      "checkDelaunay: array defining Delaunay triangulation must be a TRI-array.");
      return NULL;
    }
  }
  else
  {
    if (resd == 1) RELEASESHAREDS(arrayd, fd);
    RELEASESHAREDU(arrayc, fc, cnc);
    PyErr_SetString(PyExc_ValueError,
                    "checkDelaunay: array defining Delaunay triangulation must be a TRI-array.");
    return NULL;
  }

  // Determination des coordonnees
  E_Int posxc = K_ARRAY::isCoordinateXPresent(varStringc);
  E_Int posyc = K_ARRAY::isCoordinateYPresent(varStringc);
  E_Int poszc = K_ARRAY::isCoordinateZPresent(varStringc);
  E_Int posxd = K_ARRAY::isCoordinateXPresent(varStringd);
  E_Int posyd = K_ARRAY::isCoordinateYPresent(varStringd);
  E_Int poszd = K_ARRAY::isCoordinateZPresent(varStringd);
  if (posxc == -1 || posyc == -1 || poszc == -1 ||
      posxd == -1 || posyd == -1 || poszd == -1)
  {
    RELEASESHAREDU(arrayc, fc, cnc);
    RELEASESHAREDU(arrayd, fd, cnd);
    PyErr_SetString(PyExc_ValueError,
                    "checkDelaunay: coordinates not found in arrays.");
    return NULL;
  }
  posxc++; posyc++; poszc++; posxd++; posyd++; poszd++;

  // parcours des bars : pour une BAR si 2 points ne sont pas dans un 
  // meme triangle, insertion du milieu de la BAR
  E_Int nbars = cnc->getSize();
  E_Int ntris = cnd->getSize();
  E_Int nfld = fc->getNfld();
  E_Int npts = fc->getSize(); 
  E_Int api = fc->getApi();
  FldArrayI& cm = *(cnc->getConnect(0));
  FldArrayI& cmd = *(cnd->getConnect(0));

  E_Float* xc = fc->begin(posxc); E_Float* xd = fd->begin(posxd);
  E_Float* yc = fc->begin(posyc); E_Float* yd = fd->begin(posyd);
  E_Float* zc = fc->begin(poszc); E_Float* zd = fd->begin(poszd);

  for (E_Int tc = 0; tc < nbars; tc++) // parcours des BARS
  {
    E_Int ic1 = cm(tc, 1)-1; E_Int ic2 = cm(tc, 2)-1;
    
    E_Float xc1 = xc[ic1]; E_Float yc1 = yc[ic1]; E_Float zc1 = zc[ic1];
    E_Float xc2 = xc[ic2]; E_Float yc2 = yc[ic2]; E_Float zc2 = zc[ic2];
    short found = 0;
   
    for (E_Int td = 0; td < ntris; td++)
    {
      E_Int it1 = cmd(td, 1) - 1; 
      E_Int it2 = cmd(td, 2) - 1; 
      E_Int it3 = cmd(td, 3) - 1;
      E_Float xd1 = xd[it1]; E_Float yd1 = yd[it1]; E_Float zd1 = zd[it1];
      E_Float xd2 = xd[it2]; E_Float yd2 = yd[it2]; E_Float zd2 = zd[it2];
      E_Float xd3 = xd[it3]; E_Float yd3 = yd[it3]; E_Float zd3 = zd[it3];

      // test coincidence 1-1
      E_Float dist11 = (xd1-xc1)*(xd1-xc1) + (yd1-yc1)*(yd1-yc1) + (zd1-zc1)*(zd1-zc1);
      E_Float dist22 = (xd2-xc2)*(xd2-xc2) + (yd2-yc2)*(yd2-yc2) + (zd2-zc2)*(zd2-zc2);
      E_Float dist23 = (xd3-xc2)*(xd3-xc2) + (yd3-yc2)*(yd3-yc2) + (zd3-zc2)*(zd3-zc2);     
      if (K_FUNC::fEqualZero(dist11, eps))
      {
        if (K_FUNC::fEqualZero(dist22, eps) || 
            K_FUNC::fEqualZero(dist23, eps))
        {found = 1; break;} //coincidence des 2 pts
      }

      // test coincidence 1-2
      E_Float dist12 = (xd2-xc1)*(xd2-xc1) + (yd2-yc1)*(yd2-yc1) + (zd2-zc1)*(zd2-zc1);
      E_Float dist21 = (xd1-xc2)*(xd1-xc2) + (yd1-yc2)*(yd1-yc2) + (zd1-zc2)*(zd1-zc2);

      if (K_FUNC::fEqualZero(dist12, eps))
      {
        if (K_FUNC::fEqualZero(dist21, eps) || 
            K_FUNC::fEqualZero(dist23, eps))
        {found = 1; break;} //coincidence des 2 points
      }

      // test coincidence 1-3
      E_Float dist13 = (xd3-xc1)*(xd3-xc1) + (yd3-yc1)*(yd3-yc1) + (zd3-zc1)*(zd3-zc1);
      if (K_FUNC::fEqualZero(dist13, eps))
      {
        if (K_FUNC::fEqualZero(dist21, eps) || 
            K_FUNC::fEqualZero(dist22, eps))
        {found = 1; break;} //coincidence des 2 pts
      }

    }//parcours des triangles
    
    if ( found == 0 ) // pas de coincidence entre une arete de tri et BAR
    { 
      PyObject* tpl = K_ARRAY::buildArray3(nfld, varStringc, npts+1, nbars+1, eltTypec, 0, api);
      FldArrayF* fc2; FldArrayI* cnc2;
      K_ARRAY::getFromArray3(tpl, fc2, cnc2);
      FldArrayI& cm2 = *(cnc2->getConnect(0));

      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        E_Float* fc2p0 = fc2->begin(eq);
        E_Float* fcp = fc->begin(eq);
        for (E_Int i = 0; i < npts; i++) fc2p0[i] = fcp[i];
      }

      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        E_Float* fc2p0 = fc2->begin(eq);
        E_Float* fcp = fc->begin(eq);
        for (E_Int i = 0; i < npts; i++) fc2p0[i] = fcp[i];
        fc2p0[npts] = 0.5*(fcp[ic1]+fcp[ic2]); 
      }

      // avant la bar...
      for (E_Int tc2 = 0; tc2 < tc; tc2++)
      {
        E_Int ind1 = cm(tc2,1); E_Int ind2 = cm(tc2,2);
        cm2(tc2,1) = ind1; cm2(tc2,2) = ind2;
      }
      // la bar...
      cm2(tc,1) = cm(tc,1); cm2(tc,2) = npts+1;
      cm2(tc+1,1) = npts+1; cm2(tc+1,2) = cm(tc,2);
      // apres la bar...
      for (E_Int tc2 = tc+1; tc2 < nbars; tc2++)
      {
        E_Int tc2p = tc2+1;
        E_Int ind1 = cm(tc2,1); E_Int ind2 = cm(tc2,2);
        cm2(tc2p,1) = ind1; cm2(tc2p,2) = ind2;
      }
      RELEASESHAREDU(tpl, fc2, cnc2);
      RELEASESHAREDU(arrayc, fc, cnc);
      RELEASESHAREDU(arrayd, fd, cnd);
      return tpl;
    }
  } //parcours des bars

  PyObject* tpl = K_ARRAY::buildArray3(*fc, varStringc, *cnc, eltTypec, api);
  RELEASESHAREDU(arrayc, fc, cnc);
  RELEASESHAREDU(arrayd, fd, cnd);
  return tpl;
}
// ===================== Generator/checkDelaunay.cpp =========================
