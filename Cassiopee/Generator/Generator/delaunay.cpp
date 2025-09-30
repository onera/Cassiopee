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

// Triangulation de Delaunay a partir d'une liste de points

# include "generator.h"

using namespace std;

//=========================================================================
/* Triangulation de Delaunay a partir d'un ensemble de points */
//=========================================================================
PyObject* K_GENERATOR::delaunay(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int keepBB;
  if (!PYPARSETUPLE_(args, O_ I_, &array, &keepBB)) return NULL;
  
  // Check array
  E_Int ni, nj, nk;
  K_FLD::FldArrayF* f;
  char* varString; char* eltType;
  K_FLD::FldArrayI* cn;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, 
                                     cn, eltType);
  if (res != 1 && res != 2) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "delaunay: invalid array.");
    return NULL;
  }
  
  if (res == 1) 
  {
    if (ni != 1 && nj != 1 && nk != 1)
    {
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError,
                      "delaunay: array must define a plane.");
      return NULL;
    }
  }
  if (res == 2)
  {
    if (strcmp(eltType, "TRI") != 0 && 
        strcmp(eltType, "QUAD") != 0 &&
        strcmp(eltType, "NODE") != 0 &&
        strcmp(eltType, "BAR") != 0)
    {
      RELEASESHAREDU(array, f, cn);
      PyErr_SetString(PyExc_TypeError,
                      "delaunay: array must define a plane.");
      return NULL;
    }
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);      
    PyErr_SetString(PyExc_TypeError,
                    "delaunay: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  E_Int sizeIni = f->getSize();

  if (sizeIni < 3) 
  {
    RELEASESHAREDB(res, array, f, cn);      
    PyErr_SetString(PyExc_TypeError,
                    "delaunay: at least 3 points must be defined.");
    return NULL; 
  }

  /* Determination de l'eq. de plan */
  E_Float coefa, coefb, coefc, coefd;
  E_Int isok = compPlaneCoefficients(posx, posy, posz, *f, 
                                     coefa, coefb, coefc, coefd);
  if (isok == 0) // pas de plan
  { 
    RELEASESHAREDB(res, array, f, cn);      
    PyErr_SetString(PyExc_TypeError,
                    "delaunay: input array must be planar.");
    return NULL;
  }

  /* coordonnees */
  K_FLD::FldArrayF& field = *f;
  E_Int api = field.getApi();
  E_Int nfld = field.getNfld();
  K_FLD::FldArrayF* coordp = new K_FLD::FldArrayF(sizeIni, nfld);
  K_FLD::FldArrayF& coord = *coordp; //coord
  coord.setOneField(field, posx, 1);
  coord.setOneField(field, posy, 2);
  coord.setOneField(field, posz, 3);
  E_Int eq2 = 4;
  for (E_Int eq = 1; eq <= nfld; eq++)
  {
    if (eq != posx && eq != posy && eq != posz)
    {
      coord.setOneField(field, eq, eq2);
      eq2++;
    }
  }

  K_FLD::FldArrayI* cn2 = new K_FLD::FldArrayI();
  K_COMPGEOM::delaunay(coefa, coefb, coefc, coefd, coord, *cn2, keepBB);
  PyObject* tpl = K_ARRAY::buildArray3(coord, varString, *cn2, "TRI", api);
  delete &coord; delete cn2; 
  RELEASESHAREDB(res, array, f, cn);      
  return tpl;
}
//===========================================================================
/* Calcul des coeffs de l eq. de plan */
//===========================================================================
E_Int 
K_GENERATOR::compPlaneCoefficients(E_Int posx, E_Int posy, E_Int posz, 
                                   K_FLD::FldArrayF& f, 
                                   E_Float& alpha, E_Float& beta, 
                                   E_Float& gamma, E_Float& delta)
{
  E_Float* xt = f.begin(1);
  E_Float* yt = f.begin(2);
  E_Float* zt = f.begin(3);

  E_Int npts = f.getSize();
  E_Float eps = 1.e-12;
  if ( npts < 3 ) return 0;
  for (E_Int i = 0; i < npts; i++)
    for (E_Int j = 0; j < npts; j++)
      for (E_Int k = 0; k < npts; k++)
      {
        if ( i != j && i != k && j != k )
        { 
          E_Float xab = xt[j]-xt[i];
          E_Float yab = yt[j]-yt[i];
          E_Float zab = zt[j]-zt[i];
          E_Float xac = xt[k]-xt[i];
          E_Float yac = yt[k]-yt[i];
          E_Float zac = zt[k]-zt[i]; 

          alpha = yab * zac - yac * zab;
          beta  = xac * zab - xab * zac;
          gamma = xab * yac - xac * yab;
          
          if ( K_FUNC::fEqualZero(alpha,eps) == false || 
               K_FUNC::fEqualZero(beta, eps) == false || 
               K_FUNC::fEqualZero(gamma,eps) == false ) // ok
          {
            delta = -(alpha*xt[0]+beta*yt[0]+gamma*zt[0]);  
            return 1;
          } 
        }
      }
  // pas de points non colineaires trouves
  return 0;
}

//========================  Generator/delaunay.cpp ========================
