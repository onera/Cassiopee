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
# include "geom.h"

using namespace std;
using namespace K_FLD;

extern "C"
{
  void k6naca2_(const E_Float& e, E_Int& N,
                E_Float* x, E_Float* y, E_Float* z);
  void k6nacas4g_(E_Int& im, E_Int& ip, E_Int& it, E_Int& sharpte, 
                  E_Int& npt, E_Float* x, E_Float* y, E_Float* z, E_Float* xl);
  void k6nacas5g_(E_Int& il, E_Int& ip, E_Int& iq, E_Int& it, E_Int& sharpte, 
                  E_Int& npt, E_Float* x, E_Float* y, E_Float* z, E_Float* xl);
  void k6nacas4m_(E_Int& im, E_Int& ip, E_Int& ith, E_Int& it, E_Int& ii, E_Int& sharpte, 
                  E_Int& npt, E_Float* x, E_Float* y, E_Float* z, E_Float* xl);
}

// ============================================================================
/* Create a naca profile line of N points */
// ============================================================================
PyObject* K_GEOM::nacaMesh(PyObject* self, PyObject* args)
{
  E_Int N; E_Float e;
  E_Int im, ip, it, ith, iq, sharpte;
  if (!PYPARSETUPLE(args,
                    "dlllllll", "diiiiiii",
                    "flllllll", "fiiiiiii",
                    &e, &N,
                    &im, &ip, &it, &ith, &iq, &sharpte))
  {
      return NULL;
  }

  PyObject* tpl = NULL;
  if (im == -1) // input avec epaisseur
  {
    // Data check
    if (e < 0)
    {
      PyErr_SetString(PyExc_ValueError, "naca: thickness must be positive.");
      return NULL;
    }

    if ((N/2)*2-N == 0)
    {
      printf("Warning: naca: number of points must be odd.\n");
      printf("Warning: naca: number of points set to %d.\n", N+1);
      N = N+1;
    }
  
    E_Int n = E_Int(N);
    FldArrayF coord(N, 3);
    coord.setAllValuesAtNull();
    k6naca2_(e, n, coord.begin(1), coord.begin(2), coord.begin(3));
    coord.reAllocMat(n, 3);
    tpl = K_ARRAY::buildArray(coord, "x,y,z", n, 1, 1);
  }
  else // input avec digits
  {
    if (sharpte == 1 && (N/2)*2-N == 0)
    {
      // sharp : 2N+1
      printf("Warning: naca: number of points must be odd.\n");
      printf("Warning: naca: number of points set to %d.\n", N+1);
      N = N+1;
    }
    if (sharpte == 0 && (N/2)*2-N == 1)
    {
      // pas sharp : 2N
      printf("Warning: naca: number of points must be even.\n");
      printf("Warning: naca: number of points set to %d.\n", N+1);
      N = N+1;
    }

    E_Int npt = N/2; 
    if (sharpte == 1) npt++;

    FldArrayF xl(npt);
    FldArrayF coord(N, 3);
    //printf("sharpte=%d, size N=%d, npt=%d\n",sharpte, N, npt);
    //printf("type = %d %d %d %d %d\n", im, ip, iq, it, ith);

    // determination de la serie
    if (im > -0.5 && ip > -0.5 && ith > -0.5 && it > -0.5 && iq > -0.5)
    {
      // iq used as ii
      k6nacas4m_(im, ip, ith, it, iq, sharpte, 
                 npt, coord.begin(1), coord.begin(2), coord.begin(3), xl.begin());
    }
    else if (im > -0.5 && ip > -0.5 && iq > -0.5 && it > -0.5)
    {
      // im used as il
      k6nacas5g_(im, ip, iq, it, sharpte,
                 npt, coord.begin(1), coord.begin(2), coord.begin(3), xl.begin());
    }
    else if (im > -0.5 && ip > -0.5 && it > -0.5)
    {
      k6nacas4g_(im, ip, it, sharpte,
                 npt, coord.begin(1), coord.begin(2), coord.begin(3), xl.begin());
    }
    
    tpl = K_ARRAY::buildArray(coord, "x,y,z", N, 1, 1);
  }

  return tpl;
}
