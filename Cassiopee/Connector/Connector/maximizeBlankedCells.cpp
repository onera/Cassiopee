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

// grid connectivity: maximize blanked cells region

# include "connector.h"
using namespace K_FLD;

// ============================================================================
/* Maximize the blanked cells region */
// ============================================================================
PyObject* K_CONNECTOR::maximizeBlankedCells(PyObject* self, PyObject* args)
{
  PyObject* array; E_Int depth; E_Int dir;
  char* varCellN;
  if (!PYPARSETUPLE_(args, O_ II_ S_,
                    &array, &depth, &dir, &varCellN))
  {
      return NULL;
  }
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);
  
  if (res != 1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "maximizeBlankedCells: only for structured arrays.");
    return NULL;
  }
  E_Int poscelln = K_ARRAY::isNamePresent(varCellN,varString);
  if (poscelln == -1)
  {
     PyErr_SetString(PyExc_ValueError,
                    "maximizeBlankedCells: cell nature field must be present.");
     return NULL;
  }
  
  E_Int api = f->getApi();
  E_Float* cellnp = f->begin(poscelln+1);
  
  E_Int c = 1;
  while (c != 0) c = getRidOfInterpPoints(cellnp, im, jm, km, depth, dir);

  PyObject* tpl = K_ARRAY::buildArray3(*f, varString, im, jm, km, api);
  RELEASESHAREDS(array, f);
  return tpl;
}

//=============================================================================
// Fait passer a cellN=0 les points interpoles (cellN=2) entoures par 
// des points cellN!=1
// IN: im,jm,km: taille du champ cellN
// IN: depth: nbre de layers
// IN: dir: 0 (molecule directionnelle), 1 (molecule complete) 
//==============================================================================
E_Int K_CONNECTOR::getRidOfInterpPoints(
  E_Float* celln, 
  E_Int im, E_Int jm, E_Int km, E_Int depth, E_Int dir)
{
  E_Int imjm = im*jm;
  E_Int compt = 0;
  E_Int i, j, k, ip, jp, kp, ip1, jp1, kp1;
  E_Int ind, ind1;

  if (dir == 0)
  {
    for (k = 0; k < km; k++)
      for (j = 0; j < jm; j++)
        for (i = 0; i < im; i++)
        {
          ind = i + j*im + k*imjm;
          if (celln[ind] == 2.)
          {
            for (ip = -depth; ip <= depth; ip++)
            {
              ip1 = i+ip;
              ip1 = K_FUNC::E_min(ip1, im-1);
              ip1 = K_FUNC::E_max(ip1, 0);
              ind1 = ip1 + j*im + k*imjm;
              if (celln[ind1] == 1.) {goto fin;}
            }
            for (jp = -depth; jp <= depth; jp++)
            {
              jp1 = j+jp;
              jp1 = K_FUNC::E_min(jp1, jm-1);
              jp1 = K_FUNC::E_max(jp1, 0);
              ind1 = i + jp1*im + k*imjm;
              if (celln[ind1] == 1.) {goto fin;}
            }
            for (kp = -depth; kp <= depth; kp++)
            {
              kp1 = k+kp;
              kp1 = K_FUNC::E_min(kp1, km-1);
              kp1 = K_FUNC::E_max(kp1, 0);
              ind1 = i + j*im + kp1*imjm;
              if (celln[ind1] == 1.) {goto fin;}
            }
            celln[ind] = 0.; compt++;
            fin:;
          }
        }
  }
  else // dir == 1
  {
    for (k = 0; k < km; k++)
      for (j = 0; j < jm; j++)
        for (i = 0; i < im; i++)
        {
          ind = i + j*im + k*imjm;
          if (celln[ind] == 2.)
          {
            for (ip = -depth; ip <= depth; ip++)
              for (jp = -depth; jp <= depth; jp++)
                for (kp = -depth; kp <= depth; kp++)
                {
                  ip1 = i+ip;
                  ip1 = K_FUNC::E_min(ip1, im-1);
                  ip1 = K_FUNC::E_max(ip1, 0);
                  jp1 = j+jp;
                  jp1 = K_FUNC::E_min(jp1, jm-1);
                  jp1 = K_FUNC::E_max(jp1, 0);
                  kp1 = k+kp;
                  kp1 = K_FUNC::E_min(kp1, km-1);
                  kp1 = K_FUNC::E_max(kp1, 0);
                  ind1 = ip1 + jp1*im + kp1*imjm;
                  if (celln[ind1] == 1.) {goto fin2;}
                }
            celln[ind] = 0.; compt++;
            fin2:;
          }
        }
  }
  return compt;
}
