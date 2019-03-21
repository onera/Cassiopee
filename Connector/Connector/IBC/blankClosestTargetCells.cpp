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
#include "connector.h"
using namespace K_FLD;
using namespace std;

//=============================================================================
/*  Masque un point cible potentiel s il est "trop proche" de la paroi */
//=============================================================================
PyObject* K_CONNECTOR::_blankClosestTargetCells(PyObject* self, PyObject* args)
{
  PyObject* zone;
  char *GridCoordinates, *FlowSolutionNodes, *FlowSolutionCenters;
  char* cellNName;
  E_Int depth;
  if (!PYPARSETUPLEI(args, "Olssss","Oissss", &zone, &depth, &cellNName, 
                     &GridCoordinates, &FlowSolutionNodes, &FlowSolutionCenters))
  {
    PyErr_SetString(PyExc_TypeError,
                    "blankClosestTargetCellsZ: wrong arguments.");
    return NULL;
  }
  E_Int xyz = 1; E_Int locI = 1; // tjs aux centres
  E_Int ni, nj, nk, cnSize, cnNfld;
  char* varString; char* eltType;
  vector<E_Float*> fields; vector<E_Int> locs;
  vector<E_Int*> cn;
  vector<PyArrayObject*> hook;
  
  E_Int zoneType = K_PYTREE::getFromZone(zone, xyz, locI, varString, fields, locs, ni, nj, nk, 
                                         cn, cnSize, cnNfld, eltType, hook, GridCoordinates, 
                                         FlowSolutionNodes, FlowSolutionCenters);
  if ( zoneType == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "blankClosestTargetCellsZ: invalid zone.");
    return NULL;
  }   
  else if (zoneType==2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "blankClosestTargetCellsZ: zone must be structured.");
    RELEASESHAREDZ(hook, varString, eltType);
    return NULL;
  }
  if (locs.size() < 5)
  {
    PyErr_SetString(PyExc_TypeError,
                    "blankClosestTargetCellsZ: zone must contain at least 2 variables.");
    RELEASESHAREDZ(hook, varString, eltType);
    return NULL;
  }   
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  if ( posx == -1 )
  {
    PyErr_SetString(PyExc_TypeError,
                    "blankClosesTargetCellsZ: CoordinateX cannot be extracted from zone.");
    RELEASESHAREDZ(hook, varString, eltType);
    return NULL;
  }
  
  E_Int posc = K_ARRAY::isNamePresent(cellNName, varString);
  if ( posc == -1) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "blankClosestTargetCellsZ: cellN cannot be extracted from zone.");
    RELEASESHAREDZ(hook, varString, eltType);
    return NULL;
  }    
  E_Int posd = K_ARRAY::isNamePresent("TurbulentDistance",varString);
  if ( posd == -1) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "blankClosestTargetCellsZ: TurbulentDistance cannot be extracted from zone.");
    RELEASESHAREDZ(hook, varString, eltType);
    return NULL;
  }    
  E_Float* xp = fields[posx];
  E_Float DX = xp[1]-xp[0];//resolution du maillage cartesien
  E_Float distmax = DX*2.505;//
  //E_Float distmin = DX*0.495;//

  E_Float* cellNp = fields[posc];
  E_Float* distp = fields[posd];
  
  E_Int imc = ni-1; E_Int jmc = nj-1; E_Int kmc = nk-1;
  E_Int imcjmc = imc*jmc;
  //E_Int imcjmckmc = imcjmc*kmc;
  E_Int ncells = imc*jmc*kmc;
  FldArrayI tag(ncells); tag.setAllValuesAtNull();
  E_Int* tagp = tag.begin();
  if ( kmc == 1)//2D
  {
    FldArrayI indicesV(8); E_Int* indices = indicesV.begin();
    for (E_Int ind = 0; ind < ncells; ind++)
    {
      if (cellNp[ind]==0. && tagp[ind]==0)
      {   
        E_Int j = ind/imc;
        E_Int i = ind-j*imc;
        E_Int im1 = i-depth; E_Int ip1 = i+depth;
        E_Int jm1 = j-depth; E_Int jp1 = j+depth;
        indices[0] = im1 + jm1*imc;
        indices[1] = i + jm1*imc;
        indices[2] = ip1 + jm1*imc;

        indices[3] = im1 + j*imc;
        indices[4] = ip1 + j*imc;

        indices[5] = im1 + jp1*imc;
        indices[6] = i  +  jp1*imc;
        indices[7] = ip1 + jp1*imc;

        for (E_Int noi = 0; noi < 8; noi++)
        {
          E_Int ind2 = indices[noi];
          if ( ind2 >=0 && ind2 < ncells)
          {
            if (cellNp[ind2]==2. && tagp[ind2]==0)
            {
              if ( distp[ind2] > distmax) 
                { tagp[ind2]= 1; cellNp[ind2]=1.;}
            }
          }
        }
      }
    }
  }// fin 2D
  else// 3D
  {
    E_Int nindices=26;
    FldArrayI indicesV(nindices); E_Int* indices = indicesV.begin();
    for (E_Int ind = 0; ind < ncells; ind++)
    {
      if ( cellNp[ind]==0. && tagp[ind]==0)
      {
        E_Int k = ind/imcjmc;
        E_Int j = ( ind-k*imcjmc )/imc;
        E_Int i = ind-k*imcjmc-j*imc;
        E_Int im1 = i-depth; E_Int ip1 = i+depth;
        E_Int jm1 = j-depth; E_Int jp1 = j+depth;
        E_Int km1 = k-depth; E_Int kp1 = k+depth;

        E_Int km1imjmc = km1*imcjmc;
        E_Int kp1imjmc = kp1*imcjmc;
        E_Int kimjmc   = k*imcjmc;

        indices[0] = im1 + jm1*imc + km1imjmc;
        indices[1] = i   + jm1*imc + km1imjmc;
        indices[2] = ip1 + jm1*imc + km1imjmc;

        indices[3] = im1 + j*imc + km1imjmc;
        indices[4] = i   + j*imc + km1imjmc;
        indices[5] = ip1 + j*imc + km1imjmc;

        indices[6] = im1 + jp1*imc + km1imjmc;
        indices[7] = i  +  jp1*imc + km1imjmc;
        indices[8] = ip1 + jp1*imc + km1imjmc;

        indices[9]  = im1 + jm1*imc + kimjmc;
        indices[10] = i   + jm1*imc + kimjmc;
        indices[11] = ip1 + jm1*imc + kimjmc;

        indices[12] = im1 + j*imc + kimjmc;
        indices[13] = ip1 + j*imc + kimjmc;

        indices[14] = im1 + jp1*imc + kimjmc;
        indices[15] = i  +  jp1*imc + kimjmc;
        indices[16] = ip1 + jp1*imc + kimjmc;

        indices[17] = im1 + jm1*imc + kp1imjmc;
        indices[18] = i   + jm1*imc + kp1imjmc;
        indices[19] = ip1 + jm1*imc + kp1imjmc;

        indices[20] = im1 + j*imc + kp1imjmc;
        indices[21] = i   + j*imc + kp1imjmc;
        indices[22] = ip1 + j*imc + kp1imjmc;

        indices[23] = im1 + jp1*imc + kp1imjmc;
        indices[24] = i  +  jp1*imc + kp1imjmc;
        indices[25] = ip1 + jp1*imc + kp1imjmc;

        for (E_Int noi = 0; noi < nindices; noi++)
        {
          E_Int ind2 = indices[noi];
          if ( ind2 >=0 && ind2 < ncells)
          {
            if (cellNp[ind2]==2. && tagp[ind2]==0)
            {
              if ( distp[ind2] > distmax) 
                { tagp[ind2]= 1; cellNp[ind2]=1.;}
            }
          }
        }
      }
    }
  } // fin 3D

  RELEASESHAREDZ(hook, varString, eltType);
  return Py_None;
}
