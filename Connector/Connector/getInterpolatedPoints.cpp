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
/* Search for the fringe of interpolated nodes near blanked points; depth is 
   the number of layers of interpolated nodes. 
   IN: blankedCells: -1, point masque, 0: point interpole, 1, point normal.
   IN/OUT: cellN: -1, point masque, 0, point interpole, 1, point normal.*/
//=============================================================================
void K_CONNECTOR::searchMaskInterpolatedNodesUnstr(
  E_Int depth,  FldArrayI& cnEV,
  FldArrayI& blankedCells,
  FldArrayI& cellN)
{
  E_Int nvert = blankedCells.getSize();
  std::vector< std::vector<E_Int> > cVN(nvert);
  K_CONNECT::connectEV2VNbrs(cnEV, cVN);
  E_Int nvoisins;

  for (E_Int ind = 0; ind < nvert; ind++)
  {
    if (blankedCells[ind]  == -1) 
    {
      std::vector<E_Int>& voisins = cVN[ind];
      nvoisins = voisins.size();
      for (E_Int nov = 0; nov < nvoisins; nov++)
      {
        E_Int indv = voisins[nov]-1;
        cellN[indv] = K_FUNC::E_min(0,cellN[indv]);
      }
    }
  }
  FldArrayIS tag(nvert); tag.setAllValuesAtNull();
  for (E_Int d = 2; d<= depth; d++)
  {
    for (E_Int ind = 0; ind<nvert; ind++)
    {
      if (cellN[ind] == 0)// pt interpole 
      {
        std::vector<E_Int>& voisins = cVN[ind];
        nvoisins = voisins.size();
        for (E_Int nov = 0; nov < nvoisins; nov++)
        {
          E_Int indv = voisins[nov]-1;
          if (cellN[indv] == 1) tag[indv] = 1;
        }
      }
    }
  }
  for (E_Int ind = 0; ind < nvert; ind++)
  {if (tag[ind] == 1) cellN[ind] = 0;}
  return;
}

//=============================================================================
/* Search for the fringe of interpolated cells near blanked points; depth is 
   the number of layers of interpolated cells. 
   IN: blankedCells: -1, point masque, 0 : point interpole, 1, point normal.
   IN/OUT: cellN: -1, point masque, 0, point interpole, 1, point normal.*/
//=============================================================================
void K_CONNECTOR::searchMaskInterpolatedCellsNGON(E_Int depth, FldArrayI& cNG,
                                                  FldArrayI& blankedCells,
                                                  FldArrayI& cellN)
{
  FldArrayI cFE;
  E_Int* cnp = cNG.begin();       
  E_Int sizeFN = cnp[1];         // taille de la connectivite face/noeuds
  E_Int nelts = cnp[sizeFN+2];         // nombre d elements       
  std::vector< std::vector<E_Int> > cEEN(nelts);
  K_CONNECT::connectNG2FE(cNG, cFE);
  K_CONNECT::connectFE2EENbrs(cFE, cEEN);
  E_Int nvoisins;

  //1st layer, depth = 1
  for (E_Int et = 0; et < nelts; et++)
  {
    if (blankedCells[et] == -1)// pt masque 
    {
      std::vector<E_Int>& voisins = cEEN[et];
      nvoisins = voisins.size();
      for (E_Int noev = 0; noev < nvoisins; noev++)
      {
        E_Int et2 = voisins[noev];
        cellN[et2] = K_FUNC::E_min(0,cellN[et2]);
      }
    }
  }

  FldArrayIS tag(nelts); tag.setAllValuesAtNull();
  for (E_Int d = 2; d<= depth; d++)
  {
    for (E_Int et = 0; et < nelts; et++)
    {
      if (cellN[et] == 0)// pt interpole 
      {
        std::vector<E_Int>& voisins = cEEN[et];
        nvoisins = voisins.size();
        for (E_Int noev = 0; noev < nvoisins; noev++)
        {
          E_Int et2 = voisins[noev];
          if (cellN[et2] == 1) tag[et2] = 1;
        }
      }
    }
  }
  for (E_Int et = 0; et < nelts; et++)
  { if ( tag[et] == 1) cellN[et] = 0; }
}
//=============================================================================
/* Search for the fringe of interpolated cells near blanked points; depth is 
   the number of layers of interpolated cells. 
   IN: blankedCells: -1, point masque, 0 : point interpole, 1, point normal.
   IN/OUT: cellN: -1, point masque, 0, point interpole, 1, point normal.*/
//=============================================================================
void K_CONNECTOR::searchMaskInterpolatedCellsUnstr(char* eltType, 
                                                   E_Int depth, FldArrayI& cnEV,
                                                   FldArrayI& blankedCells,
                                                   FldArrayI& cellN)
{
  E_Int nelts = cnEV.getSize();
  E_Int nvert = nelts*cnEV.getNfld();
  std::vector< std::vector<E_Int> > cEEN(nelts);
  K_CONNECT::connectEV2EENbrs(eltType, nvert, cnEV, cEEN); 
                       
  E_Int nvoisins;

  //1st layer, depth = 1
  for (E_Int et = 0; et < nelts; et++)
  {
    if (blankedCells[et] == -1)// pt masque 
    {
      std::vector<E_Int>& voisins = cEEN[et];
      nvoisins = voisins.size();
      for (E_Int noev = 0; noev < nvoisins; noev++)
      {
        E_Int et2 = voisins[noev];
        cellN[et2] = K_FUNC::E_min(0,cellN[et2]);
      }
    }
  }

  FldArrayIS tag(nelts); tag.setAllValuesAtNull();
  for (E_Int d = 2; d<= depth; d++)
  {
    for (E_Int et = 0; et < nelts; et++)
    {
      if (cellN[et] == 0)// pt interpole 
      {
        std::vector<E_Int>& voisins = cEEN[et];
        nvoisins = voisins.size();
        for (E_Int noev = 0; noev < nvoisins; noev++)
        {
          E_Int et2 = voisins[noev];
          if (cellN[et2] == 1) tag[et2] = 1;
        }
      }
    }
  }
  for (E_Int et = 0; et < nelts; et++)
  { if ( tag[et] == 1) cellN[et] = 0; }
}

//=============================================================================
void K_CONNECTOR::searchMaskInterpolatedCellsStruct(E_Int imc, E_Int jmc, E_Int kmc, E_Int depth, E_Int dir,
                                                    FldArrayI& blankedCells, FldArrayI& cellN)
{
  E_Int imjmc = imc*jmc;
  E_Int imjmkmc = imjmc*kmc;
  E_Int i, j, k, sensor, unmsensor, ind2;
  E_Int im1, ip1, jm1, jp1, km1, kp1;
  E_Int km1imjmc, kimjmc, kp1imjmc;
  E_Int nindices;
  //On n etend que les points masques (blankedcells = -1)
  if (dir == 0) //directionnel
  {
    if (kmc == 1) 
    {
      nindices = 4;      
      vector<E_Int> indices(nindices);
      for (E_Int d = 1; d <= depth; d++)
      {
        for (E_Int ind = 0; ind < imjmc; ind++)
        {
          j = ind/imc;
          i = ind-j*imc;
          sensor = (2+blankedCells[ind])/2;
          unmsensor = 1-sensor;
              
          im1 = K_FUNC::E_max(0,i-d); ip1 = K_FUNC::E_min(i+d,imc-1);
          jm1 = K_FUNC::E_max(0,j-d); jp1 = K_FUNC::E_min(j+d,jmc-1);
          indices[0] = im1 + j*imc;
          indices[1] = ip1 + j*imc;
          indices[2] = i + jm1*imc;
          indices[3] = i + jp1*imc;      
          
          for (E_Int noi = 0; noi < nindices; noi++)
          {
            ind2 = indices[noi];
            cellN[ind2] = sensor*cellN[ind2] + unmsensor*K_FUNC::E_min(cellN[ind2],0);
          }
        }        
      }
    }// fin 2D
    else 
    {
      nindices = 6;      
      vector<E_Int> indices(nindices);
      for (E_Int d = 1; d <= depth; d++)
      {
        for (E_Int ind = 0; ind < imjmkmc; ind++)
        {
          k = ind/imjmc;
          j = ( ind-k*imjmc )/imc;
          i = ind-k*imjmc-j*imc;  
          sensor = (2+blankedCells[ind])/2;
          unmsensor = 1-sensor;
              
          im1 = K_FUNC::E_max(0,i-d); ip1 = K_FUNC::E_min(i+d,imc-1);
          jm1 = K_FUNC::E_max(0,j-d); jp1 = K_FUNC::E_min(j+d,jmc-1);
          km1 = K_FUNC::E_max(0,k-d); kp1 = K_FUNC::E_min(k+d,kmc-1);

          indices[0] = im1 + j*imc + k*imjmc;
          indices[1] = ip1 + j*imc + k*imjmc;
          indices[2] = i + jm1*imc + k*imjmc;
          indices[3] = i + jp1*imc + k*imjmc;      
          indices[4] = i + j*imc + km1*imjmc;      
          indices[5] = i + j*imc + kp1*imjmc;

          for (E_Int noi = 0; noi < nindices; noi++)
          {
            ind2 = indices[noi];
            cellN[ind2] = sensor*cellN[ind2] + unmsensor*K_FUNC::E_min(cellN[ind2],0);
          }
        }        
      }
    }//fin 3D dir = 0
  }//dir = 0
  else 
  {
    if (kmc == 1) 
    {
      nindices = 8;      
      vector<E_Int> indices(nindices);
      for (E_Int d = 1; d <= depth; d++)
      {
        for (E_Int ind = 0; ind < imjmc; ind++)
        {
          j = ind/imc;
          i = ind-j*imc;
          sensor = (2+blankedCells[ind])/2;
          unmsensor = 1-sensor;
              
          im1 = K_FUNC::E_max(0,i-d); ip1 = K_FUNC::E_min(i+d,imc-1);
          jm1 = K_FUNC::E_max(0,j-d); jp1 = K_FUNC::E_min(j+d,jmc-1);
          indices[0] = im1 + jm1*imc;
          indices[1] = i + jm1*imc;
          indices[2] = ip1 + jm1*imc;
          
          indices[3] = im1 + j*imc;
          indices[4] = ip1 + j*imc;

          indices[5] = im1 + jp1*imc;
          indices[6] = i  +  jp1*imc;
          indices[7] = ip1 + jp1*imc;
          
          for (E_Int noi = 0; noi < nindices; noi++)
          {
            ind2 = indices[noi];
            cellN[ind2] = sensor*cellN[ind2] + unmsensor*K_FUNC::E_min(cellN[ind2],0);
          }
        }        
      }
    }// 2D dir = 1
    else // 3D 
    {
      nindices = 26;
      vector<E_Int> indices(nindices);
      for (E_Int d = 1; d <= depth; d++)
      {
        for (E_Int ind = 0; ind < imjmkmc; ind++)
        {
          k = ind/imjmc;
          j = ( ind-k*imjmc )/imc;
          i = ind-k*imjmc-j*imc;  
          sensor = (2+blankedCells[ind])/2;
          unmsensor = 1-sensor;
          
          im1 = K_FUNC::E_max(0,i-d); ip1 = K_FUNC::E_min(i+d,imc-1);
          jm1 = K_FUNC::E_max(0,j-d); jp1 = K_FUNC::E_min(j+d,jmc-1);
          km1 = K_FUNC::E_max(0,k-d); kp1 = K_FUNC::E_min(k+d,kmc-1);

          km1imjmc= km1*imjmc;
          kp1imjmc= kp1*imjmc;
          kimjmc= k*imjmc;
          
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
            ind2 = indices[noi];
            cellN[ind2] = sensor*cellN[ind2] + unmsensor*K_FUNC::E_min(cellN[ind2],0);
          }
        }
      }
    }
  }//dir = 1
}

//=============================================================================
/* Determine les noeuds interpoles a partir du cellN en noeuds
   Si le celln contient des pts masques, alors les points interpoles autour 
   sont construits */
//=============================================================================
PyObject* K_CONNECTOR::getOversetHolesInterpNodes(PyObject* self, PyObject* args)
{
  PyObject *array;
  E_Int depth; E_Int dir;
  char* cellNName;
  if (!PYPARSETUPLEI(args,"Olls", "Oiis", &array, &depth, &dir, &cellNName))
  {
      return NULL;
  }
  if (dir != 0 && dir != 1) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "getOversetHolesInterpNodes: dir must be 0 or 1.");
    return NULL;
  }
  /*--------------------------------------------------*/
  /* Extraction des infos sur le domaine a interpoler */
  /*--------------------------------------------------*/
  E_Int im, jm, km;
  FldArrayF* field; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, 
                                     field, im, jm, km, cn, eltType); 
  if (res != 1 && res != 2)
  {    
    PyErr_SetString(PyExc_TypeError,
                    "getOversetHolesInterpNodes: first argument is not recognized");
    return NULL;
  }

  E_Int posc;
  if (strcmp(cellNName, "cellN") == 0)
    posc = K_ARRAY::isCellNatureField2Present(varString);
  else posc = K_ARRAY::isNamePresent(cellNName, varString);

  if (posc == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getOversetHolesInterpNodes: array must contain cellN variable.");
    delete field; if (res == 2) delete cn; return NULL;
  }
  posc++;

  E_Float* cellNp = field->begin(posc);
  /* Fin des verifs */
  E_Int npts = field->getSize();
  FldArrayI blankedCells(npts); blankedCells.setAllValuesAt(1);
  FldArrayI cellNatFld(npts); cellNatFld.setAllValuesAt(1);
  for (E_Int ind = 0; ind < npts; ind++)
  {
    if (cellNp[ind] == 2.){ blankedCells[ind] = 0; cellNatFld[ind] = 0;}
    else if (cellNp[ind] == 0.){ blankedCells[ind] = -1; cellNatFld[ind] = -1;}
  }
  if (res == 1) 
  {
    searchMaskInterpolatedCellsStruct(im, jm, km, depth, dir, blankedCells, cellNatFld);
    for (E_Int ind = 0; ind < npts; ind++)
    {
      if (cellNatFld[ind] == 0) cellNp[ind] = 2.; 
      else if (cellNatFld[ind] == -1) cellNp[ind] = 0.; 
    }
    
    PyObject* tpl =  K_ARRAY::buildArray(*field, varString, im, jm, km);
    delete field; return tpl;
  }
  else 
  {
    if ( K_STRING::cmp(eltType,"NGON")==0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "getOversetHolesInterpNodes: not implemented for NGON zones.");
      delete field; delete cn; return NULL;
    }
    searchMaskInterpolatedNodesUnstr(depth, *cn, blankedCells, cellNatFld);
    for (E_Int ind = 0; ind < npts; ind++)
    {
      if (cellNatFld[ind] == 0) cellNp[ind] = 2.; 
      else if (cellNatFld[ind] == -1) cellNp[ind] = 0.; 
    }
    
    PyObject* tpl =  K_ARRAY::buildArray(*field, varString, *cn, -1, eltType);
    delete field; delete cn; return tpl;
  }
}
//=============================================================================
/* Determine les noeuds interpoles a partir du cellN en noeuds
   Si le celln contient des pts masques, alors les points interpoles autour 
   sont construits */
//=============================================================================
PyObject* K_CONNECTOR::_getOversetHolesInterpNodes(PyObject* self, PyObject* args)
{
  PyObject *array;
  E_Int depth; E_Int dir;
  char* cellNName;
  if (!PYPARSETUPLEI(args,"Olls", "Oiis", &array, &depth, &dir, &cellNName))
  {
      return NULL;
  }
  if (dir != 0 && dir != 1) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "_getOversetHolesInterpNodes: dir must be 0 or 1.");
    return NULL;
  }
  /*--------------------------------------------------*/
  /* Extraction des infos sur le domaine a interpoler */
  /*--------------------------------------------------*/
  E_Int im, jm, km;
  FldArrayF* field; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray2(array, varString, 
                                     field, im, jm, km, cn, eltType); 
  if (res != 1)
  {    
    if ( res == 2)
    {
      PyErr_SetString(PyExc_TypeError,
                      "_getOversetHolesInterpNodes: not yet implemented for unstructured zones.");
      RELEASESHAREDU(array, field, cn);
    }
    else
      PyErr_SetString(PyExc_TypeError,
                      "getOversetHolesInterpNodes: first argument is not recognized");
    return NULL;
  }

  E_Int posc;
  if (strcmp(cellNName, "cellN") == 0)
    posc = K_ARRAY::isCellNatureField2Present(varString);
  else posc = K_ARRAY::isNamePresent(cellNName, varString);

  if (posc == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "_getOversetHolesInterpNodes: array must contain cellN variable.");
    RELEASESHAREDS(array, field);return NULL;
  }
  posc++;

  E_Float* cellNp = field->begin(posc);
  /* Fin des verifs */
  E_Int npts = field->getSize();
  FldArrayI blankedCells(npts); blankedCells.setAllValuesAt(1);
  FldArrayI cellNatFld(npts); cellNatFld.setAllValuesAt(1);
  for (E_Int ind = 0; ind < npts; ind++)
  {
    if (cellNp[ind] == 2.){ blankedCells[ind] = 0; cellNatFld[ind] = 0;}
    else if (cellNp[ind] == 0.){ blankedCells[ind] = -1; cellNatFld[ind] = -1;}
  }

  searchMaskInterpolatedCellsStruct(im, jm, km, depth, dir, blankedCells, cellNatFld);

#pragma omp parallel shared(cellNp, cellNatFld)
  {
# pragma omp for
  for (E_Int ind = 0; ind < npts; ind++)
  {
    if (cellNatFld[ind] == 0) cellNp[ind] = 2.; 
    else if (cellNatFld[ind] == -1) cellNp[ind] = 0.; 
  }
  }
  RELEASESHAREDS(array, field);
  Py_INCREF(Py_None);
  return Py_None;
}
//=============================================================================
/* Determine les centres interpoles a partir du cellN 
   Si le celln contient des pts masques, alors les points interpoles autour 
   sont construits */
//=============================================================================
PyObject* K_CONNECTOR::_getOversetHolesInterpCellCenters(PyObject* self, PyObject* args)
{
  PyObject *centersArray;
  E_Int depth; E_Int dir;
  char* cellNName;
  if (!PYPARSETUPLEI(args,"Olls", "Oiis",
                     &centersArray, &depth, &dir, &cellNName))
  {
      return NULL;
  }

  if (dir != 0 && dir != 1) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "getOversetHolesInterpNodes: dir must be 0 or 1.");
    return NULL;
  }
  /*--------------------------------------------------*/
  /* Extraction des infos sur le domaine a interpoler */
  /*--------------------------------------------------*/
  E_Int im, jm, km;
  FldArrayF* field; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray2(centersArray, varString, 
                                     field, im, jm, km, cn, eltType); 
  if (res != 1)
  {    
    if ( res == 2 )
    {
      PyErr_SetString(PyExc_TypeError,
                      "_getOversetHolesInterpNodes: not yet implemented for unstructured zones.");
      RELEASESHAREDU(centersArray, field, cn);
    }
    else 
      PyErr_SetString(PyExc_TypeError,
                      "_getOversetHolesInterpCellCenters: first argument is not recognized");
    return NULL;
  }

  E_Int posc;
  if (strcmp(cellNName, "cellN") == 0)
    posc = K_ARRAY::isCellNatureField2Present(varString);
  else posc = K_ARRAY::isNamePresent(cellNName, varString);
  if (posc == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "_getOversetHolesInterpCellCenters: array must contain cellN variable.");
    RELEASESHAREDS(centersArray, field);return NULL;
  }
  posc++;
  E_Float* cellNp = field->begin(posc);
  /* Fin des verifs */
  E_Int ncells = field->getSize();
  FldArrayI blankedCells(ncells); blankedCells.setAllValuesAt(1);
  FldArrayI cellNatFld(ncells); cellNatFld.setAllValuesAt(1);
  for (E_Int ind = 0; ind < ncells; ind++)
  {
    if (cellNp[ind] == 2.){ blankedCells[ind] = 0; cellNatFld[ind] = 0;}
    else if (cellNp[ind] == 0.){ blankedCells[ind] = -1; cellNatFld[ind] = -1;}
  }

  searchMaskInterpolatedCellsStruct(im, jm, km, depth, dir, blankedCells, cellNatFld);

#pragma omp parallel shared(cellNp, cellNatFld, ncells)
  {
# pragma omp for
  for (E_Int ind = 0; ind < ncells; ind++)
  {
    if (cellNatFld[ind] == 0) cellNp[ind] = 2.; 
    else if (cellNatFld[ind] == -1) cellNp[ind] = 0.; 
  }
  }
  RELEASESHAREDS(centersArray, field);
  Py_INCREF(Py_None);
  return Py_None;
}
//=============================================================================
/* Determine les centres interpoles a partir du cellN 
   Si le celln contient des pts masques, alors les points interpoles autour 
   sont construits */
//=============================================================================
PyObject* K_CONNECTOR::getOversetHolesInterpCellCenters(PyObject* self, PyObject* args)
{
  PyObject *centersArray;
  E_Int depth; E_Int dir;
  char* cellNName;
  if (!PYPARSETUPLEI(args,
                    "Olls", "Oiis",
                    &centersArray, &depth, &dir, &cellNName))
  {
      return NULL;
  }

  if (dir != 0 && dir != 1) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "getOversetHolesInterpNodes: dir must be 0 or 1.");
    return NULL;
  }
  /*--------------------------------------------------*/
  /* Extraction des infos sur le domaine a interpoler */
  /*--------------------------------------------------*/
  E_Int im, jm, km;
  FldArrayF* field; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(centersArray, varString, 
                                    field, im, jm, km, cn, eltType); 
  if (res != 1 && res != 2)
  {    
    PyErr_SetString(PyExc_TypeError,
                    "getOversetHolesInterpCellCenters:  first argument is not recognized");
    return NULL;
  }

  E_Int posc;
  if (strcmp(cellNName, "cellN") == 0)
    posc = K_ARRAY::isCellNatureField2Present(varString);
  else posc = K_ARRAY::isNamePresent(cellNName, varString);
  if (posc == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getOversetHolesInterpCellCenters: array must contain cellN variable.");
    delete field; if ( res == 2) delete cn; return NULL;
  }
  posc++;
  E_Float* cellNp = field->begin(posc);
  /* Fin des verifs */
  E_Int ncells = field->getSize();
  FldArrayI blankedCells(ncells); blankedCells.setAllValuesAt(1);
  FldArrayI cellNatFld(ncells); cellNatFld.setAllValuesAt(1);
  for (E_Int ind = 0; ind < ncells; ind++)
  {
    if (cellNp[ind] == 2.){ blankedCells[ind] = 0; cellNatFld[ind] = 0;}
    else if (cellNp[ind] == 0.){ blankedCells[ind] = -1; cellNatFld[ind] = -1;}
  }
  if (res == 1) 
  {
    searchMaskInterpolatedCellsStruct(im, jm, km, depth, dir, blankedCells, cellNatFld);
    for (E_Int ind = 0; ind < ncells; ind++)
    {
      if (cellNatFld[ind] == 0) cellNp[ind] = 2.; 
      else if (cellNatFld[ind] == -1) cellNp[ind] = 0.; 
    }
    
    PyObject* tpl =  K_ARRAY::buildArray(*field, varString, im, jm, km);
    delete field; return tpl;
  }
  else 
  {
    if ( K_STRING::cmp(eltType,"NGON*")==0)
      searchMaskInterpolatedCellsNGON(depth, *cn, blankedCells, cellNatFld);
    else
      searchMaskInterpolatedCellsUnstr(eltType, depth, *cn, blankedCells, cellNatFld);
    for (E_Int ind = 0; ind < ncells; ind++)
    {
      if (cellNatFld[ind] == 0) cellNp[ind] = 2.; 
      else if (cellNatFld[ind] == -1) cellNp[ind] = 0.; 
    }
    
    PyObject* tpl =  K_ARRAY::buildArray(*field, varString, *cn, -1, eltType);
    delete field; delete cn; return tpl;
  }
}

//===============================================================================
/* Retourne le numpy des indices des pts cellN=2 et les numpys des coordonnees 
  la zone en entree est le maillage des centres 
  car le cellN est localise aux noeuds pour plus d efficacite */
//===============================================================================
PyObject* K_CONNECTOR::getInterpolatedPointsZ(PyObject* self, PyObject* args)
{
  PyObject* zone;
  char *GridCoordinates, *FlowSolutionNodes, *FlowSolutionCenters;
  char* cellNName;
  if (!PyArg_ParseTuple(args, "Ossss", &zone, &cellNName, &GridCoordinates, &FlowSolutionNodes, &FlowSolutionCenters))
  {
    PyErr_SetString(PyExc_TypeError,
                    "getInterpolatedPointsZ: wrong arguments.");
    return NULL;
  }
  E_Int xyz = 1; E_Int locI = 0; // tjs aux noeuds
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
                    "getInterpolatedPointsZ: invalid zone.");
    return NULL;
  }
  if (locs.size() < 4)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getInterpolatedPointsZ: one variable missing in zone.");
    RELEASESHAREDZ(hook, varString, eltType);
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if ( posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getInterpolatedPointsZ: coordinates cannot be extracted from zone.");
    RELEASESHAREDZ(hook, varString, eltType);
    return NULL;
  }
  E_Int posc = K_ARRAY::isNamePresent(cellNName, varString);
  if ( posc == -1) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "getInterpolatedPointsZ: cellN cannot be extracted from zone.");
    RELEASESHAREDZ(hook, varString, eltType);
    return NULL;
  }
  if ( locs[posc] != 0 )
  {
    PyErr_SetString(PyExc_TypeError,
                    "getInterpolatedPointsZ: cellN must be located at nodes in input zone.");
    RELEASESHAREDZ(hook, varString, eltType);
    return NULL;
  }
  E_Int nptsTot;
  if ( zoneType == 1) nptsTot = ni*nj*nk;
  else nptsTot = ni;
  E_Float* cellNp = fields[posc];
  FldArrayI indicesInterp(nptsTot);
  FldArrayF coordX(nptsTot); FldArrayF coordY(nptsTot); FldArrayF coordZ(nptsTot);

  E_Int noi = 0;
  E_Float* xp = fields[posx];
  E_Float* yp = fields[posy];
  E_Float* zp = fields[posz];

  for (E_Int i = 0; i < nptsTot; i++)
  {
    if ( cellNp[i]==2. )
    {
      indicesInterp[noi] = i;
      coordX[noi] = xp[i];
      coordY[noi] = yp[i];
      coordZ[noi] = zp[i];      
      noi+=1;
    }
  }
  if ( noi == 0) 
  {
    RELEASESHAREDZ(hook, varString, eltType);
    Py_INCREF(Py_None);
    return Py_None;
  }
  indicesInterp.resize(noi);
  coordX.resize(noi);
  coordY.resize(noi);
  coordZ.resize(noi);
  RELEASESHAREDZ(hook, varString, eltType);

  PyObject* PyIndices = K_NUMPY::buildNumpyArray(indicesInterp,1);
  PyObject* PyCoordX = K_NUMPY::buildNumpyArray(coordX,1);
  PyObject* PyCoordY = K_NUMPY::buildNumpyArray(coordY,1);
  PyObject* PyCoordZ = K_NUMPY::buildNumpyArray(coordZ,1);
  PyObject* tpl = Py_BuildValue("[OOOO]", PyIndices, PyCoordX, PyCoordY, PyCoordZ);
  
  Py_DECREF(PyIndices); Py_DECREF(PyCoordX); Py_DECREF(PyCoordY); Py_DECREF(PyCoordZ); 
  return tpl;
}
//=============================================================================
/**/
//=============================================================================
PyObject* K_CONNECTOR::getInterpolatedPoints(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array))
  {
    PyErr_SetString(PyExc_TypeError,
                    "getInterpolatedPoints: wrong arguments.");
    return NULL;
  }
  // Check: 
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType, true); 
  if ( res != 1 && res != 2) 
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError, 
                    "getInterpolatedPoints: invalid array.");
    return NULL;   
  }
  E_Int posc = K_ARRAY::isCellNatureField2Present(varString);
  if (posc == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getInterpolatedPoints: array must contain cellN.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  posc++;

  /*fin verifs*/
  E_Int nfld = f->getNfld();
  E_Int npts = f->getSize();
  char varStringOut[K_ARRAY::VARSTRINGLENGTH]; varStringOut[0] = '\0';
  E_Int nfldOut = nfld+1;
  strcpy(varStringOut,varString); strcat(varStringOut,",indcell");
  E_Float* cellnp = f->begin(posc);
  FldArrayF* fout = new FldArrayF(npts,nfldOut);
  E_Int c=0;
  for (E_Int ind=0; ind < npts; ind++)
  {
    if (cellnp[ind] == 2.)
    { 
      for (E_Int eq = 1; eq <= nfld; eq++)
        (*fout)(c,eq) = (*f)(ind,eq);
      (*fout)(c,nfldOut) = E_Float(ind);
      c++;
    }
  }
  fout->reAllocMat(c,nfldOut);
  RELEASESHAREDB(res, array, f, cn);
  FldArrayI* cnl = new FldArrayI(0);
  PyObject* tpl = K_ARRAY::buildArray(*fout, varStringOut, *cnl, -1, 
                                      "NODE", false);
  delete fout; delete cnl;
  return tpl;
}
