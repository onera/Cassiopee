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
# include "connector.h"
using namespace std;
using namespace K_FLD;

//=============================================================================
/* Computes the sum of absolute values of coefficients for extrapolated cells*/
//=============================================================================
PyObject* K_CONNECTOR::getExtrapAbsCoefs(PyObject* self, PyObject* args)
{
  PyObject *pyIndRcv;
  PyObject *pyIndExtrap;
  PyObject *pyArrayTypes;
  PyObject *pyArrayCoefs;
  if (!PYPARSETUPLE_(args, OOO_ O_,
                      &pyIndRcv, &pyIndExtrap, &pyArrayTypes, &pyArrayCoefs))
  {
    return NULL;
  }

  /*--------------------------------------*/
  /* Extraction des indices des receveurs */
  /*--------------------------------------*/
  FldArrayI* rcvPtsI;
  E_Int res = K_NUMPY::getFromNumpyArray(pyIndRcv, rcvPtsI);
  if (res == 0) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "getExtrapCoefs: 1st arg must be a numpy of integers.");
    return NULL;
  }
  /*---------------------------------------*/
  /* Extraction des indices des extrapoles */
  /*---------------------------------------*/
  FldArrayI* extrapPtsI;
  res = K_NUMPY::getFromNumpyArray(pyIndExtrap, extrapPtsI);
  if (res == 0) 
  {
    RELEASESHAREDN(pyIndRcv, rcvPtsI);
    PyErr_SetString(PyExc_TypeError, 
                    "getExtrapCoefs: 1st arg must be a numpy of integers.");
    return NULL;
  }
  
  /*----------------------*/
  /* Extraction des types */
  /*----------------------*/
  FldArrayI* typesI;
  res = K_NUMPY::getFromNumpyArray(pyArrayTypes, typesI);
  if (res == 0) 
  {
    RELEASESHAREDN(pyIndRcv, rcvPtsI);
    RELEASESHAREDN(pyIndExtrap, extrapPtsI);
    PyErr_SetString(PyExc_TypeError, 
                    "getExtrapCoefs: 3rd arg must be a numpy of integers.");
    return NULL;
  }

  /*-----------------------*/
  /* Extraction des coefs  */
  /*-----------------------*/
  FldArrayF* donorCoefsF;
  res = K_NUMPY::getFromNumpyArray(pyArrayCoefs, donorCoefsF);
  if (res == 0) 
  {
    RELEASESHAREDN(pyIndRcv, rcvPtsI);
    RELEASESHAREDN(pyIndExtrap, extrapPtsI);
    RELEASESHAREDN(pyArrayTypes, typesI);
    PyErr_SetString(PyExc_TypeError, 
                    "getExtrapCoefs: 4th arg must be a numpy of floats.");
    return NULL;
  }
  E_Int sizecoefs = 0;
  // Types valides: 2, 3, 4, 5, 102
  E_Float* ptrCoefs = donorCoefsF->begin();
  E_Int* rcvPts = rcvPtsI->begin();
  E_Int* extrapPts = extrapPtsI->begin();
  E_Int* types = typesI->begin();
  E_Int nbExtrapPts = extrapPtsI->getSize();
  E_Int nbRcvPts = rcvPtsI->getSize();
  //E_Float* ptrCF = donorCoefsF->begin();
  
  vector<E_Int> indices;
  for (E_Int noinde = 0; noinde < nbExtrapPts; noinde++)
    for (E_Int noind = 0; noind < nbRcvPts; noind++)
    {
      if (extrapPts[noinde] == rcvPts[noind])  
      { indices.push_back(noind); break; }
    }
  E_Int nindices = indices.size();  
    
  PyObject* tpl = K_ARRAY::buildArray3(1, "extrapolated", nindices, 1, 1, 1);
  FldArrayF* sumCf;
  K_ARRAY::getFromArray3(tpl, sumCf);
  E_Float* sumCfp = sumCf->begin();

  for (E_Int noe = 0; noe < nindices; noe++)
  {        
    E_Int noind = indices[noe];
    E_Int type = types[noind];
    E_Float sum = 0.;
    switch (type)
    {
      case 22: // Structure Lineaire O2 par tetra
        sizecoefs = 4; sum = 0.; 
        for (E_Int ii=0; ii<4; ii++)
        {sum += K_FUNC::E_abs(ptrCoefs[ii]);}    
        sumCfp[noe] = sum;
        break;

      case 1: // 1 seul coeff
        sumCfp[noe] = K_FUNC::E_abs(ptrCoefs[0]);
        break;

      case 100: //setInterpolations
      case 102: //setInterpolations
      case 103: //setInterpolations
        sizecoefs = 8; sum = 0.;  
        for (E_Int ii=1; ii<=8; ii++)
        {sum += K_FUNC::E_abs((*donorCoefsF)(noind,ii));}    
        sumCfp[noe] = sum;
        break;

  
      case 2: // Structure Lineaire O2 par tetra
        sizecoefs = 8; sum = 0.; 
        for (E_Int ii=0; ii<8; ii++)
        {sum += K_FUNC::E_abs(ptrCoefs[ii]);}    
        sumCfp[noe] = sum;
        break;
        
      case 3: // Lagrange O3 
        sum = 0.; sizecoefs = 9;
        for (E_Int kk=0; kk<3; kk++)
          for (E_Int jj=0; jj<3; jj++)
            for (E_Int ii=0; ii<3; ii++)
            {sum += K_FUNC::E_abs(ptrCoefs[ii]*ptrCoefs[jj+3]*ptrCoefs[kk+6]);}        
        sumCfp[noe] = sum;
        break;
        
      case 4: // Tetra O2
        sum = 0.; sizecoefs = 4;
        for (E_Int nov = 0; nov < 4; nov++)
          sum += K_FUNC::E_abs(ptrCoefs[nov]);                   
        sumCfp[noe] = sum;
        break;
      
      case 5: // Lagrange O5
        sum = 0.;  sizecoefs = 15;
        for (E_Int kk=0; kk<5; kk++)
          for (E_Int jj=0; jj<5; jj++)
            for (E_Int ii=0; ii<5; ii++)
            {sum += K_FUNC::E_abs(ptrCoefs[ii]*ptrCoefs[jj+5]*ptrCoefs[kk+10]);}
        sumCfp[noe] = sum;
        break;
      
      default:
        RELEASESHAREDN(pyIndRcv, rcvPtsI);
        RELEASESHAREDN(pyIndExtrap, extrapPtsI);
        RELEASESHAREDN(pyArrayTypes, typesI);
        RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
        //printf("Interpolation type = %d \n", type);
        PyErr_SetString(PyExc_TypeError, 
                        "getExtrapCoefs: not a valid interpolation type.");
        return NULL;
    }
    ptrCoefs += sizecoefs;
  }     

  // sortie
  RELEASESHAREDS(tpl, sumCf);
  RELEASESHAREDN(pyIndRcv, rcvPtsI);
  RELEASESHAREDN(pyIndExtrap, extrapPtsI);
  RELEASESHAREDN(pyArrayTypes, typesI);
  RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
  return tpl;
}
