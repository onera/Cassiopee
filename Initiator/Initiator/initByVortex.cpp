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

// Initialization by vortices

# include "initiator.h"

using namespace K_FLD;
using namespace std;

extern "C"
{
  void k6lamb_(const E_Float& x0, const E_Float& y0,
               const E_Float& Gamma, const E_Float& MInf,
               const E_Int& npts,      
               const E_Float* xc, const E_Float* yc, const E_Float* zc,
               E_Float* field);

  void k6scully_(const E_Float& x0, const E_Float& y0,
                 const E_Float& Gamma, const E_Float& a, const E_Float& MInf,
                 const E_Int& npts,      
                 const E_Float* xc, const E_Float* yc, const E_Float* zc,
                 E_Float* field);

  void k6visbal_(const E_Float& x0, const E_Float& y0,
                 const E_Float& Gamma, const E_Float& MInf,
                 const E_Int& npts,
                 const E_Float* xc, const E_Float* yc, const E_Float* zc,
                 E_Float* field);

  void k6yee_(const E_Float& x0, const E_Float& y0,
              const E_Float& Gamma, const E_Float& Minf,
              const E_Int& npts,
              const E_Float* xc, const E_Float* yc, const E_Float* zc,
              E_Float* field);

  void k6scully2_(const E_Float& x0, const E_Float& y0,
                  const E_Float& Gamma, const E_Float& a, const E_Float& MInf,
                  const E_Int& npts,      
                  const E_Float* xc, const E_Float* yc, const E_Float* zc,
                  E_Float* field);
}
// ============================================================================
/* Init by a lamb vortex */
// ============================================================================
PyObject* K_INITIATOR::initLamb(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float x0, y0, Gam, MInf;

#ifdef E_DOUBLEREAL
  if (!PyArg_ParseTuple(args, "O(dd)dd", &array,
                        &x0, &y0, &Gam, &MInf) )
#else
    if (!PyArg_ParseTuple(args, "O(ff)ff", &array,  
                          &x0, &y0, &Gam, &MInf)) 
#endif
    {
      return NULL;
    }
  
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  FldArrayF cfdField, fieldTot;
  E_Int res = 
    K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType);

  char varStringOut[K_ARRAY::VARSTRINGLENGTH];

  if (res == 1 || res == 2)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      delete f;
      PyErr_SetString(PyExc_ValueError,
                      "initLamb: can't find coordinates in array.");
      return NULL;        
    }
    posx++; posy++; posz++;
    
    E_Int defVars = 0; //says if output array is built using default variables or not
    //E_Int nCfdFld;
    E_Int neqTot = f->getNfld();
    E_Int posro = K_ARRAY::isDensityPresent(varString);
    E_Int posrou = K_ARRAY::isMomentumXPresent(varString);
    E_Int posrov = K_ARRAY::isMomentumYPresent(varString);
    E_Int posrow = K_ARRAY::isMomentumZPresent(varString); 
    E_Int posroe = K_ARRAY::isEnergyStagnationDensityPresent(varString); 
    vector<E_Int> posTot;
    if (posro == -1 || posrou == -1 || posrov == -1 || posrow == -1 || posroe == -1)
    {
      defVars = 1;
      //nCfdFld = 5;
      neqTot = 8;
      strcpy(varStringOut, "x,y,z,ro,rou,rov,row,roE");
    }
    else 
    {
      posro++; posrou++; posrov++; posrow++; posroe++;
      posTot.push_back(posro);
      posTot.push_back(posrou);
      posTot.push_back(posrov);
      posTot.push_back(posrow);
      posTot.push_back(posroe);
      //nCfdFld = neqTot-3;
      strcpy(varStringOut, varString);
    }  
    E_Int npts = f->getSize();

    FldArrayF cfdField(npts, 5);
    cfdField.setAllValuesAtNull();

    k6lamb_(x0, y0, Gam, MInf, npts, 
            f->begin(posx), f->begin(posy), f->begin(posz),
            cfdField.begin());
    
    fieldTot.malloc(npts, neqTot);
    E_Float* xt = f->begin(posx);
    E_Float* yt = f->begin(posy);
    E_Float* zt = f->begin(posz);
    E_Float* xt2 = fieldTot.begin(1);
    E_Float* yt2 = fieldTot.begin(2);
    E_Float* zt2 = fieldTot.begin(3);
    E_Int posTotsize = posTot.size();
    switch (defVars)
    {
      case 0 :
        for (E_Int i = 0; i < npts; i++)
        {
          for (E_Int eq = 1; eq <= neqTot; eq++)
            fieldTot(i,eq) = (*f)(i,eq);
        
          for (E_Int c = 0; c < posTotsize; c++)
          {
            E_Int pos0 = posTot[c];
            fieldTot(i,pos0) = cfdField(i,c+1);
          }
        }
        break;
        
      default:
        for (E_Int i = 0; i < npts; i++)
        {
          xt2[i] = xt[i]; yt2[i] = yt[i]; zt2[i] = zt[i];
          for (E_Int eq = 1; eq <= 5; eq++)
            fieldTot(i,eq+3) = cfdField(i,eq);
        }
        break;
    }
    
    // build array
    delete f;
    FldArrayF* an = new FldArrayF(fieldTot);
    PyObject* tpl;
    if (res == 1)
    {
      tpl = K_ARRAY::buildArray(*an, varStringOut, im, jm, km);
      delete an;
    }
    else
    {
      tpl = K_ARRAY::buildArray(*an, varStringOut, *cn, -1, eltType);
      delete an; delete cn;
    }
    return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "initLamb: invalid type of array.");
    return NULL;
  }
}

// ===========================================================================
/* Init by a Visbal vortex */
// ============================================================================
PyObject* K_INITIATOR::initVisbal(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float x0, y0, Gam, MInf;

#ifdef E_DOUBLEREAL
  if (!PyArg_ParseTuple(args, "O(dd)dd", &array, 
                        &x0, &y0, &Gam, &MInf))
#else
    if (!PyArg_ParseTuple(args, "O(ff)ff", &array,  
                          &x0, &y0, &Gam, &MInf)) 
#endif
    {
      return NULL;
    }
  
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  FldArrayF cfdField, fieldTot;
  E_Int res = 
    K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType);
  char varStringOut[K_ARRAY::VARSTRINGLENGTH];

  if (res == 1 || res == 2)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      delete f;
      PyErr_SetString(PyExc_ValueError,
                      "initVisbal: can't find coordinates in array.");
      return NULL;        
    }
    posx++; posy++; posz++;
    
    E_Int defVars = 0;//says if output array is built using default variables or not
    //E_Int nCfdFld;
    E_Int neqTot = f->getNfld();

    E_Int posro = K_ARRAY::isDensityPresent(varString);
    E_Int posrou = K_ARRAY::isMomentumXPresent(varString);
    E_Int posrov = K_ARRAY::isMomentumYPresent(varString);
    E_Int posrow = K_ARRAY::isMomentumZPresent(varString); 
    E_Int posroe = K_ARRAY::isEnergyStagnationDensityPresent(varString); 
    vector<E_Int> posTot;

    if (posro == -1 || posrou == -1 || posrov == -1 || posrow == -1 || posroe == -1)
    {
      defVars = 1;
      //nCfdFld = 5;
      neqTot = 8;
      strcpy(varStringOut, "x,y,z,ro,rou,rov,row,roE");
    }
    else 
    {
      posro++; posrou++; posrov++; posrow++; posroe++;
      posTot.push_back(posro);
      posTot.push_back(posrou);
      posTot.push_back(posrov);
      posTot.push_back(posrow);
      posTot.push_back(posroe);
      //nCfdFld = neqTot-3;
      strcpy(varStringOut, varString);
    }  
    E_Int npts = f->getSize();

    FldArrayF cfdField(npts, 5);
    cfdField.setAllValuesAtNull();

    k6visbal_(x0, y0, Gam, MInf, npts, 
              f->begin(posx), f->begin(posy), f->begin(posz),
              cfdField.begin());
    
    fieldTot.malloc(npts, neqTot);
    E_Float* xt = f->begin(posx);
    E_Float* yt = f->begin(posy);
    E_Float* zt = f->begin(posz);
    E_Float* xt2 = fieldTot.begin(1);
    E_Float* yt2 = fieldTot.begin(2);
    E_Float* zt2 = fieldTot.begin(3);
    E_Int posTotsize = posTot.size();

    switch (defVars)
    {
      case 0 :   
        for (E_Int i = 0; i < npts; i++)
        {
          for (E_Int eq = 1; eq <= neqTot; eq++)
            fieldTot(i,eq) = (*f)(i,eq);
          
          for (E_Int c = 0; c < posTotsize; c++)
          {
            E_Int pos0 = posTot[c];
            fieldTot(i,pos0) = cfdField(i,c+1);
          }
        }
        break;
        
      default:
        for (E_Int i = 0; i < npts; i++)
        {
          xt2[i] = xt[i]; yt2[i] = yt[i]; zt2[i] = zt[i];

          for (E_Int eq = 1; eq <= 5; eq++)
            fieldTot(i,eq+3) = cfdField(i,eq);
        }
        break;
    }  
    // build array
    delete f;
    FldArrayF* an = new FldArrayF(fieldTot);
    PyObject* tpl;
    if (res == 1)
    {
      tpl = K_ARRAY::buildArray(*an, varStringOut, im, jm, km);
      delete an;
    }
    else
    {
      tpl = K_ARRAY::buildArray(*an, varStringOut, *cn, -1, eltType);
      delete an; delete cn;
    }
    return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "initVisbal : invalid type of array.");
    return NULL;
  }
}

// ============================================================================
/* Init by a scully vortex */
// ============================================================================
PyObject* K_INITIATOR::initScully(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float x0, y0, Gam, coreRadius, MInf;
  E_Int model=0;
  
#if defined E_DOUBLEREAL && defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "O(dd)dddl", &array,
                        &x0, &y0, &Gam, &coreRadius, &MInf, &model)) return NULL;
#elif defined E_DOUBLEREAL && !defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "O(dd)dddi", &array,
                        &x0, &y0, &Gam, &coreRadius, &MInf, &model)) return NULL;
#elif !defined E_DOUBLEREAL && defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "O(ff)fffl", &array,
                        &x0, &y0, &Gam, &coreRadius, &MInf, &model)) return NULL;
#else
  if (!PyArg_ParseTuple(args, "O(ff)fffi", &array,
                        &x0, &y0, &Gam, &coreRadius, &MInf, &model)) return NULL;
#endif
  
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  FldArrayF cfdField, fieldTot;
  
  E_Int res = 
    K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType);
  char varStringOut[K_ARRAY::VARSTRINGLENGTH];

  if (res == 1 || res == 2)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      delete f;
      PyErr_SetString(PyExc_ValueError,
                      "initScully: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;

    E_Int defVars = 0;//says if output array is built using default variables or not
    //E_Int nCfdFld;
    E_Int neqTot = f->getNfld();

    E_Int posro = K_ARRAY::isDensityPresent(varString);
    E_Int posrou = K_ARRAY::isMomentumXPresent(varString);
    E_Int posrov = K_ARRAY::isMomentumYPresent(varString);
    E_Int posrow = K_ARRAY::isMomentumZPresent(varString); 
    E_Int posroe = K_ARRAY::isEnergyStagnationDensityPresent(varString); 
    vector<E_Int> posTot;

    if (posro == -1 || posrou == -1 || posrov == -1 || posrow == -1 || posroe == -1)
    {
      defVars = 1;
      //nCfdFld = 5;
      neqTot = 8;
      strcpy(varStringOut, "x,y,z,ro,rou,rov,row,roE");
    }
    else 
    {
      posro++; posrou++; posrov++; posrow++; posroe++;
      posTot.push_back(posro);
      posTot.push_back(posrou);
      posTot.push_back(posrov);
      posTot.push_back(posrow);
      posTot.push_back(posroe);
      //nCfdFld = neqTot-3;
      strcpy(varStringOut, varString);
    }  
    E_Int npts = f->getSize();

    FldArrayF cfdField(npts, 5);
    cfdField.setAllValuesAtNull();
    
    switch (model)
    {
      case 0:
        k6scully_(x0, y0, Gam, coreRadius, MInf, npts, 
                  f->begin(posx), f->begin(posy), f->begin(posz),
                  cfdField.begin());
      default :
        k6scully2_(x0, y0, Gam, coreRadius, MInf, npts,
                   f->begin(posx), f->begin(posy), f->begin(posz),
                   cfdField.begin());
    }
    
    fieldTot.malloc(npts, neqTot);
    E_Float* xt = f->begin(posx);
    E_Float* yt = f->begin(posy);
    E_Float* zt = f->begin(posz);
    E_Float* xt2 = fieldTot.begin(1);
    E_Float* yt2 = fieldTot.begin(2);
    E_Float* zt2 = fieldTot.begin(3);
    E_Int posTotsize = posTot.size();

    switch (defVars)
    {
      case 0 :   
        for (E_Int i = 0; i < npts; i++)
        {
          for (E_Int eq = 1; eq <= neqTot; eq++)
            fieldTot(i,eq) = (*f)(i,eq);
          
          for (E_Int c = 0; c < posTotsize; c++)
          {
            E_Int pos0 = posTot[c];
            fieldTot(i,pos0) = cfdField(i,c+1);
          }
        }
        break;
        
      default:
        for (E_Int i = 0; i < npts; i++)
        {
          xt2[i] = xt[i]; yt2[i] = yt[i]; zt2[i] = zt[i];
          for (E_Int eq = 1; eq <= 5; eq++)
            fieldTot(i,eq+3) = cfdField(i,eq);
        }
        break;
    }
    // build array
    delete f;
    FldArrayF* an = new FldArrayF(fieldTot);
    PyObject* tpl;
    if (res == 1)
    {
      tpl = K_ARRAY::buildArray(*an, varStringOut, im, jm, km);
      delete an;
    }
    else
    {
      tpl = K_ARRAY::buildArray(*an, varStringOut, *cn, -1, eltType);
      delete an; delete cn;
    }
    return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "initScully: unknown type of array.");
    return NULL;
  }
}

// ===========================================================================
/* Init by a Yee vortex */
// ============================================================================
PyObject*
K_INITIATOR::initYee(PyObject* self,
                     PyObject* args)
{
  PyObject* array;
  E_Float x0, y0, Gam, Minf;

#ifdef E_DOUBLEREAL
  if (!PyArg_ParseTuple(args, "O(dd)dd", &array, 
                        &x0, &y0, &Gam, &Minf))
#else
    if (!PyArg_ParseTuple(args, "O(ff)ff", &array,  
                          &x0, &y0, &Gam, &Minf)) 
#endif
    {
      return NULL;
    }
  
  // Check array
  E_Int im, jm, km;
  FldArrayF* f;
  char* varString; char* eltType;
  FldArrayI* cn;
  FldArrayF cfdField, fieldTot;
  E_Int res = 
    K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType);
  char varStringOut[K_ARRAY::VARSTRINGLENGTH];

  if (res == 1 || res == 2)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      delete f;
      PyErr_SetString(PyExc_ValueError,
                      "initYee: can't find coordinates in array.");
      return NULL;        
    }
    posx++; posy++; posz++;
    
    E_Int defVars = 0;//says if output array is built using default variables or not
    //E_Int nCfdFld;
    E_Int neqTot = f->getNfld();

    E_Int posro = K_ARRAY::isDensityPresent(varString);
    E_Int posrou = K_ARRAY::isMomentumXPresent(varString);
    E_Int posrov = K_ARRAY::isMomentumYPresent(varString);
    E_Int posrow = K_ARRAY::isMomentumZPresent(varString); 
    E_Int posroe = K_ARRAY::isEnergyStagnationDensityPresent(varString); 
    vector<E_Int> posTot;

    if (posro == -1 || posrou == -1 || posrov == -1 || posrow == -1 || posroe == -1)
    {
      defVars = 1;
      //nCfdFld = 5;
      neqTot = 8;
      strcpy(varStringOut, "x,y,z,ro,rou,rov,row,roE");
    }
    else 
    {
      posro++; posrou++; posrov++; posrow++; posroe++;
      posTot.push_back(posro);
      posTot.push_back(posrou);
      posTot.push_back(posrov);
      posTot.push_back(posrow);
      posTot.push_back(posroe);
      //nCfdFld = neqTot-3;
      strcpy(varStringOut, varString);
    }  
    E_Int npts = f->getSize();

    FldArrayF cfdField(npts, 5);
    cfdField.setAllValuesAtNull();

    k6yee_(x0, y0, Gam, Minf, npts, 
              f->begin(posx), f->begin(posy), f->begin(posz),
              cfdField.begin());
    
    fieldTot.malloc(npts, neqTot);
    E_Float* xt = f->begin(posx);
    E_Float* yt = f->begin(posy);
    E_Float* zt = f->begin(posz);
    E_Float* xt2 = fieldTot.begin(1);
    E_Float* yt2 = fieldTot.begin(2);
    E_Float* zt2 = fieldTot.begin(3);
    E_Int posTotsize = posTot.size();

    switch (defVars)
    {
      case 0 :   
        for (E_Int i = 0; i < npts; i++)
        {
          for (E_Int eq = 1; eq <= neqTot; eq++)
            fieldTot(i,eq) = (*f)(i,eq);
          
          for (E_Int c = 0; c < posTotsize; c++)
          {
            E_Int pos0 = posTot[c];
            fieldTot(i,pos0) = cfdField(i,c+1);
          }
        }
        break;
        
      default:
        for (E_Int i = 0; i < npts; i++)
        {
          xt2[i] = xt[i]; yt2[i] = yt[i]; zt2[i] = zt[i];

          for (E_Int eq = 1; eq <= 5; eq++)
            fieldTot(i,eq+3) = cfdField(i,eq);
        }
        break;
    }  
    // build array
    delete f;
    FldArrayF* an = new FldArrayF(fieldTot);
    PyObject* tpl;
    if (res == 1)
    {
      tpl = K_ARRAY::buildArray(*an, varStringOut, im, jm, km);
      delete an;
    }
    else
    {
      tpl = K_ARRAY::buildArray(*an, varStringOut, *cn, -1, eltType);
      delete an; delete cn;
    }
    return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "initYee : invalid type of array.");
    return NULL;
  }
}

