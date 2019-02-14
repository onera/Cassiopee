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
# include <stdio.h>
# include <string.h>

# include "post.h"
# include <vector>

using namespace K_FUNC;
using namespace K_FLD;
using namespace std;

extern "C"
{
  void usurp_(const E_Int& nzone, const E_Int* nit, const E_Int* njt,
              const E_Int* nkt, const E_Int* iblank, const E_Int& ncellmax,
              const E_Float* coord, const E_Int& nptsmax, 
              E_Float* ratio);
}
// ============================================================================
/* USURP */
// ============================================================================
PyObject* K_POST::usurpF(PyObject* self, PyObject* args)
{ 
  PyObject* blkArrays;
  PyObject* ibArrays;
  
  if (!PyArg_ParseTuple(args, "OO",  &blkArrays, &ibArrays)) return NULL;
  
  // Check every arrays  
  if (PyList_Check(blkArrays) == 0)
  {
    PyErr_SetString(
      PyExc_TypeError, 
      "usurp: first argument must be a list.");
    return NULL;
  }
  if (PyList_Check(ibArrays) == 0)
  {
    PyErr_SetString(
      PyExc_TypeError, 
      "usurp: second argument must be a list.");
    return NULL;
  }
  
  E_Int nblkArrays = PyList_Size(blkArrays);
  E_Int nibArrays = PyList_Size(ibArrays);
  
  if (nibArrays != nblkArrays)
  {
    PyErr_SetString(
      PyExc_ValueError, 
      "usurp: number of zones in 1st and 2nd arguments must be equal.");
    return NULL;    
  }

  PyObject* tpl;  
  E_Int nil, njl, nkl;
  vector<E_Int> ni; vector<E_Int> nj; vector<E_Int> nk;
  vector<FldArrayF*> fieldc;
  char* varString;
  vector<E_Int> elt;
  FldArrayF* f; FldArrayI* cn;
  E_Int res, npts;
  char* eltType;
  vector<FldArrayI*> listOfCellNF;
  E_Int posx, posy, posz;
  vector<E_Int> posxt; vector<E_Int> posyt; vector<E_Int> poszt;
  
  // Recuperation des coordonnees et de la solution a partir de blkArrays
  for (int i = 0; i < nblkArrays; i++)
  {
    tpl = PyList_GetItem(blkArrays, i);
    res = K_ARRAY::getFromArray(tpl, varString, f, nil, njl, nkl, cn, eltType);
    
    if (res == 1) // structured
    {
      // check if coordinates exist
      posx = K_ARRAY::isCoordinateXPresent(varString); 
      posy = K_ARRAY::isCoordinateYPresent(varString); 
      posz = K_ARRAY::isCoordinateZPresent(varString);
      if (posx == -1 || posy == -1 || posz == -1)
      {
        delete f;
        printf("Warning: usurp: coordinates not found in array. Skipped...\n");
        goto end;
      }
      
      posx++; posy++; posz++;
      posxt.push_back(posx);
      posyt.push_back(posy);
      poszt.push_back(posz);

      if (nil*njl*nkl > 0)
      {        
        ni.push_back(nil); nj.push_back(njl); nk.push_back(nkl);
        fieldc.push_back(f);
        
      }
      else 
      {
        delete f;
        printf("Warning: usurp: one array is empty.\n");
        goto end;
      }
    }
    else
    {
      printf("Warning: usurp: one array is invalid.\n");
      goto end;
    }
    end:;
  }
  
  // Nfld
  if (fieldc.size() == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "usurp: nothing to write.");
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  E_Int nzone = ni.size();
  FldArrayI nit(nzone);
  FldArrayI njt(nzone);
  FldArrayI nkt(nzone);
  FldArrayIS indir(nzone); // number of nodes in each block
  indir.setAllValuesAtNull();

  for (E_Int n = 0; n < nzone; n++)
  {
    nit[n] = ni[n];
    njt[n] = nj[n];
    nkt[n] = nk[n];
    indir[n] = ni[n]*nj[n]*nk[n];
  }
  
  // Recuperation de iblank a partir de ibArrays
  E_Int pos = -1;
  E_Int nin, njn, nkn;
  for (int i = 0; i < nibArrays; i++)
  {
    tpl = PyList_GetItem(ibArrays, i);
    res = K_ARRAY::getFromArray(
      tpl, varString, f, nil, njl, nkl, cn, eltType);
    
    if (res == 1)
    {
      npts = nil*njl*nkl;
      if (npts > 0)
      {   
        // check if celln in centers correspond to the field in nodes
        nin = E_max(1, nit[i]-1);
        njn = E_max(1, njt[i]-1);
        nkn = E_max(1, nkt[i]-1);
        if (nil != nin || njl != njn || nkl != nkn)
        {
          
          PyErr_SetString(PyExc_TypeError,
                          "usurp: first argument (defined at nodes) does not correspond to the second argument (defined in centers).");
          return NULL;
        }
        pos = K_ARRAY::isCellNatureField2Present(varString);   

        FldArrayI* cellNF = new FldArrayI(npts);
        if (pos != -1)
        {
          pos++;
          for (E_Int ind = 0; ind < npts; ind++)
            (*cellNF)[ind] = E_Int((*f)(ind, pos));
        }
        else cellNF->setAllValuesAt(1);
        
        listOfCellNF.push_back(cellNF);          
      }
      else 
      {
        printf("Warning: usurp: one array is empty.\n");
        goto end2;
      }
    }
    else
    {
      printf("Warning: usurp: cellN is not defined for one zone. Set to 1.\n");
      nin = E_max(1,nit[i]-1); njn = E_max(1,njt[i]-1); nkn = E_max(1, nkt[i]-1);
      FldArrayI* cellNF = new FldArrayI(nin*njn*nkn); cellNF->setAllValuesAt(1);
      listOfCellNF.push_back(cellNF); 
    }
    end2:;
  }
  // check if cellNF list exists
  nzone = listOfCellNF.size();
  if (nzone == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "usurp: nothing to write.");
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  // Computes the number tot of cells and nodes for FORTRAN allocation
  E_Int nptsmax = 0;
  E_Int ncellmax = 0;
  for (E_Int n = 0; n < nzone; n++)
  {
    ncellmax = ncellmax + listOfCellNF[n]->getSize();
    nptsmax = nptsmax + fieldc[n]->getSize();    
  }
  //------------------------------------------------
  // Convert the cellN to iblank fortran input array
  //------------------------------------------------
  FldArrayI iblank(ncellmax);
  E_Int nelts = 0;
  E_Int celln;
  for (E_Int n = 0; n < nzone; n++)
  {
    FldArrayI& cellNF = *listOfCellNF[n];
    for (E_Int i = 0; i < cellNF.getSize(); i++)
    {
      celln = cellNF[i];
      switch (celln) 
      {
        case 0: // blanked
          iblank[i+nelts] = 0;
          break;
        case 1: // computed
          iblank[i+nelts] = 1;
           break;
        case -1000: // orphan
          iblank[i+nelts] = 101;
          break;
        default:// interpolated
          iblank[i+nelts] = -1;
          break;
      }
    }
    nelts = nelts + cellNF.getSize();
  }

  //----------------------------------------------
  // Get the grid connectivity
  //----------------------------------------------
  nelts = 0;
  FldArrayF coord(nptsmax,3);  
  
  for (E_Int n = 0; n < nzone; n++)
  {
    FldArrayF& field = *fieldc[n];
    posx = posxt[n];
    posy = posyt[n];
    posz = poszt[n];
    E_Float* fx = field.begin(posx);
    E_Float* fy = field.begin(posy);
    E_Float* fz = field.begin(posz);
    
    E_Float* coordx = coord.begin(1);
    E_Float* coordy = coord.begin(2);
    E_Float* coordz = coord.begin(3);
    E_Int npts = field.getSize();
    for (E_Int i = 0; i < npts; i++)
    {
      coordx[i+nelts] = fx[i];
      coordy[i+nelts] = fy[i];
      coordz[i+nelts] = fz[i];
    } 
    nelts = nelts + npts;
  }
   
  // Deleting fields
  E_Int fieldcSize = fieldc.size();
  E_Int listOfCellNFSize = listOfCellNF.size();
  for (E_Int i = 0; i < fieldcSize; i++)
    delete fieldc[i];
  
  for (E_Int i = 0; i < listOfCellNFSize; i++)
    delete listOfCellNF[i];

  ni.clear();
  nj.clear();
  nk.clear();
  fieldc.clear();
  listOfCellNF.clear();

  //----------------------------------------------
  // usurp
  //----------------------------------------------
  FldArrayF ratio(ncellmax);

  usurp_(nzone, nit.begin(), njt.begin(), nkt.begin(), 
         iblank.begin(), ncellmax, coord.begin(), nptsmax,
         ratio.begin());

  //------------------------
  // convert 2 arrays 
  //------------------------
  // Building numpy arrays
  E_Int istart = 0;

  vector<FldArrayF*> vectOfRatios;
  vector<E_Int> nis; vector<E_Int> njs; vector<E_Int> nks;
  
  for (E_Int n = 0; n < nzone; n++)
  {
    E_Int ni1 = ni[n]-1;
    E_Int nj1 = nj[n]-1;
    E_Int nk1 = nk[n]-1;
    if (ni1 == 0) ni1 = 1;
    if (nj1 == 0) nj1 = 1;
    if (nk1 == 0) nk1 = 1;
    E_Int ncells = ni1 * nj1 * nk1;
    
    FldArrayF* fieldr = new FldArrayF(ncells, 1); // iblank, ratio aux centres
    E_Float* ration = fieldr->begin(1);

    for (E_Int i = 0; i < ncells; i++)
    {
      E_Int ind = i+istart;
      ration[i] = ratio[ind];

      /* si le point est masque mettre un ratio nul */
      if ( iblank[ind] == 0 ) 
        ration[i] = 0.;
    }

    istart = istart + ncells;
    vectOfRatios.push_back(fieldr);
    nis.push_back(ni1);
    njs.push_back(nj1);
    nks.push_back(nk1);
  }

  /*--------------*/
  /* build arrays */
  /*--------------*/
  PyObject* l = PyList_New(0);
  for (E_Int i = 0; i < nzone; i++)
  {
    tpl = K_ARRAY::buildArray(*vectOfRatios[i], "ratio",
                              nis[i], njs[i], nks[i]);
    delete vectOfRatios[i];
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
  return l;
}

