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
# include "converter.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Add ghost cells to a NGON array */
//=============================================================================
PyObject* K_CONVERTER::addGhostCellsNGonNodes(PyObject* self, PyObject* args)
{
  PyObject *arrayN;
  E_Int depth;
  if (!PYPARSETUPLE_(args, O_ I_, &arrayN, &depth)) return NULL;

  if (depth < 1) 
  {
    PyErr_SetString(PyExc_ValueError, 
                    "addGhostCells: depth must be > 0.");
    return NULL;
  }
  E_Int ni, nj, nk;
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int res = K_ARRAY::getFromArray3(arrayN, varString, f, ni, nj, nk, 
                                     cn, eltType);
  if (res != 2)
  {
    if (res == 1) RELEASESHAREDS(arrayN,f);
    PyErr_SetString(PyExc_TypeError, 
                    "addGhostCells: array is invalid.");
    return NULL;
  }
  if (strcmp(eltType,"NGON") != 0 && strcmp(eltType,"NGON*") != 0 ) 
  {
    RELEASESHAREDU(arrayN, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "addGhostCells: array must be a NGON.");
    return NULL;
  }
  
  /* Determination des faces externes du maillage reel */
  E_Int sizeFN2 = 0;
  E_Int sizeEF2 = 0;
  // Connectivite Faces/Elts reelle
  FldArrayI cFEp; K_CONNECT::connectNG2FE(*cn, cFEp);
  vector<E_Int> facesExt;//faces externes du maillage reel: demarre a 0
  E_Int* cFEp1 = cFEp.begin(1);
  E_Int* cFEp2 = cFEp.begin(2);
  E_Int* cNGp = cn->begin();
  E_Int nfacesp = cNGp[0];// nb de faces ds le maillage reel
  E_Int sizeFN = cNGp[1];
  E_Int sizeEF = cNGp[3+sizeFN];
  E_Int etg, etd, nvertp;
  E_Int nfacesExt, nfacesAdd, neltsAdd, sizeOut, sizefout;
  E_Int dimNGON = 3;
  E_Int api = f->getApi();
  E_Int npts = f->getSize();
  E_Int nfld = f->getNfld();
  E_Int nptsAdd = 0;

  if (cNGp[2] == 2) dimNGON = 2;
  else if (cNGp[2] > 2) dimNGON = 3;
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "addGhostCells: array must be 2D or 3D.");
    RELEASESHAREDU(arrayN,f,cn); return NULL;
  }
  E_Int* ptr = cNGp+2;
  if (dimNGON == 3)
  {
    for (E_Int nof = 0; nof < nfacesp; nof++)
    {
      nvertp = ptr[0];
      etg = cFEp1[nof]; etd = cFEp2[nof];
      if (etd == 0 || etg == 0) //face exterieure
      {
        facesExt.push_back(nof); nptsAdd += nvertp;
        sizeFN2 += (5*nvertp)+(nvertp+1); // 4 vertex par face laterale + nvert sur la face opposee +1 pour dimensionner
        sizeEF2 += (2+nvertp+1);// nvertp faces laterales + 2 faces + 1 pour dimensionner
      }         
      ptr += nvertp+1;
    }
    sizeFN2 *= depth;
    sizeEF2 *= depth;
    nfacesExt = facesExt.size();
    neltsAdd = nfacesExt*depth; //nb d elts ajoutes
    nptsAdd *= depth;
  }
  else 
  {
    for (E_Int nof = 0; nof < nfacesp; nof++)
    {
      etg = cFEp1[nof]; etd = cFEp2[nof];
      if (etd == 0 || etg == 0) //face exterieure
      { facesExt.push_back(nof); nptsAdd += 2;}
    }
    nfacesExt = facesExt.size();
    nfacesAdd = 3*depth;// 1 face ajoutee identique a la face externe + 2 sur les cotes
    sizeFN2 = nfacesAdd*nfacesExt*(2+1);// 2 vertex par face * nb de faces ajoutees
    neltsAdd = nfacesExt*depth;
    nptsAdd *= depth;
    sizeEF2 = neltsAdd*5;// nb Elts ajoutes*4 faces par cellule fictive+1 pour dimensionner  
  }
  
  sizeOut = sizeFN+sizeEF+sizeFN2+sizeEF2+4;
  sizefout = npts+nptsAdd;
  FldArrayI* cnout = new FldArrayI(sizeOut); cnout->setAllValuesAtNull();
  FldArrayF* fout = new FldArrayF(sizefout,nfld); fout->setAllValuesAtNull();

  if (dimNGON == 2) 
  {
    if (nfacesExt > 0)
      addGhostCellsNGon2D(depth, *f, *cn, neltsAdd, sizeFN2, sizeEF2, facesExt, fout, cnout);
    else {*cnout = *cn; *fout = *f;} 
  }
  else addGhostCellsNGon3D(depth, *f, *cn, neltsAdd, sizeFN2, sizeEF2, facesExt, fout, cnout);

  // Sortie
  cnout->setNGonType(cn->getNGonType());
  PyObject* tpl = K_ARRAY::buildArray3(*fout, varString, *cnout, eltType, api);
  delete cnout; delete fout;
  RELEASESHAREDU(arrayN, f, cn);
  return tpl;
}  
//=============================================================================
/* Add ghost cells to a NGON array - input array defines a NGON* array
   Fields located at centers are modified */
//=============================================================================
PyObject* K_CONVERTER::addGhostCellsNGonCenters(PyObject* self, PyObject* args)
{
  PyObject *arrayC;
  E_Int depth;
  if (!PYPARSETUPLE_(args, O_ I_, &arrayC, &depth)) return NULL;

  if (depth < 1)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "addGhostCellsNGon: depth must be > 0.");
    return NULL;
  }
 
  E_Int nic, njc, nkc;
  char* varStringC; char* eltTypec;
  FldArrayF* fc; FldArrayI* cnc;
  E_Int resc = K_ARRAY::getFromArray3(arrayC, varStringC, fc, nic, njc, nkc, 
                                      cnc, eltTypec);
  if (resc != 2)
  {
    if (resc == 1) RELEASESHAREDS(arrayC,fc);
    PyErr_SetString(PyExc_TypeError, 
                    "addGhostCellsNGon: array is invalid.");
    return NULL;
  }
  if (strcmp(eltTypec, "NGON*") != 0) 
  {
    
    PyErr_SetString(PyExc_TypeError,
                    "addGhostCellsNGon: array must be NGON*.");
    RELEASESHAREDU(arrayC,fc,cnc); return NULL;
  }

  /* Determination des faces externes du maillage reel */
  E_Int sizeFN2 = 0;
  E_Int sizeEF2 = 0;
  // Connectivite Faces/Elts reelle
  FldArrayI cFEp; K_CONNECT::connectNG2FE(*cnc, cFEp);
  vector<E_Int> facesExt;//faces externes du maillage reel: demarre a 0
  E_Int* cFEp1 = cFEp.begin(1);
  E_Int* cFEp2 = cFEp.begin(2);
  E_Int* cNGp = cnc->begin();
  E_Int nfacesp = cNGp[0];// nb de faces ds le maillage reel
  E_Int sizeFN = cNGp[1];
  E_Int sizeEF = cNGp[3+sizeFN];
  E_Int indface, etl, etg, etd, ce, nvertp;
  E_Int dimNGON, neltsAdd, nfacesAdd, nfacesExt, sizeOut, sizefout;
  E_Int nptsAdd = 0;
  
  if ((*cnc)[2] == 2) dimNGON = 2;
  else if ((*cnc)[2] > 2) dimNGON = 3;
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "addGhostCells: array must be 2D or 3D.");
    RELEASESHAREDU(arrayC,fc,cnc);
    return NULL;
  }
  PyObject* tpl;
  E_Int api = fc->getApi();
  E_Int nfldc = fc->getNfld();
  E_Int nelts = fc->getSize();
  E_Int* ptr = cNGp+2;
  FldArrayI* cnout; FldArrayF *fcout, *fout;
  E_Int npts = 0; E_Int nfld = 1;
  FldArrayF* f = new FldArrayF(0,1);
  // Sortie
  //PyObject* l = PyList_New(0);
  if (dimNGON == 3)
  {
    for (E_Int nof = 0; nof < nfacesp; nof++)
    {
      nvertp = ptr[0];
      etg = cFEp1[nof]; etd = cFEp2[nof];
      if (etd == 0 || etg == 0) // face exterieure
      {
        facesExt.push_back(nof); nptsAdd += nvertp;
        sizeFN2 += (5*nvertp)+(nvertp+1); // 4 vertex par face laterale + nvert sur la face opposee +1 pour dimensionner
        sizeEF2 += (2+nvertp+1);// nvertp faces laterales + 2 faces + 1 pour dimensionner
      }         
      ptr += nvertp+1;
    }
    sizeFN2 *= depth;
    sizeEF2 *= depth;
    nfacesExt = facesExt.size();
    neltsAdd = nfacesExt*depth; //nb d elts ajoutes
    nptsAdd *= depth;
    sizeOut = sizeFN+sizeEF+sizeFN2+sizeEF2+4;
    sizefout = npts+nptsAdd;
    cnout = new FldArrayI(sizeOut); cnout->setAllValuesAtNull();
    fout = new FldArrayF(sizefout,nfld); fout->setAllValuesAtNull();
    addGhostCellsNGon3D(depth, *f, *cnc, neltsAdd, sizeFN2, sizeEF2, facesExt, fout, cnout);
    delete fout; delete f;
    fcout = new FldArrayF(nelts+neltsAdd,nfldc); fcout->setAllValuesAtNull();
    for (E_Int eq = 1; eq <= nfldc; eq++)
    {
      E_Float* fcp = fc->begin(eq);
      E_Float* fcpo = fcout->begin(eq);
      for (E_Int et = 0; et < nelts; et++)
        fcpo[et] = fcp[et];
      
      ce = nelts;
      for (E_Int nof = 0; nof < nfacesExt; nof++)
      {
        indface = facesExt[nof];
        etg = cFEp1[indface]; etd = cFEp2[indface];
        if (etg != 0) 
        {
          etl = etg-1;
          for (E_Int d = 1; d <= depth; d++)
            fcpo[ce+(d-1)*nfacesExt] = fcp[etl]; 
          ce++;
        }
        else if (etd != 0) 
        {
          etl = etd-1;
          for (E_Int d = 1; d <= depth; d++)
            fcpo[ce+(d-1)*nfacesExt] = fcp[etl]; 
          //fcpo[ce] = fcp[etl+(d-1)*nfacesExt]; 
          ce++;
        }
      }
    }
    // array en noeuds
    tpl = K_ARRAY::buildArray3(*fcout, varStringC, *cnout, "NGON", api);
    delete cnout; delete fcout;    
  }
  else // dimNGON = 2 
  {
    for (E_Int nof = 0; nof < nfacesp; nof++)
    {
      etg = cFEp1[nof]; etd = cFEp2[nof];
      if (etd == 0 || etg == 0) //face exterieure
      {facesExt.push_back(nof); nptsAdd += 2;}
    }
    nfacesExt = facesExt.size();
    nfacesAdd = 3*depth;// 1 face ajoutee identique a la face externe + 2 sur les cotes
    sizeFN2 = nfacesAdd*nfacesExt*(2+1);// 2 vertex par face * nb de faces ajoutees
    neltsAdd = nfacesExt*depth;
    sizeEF2 = neltsAdd*5;// nb Elts ajoutes*4 faces par cellule fictive+1 pour dimensionner  
    sizeOut = sizeFN+sizeEF+sizeFN2+sizeEF2+4;
    nptsAdd *= depth;
    sizefout = nptsAdd;

    if (nfacesExt > 0)
    {
      cnout = new FldArrayI(sizeOut); cnout->setAllValuesAtNull();
      fout = new FldArrayF(sizefout,nfld); fout->setAllValuesAtNull();
      addGhostCellsNGon2D(depth, *f, *cnc, neltsAdd, sizeFN2, sizeEF2, facesExt, fout, cnout);      
      delete fout;delete f;
      E_Int nfldc = fc->getNfld();
      E_Int nelts = fc->getSize();
      fcout = new FldArrayF(nelts+neltsAdd,nfldc); fcout->setAllValuesAtNull();
      for (E_Int eq = 1; eq <= nfldc; eq++)
      {
        E_Float* fcp = fc->begin(eq);
        E_Float* fcpo = fcout->begin(eq);
        for (E_Int et = 0; et < nelts; et++)
          fcpo[et] = fcp[et];
        
        ce = nelts;
        for (E_Int nof = 0; nof < nfacesExt; nof++)
        {
          indface = facesExt[nof];
          etg = cFEp1[indface]; etd = cFEp2[indface];
          
          if (etg != 0) 
          {
            etl = etg-1;
            for (E_Int d = 1; d <= depth; d++)
            { fcpo[ce] = fcp[etl]; ce++;}
          }
          else if (etd != 0) 
          {
            etl = etd-1;
            for (E_Int d = 1; d <= depth; d++)
            { fcpo[ce] = fcp[etl]; ce++;}
          }    
        }
      }
    }
    else
    {
      cnout = new FldArrayI(cnc->getSize()*cnc->getNfld());
      *cnout = *cnc; 
      
      fcout = new FldArrayF(nelts,nfldc);
      *fcout = *fc;//copie
    }
    // array en centres
    tpl = K_ARRAY::buildArray3(*fcout, varStringC, *cnout, "NGON", api);
    delete cnout; delete fcout;
  }

  RELEASESHAREDU(arrayC, fc, cnc);
  return tpl;
}  
//=============================================================================
/* Add ghost cells to a NGON array (both nodes and centers are modified) 
   Fields located at nodes are not modified
*/
//=============================================================================
PyObject* K_CONVERTER::addGhostCellsNGonBoth(PyObject* self, PyObject* args)
{
  PyObject *arrayN, *arrayC;
  E_Int depth;
  if (!PYPARSETUPLE_(args, OO_ I_, &arrayN, &arrayC, &depth)) return NULL;

  if (depth < 1) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "addGhostCells: depth must be > 0.");
    return NULL;
  }
  E_Int ni, nj, nk;
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int res = K_ARRAY::getFromArray3(arrayN, varString, f, ni, nj, nk, 
                                     cn, eltType);
  if (res != 2)
  {
    if (res == 1) RELEASESHAREDS(arrayN, f);
    PyErr_SetString(PyExc_TypeError, 
                    "addGhostCells: array is invalid.");
    return NULL;
  }
  if (strcmp(eltType, "NGON") != 0 ) 
  {
    RELEASESHAREDU(arrayN, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "addGhostCells: array must be NGON.");
    return NULL;
  }
  
  E_Int nic, njc, nkc;
  char* varStringC; char* eltTypec;
  FldArrayF* fc; FldArrayI* cnc;
  E_Int resc = K_ARRAY::getFromArray3(arrayC, varStringC, fc, nic, njc, nkc, 
                                      cnc, eltTypec);
  if (resc != 2)
  {
    if (resc == 1) RELEASESHAREDS(arrayC, fc);
    PyErr_SetString(PyExc_TypeError,
                    "addGhostCells: array is invalid.");
    RELEASESHAREDU(arrayN, f, cn); return NULL;
  }
  if (strcmp(eltTypec, "NGON*") != 0 ) 
  {
    RELEASESHAREDU(arrayN, f, cn);
    RELEASESHAREDU(arrayC, fc, cnc);
    PyErr_SetString(PyExc_TypeError,
                    "addGhostCells: array must be NGON*.");
    return NULL;
  }

  /* Determination des faces externes du maillage reel */
  E_Int sizeFN2 = 0;
  E_Int sizeEF2 = 0;
  // Connectivite Faces/Elts reelle
  FldArrayI cFEp; K_CONNECT::connectNG2FE(*cn, cFEp);
  vector<E_Int> facesExt;//faces externes du maillage reel: demarre a 0
  E_Int* cFEp1 = cFEp.begin(1);
  E_Int* cFEp2 = cFEp.begin(2);
  E_Int* cNGp = cn->begin();
  E_Int nfacesp = cNGp[0];// nb de faces ds le maillage reel
  E_Int sizeFN = cNGp[1];
  E_Int sizeEF = cNGp[3+sizeFN];
  E_Int indface, etl, etg, etd, ce, nvertp;
  E_Int dimNGON, neltsAdd, nfacesAdd, nfacesExt, sizeOut, sizefout;

  if (cNGp[2] == 2) dimNGON = 2;
  else if (cNGp[2] > 2) dimNGON = 3;
  else 
  {
    PyErr_SetString(PyExc_TypeError,
                    "addGhostCells: array must be 2D or 3D.");
    RELEASESHAREDU(arrayN,f,cn); 
    RELEASESHAREDU(arrayC,fc,cnc);
    return NULL;
  }

  E_Int api = f->getApi();
  E_Int npts = f->getSize();
  E_Int nfld = f->getNfld();
  E_Int nfldc = fc->getNfld();
  E_Int nelts = fc->getSize();
  E_Int* ptr = cNGp+2;
  FldArrayI* cnout; FldArrayF* fcout; FldArrayF* fout; 
  E_Int nptsAdd = 0;

  // Sortie
  PyObject* l = PyList_New(0);
  if ( dimNGON == 3)
  {
    for (E_Int nof = 0; nof < nfacesp; nof++)
    {
      nvertp = ptr[0];
      etg = cFEp1[nof]; etd = cFEp2[nof];
      if (etd == 0 || etg == 0) //face exterieure
      {
        facesExt.push_back(nof);  nptsAdd += nvertp;
        sizeFN2 += (5*nvertp)+(nvertp+1); // 4 vertex par face laterale + nvert sur la face opposee +1 pour dimensionner
        sizeEF2 += (2+nvertp+1);// nvertp faces laterales + 2 faces + 1 pour dimensionner
      }         
      ptr += nvertp+1;
    }
    sizeFN2 *= depth;
    sizeEF2 *= depth;
    nfacesExt = facesExt.size();
    neltsAdd = nfacesExt*depth; //nb d elts ajoutes
    nptsAdd *= depth;
    sizeOut = sizeFN+sizeEF+sizeFN2+sizeEF2+4;
    sizefout = npts+nptsAdd;
    cnout = new FldArrayI(sizeOut);cnout->setAllValuesAtNull();
    fout = new FldArrayF(sizefout,nfld); fout->setAllValuesAtNull();
    addGhostCellsNGon3D(depth, *f, *cn, neltsAdd, sizeFN2, sizeEF2, facesExt, fout, cnout);

    fcout = new FldArrayF(nelts+neltsAdd,nfldc);fcout->setAllValuesAtNull();
    for (E_Int eq = 1; eq <= nfldc; eq++)
    {
      E_Float* fcp = fc->begin(eq);
      E_Float* fcpo = fcout->begin(eq);
      for (E_Int et = 0; et < nelts; et++)
        fcpo[et] = fcp[et];
      
      ce = nelts;
      for (E_Int nof = 0; nof < nfacesExt; nof++)
      {
        indface = facesExt[nof];
        etg = cFEp1[indface]; etd = cFEp2[indface];
        
        if (etg != 0) 
        { 
          etl = etg-1;
          for (E_Int d = 1; d <= depth; d++)
            fcpo[ce+(d-1)*nfacesExt] = fcp[etl]; 
          ce++;
        }
        else if (etd != 0) 
        {
          etl = etd-1;
          for (E_Int d = 1; d <= depth; d++)
            fcpo[ce+(d-1)*nfacesExt] = fcp[etl]; 
          ce++;
        }
      }
    }
    // array en noeuds
    PyObject* tpl = K_ARRAY::buildArray3(*fout, varString, *cnout, eltType, api);
    PyList_Append(l, tpl); Py_DECREF(tpl); delete fout;
    // array en centres
    tpl = K_ARRAY::buildArray3(*fcout, varStringC, *cnout, eltType, api);
    PyList_Append(l, tpl); Py_DECREF(tpl); delete cnout; delete fcout;    
  }
  else // dimNGON = 2 
  {
    for (E_Int nof = 0; nof < nfacesp; nof++)
    {
      etg = cFEp1[nof]; etd = cFEp2[nof];
      if (etd == 0 || etg == 0) //face exterieure
      {facesExt.push_back(nof); nptsAdd +=2;}
    }
    nfacesExt = facesExt.size();
    nfacesAdd = 3*depth;// 1 face ajoutee identique a la face externe + 2 sur les cotes
    sizeFN2 = nfacesAdd*nfacesExt*(2+1);// 2 vertex par face * nb de faces ajoutees
    neltsAdd = nfacesExt*depth;
    sizeEF2 = neltsAdd*5;// nb Elts ajoutes*4 faces par cellule fictive+1 pour dimensionner  
    sizeOut = sizeFN+sizeEF+sizeFN2+sizeEF2+4;
    nptsAdd *= depth;
    sizefout = npts+nptsAdd;

    if (nfacesExt > 0)
    {
      cnout = new FldArrayI(sizeOut);cnout->setAllValuesAtNull();
      fout = new FldArrayF(sizefout,nfld);fout->setAllValuesAtNull();
      addGhostCellsNGon2D(depth, *f, *cn, neltsAdd, sizeFN2, sizeEF2, facesExt, fout, cnout);      
      E_Int nfldc = fc->getNfld();
      E_Int nelts = fc->getSize();
      fcout = new FldArrayF(nelts+neltsAdd,nfldc); fcout->setAllValuesAtNull();
      for (E_Int eq = 1; eq <= nfldc; eq++)
      {
        E_Float* fcp = fc->begin(eq);
        E_Float* fcpo = fcout->begin(eq);
        for (E_Int et = 0; et < nelts; et++)
          fcpo[et] = fcp[et];
        
        ce = nelts;
        for (E_Int nof = 0; nof < nfacesExt; nof++)
        {
          indface = facesExt[nof];
          etg = cFEp1[indface]; etd = cFEp2[indface];

          if (etg != 0) 
          {
            etl = etg-1;
            for (E_Int d = 1; d <= depth; d++)
            { fcpo[ce+(d-1)*nfacesExt] = fcp[etl];} 
            ce++;
          }
          else if (etd != 0) 
          {
            etl = etd-1;
            for (E_Int d = 1; d <= depth; d++)
            { fcpo[ce+(d-1)*nfacesExt] = fcp[etl]; }
            ce++;
          }         
        }
      }
    }
    else
    {
      cnout = new FldArrayI(cn->getSize()*cn->getNfld());
      *cnout = *cn; 
      
      fcout = new FldArrayF(nelts,nfldc);
      *fcout = *fc;//copie
      
      fout = new FldArrayF(npts, nfld);
      *fout = *f;
    }

    // array en noeuds
    PyObject* tpl = K_ARRAY::buildArray3(*fout, varString, *cnout, eltType, api);
    PyList_Append(l, tpl); Py_DECREF(tpl); delete fout;
    //array en centres
    tpl = K_ARRAY::buildArray3(*fcout, varStringC, *cnout, eltType, api);
    PyList_Append(l, tpl); Py_DECREF(tpl); delete cnout; delete fcout;
  }

  RELEASESHAREDU(arrayN, f, cn); RELEASESHAREDU(arrayC, fc, cnc);
  return l;
}  

//=============================================================================
/* Add ghost cells for a 2D NGON - Redimensionne cno et fo */
//=============================================================================
void K_CONVERTER::addGhostCellsNGon2D(E_Int depth,  
                                      FldArrayF& f, FldArrayI& cn, 
                                      E_Int neltsAdd, 
                                      E_Int sizeFN2, E_Int sizeEF2, vector<E_Int>& facesExt,
                                      FldArrayF*& fo, FldArrayI*& cno)
{
  // Construction de fo :
  E_Int npts = f.getSize(); E_Int nfld = f.getNfld();
  for (E_Int eq = 1; eq <= nfld; eq++)
  {
    E_Float* fp = f.begin(eq);
    E_Float* fpo = fo->begin(eq);
    // 1/ copie de f pour les npts premiers pts
    for (E_Int ind = 0; ind < npts; ind++)
      fpo[ind] = fp[ind];
  }
  FldArrayI posFace; K_CONNECT::getPosFaces(cn, posFace);
  E_Int* cNGp = cn.begin();
  
  E_Int nfacesp = cNGp[0];// nb de faces ds le maillage reel
  E_Int sizeFNp = cNGp[1]; 
  E_Int neltsp = cNGp[2+sizeFNp];
  E_Int sizeEFp = cNGp[3+sizeFNp];
  
  E_Int nfacesExt = facesExt.size();

  //----------
  //Construction des elts degeneres en QUAD
  FldArrayI cFNp2(sizeFN2); FldArrayI cEFp2(sizeEF2);

  E_Int* ptr1;
  E_Int* ptrFN2 = cFNp2.begin(); 
  E_Int* ptrEF2 = cEFp2.begin();
  E_Int cf = 0;//compteur sur les faces 
  E_Int ce = 0; // compteur sur les elts  
  E_Int indv, indvd, indface, posface;
  FldArrayI tagV(2*nfacesp*depth); 

  /* 1. creation de l elt QUAD contenant ces faces degenerees
     1.1. on cree les faces opposees a la face externe 
     1.2. on cree les faces laterales  */
  E_Int sizeFNP = 0; E_Int sizeEFP = 0;
  FldArrayI facesOpp(nfacesExt); 
  for (E_Int nof = 0; nof < nfacesExt; nof++)
  {facesOpp[nof] = facesExt[nof]+1;}
  FldArrayI indirN(npts); FldArrayI indirP(npts);
  for (E_Int i = 0; i < npts; i++) indirP[i] = i;
  E_Int indstart = npts;
  E_Int indorig;
  for (E_Int d = 1; d <= depth; d++)
  {
    indirN.setAllValuesAtNull();
    tagV.setAllValuesAtNull();
    for (E_Int nof = 0; nof < nfacesExt; nof++)
    {
      indface = facesExt[nof];
      posface = posFace[indface];
      ptr1 = cNGp+posface;// connectivite Faces/Noeuds du maillage reel
    
      //Construction des 4 faces pour l element
      ptrEF2[0] = 4;
      ptrEF2[1] = facesOpp[nof];// face externe
    
      // Creation de la face opposee de 1er niveau
      ptrFN2[0] = 2;
      for (E_Int nov = 1; nov <= 2; nov++) 
      {
        indorig = ptr1[nov]-1;
        if ( indirN[indorig] == 0 )//pt pas encore insere
        {
          indirN[indorig] = indstart;
          //on ajoute les n champs
          for (E_Int eq = 1; eq <= nfld; eq++)
            (*fo)(indstart,eq) = f(indorig,eq);
          indstart++;
          ptrFN2[nov] = indstart;
        }
        else
        {
          ptrFN2[nov] = indirN[indorig]+1;
        }
      }
      ptrFN2 += 3; sizeFNP+=3;
      cf++;  

      ptrEF2[3] = nfacesp+cf;
      facesOpp[nof] = nfacesp+cf;

      // Creation des faces laterales
      for (E_Int nov = 1; nov <= 2; nov++)
      {
        indv = ptr1[nov];
        indvd = indv-1;
        if (tagV[indvd] == 0)
        {       
          ptrFN2[0] = 2; 
          indorig = indv-1;
          // ptrFN2[1] = indv; ptrFN2[2] = indv;
          ptrFN2[1] = indirP[indorig]+1;
          ptrFN2[2] = indirN[indorig]+1;
          
          ptrFN2 += 3; sizeFNP+=3;
          cf++;

          ptrEF2[2*nov] = nfacesp+cf;
          tagV[indvd] = ptrEF2[2+nov]; // pour retrouver le numero de la face 
        }
        else 
          ptrEF2[2*nov] = tagV[indvd];      
      }
      ptrEF2 += 5; sizeEFP+= 5; ce++;        
    }
    indirP = indirN;
  }  
  // Creation de la nouvelle connectivite avec degenerescences 
  nfacesp += cf; neltsp += ce;

  E_Int sizeFN = sizeFNp+sizeFNP;
  E_Int sizeEF = sizeEFp+sizeEFP;
  cno->resize(sizeEF+sizeFN+4);
  E_Int* cnop = cno->begin();
  
  cnop[0] = nfacesp; cnop[1] = sizeFN;

  E_Int* ptrFNp = cn.begin()+2;
  E_Int* ptrEFp = cn.begin()+4+sizeFNp;

  cFNp2.resize(sizeFNP); ptrFN2 = cFNp2.begin(); 
  cEFp2.resize(sizeEFP); ptrEF2 = cEFp2.begin();

  cnop = cno->begin()+2;
  for (E_Int i = 0; i < sizeFNp; i++)
  {cnop[0] = ptrFNp[i]; cnop++; }
  
  for (E_Int i = 0; i < sizeFNP; i++)
  {cnop[0] = ptrFN2[i];  cnop++; }
  
  cnop = cno->begin()+2+sizeFN;
  cnop[0] = neltsp; //cEF initial
  cnop[1] = sizeEF; 
  cnop+=2;
  
  for (E_Int i = 0; i < sizeEFp; i++)
  {cnop[0] = ptrEFp[i]; cnop++; }
  
  for (E_Int i = 0; i < sizeEFP; i++)
  {cnop[0] = ptrEF2[i]; cnop++;}
  
  fo->reAllocMat(indstart,nfld);
}
//=============================================================================
/* Add ghost cells for a 3D NGON       
   redimensionne cno */
//=============================================================================
void K_CONVERTER::addGhostCellsNGon3D(E_Int depth,
                                      FldArrayF& f, FldArrayI& cn, 
                                      E_Int neltsAdd,
                                      E_Int sizeFN2, E_Int sizeEF2, vector<E_Int>& facesExt,
                                      FldArrayF*& fo, FldArrayI*& cno)
{
  // Construction de fo :
  E_Int npts = f.getSize(); E_Int nfld = f.getNfld();
  for (E_Int eq = 1; eq <= nfld; eq++)
  {
    E_Float* fp = f.begin(eq);
    E_Float* fpo = fo->begin(eq);
    // 1/ copie de f pour les npts premiers pts
    for (E_Int ind = 0; ind < npts; ind++)
      fpo[ind] = fp[ind];
  }
  FldArrayI posFace; K_CONNECT::getPosFaces(cn, posFace);
  E_Int* cNGp = cn.begin();
  
  E_Int nfacesp = cNGp[0];// nb de faces ds le maillage reel
  E_Int sizeFNp = cNGp[1]; 
  E_Int neltsp = cNGp[2+sizeFNp];
  E_Int sizeEFp = cNGp[3+sizeFNp];
  
  E_Int nfacesExt = facesExt.size();
  //----------  
  FldArrayI cFNp2(sizeFN2); FldArrayI cEFp2(sizeEF2);

  E_Int* ptr1;
  E_Int* ptrFN2 = cFNp2.begin(); 
  E_Int* ptrEF2 = cEFp2.begin();
  E_Int nvertp, indface, posface, indv1, indv2;
 
  /* 1. Creation du NGON 2D des faces exterieures facesExt
     Les edges 3D sont des faces du NGON2D 
     Les faces sont uniques - 
     Le numero de l elt NGON2D correspond au numero de la face ext */
  E_Int nvertpmax = 0;
  for (E_Int nof = 0; nof < nfacesExt; nof++)
  {
    indface = facesExt[nof];
    posface = posFace[indface];
    ptr1 = cNGp+posface;// connectivite Faces/Noeuds du maillage reel
    nvertp = ptr1[0];
    nvertpmax = K_FUNC::E_max(nvertp,nvertpmax);
  }
  E_Int nedgesmax = nvertpmax*nfacesExt;
  vector< vector<E_Int> > edgesOfFaces(nfacesExt); // liste des edges par face externe
  FldArrayI vertexOfEdges(nedgesmax,2); vertexOfEdges.setAllValuesAt(-1);// numero des sommets de l edge
  E_Int nedges = 0;
  E_Int found = 0;
  E_Int indvmin, indvmax;
  E_Int* vertexOfEdge1 = vertexOfEdges.begin(1);
  E_Int* vertexOfEdge2 = vertexOfEdges.begin(2);
  for (E_Int nof = 0; nof < nfacesExt; nof++)
  {
    indface = facesExt[nof];
    posface = posFace[indface];
    ptr1 = cNGp+posface;// connectivite Faces/Noeuds du maillage reel
    nvertp = ptr1[0];
    vector<E_Int>& edgesOfFaceLoc = edgesOfFaces[nof];
    
    for (E_Int nov = 1; nov < nvertp; nov++)
    {
      indv1 = ptr1[nov]; indv2 = ptr1[nov+1];
      indvmin = K_FUNC::E_min(indv1, indv2);
      indvmax = K_FUNC::E_max(indv1, indv2);
      // recherche si l edge est deja dans la liste des edges
      // on pourrait ameliorer encore en regardant si une face voisine a deja ete ajoutee
      found = 0;
      for (E_Int v1 = 0; v1 < nedges; v1++)
      {
        if ( vertexOfEdge1[v1] == indvmin && vertexOfEdge2[v1] == indvmax) 
        { 
          found = 1; 
          edgesOfFaceLoc.push_back(v1+1);
          break;
        }        
      }
      if ( found == 0) 
      {
        vertexOfEdge1[nedges] = indvmin; 
        vertexOfEdge2[nedges] = indvmax;
        edgesOfFaceLoc.push_back(nedges+1);
        nedges++;
      }
    }
    // pour cycler
    indv1 = ptr1[nvertp]; indv2 = ptr1[1];
    indvmin = K_FUNC::E_min(indv1, indv2);
    indvmax = K_FUNC::E_max(indv1, indv2);
    // recherche si l edge est deja dans la liste des edges
    // on pourrait ameliorer encore en regardant si une face voisine a deja ete ajoutee
    found = 0;
    for (E_Int v1 = 0; v1 < nedges; v1++)
    {
      if ( vertexOfEdge1[v1] == indvmin && vertexOfEdge2[v1] == indvmax) 
      { 
        found = 1; 
        edgesOfFaceLoc.push_back(v1+1);
        break;
      }        
    }
    if ( found == 0) 
    {
      vertexOfEdge1[nedges] = indvmin; 
      vertexOfEdge2[nedges] = indvmax;
      edgesOfFaceLoc.push_back(nedges+1);
      nedges++;
    }
    // Fin ajout de l edge
  }
 
  vertexOfEdges.reAllocMat(nedges,2);
  vertexOfEdge1 = vertexOfEdges.begin(1);
  vertexOfEdge2 = vertexOfEdges.begin(2);
  FldArrayI tagE(nedges); //stocke le numero de la face quad creee a partir de l edge
  E_Int nedgesOfFace, indedge, countEF;
  /* 2. creation de l elt NGON contenant ces faces degenerees
     2.1. on cree les faces opposees a la face externe 
     2.2. on cree les faces laterales  */
  E_Int sizeFNP = 0; E_Int sizeEFP = 0;
  FldArrayI facesOpp(nfacesExt); 
  // init faces opposees
  for (E_Int nof = 0; nof < nfacesExt; nof++)
  {facesOpp[nof] = facesExt[nof]+1;}
  FldArrayI indirN(npts); FldArrayI indirP(npts);
  for (E_Int i = 0; i < npts; i++) indirP[i] = i;
  E_Int indstart = npts;
  E_Int indorig, indorig1, indorig2;

  E_Int cf = 0;//compteur sur les faces 
  E_Int ce = 0; // compteur sur les elts 
  for (E_Int d = 1; d <= depth; d++)
  {
    tagE.setAllValuesAtNull(); 
    indirN.setAllValuesAtNull();

    for (E_Int nof = 0; nof < nfacesExt; nof++)
    { 
      indface = facesExt[nof];
      posface = posFace[indface];
      ptr1 = cNGp+posface; // connectivite Faces/Noeuds du maillage reel
      nvertp = ptr1[0];    
      ptrEF2[0] = nvertp+2; //nb de faces ds l elt + 1            
      // (d-1)-ieme face opposee a la face externe
      ptrEF2[1] = facesOpp[nof];
      // Creation de la d-ieme face opposee a la face externe 
      ptrFN2[0] = nvertp;
      // for (E_Int nov = 1; nov <= nvertp; nov++) ptrFN2[nov] = ptr1[nov];
      for (E_Int nov = 1; nov <= nvertp; nov++) 
      {
        indorig = ptr1[nov]-1;
        if (indirN[indorig] == 0) //pt pas encore traite
        {
          indirN[indorig] = indstart;
          //on ajoute les n champs
          for (E_Int eq = 1; eq <= nfld; eq++)
            (*fo)(indstart,eq) = f(indorig,eq);
          indstart++;
          ptrFN2[nov] = indstart;
        }
        else
        {
          ptrFN2[nov] = indirN[indorig]+1;
        }
      }      
      ptrFN2 += nvertp+1; sizeFNP+=nvertp+1;
      cf++;
        
      ptrEF2[2] = nfacesp+cf; facesOpp[nof] = nfacesp+cf;
      countEF = 3;

      //Creation des faces laterales QUAD (degenerees en BAR)
      vector<E_Int>& edgesOfFaceLoc = edgesOfFaces[nof];
      nedgesOfFace = edgesOfFaceLoc.size();
      for (E_Int noedge = 0; noedge < nedgesOfFace; noedge++)
      {
        indedge = edgesOfFaceLoc[noedge]-1;      
        if (tagE[indedge] == 0) 
        {
          //on cree la face quad
          indv1 = vertexOfEdge1[indedge];
          indv2 = vertexOfEdge2[indedge];
          ptrFN2[0] = 4;
          indorig1 = indv1-1; indorig2 = indv2-1;
          // ptrFN2[1] = indv1; ptrFN2[2] = indv2; ptrFN2[3] = indv2; ptrFN2[4] = indv1;
          ptrFN2[1] = indirP[indorig1]+1; ptrFN2[2] = indirP[indorig2]+1; 
          ptrFN2[3] = indirN[indorig2]+1; ptrFN2[4] = indirN[indorig1]+1;
          ptrFN2+= 5; sizeFNP+= 5;
          cf++;

          ptrEF2[countEF] = nfacesp+cf;
          tagE[indedge] = nfacesp+cf;// numero de la face quad creee
          countEF++;
        }
        else // on recupere la face quad deja creee 
        {
          ptrEF2[countEF] = tagE[indedge];
          countEF++;
        }
      }
      ptrEF2 += nvertp+3; sizeEFP+= nvertp+3; ce++;    
    }
    indirP = indirN;
  }
  // Creation de la nouvelle connectivite avec degenerescences 
  nfacesp += cf; neltsp += ce;
  E_Int sizeFN = sizeFNp+sizeFNP;
  E_Int sizeEF = sizeEFp+sizeEFP;
  cno->resize(sizeEF+sizeFN+4);
  E_Int* cnop = cno->begin();
  
  cnop[0] = nfacesp; cnop[1] = sizeFN;

  E_Int* ptrFNp = cn.begin()+2;
  E_Int* ptrEFp = cn.begin()+4+sizeFNp;

  cFNp2.resize(sizeFNP); ptrFN2 = cFNp2.begin(); 
  cEFp2.resize(sizeEFP); ptrEF2 = cEFp2.begin();

  cnop = cno->begin()+2;
  for (E_Int i = 0; i < sizeFNp; i++)
  {cnop[0] = ptrFNp[i]; cnop++; }
  
  for (E_Int i = 0; i < sizeFNP; i++)
  {cnop[0] = ptrFN2[i];  cnop++; }
  
  cnop = cno->begin()+2+sizeFN;
  cnop[0] = neltsp; //cEF initial
  cnop[1] = sizeEF; 
  cnop += 2;
  
  for (E_Int i = 0; i < sizeEFp; i++)
  {cnop[0] = ptrEFp[i]; cnop++;}
  
  for (E_Int i = 0; i < sizeEFP; i++)
  {cnop[0] = ptrEF2[i]; cnop++;}
  fo->reAllocMat(indstart,nfld);
}
