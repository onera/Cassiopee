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
# include "converter.h"

using namespace K_FLD;

//=============================================================================
/* Convert a ElsaHybrid node
   method 0: TRI-int,QUAD-int,TRI-ext,QUAD-ext
   method 1: int, ext */
//=============================================================================
PyObject* K_CONVERTER::createElsaHybrid(PyObject* self, PyObject* args)
{
  PyObject* NGon; PyObject* PE; PyObject* ict; PyObject* bcct;
  PyObject* inct;
  PyObject* ox; PyObject* oy; PyObject* oz;
  E_Int axe2D, method;
  if (!PYPARSETUPLEI(args, "OOOOOllOOO", "OOOOOiiOOO", 
                     &NGon, &PE, &ict, &bcct, &inct, 
                     &method, &axe2D, &ox, &oy, &oz))
    return NULL;
  // Check numpy NGon
  FldArrayI* cNGon;
  E_Int res = K_NUMPY::getFromNumpyArray(NGon, cNGon, true);
  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "createElsaHybrid: numpy is invalid.");
    return NULL;
  }

  // Check numpy PE (ParentElements)
  FldArrayI* cPE;
  res = K_NUMPY::getFromNumpyArray(PE, cPE, true);
  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "createElsaHybrid: numpy is invalid.");
    return NULL;
  }

  // Check numpy ICT (InverseCrossTable)
  FldArrayI* cICT;
  res = K_NUMPY::getFromNumpyArray(ict, cICT, true);
  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "createElsaHybrid: numpy is invalid.");
    return NULL;
  }

  // Check numpy BCCT (BCCrossTable)
  FldArrayI* cBCCT;
  res = K_NUMPY::getFromNumpyArray(bcct, cBCCT, true);
  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "createElsaHybrid: numpy is invalid.");
    return NULL;
  }

  // Check numpy INCT (IndexNGONCrossTable)
  FldArrayI* cINCT;
  res = K_NUMPY::getFromNumpyArray(inct, cINCT, true);
  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "createElsaHybrid: numpy is invalid.");
    return NULL;
  }

  E_Float* x=NULL; E_Float* y=NULL; E_Float* z=NULL;
  if (axe2D > 0)
  {
    E_Int np; E_Int nfld;
    res = K_NUMPY::getFromNumpyArray(ox, x, np, nfld, true);
    res = K_NUMPY::getFromNumpyArray(oy, y, np, nfld, true);
    res = K_NUMPY::getFromNumpyArray(oz, z, np, nfld, true);
  }

  E_Int nfaces = cICT->getSize();
  E_Int* ptNGon = cNGon->begin();
  FldArrayI posFaces(nfaces);
  K_CONNECT::getPosFacets(ptNGon, 0, nfaces, posFaces);
  E_Int* posFacesp = posFaces.begin();

  //printf("here nfaces=%d\n", nfaces);
  FldArrayI CT(nfaces);
  E_Int iTRI=0, iQUAD=0, eTRI=0, eQUAD=0;
  E_Int ints=0, exts=0;

  E_Int* ptr;
  E_Int np;

  if (method == 0)
  {
    // Compte
    E_Int* PEG = cPE->begin(1);
    E_Int* PED = cPE->begin(2);

    for (E_Int i = 0; i < nfaces; i++)
    {
      ptr = ptNGon+posFacesp[i];
      np = ptr[0];
      if (np == 3)
      {
        if (PED[i]*PEG[i] == 0) eTRI++;
        else iTRI++;
      }
      else if (np == 4)
      {
        if (PED[i]*PEG[i] == 0) eQUAD++;
        else iQUAD++;
      }
    }

    //printf("here %d %d %d %d\n", iTRI, iQUAD, eTRI, eQUAD);

    // Remplit CT
    E_Int* piTRI = CT.begin();
    E_Int* piQUAD = CT.begin()+iTRI;
    E_Int* peTRI = CT.begin()+iTRI+iQUAD;
    E_Int* peQUAD = CT.begin()+iTRI+iQUAD+eTRI;

    iTRI=0; iQUAD=0; eTRI=0; eQUAD=0;

    for (E_Int i = 0; i < nfaces; i++)
    {
      ptr = ptNGon+posFacesp[i];
      np = ptr[0];
      if (np == 3)
      {
        if (PED[i]*PEG[i] == 0) { peTRI[eTRI] = i+1; eTRI++; }
        else { piTRI[iTRI] = i+1; iTRI++; }
      }
      else if (np == 4)
      {
        if (PED[i]*PEG[i] == 0) { peQUAD[eQUAD] = i+1; eQUAD++; }
        else { piQUAD[iQUAD] = i+1; iQUAD++; }
      }
      else printf("Face no=%d non TRI et non QUAD (%d)\n",i,np);
    }
    ints = iTRI+iQUAD; exts = eTRI+eQUAD;
  }
  else // method 1
  {
    // Compte les faces int et ext
    E_Int* PEG = cPE->begin(1);
    E_Int* PED = cPE->begin(2);

    for (E_Int i = 0; i < nfaces; i++)
    {
      ptr = ptNGon+posFacesp[i];
      if (PED[i]*PEG[i] == 0) exts++;
      else ints++;
    }
    //printf("here %d %d - %d %d\n", ints, exts, nfaces, ints+exts);

    // Remplit CT
    E_Int* pints = CT.begin();
    E_Int* pexts = CT.begin()+ints;
    ints=0; exts=0;

    for (E_Int i = 0; i < nfaces; i++)
    {
      ptr = ptNGon+posFacesp[i];
      if (PED[i]*PEG[i] == 0) { pexts[exts] = i+1; exts++; }
      else { pints[ints] = i+1; ints++; }
    }
    iTRI=ints; iQUAD=ints; eTRI=exts; eQUAD=exts; // only for output
  }

  // inverse CT
  E_Int* pICT = cICT->begin();
  E_Int* pCT = CT.begin();
  for (E_Int i = 0; i < nfaces; i++)
  {
    //if (pCT[i]-1 <= 0 || pCT[i]-1 > cICT->getSize())
    //{ printf("%d %d %d\n",i,pCT[i],cICT->getSize());}
    pICT[pCT[i]-1] = i+1;
  }

  // Remplit BCCT
  E_Int* pt = cBCCT->begin();

  // pour les faces interieures (indirection directe)
  for (E_Int i = 0; i < ints; i++)
  {
    pt[i] = i+1; // direct
  }

  if (axe2D == 1) // x
  {
    // detecte les faces x=cte
    E_Float xmin, xmax; E_Int ind;
    E_Int cpt = ints+1;
    for (E_Int i = ints; i < exts+ints; i++)
    {
      ind = CT[i]-1;
      ptr = ptNGon+posFacesp[ind];
      np = ptr[0];
      xmin = 1.e6; xmax = -1.e6;
      for (E_Int j = 0; j < np; j++)
      {
        ind = ptr[j+1]-1;
        xmin = K_FUNC::E_min(x[ind], xmin);
        xmax = K_FUNC::E_max(x[ind], xmax);
      }
      if (xmax-xmin < 1.e-11) pt[i] = -1;
      else { pt[i] = cpt; cpt++; }
    }
  }
  else if (axe2D == 2) // y
  {
    // detecte les faces y=cte
    E_Float ymin, ymax; E_Int ind;
    E_Int cpt = ints+1;
    for (E_Int i = ints; i < ints+exts; i++)
    {
      ind = CT[i]-1;
      ptr = ptNGon+posFacesp[ind];
      np = ptr[0];
      ymin = 1.e6; ymax = -1.e6;
      for (E_Int j = 0; j < np; j++)
      {
        ind = ptr[j+1]-1;
        ymin = K_FUNC::E_min(y[ind], ymin);
        ymax = K_FUNC::E_max(y[ind], ymax);
      }
      if (ymax-ymin < 1.e-11) pt[i] = -1;
      else { pt[i] = cpt; cpt++; }
    }
  }
  else if (axe2D == 3) // z
  {
    // detecte les faces z=cte
    E_Float zmin, zmax; E_Int ind;
    E_Int cpt = ints+1;
    for (E_Int i = ints; i < ints+exts; i++)
    {
      ind = CT[i]-1;
      ptr = ptNGon+posFacesp[ind];
      np = ptr[0];
      zmin = 1.e6; zmax = -1.e6;
      for (E_Int j = 0; j < np; j++)
      {
        ind = ptr[j+1]-1;
        zmin = K_FUNC::E_min(z[ind], zmin);
        zmax = K_FUNC::E_max(z[ind], zmax);
      }
      if (zmax-zmin < 1.e-11) pt[i] = -1;
      else { pt[i] = cpt; cpt++; }
    }
  }
  else
  { // default
    for (E_Int i = ints; i < ints+exts; i++)
    {
      pt[i] = i+1; // direct
    }
  }

  // IndexNGONCrossTable (c'est le posFace)
  E_Int* pINCT = cINCT->begin();
  for (E_Int i = 0; i < nfaces; i++) pINCT[i]=posFacesp[i];

  RELEASESHAREDN(NGon, cNGon);
  RELEASESHAREDN(PE, cPE);
  RELEASESHAREDN(ict, cICT);
  RELEASESHAREDN(bcct, cBCCT);
  return Py_BuildValue("(llll)", iTRI, iQUAD, eTRI, eQUAD);
}
