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

// stack the planes of two meshes (with same nixnj) into a single mesh

#include "generator.h"
using namespace K_FLD;
using namespace std;

//=============================================================================
// Conactenate mesh layers
//=============================================================================
PyObject* K_GENERATOR::stackMesh(PyObject* self, PyObject* args)
{
  PyObject* arrays;
  if (!PYPARSETUPLE_(args, O_, &arrays)) return NULL;

  vector<E_Int> res;
  vector<char*> structVarString; vector<char*> unstructVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstructF;
  vector<E_Int> ni; vector<E_Int> nj; vector<E_Int> nk;
  vector<FldArrayI*> cn; vector<char*> eltType;
  vector<PyObject*> objs; vector<PyObject*> obju;
  K_ARRAY::getFromArrays(arrays,
                         res, structVarString, unstructVarString,
                         structF, unstructF, ni, nj, nk, cn,
                         eltType, objs, obju,
                         true, false, false, false, true);

  E_Int ns = structF.size();
  E_Int nu = unstructF.size();
  PyObject* tpl=NULL;

  if (ns >= 2) // concatenate all structs, il doivent etre tous ni*nj et memes champs
  {
    E_Int nfld = structF[0]->getNfld();
    E_Int api = structF[0]->getApi();
    E_Int ni1 = ni[0]; E_Int nj1 = nj[0];
    E_Int ni1nj1 = ni1*nj1;
    E_Int nk1 = 0;
    for (E_Int p = 0; p < ns; p++) nk1 += nk[p];
    tpl = K_ARRAY::buildArray3(nfld, structVarString[0], ni1, nj1, nk1, api);
    FldArrayF* coord;
    K_ARRAY::getFromArray3(tpl, coord);

    nk1 = 0;
    for (E_Int p = 0; p < ns; p++)
    {
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* f1p = structF[p]->begin(n);
        E_Float* coordp = coord->begin(n);
        for (E_Int k = 0; k < nk[p]; k++)
          for (E_Int i = 0; i < ni1nj1; i++)
          {
            coordp[i+(k+nk1)*ni1nj1] = f1p[i+k*ni1nj1];
          }
      }
      nk1 += nk[p];
    }
    RELEASESHAREDS(tpl, coord);
  }
  else if (nu >= 2) // concatenate all unstructs
  {
    char outEltType[10];
    if (strcmp(eltType[0], "BAR") == 0) { strcpy(outEltType, "QUAD"); }
    if (strcmp(eltType[0], "QUAD") == 0) { strcpy(outEltType, "HEXA"); }
    if (strcmp(eltType[0], "TRI") == 0) { strcpy(outEltType, "PENTA"); }
  
    E_Int nfld = unstructF[0]->getNfld();
    E_Int api = unstructF[0]->getApi();
    E_Int ne = cn[0]->getSize();
    E_Int nv1 = 0;
    for (E_Int p = 0; p < nu; p++) nv1 += unstructF[p]->getSize();
    
    E_Int nv0 = unstructF[0]->getSize();
    E_Int ne1 = (nu-1)*cn[0]->getSize();

    tpl = K_ARRAY::buildArray3(nfld, unstructVarString[0], nv1, ne1,
                               outEltType, false, api);
    FldArrayF* coord; FldArrayI* cno;
    K_ARRAY::getFromArray3(tpl, coord, cno);

    for (E_Int p = 0; p < nu; p++)
    {
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* f1p = unstructF[p]->begin(n);
        E_Float* coordp = coord->begin(n);
        for (E_Int i = 0; i < nv0; i++)
        {
          coordp[i+nv0*p] = f1p[i];
        }
      }
    }
    
    for (E_Int p = 0; p < nu-1; p++)
    {
      if (strcmp(eltType[0], "BAR") == 0)
      {
        FldArrayI& cn0 = *cn[p];
        FldArrayI& cn1 = *cn[p+1];
        FldArrayI& cmo = *(cno->getConnect(0));

        for (E_Int i = 0; i < ne; i++)
        {
          cmo(i + ne*p, 1) = cn0(i, 1) + nv0*p;
          cmo(i + ne*p, 2) = cn0(i, 2) + nv0*p;
          cmo(i + ne*p, 3) = cn1(i, 2) + nv0*(p + 1);
          cmo(i + ne*p, 4) = cn1(i, 1) + nv0*(p + 1);
        }
      }
      else if (strcmp(eltType[0], "QUAD") == 0)
      {
        FldArrayI& cn0 = *cn[p];
        FldArrayI& cn1 = *cn[p+1];
        FldArrayI& cmo = *(cno->getConnect(0));

        for (E_Int i = 0; i < ne; i++)
        {
          cmo(i + ne*p, 1) = cn0(i, 1) + nv0*p;
          cmo(i + ne*p, 2) = cn0(i, 2) + nv0*p;
          cmo(i + ne*p, 3) = cn0(i, 3) + nv0*p;
          cmo(i + ne*p, 4) = cn0(i, 4) + nv0*p;
          cmo(i + ne*p, 5) = cn1(i, 1) + nv0*(p + 1);
          cmo(i + ne*p, 6) = cn1(i, 2) + nv0*(p + 1);
          cmo(i + ne*p, 7) = cn1(i, 3) + nv0*(p + 1);
          cmo(i + ne*p, 8) = cn1(i, 4) + nv0*(p + 1);
        }
      }
      else if (strcmp(eltType[0], "TRI") == 0)
      {
        FldArrayI& cn0 = *cn[p];
        FldArrayI& cn1 = *cn[p+1];
        FldArrayI& cmo = *(cno->getConnect(0));

        for (E_Int i = 0; i < ne; i++)
        {
          cmo(i + ne*p, 1) = cn0(i, 1) + nv0*p;
          cmo(i + ne*p, 2) = cn0(i, 2) + nv0*p;
          cmo(i + ne*p, 3) = cn0(i, 3) + nv0*p;
          cmo(i + ne*p, 4) = cn1(i, 1) + nv0*(p + 1);
          cmo(i + ne*p, 5) = cn1(i, 2) + nv0*(p + 1);
          cmo(i + ne*p, 6) = cn1(i, 3) + nv0*(p + 1);     
        }
      }
    }

    RELEASESHAREDU(tpl, coord, cno);
  }
    
  for (E_Int noz = 0; noz < ns; noz++)
    RELEASESHAREDS(objs[noz], structF[noz]);
  for (E_Int noz = 0; noz < nu; noz++)
    RELEASESHAREDU(obju[noz], unstructF[noz], cn[noz]);

  return tpl;
}
