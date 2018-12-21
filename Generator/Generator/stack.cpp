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

// stack the planes of two meshes (with same nixnj) into a single mesh

# include "generator.h"
using namespace K_FLD;
using namespace std;

//=============================================================================
// Conactenate mesh layers
//=============================================================================
PyObject* K_GENERATOR::stackMesh(PyObject* self, PyObject* args)
{
  PyObject* arrays;
  if (!PyArg_ParseTuple(args, "O", &arrays)) return NULL;

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
    E_Int ni1 = ni[0]; E_Int nj1 = nj[0];
    E_Int ni1nj1 = ni1*nj1;
    E_Int nk1 = 0;
    for (E_Int p = 0; p < ns; p++) nk1 += nk[p];
    tpl = K_ARRAY::buildArray(nfld, structVarString[0], 
                              ni1, nj1, nk1);
    E_Float* coordp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF coord(ni1nj1*nk1, nfld, coordp, true);

    for (E_Int p = 0; p < ns; p++)
    {
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* f1p = structF[p]->begin(n);
        E_Float* coordp = coord.begin(n);
        for (E_Int k = 0; k < nk[p]; k++)
          for (E_Int i = 0; i < ni1nj1; i++)
          {
            coordp[i+(k+p)*ni1nj1] = f1p[i+k*ni1nj1];
          }
      }
    }
  }
  else if (nu >= 2) // concatenate all unstructs
  {
    /* A finir
    char outEltType[10]; E_Int nf = 1;
    if (strcmp(eltType[0], "BAR") == 0) { strcpy(outEltType, "QUAD"); nf=4; }
    if (strcmp(eltType[0], "QUAD") == 0) { strcpy(outEltType, "HEXA"); nf=8; }
    if (strcmp(eltType[0], "TRI") == 0) { strcpy(outEltType, "PENTA"); nf=6; }
  
    E_Int nfld = unstructF[0]->getNfld();
    E_Int nv = unstructF[0]->getSize(); // same nvertex for all
    E_Int ne = cn[0]->getSize();
    E_Int nv1 = 0;
    E_Int ne1 = 0;
    for (E_Int p = 0; p < nu; p++) nv1 += unstructF[p]->getSize();
    ne1 = (nu-1)*cn[0]->getSize();

    tpl = K_ARRAY::buildArray(nfld, unstructVarString[0], nv1, ne1, -1, 
                              outEltType);
    E_Float* coordp = K_ARRAY::getFieldPtr(tpl);
    E_Int* cnp = K_ARRAY::getConnectPtr(tpl);
    FldArrayF coord(nv1, nfld, coordp, true);
    //FldArrayI cni();
    */
  }
  
  for (E_Int noz = 0; noz < ns; noz++)
    RELEASESHAREDS(objs[noz], structF[noz]);
  for (E_Int noz = 0; noz < nu; noz++)
    RELEASESHAREDU(obju[noz], unstructF[noz], cn[noz]);

  return tpl;
}
