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

# include "transform.h"

//=============================================================================
PyObject* K_TRANSFORM::breakMixedElements(
  K_FLD::FldArrayF& field, K_FLD::FldArrayI& ce, char* varString
)
{
  E_Int npts = field.getSize();
  E_Int nfld = field.getNfld();
  E_Int api = field.getApi();
  E_Int* cnp = ce.begin();
  E_Int size = ce.getSize();

  // Number of elements for each connectivity of the output list of BEs
  // NB1: 'tmp_' is uncompressed: all possible connectivities listed
  const E_Int nbuckets = 7;
  std::vector<E_Int> tmp_nepc2(nbuckets, 0);

  // List all possible BE element types. 'tmp_' as in uncompressed: some 
  // element types will not be found in the input array
  std::vector<const char*> tmp_eltType2;
  tmp_eltType2.push_back("BAR");
  tmp_eltType2.push_back("TRI"); tmp_eltType2.push_back("QUAD");
  tmp_eltType2.push_back("TETRA"); tmp_eltType2.push_back("PYRA");
  tmp_eltType2.push_back("PENTA"); tmp_eltType2.push_back("HEXA");

  // Map associating a MIXED element type index to the number of vertices per
  // element of this type, nvpe
  E_Int nvpe, nvpeloc, icBE;
  std::vector<E_Int> mixedType2nvpe(18, -1);
  mixedType2nvpe[3] = 2; // BAR
  mixedType2nvpe[5] = 3; // TRI
  mixedType2nvpe[7] = 4; // QUAD
  mixedType2nvpe[10] = 4; // TETRA
  mixedType2nvpe[12] = 5; // PYRA
  mixedType2nvpe[14] = 6; // PENTA
  mixedType2nvpe[17] = 8; // HEXA

  // Map associating a MIXED element type index to its corresponding BE bucket
  // index (as given by tmp_eltType2)
  std::vector<E_Int> mixed2ConnIc(18, -1); 
  mixed2ConnIc[3] = 0; // BAR
  mixed2ConnIc[5] = 1; // TRI
  mixed2ConnIc[7] = 2; // QUAD
  mixed2ConnIc[10] = 3; // TETRA
  mixed2ConnIc[12] = 4; // PYRA
  mixed2ConnIc[14] = 5; // PENTA
  mixed2ConnIc[17] = 6; // HEXA

  // First pass to get the number of elements of each type
  E_Int ps = 0; E_Int ntype;
  while (ps < size)
  {
    ntype = cnp[0];
    if (ntype < 0 || ntype >= (E_Int)mixedType2nvpe.size())
    {
      PyErr_SetString(PyExc_ValueError,
                      "breakMixedElements: unknown type of element.");
      return NULL;
    }

    nvpeloc = mixedType2nvpe[ntype];
    icBE = mixed2ConnIc[ntype];
    if (nvpeloc == -1)
    {
      PyErr_SetString(PyExc_ValueError,
                      "breakMixedElements: unknown type of element.");
      return NULL;
    }
    else
    {
      tmp_nepc2[icBE]++;
      ps += nvpeloc+1; cnp += nvpeloc+1;
    }
  }

  // Build new uncompressed eltType from connectivities that have at least
  // one element
  E_Int nc2 = 0;
  for (E_Int ic = 0; ic < nbuckets; ic++)
    if (tmp_nepc2[ic] > 0) nc2++;

  // Build new list of BE connectivities
  PyObject* l = PyList_New(0);
  PyObject* tpl;
  FldArrayF* f2; FldArrayI* cn2;
  
  for (E_Int ic = 0; ic < nbuckets; ic++)  // from BAR (0) to HEXA (6)
  {
    if (tmp_nepc2[ic] == 0) continue;  // no elements in this conn., skip
    
    tpl = K_ARRAY::buildArray3(nfld, varString, npts,
                               tmp_nepc2[ic], tmp_eltType2[ic], false, api);

    K_ARRAY::getFromArray3(tpl, f2, cn2);
    nvpe = cn2->getNfld();

    #pragma omp parallel
    {
      // Copy fields
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fp = field.begin(n);
        E_Float* f2p = f2->begin(n);
        #pragma omp for nowait
        for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
      }
    }

    // Second pass: copy connectivity
    tmp_nepc2[ic] = 0;
    ps = 0; cnp = ce.begin();

    while (ps < size)
    {
      ntype = cnp[0];
      if (ic == mixed2ConnIc[ntype])  // mixed element corresponding to this bucket
      {
        for (E_Int j = 1; j < nvpe+1; j++) (*cn2)(tmp_nepc2[ic], j) = cnp[j];
        tmp_nepc2[ic]++;
        cnp += nvpe+1; ps += nvpe+1;
      }
      else  // skip any other mixed element type
      {
        nvpeloc = mixedType2nvpe[ntype];
        cnp += nvpeloc+1; ps += nvpeloc+1;
      }
    }

    // Remove orphans when more than one type of element is found
    if (nc2 > 1)
    {
      E_Float eps = 1.e-12;
      PyObject* tplClean = K_CONNECT::V_cleanConnectivity(
        varString, *f2, *cn2, tmp_eltType2[ic], eps,
        false, true, false, false,
        false, false
      );

      RELEASESHAREDU(tpl, f2, cn2); Py_DECREF(tpl);
      PyList_Append(l, tplClean); Py_DECREF(tplClean);
    }
    else
    {
      RELEASESHAREDU(tpl, f2, cn2);
      PyList_Append(l, tpl); Py_DECREF(tpl);
    }
  }

  return l;
}
