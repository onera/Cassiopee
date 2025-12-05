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

# include "transform.h"
#include <map>

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

  // Number of elements per connectivity of the output ME
  // NB1: 'tmp_' is uncompressed: all possible connectivities listed
  const E_Int nbuckets = 7;
  std::vector<E_Int> tmp_nepc2(nbuckets, 0);

  // First pass to get the number of elements of each type
  E_Int ps = 0; E_Int ntype;
  while (ps < size)
  {
    ntype = cnp[0];
    if (ntype == 3) // BAR
    {
      tmp_nepc2[0]++;
      ps += 3; cnp += 3;
    }
    else if (ntype == 5) // TRI
    {
      tmp_nepc2[1]++;
      ps += 4; cnp += 4;
    }
    else if (ntype == 7) // QUAD
    {
      tmp_nepc2[2]++;
      ps += 5; cnp += 5;
    }
    else if (ntype == 10) // TETRA
    {
      tmp_nepc2[3]++;
      ps += 5; cnp += 5;
    }
    else if (ntype == 12) // PYRA
    {
      tmp_nepc2[4]++;
      ps += 6; cnp += 6;
    }
    else if (ntype == 14) // PENTA
    {
      tmp_nepc2[5]++;
      ps += 7; cnp += 7;
    }
    else if (ntype == 17) // HEXA
    {
      tmp_nepc2[6]++;
      ps += 9; cnp += 9;
    }
    else
    {
      printf("Warning: breakMixedElements: unknow type of element.\n");
    }
  }

  // Build new eltType from connectivities that have at least one element
  E_Int nc2 = 0;
  std::map<E_Int, E_Int> old2newIc;
  char* eltType2 = new char[K_ARRAY::VARSTRINGLENGTH];
  eltType2[0] = '\0';
  if (tmp_nepc2[0] > 0)
  {
    old2newIc.insert({0, nc2});
    strcat(eltType2, "BAR");
    nc2++;
  }
  if (tmp_nepc2[1] > 0)
  {
    old2newIc.insert({1, nc2});
    if (nc2 > 0) strcat(eltType2, ",");
    strcat(eltType2, "TRI");
    nc2++;
  }
  if (tmp_nepc2[2] > 0)
  {
    old2newIc.insert({2, nc2});
    if (nc2 > 0) strcat(eltType2, ",");
    strcat(eltType2, "QUAD");
    nc2++;
  }
  if (tmp_nepc2[3] > 0)
  {
    old2newIc.insert({3, nc2});
    if (nc2 > 0) strcat(eltType2, ",");
    strcat(eltType2, "TETRA");
    nc2++;
  }
  if (tmp_nepc2[4] > 0)
  {
    old2newIc.insert({4, nc2});
    if (nc2 > 0) strcat(eltType2, ",");
    strcat(eltType2, "PYRA");
    nc2++;
  }
  if (tmp_nepc2[5] > 0)
  {
    old2newIc.insert({5, nc2});
    if (nc2 > 0) strcat(eltType2, ",");
    strcat(eltType2, "PENTA");
    nc2++;
  }
  if (tmp_nepc2[6] > 0)
  {
    old2newIc.insert({6, nc2});
    if (nc2 > 0) strcat(eltType2, ",");
    strcat(eltType2, "HEXA");
    nc2++;
  }

  // Compress the number of elements per connectivity of the output ME, ie,
  // drop connectivities containing no elements
  std::vector<E_Int> nepc2(nc2);
  nc2 = 0;
  for (E_Int ic = 0; ic < nbuckets; ic++)  // from BAR (0) to HEXA (6)
  {
    if (tmp_nepc2[ic] > 0) { nepc2[nc2] = tmp_nepc2[ic]; nc2++; }
  }

  // Build new ME connectivity
  PyObject* l = PyList_New(0);
  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts,
                                       nepc2, eltType2, false, api);
  FldArrayF* f2; FldArrayI* cn2;
  K_ARRAY::getFromArray3(tpl, f2, cn2);

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

  for (E_Int ic = 0; ic < nbuckets; ic++) tmp_nepc2[ic] = 0;

  std::vector<FldArrayI*> cms2(nc2);
  for (E_Int ic = 0; ic < nc2; ic++) cms2[ic] = cn2->getConnect(ic);

  // Copy connectivities
  E_Int ic2;
  while (ps < size)
  {
    ntype = cnp[0];

    if (ntype == 3) // BAR
    {
      ic2 = old2newIc[0];
      for (E_Int j = 1; j <= 2; j++)
      {
        (*cms2[ic2])(tmp_nepc2[0], j) = cnp[j];
      }
      tmp_nepc2[0]++;
      cnp += 2+1; ps += 2+1;
    }
    else if (ntype == 5) // TRI
    {
      ic2 = old2newIc[1];
      for (E_Int j = 1; j <= 3; j++)
      {
        (*cms2[ic2])(tmp_nepc2[1], j) = cnp[j];
      }
      tmp_nepc2[1]++;
      cnp += 3+1; ps += 3+1;
    }
    else if (ntype == 7) // QUAD
    {
      ic2 = old2newIc[2];
      for (E_Int j = 1; j <= 4; j++)
      {
        (*cms2[ic2])(tmp_nepc2[2], j) = cnp[j];
      }
      tmp_nepc2[2]++;
      cnp += 4+1; ps += 4+1;
    }
    else if (ntype == 10) // TETRA
    {
      ic2 = old2newIc[3];
      for (E_Int j = 1; j <= 4; j++)
      {
        (*cms2[ic2])(tmp_nepc2[3], j) = cnp[j];
      }
      tmp_nepc2[3]++;
      cnp += 4+1; ps += 4+1;
    }
    else if (ntype == 12) // PYRA
    {
      ic2 = old2newIc[4];
      for (E_Int j = 1; j <= 5; j++)
      {
        (*cms2[ic2])(tmp_nepc2[4], j) = cnp[j];
      }
      tmp_nepc2[4]++;
      cnp += 5+1; ps += 5+1;
    }
    else if (ntype == 14) // PENTA
    {
      ic2 = old2newIc[5];
      for (E_Int j = 1; j <= 6; j++)
      {
        (*cms2[ic2])(tmp_nepc2[5], j) = cnp[j];
      }
      tmp_nepc2[5]++;
      cnp += 6+1; ps += 6+1;
    }
    else if (ntype == 17) // HEXA
    {
      ic2 = old2newIc[6];
      for (E_Int j = 1; j <= 8; j++)
      {
        (*cms2[ic2])(tmp_nepc2[6], j) = cnp[j];
      }
      tmp_nepc2[6]++;
      cnp += 8+1; ps += 8+1;
    }
  }

  RELEASESHAREDU(tpl, f2, cn2);
  delete [] eltType2;
  PyList_Append(l, tpl); Py_DECREF(tpl);
  return l;
}
