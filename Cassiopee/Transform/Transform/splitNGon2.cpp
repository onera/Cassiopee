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
# include "Connect/connect.h"
# include "Metis/metis.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
// Split NGon using METIS
// IN: array NGON, arrayc in center containing "part"
// IN: nparts: nbre de parties pour le decoupage de premier niveau
// IN: nparts2: nbre de parties pour le decoupage de niveau 2
// OUT: modified center field with NPart indicator
//==============================================================================
PyObject* K_TRANSFORM::splitNGon2(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* arrayc;
  E_Int nparts; E_Int nparts2; E_Int shift;
  if (!PYPARSETUPLE_(args, OO_ III_,
                    &array, &arrayc, &nparts, &nparts2, &shift))
  {
    return NULL;
  }
  // Check array nodes
  E_Int ni, nj, nk, res;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray3(array, varString,
                               f, ni, nj, nk, cn, eltType);
  if (res != 2)
  {
    if (res == 1) RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError,
                    "splitNGon: array is invalid.");
    return NULL;
  }
  if (strcmp(eltType, "TRI")   == 0 || strcmp(eltType, "QUAD") == 0 ||
      strcmp(eltType, "TETRA") == 0 || strcmp(eltType, "HEXA") == 0 ||
      strcmp(eltType, "PENTA") == 0 || strcmp(eltType, "BAR")  == 0 ||
      strcmp(eltType, "PYRA")  == 0 || strcmp(eltType, "NODE") == 0)
  { RELEASESHAREDU(array, f, cn); return array; }
  if (strcmp(eltType, "NGON") != 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "splitNGon: elt type must be NGON.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }

  // Check array centers
  E_Int nic, njc, nkc;
  FldArrayF* fc; FldArrayI* cnc;
  char* varStringc; char* eltTypec;
  res = K_ARRAY::getFromArray3(arrayc, varStringc,
                               fc, nic, njc, nkc, cnc, eltTypec);

  E_Int posc = K_ARRAY::isNamePresent("part", varStringc); posc++;
  E_Float* fp = fc->begin(posc);

  // Construit le graph
  FldArrayI cFE;
  K_CONNECT::connectNG2FE(*cn, cFE);
  E_Int* cFE1 = cFE.begin(1);
  E_Int* cFE2 = cFE.begin(2);

  E_Int* cnp = cn->begin();
  E_Int* cne = cnp+2+cnp[1];
  E_Int nelts = cne[0];
  //printf("nelts=%d\n", nelts);
  cne += 2;

  E_Int nf;
  E_Int size = 0; // size of adj
  for (E_Int i = 0; i < nelts; i++)
  {
    nf = cne[0];
    size += nf;
    cne += nf+1;
  }
  //printf("size = %d\n", size);

  E_Int e1, e2, indf;
  idx_t* adj1 = new idx_t [size];
  idx_t* adj = adj1;
  idx_t* xadj = new idx_t [nelts+1];
  cne = cnp+4+cnp[1];
  size = 0;
  for (E_Int i = 0; i < nelts; i++)
  {
    xadj[i] = size;
    nf = cne[0];

    for (E_Int n = 0; n < nf; n++)
    {
      indf = cne[n+1]-1;
      e1 = cFE1[indf];
      e2 = cFE2[indf];
      //printf("%d - %d %d\n",i+1, e1,e2);
      if (e1 > 0 && e1 != i+1) { adj[size] = e1-1; size++; }
      else if (e2 > 0 && e2 != i+1) { adj[size] = e2-1; size++; }
    }
    cne += nf+1;
  }
  xadj[nelts] = size;
  adj = adj1;
  //for (E_Int i = 0; i < nelts+1; i++) printf("%d ",xadj[i]); printf("\n\n");
  //for (E_Int i = 0; i < size; i++) printf("%d ",adj[i]);

  //cFE.malloc(0);

  E_Int ncon = 1;
  E_Int objval = 0;
  idx_t* parts = new idx_t [nelts];
  if (nparts == 1) { for (E_Int i = 0; i < nelts; i++) parts[i] = 0; }
  else
  METIS_PartGraphKway(&nelts, &ncon, xadj, adj, NULL, NULL, NULL,
                      &nparts, NULL, NULL, NULL, &objval, parts);
  //METIS_PartGraphRecursive(&nelts, &ncon, xadj, adj, NULL, NULL, NULL,
  //                         &nparts, NULL, NULL, NULL, &objval, parts);

  delete [] xadj; delete [] adj1;

  // output: sortie de la liste des elements pour chaque part
  E_Int p;
  E_Int* partSize = new E_Int [nparts]; // nbre d'elements dans chaque partition

  for (E_Int i = 0; i < nparts; i++) partSize[i] = 0;
  for (E_Int i = 0; i < nelts; i++)
  {
    p = parts[i]; partSize[p] += 1;
  }
  //for (E_Int i = 0; i < nparts; i++) printf("Info: partSize=%d\n", partSize[i]);

  // output on field
  for (E_Int i = 0; i < nelts; i++)
  {
    fp[i] = shift*parts[i];
  }

  // Deuxieme niveau de decoupage
  if (nparts2 > 0)
  {
    for (E_Int np = 0; np < nparts; np++) // pour chaque partie de niveau 1
    {
      // nbre d'elements de cette partie
      E_Int nelts2 = partSize[np];
      cne = cnp+4+cnp[1];
      size = 0; // size of adj
      E_Int c = 0;
      E_Int* indir = new E_Int [nelts]; // indirection
      for (E_Int i = 0; i < nelts; i++)
      {
        nf = cne[0];
        if (parts[i] == np)
        { size += nf; indir[i] = c; c +=1; }
        cne += nf+1;
      }
      //printf("%d %d\n", nelts2, size);

      adj1 = new idx_t [size];
      adj = adj1;
      xadj = new idx_t [nelts2+1];
      cne = cnp+4+cnp[1];
      size = 0; c = 0;
      for (E_Int i = 0; i < nelts; i++)
      {
        xadj[c] = size;
        if (parts[i] == np) c += 1;
        nf = cne[0];

        for (E_Int n = 0; n < nf; n++)
        {
          indf = cne[n+1]-1;
          e1 = cFE1[indf];
          e2 = cFE2[indf];
          if (e1 > 0 && e1 != i+1 && parts[i] == np && parts[e1-1] == np) { adj[size] = indir[e1-1]; size++; }
          else if (e2 > 0 && e2 != i+1 && parts[i] == np && parts[e2-1] == np) { adj[size] = indir[e2-1]; size++; }
        }
        cne += nf+1;
      }
      xadj[nelts2] = size;
      adj = adj1;

      ncon = 1;
      objval = 0;
      idx_t* parts2 = new idx_t [nelts2];
      if (nparts2 == 1) { for (E_Int i = 0; i < nelts2; i++) parts2[i] = 0; }
      else
      METIS_PartGraphKway(&nelts2, &ncon, xadj, adj, NULL, NULL, NULL,
                          &nparts2, NULL, NULL, NULL, &objval, parts2);
      delete [] xadj; delete [] adj1;

      c = 0;
      for (E_Int i = 0; i < nelts; i++)
      {
        if (parts[i] == np)
        {
          fp[i] += parts2[c]; c += 1;
        }
      }
      delete [] parts2;

    } // Fin parts de niveau 1
  }
  delete [] parts;
  delete [] partSize;

  RELEASESHAREDU(array, f, cn);
  Py_INCREF(Py_None);
  return Py_None;
}
