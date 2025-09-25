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

// Build an isoSurf in a volume tetra mesh by marching tetra

# include "post.h"
using namespace std;
using namespace K_FLD;

//=============================================================================
/* Construit une isoSurf dans un maillage tetra */
//=============================================================================
PyObject* K_POST::isoSurf(PyObject* self, PyObject* args)
{
  // grid: maillage volumique tetra (x,y,z+sol)
  // field: nom du field dont on cherche l'iso
  // value: valeur de l'iso
  PyObject* grid;
  char* field; E_Float value;
  if (!PYPARSETUPLE_(args, O_ S_ R_, &grid, &field, &value))
  {
    return NULL;
  }

  /*----------------------------------------------*/
  /* Extraction des donnees du maillage volumique */ 
  /*----------------------------------------------*/
  char* varString0; char* eltType0;
  FldArrayF* f; FldArrayI* cn;
  E_Int nil, njl, nkl;
  E_Int res = 
    K_ARRAY::getFromArray3(grid, varString0, f, nil, njl, nkl, 
                           cn, eltType0);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "isoSurf: input array is invalid.");
    return NULL;
  }
  if (res != 2 || strcmp(eltType0, "TETRA") != 0)
  {
    RELEASESHAREDB(res, grid, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "isoSurf: input array must be TETRA.");
    return NULL;
  }
  
  // Check size of array
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString0);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString0);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString0);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "isoSurf: coordinates not found in array.");
    RELEASESHAREDU(grid, f, cn);
    return NULL;
  }
  posx++; posy++; posz++;
  
  // position de la variable iso
  E_Int posf = K_ARRAY::isNamePresent(field, varString0);
  if (posf == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "isoSurf: variable doesn't exist in array.");
    RELEASESHAREDU(grid, f, cn);
    return NULL;
  }
  posf++;

  // Position du cellN eventuellement
  E_Int poscellN = K_ARRAY::isCellNatureField2Present(varString0)+1;

  FldArrayF fiso; FldArrayI ciso;
  doIsoSurf(*f, *cn, posf, value, poscellN, fiso, ciso);
  RELEASESHAREDU(grid, f, cn);
  E_Float tolc = 1.e-12;
  K_CONNECT::cleanConnectivity(posx, posy, posz, tolc, 
                               "TRI", fiso, ciso);

  if (fiso.getSize() == 0 || ciso.getSize() == 0)
  {
    PyErr_SetString(PyExc_ValueError,
                    "isoSurf: isosurf is empty.");
    return NULL;
  }

  E_Int api = 1; //f->getApi();
  PyObject* t = K_ARRAY::buildArray3(fiso, varString0, ciso, "TRI", api);
  return t;
}

//=============================================================================
/* 
   IN: f: le champ a iso-surfacer.
   IN: la connectivite de f
   IN: posf: la position de la variable d'isosurface dans f
   IN: value: la valeur de l'isosurface
   IN: posCellN: la position de cellN dans f
   OUT: fiso: le champ de l'iso
   OUT: ciso: la connectivite de l'iso
*/
//==============================================================================
void K_POST::doIsoSurf(FldArrayF& f, FldArrayI& cn, E_Int posf, E_Float value,
                       E_Int poscellN,
                       FldArrayF& fiso, FldArrayI& ciso)
{
  E_Int nelts = cn.getSize();
  E_Int nfld = f.getNfld();
  E_Float* fp = f.begin(posf);
  E_Int* cn1 = cn.begin(1); E_Int* cn2 = cn.begin(2);
  E_Int* cn3 = cn.begin(3); E_Int* cn4 = cn.begin(4);

  E_Int nthreads = __NUMTHREADS__;

  // Dimensionnement: npts et ntri (par thread)
  E_Int* ntris = new E_Int [nthreads];
  E_Int* npts = new E_Int [nthreads];
#pragma omp parallel default(shared)
  {
  E_Int  ithread = __CURRENT_THREAD__;
  int triindex;
  E_Float f0, f1, f2, f3;
  E_Int ind0, ind1, ind2, ind3;
  E_Int np = 0; E_Int ntri = 0;

#pragma omp for
  for (E_Int i = 0; i < nelts; i++)
  {
    ind0 = cn1[i]-1; ind1 = cn2[i]-1; ind2 = cn3[i]-1; ind3 = cn4[i]-1;
    f0 = fp[ind0]; f1 = fp[ind1]; f2 = fp[ind2]; f3 = fp[ind3];
    triindex = 0;
    if (f0 < value) triindex |= 1;
    if (f1 < value) triindex |= 2;
    if (f2 < value) triindex |= 4;
    if (f3 < value) triindex |= 8;
    switch (triindex)
    {
      case 0x00:
      case 0x0F:
        break;

      case 0x0E: // OK
        np += 3; ntri++;
        break;

      case 0x01: // OK
        np += 3; ntri++;
        break;

      case 0x0D: // OK
        np += 3; ntri++;
      break;

      case 0x02: // OK
        np += 3; ntri++;
      break;

      case 0x0C: // OK
        np += 6; ntri += 2;
        break;

      case 0x03:
        np += 6; ntri += 2;
        break;

      case 0x0B: // OK
        np += 3; ntri++;
        break;

      case 0x04: // OK
        np += 3; ntri++;
        break;

      case 0x0A:
        np += 6; ntri += 2;
        break;

      case 0x05:
        np += 6; ntri += 2;
        break;

      case 0x09: // OK
        np += 6; ntri += 2;        
        break;

      case 0x06:
        np += 6; ntri += 2;        
        break;

      case 0x07:
        np += 3; ntri++;
        break;

      case 0x08: // OK
        np += 3; ntri++;
        break;
    }
  }
  npts[ithread] = np;
  ntris[ithread] = ntri;
  }

  FldArrayI** cisos = new FldArrayI* [nthreads];
  for (E_Int i = 0; i < nthreads; i++) cisos[i] = new FldArrayI(ntris[i], 3);
  FldArrayF** fisos = new FldArrayF* [nthreads];
  for (E_Int i = 0; i < nthreads; i++) fisos[i] = new FldArrayF(npts[i],nfld);
  E_Int* prevT = new E_Int [nthreads];
  E_Int* prevF = new E_Int [nthreads];
  
#pragma omp parallel default(shared)
  {
  E_Int  ithread = __CURRENT_THREAD__;
  E_Float f0, f1, f2, f3;
  E_Int ind0, ind1, ind2, ind3;
  int triindex;

  E_Int ntri = 0; // nombre de tri dans l'iso
  E_Int npt = 0; // nombre de pts dans l'iso
  
  FldArrayI& cisop = *cisos[ithread];
  E_Int* ciso1 = cisop.begin(1);
  E_Int* ciso2 = cisop.begin(2);
  E_Int* ciso3 = cisop.begin(3);
  FldArrayF& fisop = *fisos[ithread];

#pragma omp for
  for (E_Int i = 0; i < nelts; i++)
  {
    ind0 = cn1[i]-1; ind1 = cn2[i]-1; ind2 = cn3[i]-1; ind3 = cn4[i]-1;
    f0 = fp[ind0]; f1 = fp[ind1]; f2 = fp[ind2]; f3 = fp[ind3];

    triindex = 0;
    if (f0 < value) triindex |= 1;
    if (f1 < value) triindex |= 2;
    if (f2 < value) triindex |= 4;
    if (f3 < value) triindex |= 8;

    /* Form the vertices of the triangles for each case */
    switch (triindex)
    {
      case 0x00:
      case 0x0F:
        break;

      case 0x0E: // OK
        vertexInterp(nfld, value, f, poscellN,
                     f0, f1, ind0, ind1,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f0, f2, ind0, ind2,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f0, f3, ind0, ind3,
                     fisop, npt);
        ciso1[ntri] = npt;
        ciso2[ntri] = npt-1;
        ciso3[ntri] = npt-2;
        ntri++;
        break;

      case 0x01: // OK
        vertexInterp(nfld, value, f, poscellN,
                     f0, f1, ind0, ind1,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f0, f2, ind0, ind2,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f0, f3, ind0, ind3,
                     fisop, npt);
        ciso1[ntri] = npt-2;
        ciso2[ntri] = npt-1;
        ciso3[ntri] = npt;
        ntri++;
        break;

      case 0x0D: // OK
        vertexInterp(nfld, value, f, poscellN,
                     f1, f0, ind1, ind0,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f1, f3, ind1, ind3,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f1, f2, ind1, ind2,
                     fisop, npt);
        ciso1[ntri] = npt;
        ciso2[ntri] = npt-1;
        ciso3[ntri] = npt-2;
        ntri++;
      break;

      case 0x02: // OK
        vertexInterp(nfld, value, f, poscellN,
                     f1, f0, ind1, ind0,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f1, f3, ind1, ind3,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f1, f2, ind1, ind2,
                     fisop, npt);
        ciso1[ntri] = npt-2;
        ciso2[ntri] = npt-1;
        ciso3[ntri] = npt;
        ntri++;
      break;

      case 0x0C: // OK
        vertexInterp(nfld, value, f, poscellN,
                     f0, f3, ind0, ind3,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f0, f2, ind0, ind2,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f1, f3, ind1, ind3,
                     fisop, npt);
        ciso1[ntri] = npt-2;
        ciso2[ntri] = npt-1;
        ciso3[ntri] = npt;
        ntri++;
        
        vertexInterp(nfld, value, f, poscellN,
                     f1, f3, ind1, ind3,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f1, f2, ind1, ind2,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f0, f2, ind0, ind2,
                     fisop, npt);
        ciso1[ntri] = npt;
        ciso2[ntri] = npt-1;
        ciso3[ntri] = npt-2;
        ntri++;
        break;

      case 0x03:
        vertexInterp(nfld, value, f, poscellN,
                     f0, f3, ind0, ind3,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f0, f2, ind0, ind2,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f1, f3, ind1, ind3,
                     fisop, npt);
        ciso1[ntri] = npt-2;
        ciso2[ntri] = npt-1;
        ciso3[ntri] = npt;
        ntri++;
        
        vertexInterp(nfld, value, f, poscellN,
                     f1, f3, ind1, ind3,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f1, f2, ind1, ind2,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f0, f2, ind0, ind2,
                     fisop, npt);
        ciso1[ntri] = npt-2;
        ciso2[ntri] = npt-1;
        ciso3[ntri] = npt;
        ntri++;
        break;

      case 0x0B: // OK
        vertexInterp(nfld, value, f, poscellN,
                     f2, f0, ind2, ind0,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f2, f1, ind2, ind1,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f2, f3, ind2, ind3,
                     fisop, npt);
        ciso1[ntri] = npt;
        ciso2[ntri] = npt-1;
        ciso3[ntri] = npt-2;
        ntri++;
        break;

      case 0x04: // OK
        vertexInterp(nfld, value, f, poscellN,
                     f2, f0, ind2, ind0,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f2, f1, ind2, ind1,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f2, f3, ind2, ind3,
                     fisop, npt);
        ciso1[ntri] = npt-2;
        ciso2[ntri] = npt-1;
        ciso3[ntri] = npt;
        ntri++;
        break;

      case 0x0A:
        vertexInterp(nfld, value, f, poscellN,
                     f0, f1, ind0, ind1,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f2, f3, ind2, ind3,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f0, f3, ind0, ind3,
                     fisop, npt);
        ciso1[ntri] = npt;
        ciso2[ntri] = npt-1;
        ciso3[ntri] = npt-2;
        ntri++;

        vertexInterp(nfld, value, f, poscellN,
                     f0, f1, ind0, ind1,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f1, f2, ind1, ind2,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f2, f3, ind2, ind3,
                     fisop, npt);
        ciso1[ntri] = npt;
        ciso2[ntri] = npt-1;
        ciso3[ntri] = npt-2;
        ntri++;
        break;

      case 0x05:
        vertexInterp(nfld, value, f, poscellN,
                     f0, f1, ind0, ind1,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f2, f3, ind2, ind3,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f0, f3, ind0, ind3,
                     fisop, npt);
        ciso1[ntri] = npt-2;
        ciso2[ntri] = npt-1;
        ciso3[ntri] = npt;
        ntri++;

        vertexInterp(nfld, value, f, poscellN,
                     f0, f1, ind0, ind1,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f1, f2, ind1, ind2,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f2, f3, ind2, ind3,
                     fisop, npt);
        ciso1[ntri] = npt-2;
        ciso2[ntri] = npt-1;
        ciso3[ntri] = npt;
        ntri++;
        break;

      case 0x09: // OK
        vertexInterp(nfld, value, f, poscellN,
                     f0, f1, ind0, ind1,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f1, f3, ind1, ind3,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f2, f3, ind2, ind3,
                     fisop, npt);
        ciso1[ntri] = npt; // OK
        ciso2[ntri] = npt-1;
        ciso3[ntri] = npt-2;
        ntri++;

        vertexInterp(nfld, value, f, poscellN,
                     f0, f1, ind0, ind1,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f0, f2, ind0, ind2,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f2, f3, ind2, ind3,
                     fisop, npt);
        ciso1[ntri] = npt-2; // OK
        ciso2[ntri] = npt-1;
        ciso3[ntri] = npt;
        ntri++;
      break;

      case 0x06:
        vertexInterp(nfld, value, f, poscellN,
                     f0, f1, ind0, ind1,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f1, f3, ind1, ind3,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f2, f3, ind2, ind3,
                     fisop, npt);
        ciso1[ntri] = npt-2; // OK
        ciso2[ntri] = npt-1;
        ciso3[ntri] = npt;
        ntri++;

        vertexInterp(nfld, value, f, poscellN,
                     f0, f1, ind0, ind1,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f0, f2, ind0, ind2,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f2, f3, ind2, ind3,
                     fisop, npt);
        ciso1[ntri] = npt;
        ciso2[ntri] = npt-1;
        ciso3[ntri] = npt-2;
        ntri++;
      break;

      case 0x07:
        vertexInterp(nfld, value, f, poscellN,
                     f3, f0, ind3, ind0,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f3, f2, ind3, ind2,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f3, f1, ind3, ind1,
                     fisop, npt);
        ciso1[ntri] = npt-2;
        ciso2[ntri] = npt-1;
        ciso3[ntri] = npt;
        ntri++;
        break;

      case 0x08: // OK
        vertexInterp(nfld, value, f, poscellN,
                     f3, f0, ind3, ind0,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f3, f2, ind3, ind2,
                     fisop, npt);
        vertexInterp(nfld, value, f, poscellN,
                     f3, f1, ind3, ind1,
                     fisop, npt);
        ciso1[ntri] = npt-2;
        ciso2[ntri] = npt-1;
        ciso3[ntri] = npt;
        ntri++;
        break;
    }
  }
  }

  // rebuild
  E_Int ntri = 0; E_Int npt = 0;
  for (E_Int i = 0; i < nthreads; i++) 
  { prevT[i] = ntri; ntri += ntris[i]; 
    prevF[i] = npt; npt += npts[i]; }
  //printf("%d %d\n", npt, ntri);

  fiso.malloc(npt, nfld);
  ciso.malloc(ntri, 3);
  
#pragma omp parallel default(shared)
  {
  E_Int ithread = __CURRENT_THREAD__;
  E_Int nq = ntris[ithread];
  E_Int p = prevT[ithread];
  E_Int f = prevF[ithread];
  for (E_Int n = 1; n <= 3; n++)
  {
    E_Int* cisop = ciso.begin(n);
    E_Int* cisol = cisos[ithread]->begin(n);
    for (E_Int e = 0; e < nq; e++) cisop[e+p] = cisol[e]+f;
  }
  E_Int np = npts[ithread];
  for (E_Int n = 1; n <= nfld; n++)
  {
    E_Float* fisop = fiso.begin(n);
    E_Float* fisol = fisos[ithread]->begin(n);
    for (E_Int e = 0; e < np; e++) fisop[e+f] = fisol[e];
  }
  }
  delete [] prevT; delete [] prevF;
  delete [] npts; delete [] ntris;
  for (E_Int i = 0; i < nthreads; i++) delete fisos[i];
  for (E_Int i = 0; i < nthreads; i++) delete cisos[i];
  delete [] fisos; delete [] cisos;
}
