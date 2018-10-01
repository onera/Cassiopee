/*    
    Copyright 2013-2018 Onera.

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

// Build an isoSurf in a volume NGON mesh using marching tetra

# include "post.h"
using namespace std;
using namespace K_FLD;

PyObject* K_POST::isoSurfNGon(PyObject* self, PyObject* args)
{
  PyObject* array;
  char* fieldName; E_Float value;
  if (!PYPARSETUPLEF(args,
                    "Osd", "Osf",
                    &array, &fieldName, &value)) return NULL;

  /*----------------------------------------------*/
  /* Extraction des donnees du maillage volumique */ 
  /*----------------------------------------------*/
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int nil, njl, nkl;
  E_Int res = 
    K_ARRAY::getFromArray2(array, varString, f, nil, njl, nkl, 
                          cn, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "isoSurf: input array is invalid.");
    return NULL;
  }
  if (res != 2 || strcmp(eltType, "NGON") != 0)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "isoSurf: input array must be NGON.");
    return NULL;
  }
  
  // Check size of array
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "isoSurf: coordinates not found in array.");
    RELEASESHAREDU(array, f, cn);
    return NULL;
  }
  posx++; posy++; posz++;

  // position de la variable iso
  E_Int posf = K_ARRAY::isNamePresent(fieldName, varString);
  if (posf == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "isoSurf: variable doesn't exist in array.");
    RELEASESHAREDU(array, f, cn);
    return NULL;
  }
  posf++;

  E_Int poscellN = K_ARRAY::isCellNatureField2Present(varString)+1;

  FldArrayF fiso; FldArrayI ciso;
  doIsoSurfNGon(*f, *cn, posf, value, poscellN, fiso, ciso);
  RELEASESHAREDU(array, f, cn);
  return NULL;

  /*
  E_Float tolc = 1.e-12;
  K_CONNECT::cleanConnectivity(posx, posy, posz, tolc, 
                               "TRI", fiso, ciso);

  if (fiso.getSize() == 0 || ciso.getSize() == 0)
  {
    PyErr_SetString(PyExc_ValueError,
                    "isoSurf: isosurf is empty.");
    return NULL;
  }

  PyObject* t = K_ARRAY::buildArray2(fiso, varString, ciso, -1, "TRI");
  return t;
  */
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
void K_POST::doIsoSurfNGon(FldArrayF& f, FldArrayI& cn, E_Int posf, E_Float value,
                           E_Int poscellN,
                           FldArrayF& fiso, FldArrayI& ciso)
{
  E_Int nfld = f.getNfld();
  E_Float* fp = f.begin(posf);
  
  E_Int nthreads = __NUMTHREADS__;

  E_Int nelts = cn.getNElts();
  E_Int nfaces = cn.getNFaces();
  printf("nelts=%d, nfaces=%d\n", nelts, nfaces);
  printf("api=%d, ngon=%d\n", cn.getApi(), cn.isNGon());
  printf("nfld=%d\n", cn.getNfld());
  E_Int* ptrf = cn.getNGon();
  E_Int* ptre = cn.getNFace();

  // Tableau de position des faces dans la connectivite
  E_Int* indPG = cn.getIndPG();

  /*
  printf("indPG\n");
  for (E_Int i = 0; i < nfaces; i++) printf("%d ", indPG[i]);
  printf("\n");
  */

  // Tableau de position des elements dans la connectivite
  E_Int* indPH = cn.getIndPH();
  
  /*
  printf("indPH\n");
  for (E_Int i = 0; i < nelts; i++) printf("%d ", indPH[i]);
  printf("\n");
  */

  // Dimension du NGON
  FldArrayI dimElts;
  K_CONNECT::getDimElts(cn, indPG, indPH, dimElts);
  /*
  printf("dimElts\n");
  for (E_Int i = 0; i < dimElts.getSize(); i++) printf("%d ", dimElts[i]);
  printf("\n");
  */
  E_Int dim = dimElts[0];

  if (dim == 1)
  {
    for (E_Int elt = 0; elt < nelts; elt++)
    {
      E_Int pe = indPH[elt];
      E_Int* pte = ptre+pe;
      E_Int nbFaces = pte[0];
      /*
      for (E_Int fa = 0; fa < nbFaces; fa++)
      {
        numface = pte[fa+1]-1;
        pf = indPH[numface];
        ptf = ptrf[pf];
        npoints = ptf[pos]; pos++;
        ind = cn1[pos]-1;

        // connectivite du nouvel element BAR
        newcnp[fa][indelt] = ind+1;
        // champs du nouvel element BAR
        for (E_Int p=0;p<nfld;p++) 
        {
            fnewp[p][ind] = fp[p][ind]; // premier point de l'arete de l ancien element
          }
        }
        indelt++;
        cn2 += nbFaces+1;
        
      }
      */
    }

  }
  else if (dim == 2)
  {

  }
  else if (dim == 3)
  {
    E_Int pe, pf, indFace, nbFaces, nbPts, ind;
    E_Int* ptf; E_Int* pte;
    FldArrayF fco(nfld); E_Float* fc = fco.begin();
    FldArrayF ffo(nfld); E_Float* ff = ffo.begin();
    for (E_Int elt = 0; elt < nelts; elt++)
    {
      // Construit centre de l'elements + centres des faces
      pe = indPH[elt]; pte = ptre+pe;
      nbFaces = pte[0];
      for (E_Int p = 0; p < nfld; p++) fc[p] = 0.;
      for (E_Int fa = 0; fa < nbFaces; fa++)
      {
        indFace = pte[fa+1]-1;
        pf = indPG[indFace]; ptf = ptrf+pf;
        nbPts = ptf[0];
        for (E_Int p = 0; p < nfld; p++) ff[p] = 0.;
        for (E_Int pt = 0; pt < nbPts; pt++)
        {
          ind = ptf[1+pt];
          for (E_Int p = 0; p < nfld; p++) ff[p] += f(ind,p+1);
        }
        for (E_Int p = 0; p < nfld; p++) ff[p] = ff[p]/nbPts;

        for (E_Int p = 0; p < nfld; p++) fc[p] += ff[p];
      }
      for (E_Int p = 0; p < nfld; p++) fc[p] = fc[p]/nbFaces;

      for (E_Int fa = 0; fa < nbFaces; fa++)
      {
        indFace = pte[fa+1]-1;
        pf = indPG[indFace]; ptf = ptrf+pf;
        nbPts = ptf[0];
        for (E_Int p = 0; p < nfld; p++) ff[p] = 0.;
        for (E_Int pt = 0; pt < nbPts; pt++)
        {
          ind = ptf[1+pt];
          for (E_Int p = 0; p < nfld; p++) ff[p] += f(ind,p+1);
        }
        for (E_Int p = 0; p < nfld; p++) ff[p] = ff[p]/nbPts;

        for (E_Int pt = 0; pt < nbPts-1; pt++)
        {
          ind = ptf[1+pt]; ind2 = ptrf[2+pt];
          // tetra = centre, face, pts de la face
          f0 = fc[posf-1]; f1 = ff[posf-1]; f2 = f(ind,posf); f3 = f(ind,posf);
          
        } 
      }
    }
  }
  return;



  E_Int* cn1 = cn.begin(1); E_Int* cn2 = cn.begin(2);
  E_Int* cn3 = cn.begin(3); E_Int* cn4 = cn.begin(4);

  
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
