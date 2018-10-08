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

#define VERTEXINTERP(nfld, value, f, poscellN, f0, f1, ind0, ind1, fisos, npt) \
  { E_Float alpha, alpha1, val; E_Float df = f1-f0;              \
    if (K_FUNC::fEqualZero(df) == true) alpha = 1.;              \
    else alpha = (value-f0)/df;                                  \
    alpha1 = 1.-alpha;                                           \
    for (E_Int j = 1; j <= nfld; j++)                            \
     { val = alpha1*f(ind0, j)+alpha*f(ind1, j);             \
      fisos(npt, j) = val; }                                    \
    if (poscellN != 0)                                         \
    { if (f(ind0, poscellN) == 0. || f(ind1, poscellN) == 0.)  \
     fisos(npt, poscellN) = 0.; }                    \
     npt++; }

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
  
  //E_Float tolc = 1.e-12;
  //K_CONNECT::cleanConnectivity(posx, posy, posz, tolc, 
  //                             "TRI", fiso, ciso);

  if (fiso.getSize() == 0 || ciso.getSize() == 0)
  {
    PyErr_SetString(PyExc_ValueError,
                    "isoSurf: isosurf is empty.");
    return NULL;
  }

  PyObject* t = K_ARRAY::buildArray(fiso, varString, ciso, -1, "TRI");
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
#define COMPUTEFF(fa) \
    indFace = pte[fa+1]-1; \
    pf = indPG[indFace]; ptf = ptrf+pf; \
    nbPts = ptf[0]; \
    for (E_Int p = 0; p < nfld; p++) ff[p] = 0.; \
      for (E_Int pt = 0; pt < nbPts; pt++) \
      { ind = ptf[1+pt]-1; \
        for (E_Int p = 0; p < nfld; p++) ff[p] += f(ind,p+1); } \
    for (E_Int p = 0; p < nfld; p++) ff[p] = ff[p]/nbPts;

#define COMPUTEFFS(fa) \
    indFace = pte[fa+1]-1; \
    pf = indPG[indFace]; ptf = ptrf+pf; \
    nbPts = ptf[0]; \
    ffs = 0.; \
    for (E_Int pt = 0; pt < nbPts; pt++) \
    { ind = ptf[1+pt]-1; \
      ffs += fp[ind]; } \
    ffs = ffs/nbPts;

    // Dimensionnement: npts et ntri (par thread)
    E_Int* ntris = new E_Int [nthreads];
    E_Int* npts = new E_Int [nthreads];

#pragma omp parallel default(shared)
    {
      E_Int  ithread = __CURRENT_THREAD__;
      int triindex;
      E_Int np = 0; E_Int ntri = 0;

      E_Int pe, pf, indFace, nbFaces, nbPts, ind, ind2;
      E_Int* ptf; E_Int* pte;
      E_Float f0, f1, f2, f3;
      FldArrayF fco(nfld); E_Float* fc = fco.begin();
      FldArrayF ffo(nfld); E_Float* ff = ffo.begin();
      E_Float ffs, fcs;

#pragma omp for
      for (E_Int elt = 0; elt < nelts; elt++)
      {
        // Construit centre de l'element
        pe = indPH[elt]; pte = ptre+pe;
        nbFaces = pte[0];
        fcs = 0.;
        for (E_Int fa = 0; fa < nbFaces; fa++)
        {
          COMPUTEFFS(fa);
          fcs += ffs;
        }
        fcs = fcs/nbFaces;

        // construit les tetras
        for (E_Int fa = 0; fa < nbFaces; fa++)
        {
          COMPUTEFFS(fa);
          for (E_Int pt = 0; pt < nbPts; pt++)
          {
            ind = ptf[1+pt]-1;
            if (pt == nbPts-1) ind2 = ptf[1]-1;
            else ind2 = ptf[2+pt]-1;
            // tetra = centre, face, pts de la face
            f0 = fcs; f1 = ffs; f2 = fp[ind]; f3 = fp[ind2];
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
      }
    }
    printf("dimensionnement...\n");


    // equilibrage dynamique
    /*
    E_Int* istart = new E_Int [nthreads];
    E_Int* iend = new E_Int [nthreads];
    E_Int* ntris = new E_Int [nthreads];
    E_Int* npts = new E_Int [nthreads];
    E_Float alpha = 0.;
    for (E_Int i = 0; i < 10*nthreads; i++) alpha += ntris2[i];
    alpha = alpha/nthreads;

    istart[0] = 0; ib = 0;
    for (E_Int i = 1; i < nthreads; i++)
    {
      E_Int nc = 0; E_Int np = 0;
      while (nc < alpha)
      {
        nc += ntris[ib];
        np += npts[ib]; ib++;
      }
      ntris[i] = nc; npts[i] = np;
      istart[i] = istart[i-1]+nc;
      iend[i-1] = istart[i];
    }
    iend[nthreads-1] = nelts;
    printf("nthreads=%d\n", nthreads);
    for (E_Int i = 0; i < nthreads; i++) printf("%d %d %d\n", ntris2[i], istart[i],iend[i]);
    */
    // fin equilibrage dynamique

    FldArrayI** cisos = new FldArrayI* [nthreads];
    for (E_Int i = 0; i < nthreads; i++) cisos[i] = new FldArrayI(ntris[i], 3);
    FldArrayF** fisos = new FldArrayF* [nthreads];
    for (E_Int i = 0; i < nthreads; i++) fisos[i] = new FldArrayF(npts[i],nfld);
    E_Int* prevT = new E_Int [nthreads];
    E_Int* prevF = new E_Int [nthreads];
    //for (E_Int i = 0; i < nthreads; i++) printf("%d %d\n", ntris[i], npts[i]);

#pragma omp parallel default(shared)
    {
      E_Int ithread = __CURRENT_THREAD__;
      E_Int pe, pf, indFace, nbFaces, nbPts, ind, ind2;
      E_Int ind0, ind1, ind3;
      E_Int* ptf; E_Int* pte;
      FldArrayF fco(nfld); E_Float* fc = fco.begin();
      FldArrayF ffo(nfld); E_Float* ff = ffo.begin();
      FldArrayF ffp(4,nfld);
      E_Float f0, f1, f2, f3;
      int triindex;

      E_Int ntri = 0; // nombre de tri dans l'iso
      E_Int npt = 0; // nombre de pts dans l'iso

      FldArrayI& cisop = *cisos[ithread];
      E_Int* ciso1 = cisop.begin(1);
      E_Int* ciso2 = cisop.begin(2);
      E_Int* ciso3 = cisop.begin(3);
      FldArrayF& fisop = *fisos[ithread];

#pragma omp for
      for (E_Int elt = 0; elt < nelts; elt++)
      {
        // Construit centre de l'element
        pe = indPH[elt]; pte = ptre+pe;
        nbFaces = pte[0];
        for (E_Int p = 0; p < nfld; p++) fc[p] = 0.;
        for (E_Int fa = 0; fa < nbFaces; fa++)
        {
          COMPUTEFF(fa);
          for (E_Int p = 0; p < nfld; p++) fc[p] += ff[p];
        }
        for (E_Int p = 0; p < nfld; p++) fc[p] = fc[p]/nbFaces;

        // construit les tetras
        for (E_Int fa = 0; fa < nbFaces; fa++)
        {
          COMPUTEFF(fa);
          for (E_Int pt = 0; pt < nbPts; pt++)
          {
            ind = ptf[1+pt]-1; 
            if (pt == nbPts-1) ind2 = ptf[1]-1;
            else ind2 = ptf[2+pt]-1;
            // tetra = centre, face, pts de la face
            f0 = fc[posf-1]; f1 = ff[posf-1]; f2 = fp[ind]; f3 = fp[ind2];
            for (E_Int n = 0; n < nfld; n++) ffp(0,n+1) = fc[n];
            for (E_Int n = 0; n < nfld; n++) ffp(1,n+1) = ff[n];
            for (E_Int n = 0; n < nfld; n++) ffp(2,n+1) = f(ind,n+1);
            for (E_Int n = 0; n < nfld; n++) ffp(3,n+1) = f(ind2,n+1);
            ind0 = 0; ind1 = 1; ind2 = 2; ind3 = 3;
              
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
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f1, ind0, ind1,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f2, ind0, ind2,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f3, ind0, ind3,
                           fisop, npt);
              ciso1[ntri] = npt;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt-2;
              ntri++;
              break;

              case 0x01: // OK
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f1, ind0, ind1,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f2, ind0, ind2,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f3, ind0, ind3,
                           fisop, npt);
              ciso1[ntri] = npt-2;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt;
              ntri++;
              break;

              case 0x0D: // OK
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f0, ind1, ind0,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f3, ind1, ind3,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f2, ind1, ind2,
                           fisop, npt);
              ciso1[ntri] = npt;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt-2;
              ntri++;
              break;

              case 0x02: // OK
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f0, ind1, ind0,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f3, ind1, ind3,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f2, ind1, ind2,
                           fisop, npt);
              ciso1[ntri] = npt-2;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt;
              ntri++;
              break;

              case 0x0C: // OK
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f3, ind0, ind3,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f2, ind0, ind2,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f3, ind1, ind3,
                           fisop, npt);
              ciso1[ntri] = npt-2;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt;
              ntri++;

              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f3, ind1, ind3,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f2, ind1, ind2,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f2, ind0, ind2,
                           fisop, npt);
              ciso1[ntri] = npt;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt-2;
              ntri++;
              break;

              case 0x03:
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f3, ind0, ind3,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f2, ind0, ind2,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f3, ind1, ind3,
                           fisop, npt);
              ciso1[ntri] = npt-2;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt;
              ntri++;

              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f3, ind1, ind3,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f2, ind1, ind2,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f2, ind0, ind2,
                           fisop, npt);
              ciso1[ntri] = npt-2;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt;
              ntri++;
              break;

              case 0x0B: // OK
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f0, ind2, ind0,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f1, ind2, ind1,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f3, ind2, ind3,
                           fisop, npt);
              ciso1[ntri] = npt;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt-2;
              ntri++;
              break;

              case 0x04: // OK
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f0, ind2, ind0,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f1, ind2, ind1,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f3, ind2, ind3,
                           fisop, npt);
              ciso1[ntri] = npt-2;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt;
              ntri++;
              break;

              case 0x0A:
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f1, ind0, ind1,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f3, ind2, ind3,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f3, ind0, ind3,
                           fisop, npt);
              ciso1[ntri] = npt;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt-2;
              ntri++;

              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f1, ind0, ind1,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f2, ind1, ind2,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f3, ind2, ind3,
                           fisop, npt);
              ciso1[ntri] = npt;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt-2;
              ntri++;
              break;

              case 0x05:
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f1, ind0, ind1,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f3, ind2, ind3,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f3, ind0, ind3,
                           fisop, npt);
              ciso1[ntri] = npt-2;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt;
              ntri++;

              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f1, ind0, ind1,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f2, ind1, ind2,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f3, ind2, ind3,
                           fisop, npt);
              ciso1[ntri] = npt-2;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt;
              ntri++;
              break;

              case 0x09: // OK
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f1, ind0, ind1,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f3, ind1, ind3,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f3, ind2, ind3,
                           fisop, npt);
              ciso1[ntri] = npt; // OK
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt-2;
              ntri++;

              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f1, ind0, ind1,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f2, ind0, ind2,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f3, ind2, ind3,
                           fisop, npt);
              ciso1[ntri] = npt-2; // OK
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt;
              ntri++;
              break;

              case 0x06:
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f1, ind0, ind1,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f3, ind1, ind3,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f3, ind2, ind3,
                           fisop, npt);
              ciso1[ntri] = npt-2; // OK
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt;
              ntri++;

              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f1, ind0, ind1,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f2, ind0, ind2,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f3, ind2, ind3,
                           fisop, npt);
              ciso1[ntri] = npt;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt-2;
              ntri++;
              break;

              case 0x07:
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f3, f0, ind3, ind0,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f3, f2, ind3, ind2,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f3, f1, ind3, ind1,
                           fisop, npt);
              ciso1[ntri] = npt-2;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt;
              ntri++;
              break;

              case 0x08: // OK
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f3, f0, ind3, ind0,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f3, f2, ind3, ind2,
                           fisop, npt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
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

  } // dim=3

}
