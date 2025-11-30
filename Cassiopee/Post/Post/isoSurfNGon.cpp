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

// Build an isoSurf in a volume NGON mesh using marching tetra

# include <unordered_map>
# include "post.h"
using namespace std;
using namespace K_FLD;

#define VERTEXINTERP(nfld, value, f, poscellN, f0, f1, ind0, ind1, fisos, npt) \
  { df = f1-f0;              \
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
  if (!PYPARSETUPLE_(args, O_ S_ R_,
                    &array, &fieldName, &value)) return NULL;

  /*----------------------------------------------*/
  /* Extraction des donnees du maillage volumique */ 
  /*----------------------------------------------*/
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int nil, njl, nkl;
  E_Int res = 
    K_ARRAY::getFromArray3(array, varString, f, nil, njl, nkl, 
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
  
  if (fiso.getSize() == 0 || ciso.getSize() == 0)
  {
    PyErr_SetString(PyExc_ValueError,
                    "isoSurf: isosurf is empty.");
    return NULL;
  }

  E_Int api = 1; //f->getApi();
  PyObject* t = K_ARRAY::buildArray3(fiso, varString, ciso, "TRI", api);
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
  E_Int npoints = f.getSize();
  E_Float* fp = f.begin(posf);
  
  E_Int nthreads = __NUMTHREADS__;

  E_Int nelts = cn.getNElts();
  //E_Int nfaces = cn.getNFaces();
  //printf("nelts=" SF_D_ ", nfaces=" SF_D_ "\n", nelts, nfaces);
  //printf("api=" SF_D_ ", ngon=" SF_D_ "\n", cn.getApi(), cn.getNGonType());
  //printf("nfld=" SF_D_ "\n", cn.getNfld());
  fflush(stdout);
  E_Int* ptrf = cn.getNGon();
  E_Int* ptre = cn.getNFace();

  // Tableau de position des faces dans la connectivite
  E_Int* indPG = cn.getIndPG();

  /*
  printf("indPG\n");
  for (E_Int i = 0; i < nfaces; i++) printf(SF_D_ " ", indPG[i]);
  printf("\n");
  */

  // Tableau de position des elements dans la connectivite
  E_Int* indPH = cn.getIndPH();
  
  /*
  printf("indPH\n");
  for (E_Int i = 0; i < nelts; i++) printf(SF_D_ " ", indPH[i]);
  printf("\n");
  */

  // Dimension du NGON
  //FldArrayI dimElts;
  //K_CONNECT::getDimElts(cn, dimElts);
  /*
  printf("dimElts\n");
  for (E_Int i = 0; i < dimElts.getSize(); i++) printf(SF_D_ " ", dimElts[i]);
  printf("\n");
  */
  //E_Int dim = dimElts[0];

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

    // Dimensionnement: npts et ntri (par thread*10)
    E_Int* ntris2 = new E_Int [nthreads*10];
    E_Int* npts2 = new E_Int [nthreads*10];
    
#pragma omp parallel default(shared)
    {
      E_Int ithread = __CURRENT_THREAD__;
      int triindex;
      E_Int np = 0; E_Int ntri = 0;
      E_Int pe, pf, indFace, nbFaces, nbPts, ind, ind2;
      E_Int* ptf; E_Int* pte;
      E_Float f0, f1, f2, f3;
      E_Float ffs, fcs;
      E_Float delta = (nelts*1.)/(nthreads*1.);
      E_Int ieltstart = int(ithread*delta);
      E_Int ieltend = int((ithread+1)*delta);
      //printf("ieltstart = " SF_D_ " , " SF_D_ " check=" SF_D_ "\n",ieltstart, ieltend, nelts);
      E_Float deltap = (ieltend-ieltstart)/(10.);
      E_Int elt;

      elt = ieltstart;
      //printf("borne start=" SF_D_ "\n", elt);
      for (E_Int j = 0; j < 10; j++)
      {
        np = 0; ntri = 0;
        //printf("" SF_D2_ "\n", int((j)*deltap), int((j+1)*deltap));
        for (E_Int k = 0; k < int((j+1)*deltap)-int(j*deltap); k++)
        {
          //printf("" SF_D2_ "\n", elt, nelts);
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
          }
          elt++;
        }
        npts2[j+ithread*10] = np;
        ntris2[j+ithread*10] = ntri;
      }
      //printf("borne end=" SF_D_ "\n", elt);
    }

    // Nbre de tris
    E_Float alpha = 0.;
    for (E_Int i = 0; i < 10*nthreads; i++) alpha += ntris2[i];
    alpha = alpha/nthreads;
    //printf("ntri=" SF_D_ ", ntri moyen par thread=" SF_D_ "\n", int(alpha*nthreads),int(alpha));
    //fflush(stdout);
    if (K_FUNC::fEqualZero(alpha))
    {
      delete [] npts2; delete [] ntris2;
      fiso.malloc(0,nfld);
      ciso.malloc(0,3);
      return;
    }

    //for (E_Int i = 0; i < nthreads; i++)
    //{
    //  for (E_Int j = 0; j < 10; j++) printf("thread=" SF_D_ ", j=" SF_D_ " -> " SF_D_ "," SF_D_ "\n", i,j,npts2[j+i*10],ntris2[j+i*10]);
    //}
    //printf("dimensionnement...\n");
    //fflush(stdout);

    // equilibrage dynamique
    FldArrayI iestart(10*nthreads);
    FldArrayI ieend(10*nthreads);
    E_Float delta = (nelts*1.)/(nthreads*1.);
    for (E_Int i = 0; i < nthreads; i++)
    {
      E_Int ieltstart = int(i*delta);
      E_Int ieltend = int((i+1)*delta);
      E_Float deltap = (ieltend-ieltstart)/(10.);
      for (E_Int j = 0; j < 10; j++)
      {
        iestart[j+10*i] = ieltstart+int(deltap*j);
        ieend[j+10*i] = ieltstart+int(deltap*(j+1));
      }
    }
    ieend[10*nthreads-1] = nelts;
    //for (E_Int i = 0; i < nthreads; i++)
    //  for (E_Int j = 0; j < 10; j++)
    //    printf("" SF_D2_ ": " SF_D2_ "\n",i,j,iestart[j+10*i],ieend[j+10*i]);
    //fflush(stdout);

    E_Int* istart = new E_Int [nthreads];
    E_Int* iend = new E_Int [nthreads];
    E_Int* ntris = new E_Int [nthreads];
    E_Int* npts = new E_Int [nthreads];
    
    istart[0] = 0; E_Int ibold = 0; E_Int ib = 0;
    E_Float plus = 0.;
    for (E_Int i = 0; i < nthreads; i++)
    {
      E_Int nc = 0; E_Int np = 0;
      while (ib < nthreads*10 && nc+plus*ntris2[ib] < int(alpha))
      {
        nc += ntris2[ib];
        np += npts2[ib]; ib++;
      }
      if (plus == 0.) plus = 1.;
      else plus = 0.;
      
      if (i == nthreads-1) // ajoute la fin (si necessaire)
      {
        while (ib < nthreads*10)
        {
          nc += ntris2[ib];
          np += npts2[ib]; ib++;
        }
      }
      ntris[i] = nc; npts[i] = np;
      //printf("DEBUG ib=" SF_D2_ "\n",ibold,ib);
      if (ib > ibold)
      {
        istart[i] = iestart[ibold];
        iend[i] = ieend[ib-1];
      }
      else if (ibold < 10*nthreads)
      {
        istart[i] = iestart[ibold];
        iend[i] = iestart[ibold];
      }
      else
      {
        istart[i] = ieend[10*nthreads-1];
        iend[i] = ieend[10*nthreads-1];
      }
      
      //printf("DEBUG istart=" SF_D2_ "\n",istart[i],iend[i]);
      ibold = ib;
    }
    //iend[nthreads-1] = nelts;
    //printf("reequilibrage: nthreads=" SF_D_ "\n", nthreads);
    //for (E_Int i = 0; i < nthreads; i++) printf("thread=" SF_D_ ": ntri=" SF_D_ " / start=" SF_D_ " end=" SF_D_ "\n", i, ntris[i], istart[i], iend[i]);
    //fflush(stdout);
    // fin equilibrage dynamique
    delete [] npts2; delete [] ntris2;

    FldArrayI** cisos = new FldArrayI* [nthreads];
    for (E_Int i = 0; i < nthreads; i++) cisos[i] = new FldArrayI(ntris[i], 3);
    FldArrayF** fisos = new FldArrayF* [nthreads];
    for (E_Int i = 0; i < nthreads; i++) fisos[i] = new FldArrayF(npts[i],nfld);
    E_Int* prevT = new E_Int [nthreads];
    E_Int* prevF = new E_Int [nthreads];
    
//#define KEY(p1,p2,e,f) p1+p2+2*npoints*e+(2*npoints*nelts)*f
#define KEY(p1,p2,e,f) (ETK)p1+(ETK)p2*npoints+npoints*npoints*(ETK)e+(npoints*npoints*nelts)*(ETK)f
#define ETK E_LONG
#define SETKEY12(p1,p2) \
  if (K_FUNC::fEqualZero(alpha) == true) \
    key[npt-1] = KEY(p1,-1,-1,-1); \
  else if (K_FUNC::fEqualZero(alpha1) == true) \
    key[npt-1] = KEY(p2,-1,-1,-1); \
  else if (p1 > p2) \
    key[npt-1] = KEY(p1,p2,-1,-1); \
  else key[npt-1] = KEY(p2,p1,-1,-1);
  
#define SETKEY13(p1,e) \
  if (K_FUNC::fEqualZero(alpha) == true) \
    key[npt-1] = KEY(p1,-1,-1,-1); \
  else if (K_FUNC::fEqualZero(alpha1) == true) \
    key[npt-1] = KEY(-1,-1,e,-1); \
  else \
    key[npt-1] = KEY(p1,-1,e,-1);
  
#define SETKEY31(e,p1) \
  if (K_FUNC::fEqualZero(alpha) == true) \
    key[npt-1] = KEY(-1,-1,e,-1); \
  else if (K_FUNC::fEqualZero(alpha1) == true) \
    key[npt-1] = KEY(p1,-1,-1,-1); \
  else \
    key[npt-1] = KEY(p1,-1,e,-1);
  
#define SETKEY14(p1,f) \
  if (K_FUNC::fEqualZero(alpha) == true) \
    key[npt-1] = KEY(p1,-1,-1,-1); \
  else if (K_FUNC::fEqualZero(alpha1) == true) \
    key[npt-1] = KEY(-1,-1,-1,f); \
  else \
    key[npt-1] = KEY(p1,-1,-1,f);

#define SETKEY41(f,p1) \
  if (K_FUNC::fEqualZero(alpha) == true) \
    key[npt-1] = KEY(-1,-1,-1,f); \
  else if (K_FUNC::fEqualZero(alpha1) == true) \
    key[npt-1] = KEY(p1,-1,-1,-1); \
  else \
    key[npt-1] = KEY(p1,-1,-1,f);
  
#define SETKEY34(e,f) \
  if (K_FUNC::fEqualZero(alpha) == true) \
    key[npt-1] = KEY(-1,-1,e,-1); \
  else if (K_FUNC::fEqualZero(alpha1) == true) \
    key[npt-1] = KEY(-1,-1,-1,f); \
  else \
    key[npt-1] = KEY(-1,-1,e,f);
  
#define SETKEY43(f,e) \
  if (K_FUNC::fEqualZero(alpha) == true) \
    key[npt-1] = KEY(-1,-1,-1,f); \
  else if (K_FUNC::fEqualZero(alpha1) == true) \
    key[npt-1] = KEY(-1,-1,e,-1); \
  else \
    key[npt-1] = KEY(-1,-1,e,f);

    FldArray<ETK>** keys = new FldArray<ETK>* [nthreads];
    for (E_Int i = 0; i < nthreads; i++) keys[i] = new FldArray<ETK>(npts[i]);

#pragma omp parallel default(shared)
    {
      E_Int ithread = __CURRENT_THREAD__;
      E_Int pe, pf, indFace, nbFaces, nbPts, ind, indp;
      E_Int ind0, ind1, ind2, ind3;
      E_Int* ptf; E_Int* pte;
      FldArrayF fco(nfld); E_Float* fc = fco.begin();
      FldArrayF ffo(nfld); E_Float* ff = ffo.begin();
      FldArrayF ffp(4,nfld);
      E_Float f0, f1, f2, f3, alpha, alpha1, val, df;
      int triindex;

      E_Int ntri = 0; // nombre de tri dans l'iso
      E_Int npt = 0; // nombre de pts dans l'iso

      FldArrayI& cisop = *cisos[ithread];
      E_Int* ciso1 = cisop.begin(1);
      E_Int* ciso2 = cisop.begin(2);
      E_Int* ciso3 = cisop.begin(3);
      FldArrayF& fisop = *fisos[ithread];

      ETK* key = keys[ithread]->begin();
      
      //printf("" SF_D_ ": " SF_D2_ "\n", ithread,istart[ithread],iend[ithread]); fflush(stdout);
      
      for (E_Int elt = istart[ithread]; elt < iend[ithread]; elt++)
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
            if (pt == nbPts-1) indp = ptf[1]-1;
            else indp = ptf[2+pt]-1;
            // tetra = centre, face, pts de la face
            f0 = fc[posf-1]; f1 = ff[posf-1]; f2 = fp[ind]; f3 = fp[indp];
            for (E_Int n = 0; n < nfld; n++) ffp(0,n+1) = fc[n];
            for (E_Int n = 0; n < nfld; n++) ffp(1,n+1) = ff[n];
            for (E_Int n = 0; n < nfld; n++) ffp(2,n+1) = f(ind,n+1);
            for (E_Int n = 0; n < nfld; n++) ffp(3,n+1) = f(indp,n+1);
            ind0 = 0; ind1 = 1; ind2 = 2; ind3 = 3;
            
            // tetra = ind, indp, elt, face
            /*
            E_Float tf0, tf1, tf2, tf3;
            E_Int tind0, tind1, tind2, tind3, telt, tindFace, tind, tindp;
            tf0 = f2; tf1 = f3; tf2 = f1; tf3 = f0;
            f0 = tf0; f1 = tf1; f2 = tf2; f3 = tf3;
            tind0 = ind2; tind1 = ind3; tind2 = ind1; tind3 = ind0;
            ind0 = tind0; ind1 = tind1; ind2 = tind2; ind3 = tind3;
            */
              
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
              SETKEY34(elt,indFace);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f2, ind0, ind2,
                           fisop, npt);
              SETKEY31(elt,ind);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f3, ind0, ind3,
                           fisop, npt);
              SETKEY31(elt,indp);
              ciso1[ntri] = npt;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt-2;
              ntri++;
              break;

              case 0x01: // OK
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f1, ind0, ind1,
                           fisop, npt);
              SETKEY34(elt,indFace);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f2, ind0, ind2,
                           fisop, npt);
              SETKEY31(elt,ind);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f3, ind0, ind3,
                           fisop, npt);
              SETKEY31(elt,indp);
              ciso1[ntri] = npt-2;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt;
              ntri++;
              break;

              case 0x0D: // OK
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f0, ind1, ind0,
                           fisop, npt);
              SETKEY34(elt,indFace);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f3, ind1, ind3,
                           fisop, npt);
              SETKEY41(indFace, indp);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f2, ind1, ind2,
                           fisop, npt);
              SETKEY41(indFace, ind);
              ciso1[ntri] = npt;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt-2;
              ntri++;
              break;

              case 0x02: // OK
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f0, ind1, ind0,
                           fisop, npt);
              SETKEY43(indFace, elt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f3, ind1, ind3,
                           fisop, npt);
              SETKEY41(indFace, indp);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f2, ind1, ind2,
                           fisop, npt);
              SETKEY41(indFace, ind);
              ciso1[ntri] = npt-2;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt;
              ntri++;
              break;

              case 0x0C: // OK
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f3, ind0, ind3,
                           fisop, npt);
              SETKEY31(elt, indp);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f2, ind0, ind2,
                           fisop, npt);
              SETKEY31(elt, ind);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f3, ind1, ind3,
                           fisop, npt);
              SETKEY41(indFace, indp);
              ciso1[ntri] = npt-2;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt;
              ntri++;
              
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f3, ind1, ind3,
                           fisop, npt);
              SETKEY41(indFace, indp);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f2, ind1, ind2,
                           fisop, npt);
              SETKEY41(indFace, ind);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f2, ind0, ind2,
                           fisop, npt);
              SETKEY31(elt, ind);
              ciso1[ntri] = npt;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt-2;
              ntri++;
              break;

              case 0x03:
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f3, ind0, ind3,
                           fisop, npt);
              SETKEY31(elt, indp);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f2, ind0, ind2,
                           fisop, npt);
              SETKEY31(elt, ind);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f3, ind1, ind3,
                           fisop, npt);
              SETKEY41(indFace, indp);
              ciso1[ntri] = npt-2;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt;
              ntri++;
              
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f3, ind1, ind3,
                           fisop, npt);
              SETKEY41(indFace, indp);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f2, ind1, ind2,
                           fisop, npt);
              SETKEY41(indFace, ind);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f2, ind0, ind2,
                           fisop, npt);
              SETKEY31(elt, ind);
              ciso1[ntri] = npt-2;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt;
              ntri++;
              break;

              case 0x0B: // OK
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f0, ind2, ind0,
                           fisop, npt);
              SETKEY13(ind, elt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f1, ind2, ind1,
                           fisop, npt);
              SETKEY14(ind, indFace);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f3, ind2, ind3,
                           fisop, npt);
              SETKEY12(ind, indp);
              ciso1[ntri] = npt;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt-2;
              ntri++;
              break;

              case 0x04: // OK
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f0, ind2, ind0,
                           fisop, npt);
              SETKEY13(ind, elt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f1, ind2, ind1,
                           fisop, npt);
              SETKEY14(ind, indFace);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f3, ind2, ind3,
                           fisop, npt);
              SETKEY12(ind, indp);
              ciso1[ntri] = npt-2;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt;
              ntri++;
              break;

              case 0x0A:
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f1, ind0, ind1,
                           fisop, npt);
              SETKEY34(elt, indFace);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f3, ind2, ind3,
                           fisop, npt);
              SETKEY12(ind, indp);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f3, ind0, ind3,
                           fisop, npt);
              SETKEY31(elt, indp);
              ciso1[ntri] = npt;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt-2;
              ntri++;
              
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f1, ind0, ind1,
                           fisop, npt);
              SETKEY34(elt, indFace);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f2, ind1, ind2,
                           fisop, npt);
              SETKEY41(indFace, ind);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f3, ind2, ind3,
                           fisop, npt);
              SETKEY12(ind, indp);
              ciso1[ntri] = npt;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt-2;
              ntri++;
              break;

              case 0x05:
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f1, ind0, ind1,
                           fisop, npt);
              SETKEY34(elt, indFace);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f3, ind2, ind3,
                           fisop, npt);
              SETKEY12(ind, indp);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f3, ind0, ind3,
                           fisop, npt);
              SETKEY31(elt, indp);
              ciso1[ntri] = npt-2;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt;
              ntri++;
              
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f1, ind0, ind1,
                           fisop, npt);
              SETKEY34(elt, indFace);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f2, ind1, ind2,
                           fisop, npt);
              SETKEY41(indFace, ind);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f3, ind2, ind3,
                           fisop, npt);
              SETKEY12(ind, indp);
              ciso1[ntri] = npt-2;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt;
              ntri++;
              break;

              case 0x09: // OK
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f1, ind0, ind1,
                           fisop, npt);
              SETKEY34(elt, indFace);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f3, ind1, ind3,
                           fisop, npt);
              SETKEY41(indFace, indp);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f3, ind2, ind3,
                           fisop, npt);
              SETKEY12(ind, indp);
              ciso1[ntri] = npt; // OK
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt-2;
              ntri++;
              
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f1, ind0, ind1,
                           fisop, npt);
              SETKEY34(elt, indFace);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f2, ind0, ind2,
                           fisop, npt);
              SETKEY31(elt, ind);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f3, ind2, ind3,
                           fisop, npt);
              SETKEY12(ind, indp);
              ciso1[ntri] = npt-2; // OK
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt;
              ntri++;
              break;

              case 0x06:
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f1, ind0, ind1,
                           fisop, npt);
              SETKEY34(elt, indFace);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f1, f3, ind1, ind3,
                           fisop, npt);
              SETKEY41(indFace, indp);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f3, ind2, ind3,
                           fisop, npt);
              SETKEY12(ind, indp);
              ciso1[ntri] = npt-2; // OK
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt;
              ntri++;
              
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f1, ind0, ind1,
                           fisop, npt);
              SETKEY34(elt, indFace);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f0, f2, ind0, ind2,
                           fisop, npt);
              SETKEY31(elt, ind);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f2, f3, ind2, ind3,
                           fisop, npt);
              SETKEY12(ind, indp);
              ciso1[ntri] = npt;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt-2;
              ntri++;
              break;

              case 0x07:
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f3, f0, ind3, ind0,
                           fisop, npt);
              SETKEY13(indp, elt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f3, f2, ind3, ind2,
                           fisop, npt);
              SETKEY12(indp, ind);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f3, f1, ind3, ind1,
                           fisop, npt);
              SETKEY14(indp, indFace);
              ciso1[ntri] = npt-2;
              ciso2[ntri] = npt-1;
              ciso3[ntri] = npt;
              ntri++;
              break;

              case 0x08: // OK
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f3, f0, ind3, ind0,
                           fisop, npt);
              SETKEY13(indp, elt);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f3, f2, ind3, ind2,
                           fisop, npt);
              SETKEY12(indp, ind);
              VERTEXINTERP(nfld, value, ffp, poscellN,
                           f3, f1, ind3, ind1,
                           fisop, npt);
              SETKEY14(indp, indFace);
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

  //printf("Analyse\n");
  // Analyse
  /*
  for (E_Int ithread = 0; ithread < nthreads; ithread++)
  {
    E_Float* x = fisos[ithread]->begin(1);
    E_Float* y = fisos[ithread]->begin(2);
    E_Float* z = fisos[ithread]->begin(3);
    ETK* key = keys[ithread]->begin();
    ETK* key2 = keys[ithread]->begin(2);
    for (E_Int i = 0; i < npts[ithread]; i++)
      printf("" SF_D_ ": %f %f %f (key=" SF_D_ ", source=" SF_D_ ")\n",i,x[i],y[i],z[i],key[i],key2[i]);
    fflush(stdout);
  }
  */
  
  // Nbre de pts dup + nbre de tri
  E_Int ntri = 0; E_Int npt = 0;
  for (E_Int i = 0; i < nthreads; i++) 
  { prevT[i] = ntri; ntri += ntris[i];
    prevF[i] = npt; npt += npts[i]; }
  //printf("nbre de pts dup=" SF_D_ ", nbre de tris dup=" SF_D_ "\n",npt,ntri);
  //fflush(stdout);
  
  // DEBUG: force une key unique
  /*
  for (E_Int i = 0; i < nthreads; i++)
  {
    E_Int f = prevF[i];
    ETK* key = keys[i]->begin();
    for (E_Int j = 0; j < npts[i]; j++)
      key[j] = (ETK)(j+f);
  }
  */
  
  // Construction de la map (cher) key->indDup
  //printf("construction de la map\n");
  std::unordered_map<ETK, E_Int> map;
  for (E_Int i = 0; i < nthreads; i++)
  {
    E_Int f = prevF[i];
    ETK* key = keys[i]->begin();
    for (E_Int j = 0; j < keys[i]->getSize(); j++)
      map[key[j]] = j+f;
  }
  
  // invMap: ind dup -> ind 
  FldArrayI invMap(npt);
  E_Int* invMapp = invMap.begin();
    
  // Nouveau nombre de points (non dup)
  npt = map.size();
  //printf("nbre de pts uniques=" SF_D_ "\n",npt); fflush(stdout);
  E_Int c = 0;
  for (std::pair<ETK,E_Int> elt : map)
  {
    //ETK k = elt.first;
    E_Int ind = elt.second;
    //printf("map c=" SF_D_ " inddup=" SF_D_ " key=" SF_D_ "\n",c,ind,k);
    invMapp[ind] = c;
    c++;
  }
  //fflush(stdout);
  //printf("invmap0\n");
  //for (E_Int i = 0; i < invMap.getSize(); i++) printf("invdup=" SF_D_ ": ind=" SF_D_ "\n",i,invMap[i]);
  fflush(stdout);
  
  // complete invMap
#pragma omp parallel default(shared)
  {
    E_Int ithread = __CURRENT_THREAD__;
    E_Int f = prevF[ithread];
    ETK* key = keys[ithread]->begin();
    for (E_Int i = 0; i < npts[ithread]; i++)
    { 
      ETK k = key[i];
      //printf("check f=" SF_D_ " key=" SF_D_ " inddup=" SF_D_ " [" SF_D_ "]\n",f+i,k,map[k],invMap[map[k]]);
      //if (f+i != map[k]) 
      invMapp[f+i] = invMapp[map[k]];
    }
  }
  
  // free the map
  map.clear();
  
  //printf("invmap\n");
  //for (E_Int i = 0; i < invMap.getSize(); i++) printf("invdup=" SF_D_ ": ind=" SF_D_ "\n",i,invMap[i]);
  //fflush(stdout);
  
  //printf("reconstruction fiso (" SF_D_ " points)\n", npt); fflush(stdout);
  fiso.malloc(npt, nfld);
  
#pragma omp parallel default(shared)
  {
    E_Int ithread = __CURRENT_THREAD__;
  
    E_Int f = prevF[ithread];
    E_Int np = npts[ithread];
    //printf("" SF_D2_ "\n", np, f); fflush(stdout);
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fisop = fiso.begin(n);
      E_Float* fisol = fisos[ithread]->begin(n);
      for (E_Int e = 0; e < np; e++) fisop[invMapp[e+f]] = fisol[e];
    }
  }
  //for (E_Int i = 0; i < npt; i++) printf("f " SF_D_ ": %f %f %f\n",i,fiso(i,1),fiso(i,2),fiso(i,3));
  //fflush(stdout);
  
  //printf("reconstruction ciso\n"); fflush(stdout);
  //ciso.malloc(ntri, 3);
  FldArrayI ciso2(ntri,3);
  
#pragma omp parallel default(shared)
  {
    E_Int ithread = __CURRENT_THREAD__;

    E_Int f = prevF[ithread];
    E_Int p = prevT[ithread];
    E_Int ne = ntris[ithread];
    for (E_Int n = 1; n <= 3; n++)
    {
      E_Int* cisop = ciso2.begin(n);
      E_Int* cisol = cisos[ithread]->begin(n);
      for (E_Int e = 0; e < ne; e++) cisop[e+p] = invMapp[cisol[e]+f-1]+1;
    }
  } 
  
  //for (E_Int i = 0; i < ntri; i++) printf("c " SF_D_ ": " SF_D3_ "\n",i,ciso2(i,1),ciso2(i,2),ciso2(i,3));
  //fflush(stdout);
  //printf("done\n"); fflush(stdout);
  
  // delete
  for (E_Int i = 0; i < nthreads; i++) delete keys[i];
  delete [] keys;

  delete [] prevT; delete [] prevF;
  delete [] npts; delete [] ntris;
  delete [] istart; delete [] iend;
  for (E_Int i = 0; i < nthreads; i++) delete fisos[i];
  for (E_Int i = 0; i < nthreads; i++) delete cisos[i];
  delete [] fisos; delete [] cisos;
  invMap.malloc(0);
  
  // Elimination des elements identiques (eventuels)
  //printf("Elimination des elements identiques\n");
  
#define KEY3S(c1,c2,c3) (ETK)(c1-1)+npt*(ETK)(c2-1)+npt*npt*(ETK)(c3-1)
//#define KEY3(c1,c2,c3) KEY3S(c1,c2,c3) + KEY3S(c1,c3,c2) + KEY3S(c2,c1,c3) + KEY3S(c2,c3,c1) + KEY3S(c3,c1,c2) + KEY3S(c3,c2,c1)
//#define KEY3(c1,c2,c3) KEY3S(c1,c2,c3)*KEY3S(c1,c3,c2)*KEY3S(c2,c1,c3)*KEY3S(c2,c3,c1)*KEY3S(c3,c1,c2)*KEY3S(c3,c2,c1)
//#define KEY3(c1,c2,c3) (ETK)(c1-1)+npt*(ETK)(c2-1)+npt*npt*(ETK)(c3-1)
#define KEY3(c1,c2,c3) \
  if (c1 == c2 || c1 == c3 || c2 == c3) k = -1; \
  else if (c1 < c2 && c2 < c3) k=KEY3S(c1,c2,c3); \
    else if (c1 < c3 && c3 < c2) k=KEY3S(c1,c3,c2); \
      else if (c2 < c1 && c1 < c3) k=KEY3S(c2,c1,c3); \
        else if (c2 < c3 && c3 < c1) k=KEY3S(c2,c3,c1); \
          else if (c3 < c1 && c1 < c2) k=KEY3S(c3,c1,c2); \
            else k=KEY3S(c3,c2,c1);
      
  FldArray<ETK> key(ntri); ETK* keyp = key.begin(); 
  E_Int* ct1 = ciso2.begin(1); E_Int* ct2 = ciso2.begin(2); E_Int* ct3 = ciso2.begin(3);
  
#pragma omp parallel default(shared)
  {
    ETK k;
#pragma omp for
    for (E_Int i = 0; i < ntri; i++) 
    {
      KEY3(ct1[i],ct2[i],ct3[i]);
      keyp[i] = k;
      //keyp[i] = i; // DEBUG
    }
  }  
  //for (E_Int i = 0; i < ntri; i++) printf("ct=" SF_D3_ " = " SF_D_ "\n",ct1[i],ct2[i],ct3[i],keyp[i]);
  //fflush(stdout);

  for (E_Int i = 0; i < ntri; i++) { if (keyp[i] != -1) map[keyp[i]] = i; }
    
  E_Int newTri = map.size();
  c = 0;
  ciso.malloc(newTri,3);
  E_Int* cp1 = ciso.begin(1); E_Int* cp2 = ciso.begin(2); E_Int* cp3 = ciso.begin(3);
  E_Int ind;
  //printf("elimination TRI multiples ou degeneres =" SF_D_ "\n", newTri); fflush(stdout);
  for (std::pair<ETK,E_Int> elt : map)
  {
    //ETK k = elt.first;
    ind = elt.second;
    //printf("map c=" SF_D_ " indTri=" SF_D_ " key=" SF_D_ "\n",c,ind,k);
    cp1[c] = ct1[ind];
    cp2[c] = ct2[ind];
    cp3[c] = ct3[ind];
    c++;
  }
  
}
