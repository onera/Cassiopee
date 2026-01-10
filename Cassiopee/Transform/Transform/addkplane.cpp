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

using namespace K_FLD;

// ============================================================================
/* Ajoute un plan k a un maillage */
// ============================================================================
PyObject* K_TRANSFORM::addkplane(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km,
                                     cn, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "addkplane: unknown type of array.");
    return NULL;
  }

  E_Int api = f->getApi();
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  posx++; posy++; posz++;

  // Vector add
  E_Float vx = 0.; E_Float vy = 0.; E_Float vz = 1.;

  if (res == 1)
  {
    E_Int imjm, imjmkm, km1;
    imjm = im*jm;
    imjmkm = imjm*km;
    km1 = km+1;
    E_Int nfld = f->getNfld();
    PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, im, jm, km1, api);
    FldArrayF* f2;
    K_ARRAY::getFromArray3(tpl, f2);

    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* f2p = f2->begin(n);
      E_Float* fpn = f->begin(n);
      #pragma omp parallel
      {
        E_Int ind, ind2;
        #pragma omp for nowait
        for (E_Int i = 0; i < imjmkm; i++)
        {
          f2p[i] = fpn[i];
        }
      
        #pragma omp for
        for (E_Int i = 0; i < im*jm; i++)
        {
          ind = i + (km-1)*imjm;
          ind2 = i + imjmkm;
          f2p[ind2] = fpn[ind];
        }
      }
    }

    if (posx > 0 && posy > 0 && posz > 0)
    {
      E_Float* fx = f->begin(posx); E_Float* f2x = f2->begin(posx);
      E_Float* fy = f->begin(posy); E_Float* f2y = f2->begin(posy);
      E_Float* fz = f->begin(posz); E_Float* f2z = f2->begin(posz);
      
      #pragma omp parallel
      {
        E_Int ind, ind1, ind2;
        #pragma omp for collapse(2)
        for (E_Int j = 0; j < jm; j++)
        for (E_Int i = 0; i < im; i++)
        {
          ind = i + j*im;
          ind1 = ind + (km-1)*imjm;
          ind2 = ind + imjmkm;
          f2x[ind2] = fx[ind1] + vx;
          f2y[ind2] = fy[ind1] + vy;
          f2z[ind2] = fz[ind1] + vz;
        }
      }
    }

    RELEASESHAREDS(tpl, f2);
    RELEASESHAREDS(array, f);
    return tpl;
  }
  else if (res == 2)
  {
    E_Bool isValidNGon = (
      K_STRING::cmp(eltType, 4, "NGON") == 0 && 
      (cn->getDim() == 1 || cn->getDim() == 2)
    );
    E_Bool isValidME = (
      K_STRING::cmp(eltType, 4, "NGON") != 0 &&
      (K_CONNECT::getDimME(eltType) == 1 || K_CONNECT::getDimME(eltType) == 2)
    );
                         
    if (!(isValidNGon || isValidME))
    {
      RELEASESHAREDU(array, f, cn);
      PyErr_SetString(PyExc_TypeError,
                      "addkplane: only for 1D/2D NGON or BE/ME arrays.");
      return NULL;
    }
    
    PyObject* tpl;
    FldArrayF* f2; FldArrayI* cn2;
    E_Int nfld = f->getNfld();

    if (isValidME)  // ME
    { 
      E_Int npts = f->getSize(); E_Int npts2 = 2*npts;
      E_Int nc = cn->getNConnect();
      std::vector<E_Int> nepc(nc);

      std::vector<char*> eltTypes;
      K_ARRAY::extractVars(eltType, eltTypes);
      char* eltType2 = new char[K_ARRAY::VARSTRINGLENGTH];
      eltType2[0] = '\0';
      
      for (E_Int ic = 0; ic < nc; ic++)
      {
        FldArrayI& cm = *(cn->getConnect(ic));
        nepc[ic] = cm.getSize();
        if (eltType2[0] != '\0') strcat(eltType2, ",");
        if (K_STRING::cmp(eltTypes[ic], "BAR") == 0) strcat(eltType2, "QUAD");
        else if (K_STRING::cmp(eltTypes[ic], "TRI") == 0) strcat(eltType2, "PENTA");
        else if (K_STRING::cmp(eltTypes[ic], "QUAD") == 0) strcat(eltType2, "HEXA");
      }

      for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
      
      tpl = K_ARRAY::buildArray3(nfld, varString, npts2, nepc, eltType2, false, api);
      K_ARRAY::getFromArray3(tpl, f2, cn2);

      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fp = f->begin(n);
        E_Float* f2p = f2->begin(n);
        for (E_Int i = 0; i < npts; i++)
        {
          f2p[i] = fp[i];
          f2p[i + npts] = fp[i];
        }
      }

      if (posx > 0 && posy > 0 && posz > 0)
      {
        E_Float* fx = f->begin(posx); E_Float* f2x = f2->begin(posx);
        E_Float* fy = f->begin(posy); E_Float* f2y = f2->begin(posy);
        E_Float* fz = f->begin(posz); E_Float* f2z = f2->begin(posz);
        for (E_Int i = 0; i < npts; i++)
        {
          f2x[i + npts] = fx[i] + vx;
          f2y[i + npts] = fy[i] + vy;
          f2z[i + npts] = fz[i] + vz;
        }
      }

      for (E_Int ic = 0; ic < nc; ic++)
      {
        FldArrayI& cm = *(cn->getConnect(ic));
        FldArrayI& cm2 = *(cn2->getConnect(ic));
        E_Int nvpe2 = cm2.getNfld();
      
        if (nvpe2 == 4) // to QUAD
        {
          for (E_Int i = 0; i < nepc[ic]; i++)
          {
            cm2(i, 1) = cm(i, 1);
            cm2(i, 2) = cm(i, 2);
            cm2(i, 3) = cm(i, 2) + npts;
            cm2(i, 4) = cm(i, 1) + npts;
          }
        }
        else if (nvpe2 == 6) // to PENTA
        {
          for (E_Int i = 0; i < nepc[ic]; i++)
          {
            cm2(i, 1) = cm(i, 1);
            cm2(i, 2) = cm(i, 2);
            cm2(i, 3) = cm(i, 3);
            cm2(i, 4) = cm(i, 1) + npts;
            cm2(i, 5) = cm(i, 2) + npts;
            cm2(i, 6) = cm(i, 3) + npts;
          }
        }
        else if (nvpe2 == 8) // to HEXA
        {
          for (E_Int i = 0; i < nepc[ic]; i++)
          {
            cm2(i, 1) = cm(i, 1);
            cm2(i, 2) = cm(i, 2);
            cm2(i, 3) = cm(i, 3);
            cm2(i, 4) = cm(i, 4);
            cm2(i, 5) = cm(i, 1) + npts;
            cm2(i, 6) = cm(i, 2) + npts;
            cm2(i, 7) = cm(i, 3) + npts;
            cm2(i, 8) = cm(i, 4) + npts;
          }
        }
      }

      delete[] eltType2;
      RELEASESHAREDU(tpl, f2, cn2);
    }
    else  // NGONs
    {
      E_Int* ngon = cn->getNGon(); E_Int* nface = cn->getNFace();
      E_Int *indPG = cn->getIndPG(), *indPH = cn->getIndPH();
      
      E_Int npts = f->getSize(); E_Int npts2 = 2*npts;
      E_Int nfaces = cn->getNFaces(); E_Int nelts = cn->getNElts();
      E_Int sizeFN = cn->getSizeNGon();
      E_Int sizeEF = cn->getSizeNFace();
      E_Int ngonType = cn->getNGonType();
      E_Int shift = 1; if (ngonType == 3) shift = 0;
     
      // E_Int sizeEF2 = sizeEF + nelts*2;// (nfaces+2) faces dans le volume
      // //E_Int nelts2 = nelts; // nb d elts dans le NGON volumique
      // E_Int nfaces2 = nfaces + 2*nelts;//nb de faces ds le NGON volumique
      // E_Int sumFS = 0;// dimensionnement du tableau faces/noeuds pour les faces correspondant aux elts surfaciques
      // E_Int nf;
      // for (E_Int i = 0; i < nelts; i++)
      // {
      //   cn->getElt(i, nf, nface, indPH);
      //   sumFS += nfloc + shift; // pour chq face vol: nfacesloc vertex + shift pour dimensionner
      // }
      // E_Int sizeFN2 = nfaces*(4 + shift) + 2*sumFS; // (nb de sommets + 1)
      
      // tpl = K_ARRAY::buildArray3(
      //   nfld, varString, npts2, nelts, nfaces2,
      //   "NGON", sizeFN2, sizeEF2, ngonType, false, api
      // );
      // K_ARRAY::getFromArray3(tpl, f2, cn2);

      // E_Int* ngon2 = cn2->getNGon();
      // E_Int* nface2 = cn2->getNFace();
      // E_Int *indPG2 = NULL, *indPH2 = NULL; 
      // if (api == 2 || api == 3) // array2 ou array3
      // {
      //   indPG2 = cn2->getIndPG(); indPH2 = cn2->getIndPH();
      // }

      // // duplication des champs, avec z = z+1 
      // for (E_Int n = 1; n <= nfld; n++)
      // {
      //   E_Float* fp = f->begin(n); E_Float* f2p = f2->begin(n);
      //   for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
      //   for (E_Int i = 0; i < npts; i++) f2p[i+npts] = fp[i];
      // }

      // if (posx > 0 && posy > 0 && posz > 0)
      // {
      //   E_Float* fx = f->begin(posx); E_Float* f2x = f2->begin(posx);
      //   E_Float* fy = f->begin(posy); E_Float* f2y = f2->begin(posy);
      //   E_Float* fz = f->begin(posz); E_Float* f2z = f2->begin(posz);

      //   for (E_Int i = 0; i < npts; i++)
      //   {
      //     f2x[i + npts] = fx[i] + vx;
      //     f2y[i + npts] = fy[i] + vy;
      //     f2z[i + npts] = fz[i] + vz;
      //   }
      // }

      // //=======================================================================
      // // connectivites
      // //=======================================================================
      // // a partir de chq face construction des faces laterales "quad" 
      // // extrudee a partir des faces surfaciques
      // E_Int c1, c2;
      // for (E_Int i = 0; i < nfaces; i++)
      // {
      //   c1 = i*(2 + shift); c2 = i*(4 + shift);
      //   ngon2[c2] = 4;
      //   ngon2[c2 + shift] = ngon[c1 + shift];
      //   ngon2[c2 + shift + 1] = ngon[c1 + shift + 1];
      //   ngon2[c2 + shift + 2] = ngon[c1 + shift + 1] + npts;
      //   ngon2[c2 + shift + 3] = ngon[c1 + shift] + npts;
      // }

      // // a partir des elts: recup des faces laterales: meme numerotation 
      // // qu'en surfacique
      // for (E_Int i = 0; i < nelts; i++)
      // {
      //   E_Int nfacessloc = ptrEFs[0];
      //   ptrEFv[0] = nfacessloc+2;
      //   for (E_Int nof = 1; nof <= nfacessloc; nof++)
      //     ptrEFv[nof] = ptrEFs[nof];
      //   ptrEFs += nfacessloc+1;
      //   ptrEFv += nfacessloc+3;
      // }
      
      // // construction des faces NGons
      // std::vector<E_Int> indices;
      // for (E_Int i = 0; i < nelts; i++)
      // {
      //   // les vertex surfaciques sont dans l'ordre rotatif
      //   indices.clear();
      //   K_CONNECT::getVertexIndices(cn->begin(), posFacesSurf.begin(), posEltsSurf[noe], indices);
      //   E_Int nvert = indices.size();

      //   //creation de la face correspondant a l elt surfacique
      //   ptrFNv[0] = nvert;
      //   for (E_Int i = 0; i < nvert; i++) ptrFNv[i+1] = indices[i];        
      //   ptrFNv+= nvert+1; 
      //   //creation de la face shiftee en z+1
      //   ptrFNv[0] = nvert;
      //   for (E_Int i = 0; i < nvert; i++) ptrFNv[i+1] = indices[i]+npts;        
      //   ptrFNv+= nvert+1;

      //   //modif de l'elt: on remplit les 2 derniers
      //   E_Int nfacesV = ptrEFv[0];
      //   ptrEFv[nfacesV-1] = nofv+1;
      //   ptrEFv[nfacesV]   = nofv+2;
      //   nofv += 2; ptrEFv += nfacesV+1;
      // }    

      // RELEASESHAREDU(tpl, f2, cn2);
      PyErr_SetString(PyExc_TypeError,
                      "addkplane: NGON array not supported yet.");
      RELEASESHAREDU(array, f, cn);
      return NULL;
    }
    
    RELEASESHAREDU(array, f, cn);
    return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "addkplane: unknow type of array.");
    return NULL;
  }
}

//=============================================================================
/* copie les champs de arrayC dans les centres de arrayK
   arrayC doit definir les champs en centres de l'array arrayK avant
   son addkplane(N)
   On fournit arrayK et pas seulement N car on a besoin de connaitre la
   nature de l array (a cause du passage de 2D a 3D si N>1) */
//=============================================================================
PyObject* K_TRANSFORM::addkplaneCenters(PyObject* self, PyObject* args)
{
  PyObject *arrayC, *arrayK;
  E_Int N;
  if (!PYPARSETUPLE_(args, OO_ I_,
                    &arrayC, &arrayK, &N)) return NULL;

  // Check array of centers
  E_Int imc, jmc, kmc;
  FldArrayF* fc; FldArrayI* cnc;
  char* varStringc; char* eltTypec;
  E_Int resc = K_ARRAY::getFromArray3(arrayC, varStringc, fc, imc, jmc, kmc,
                                      cnc, eltTypec);

  if (resc != 1 && resc != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "addkplane: unknown type of array for centers.");
    return NULL;
  }
  // Check array of nodes already with N kplanes
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(arrayK, varString, f, im, jm, km,
                                     cn, eltType);

  E_Bool isValidNGon = false, isValidME = false;
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "addkplane: unknown type of array for nodes.");
    return NULL;
  }
  if (res != resc)
  {
    RELEASESHAREDB(resc, arrayC, fc, cnc);
    RELEASESHAREDB(res, arrayK, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "addkplane: array of centers and nodes must be both structured or unstructured.");
    return NULL;
  }
  if (resc == 2)
  {
    isValidNGon = (
      K_STRING::cmp(eltTypec, "NGON*") == 0 && 
      (cnc->getDim() == 1 || cnc->getDim() == 2)
    );
    isValidME = (
      K_STRING::cmp(eltTypec, "NGON*") != 0 &&
      strchr(eltTypec, '*') != NULL &&
      (K_CONNECT::getDimME(eltTypec) == 1 || K_CONNECT::getDimME(eltTypec) == 2)
    );
    
    if (!(isValidNGon || isValidME))
    {
      RELEASESHAREDB(resc, arrayC, fc, cnc);
      RELEASESHAREDB(res, arrayK, f, cn);
      PyErr_SetString(PyExc_TypeError,
                      "addkplane: array of centers and nodes must be both structured or unstructured.");
      return NULL;
    }
  }

  E_Int api = fc->getApi();
  E_Int nfld = fc->getNfld();
  PyObject* tpl;
  FldArrayF* f2;

  if (resc == 1)
  {
    E_Int size = imc*jmc*kmc;
    E_Int imcjmc = imc*jmc;
    tpl = K_ARRAY::buildArray3(nfld, varStringc, imc, jmc, km-1, api);
    K_ARRAY::getFromArray3(tpl, f2);
    if (km == 2) // cas 2D
    {
      #pragma omp parallel
      {
        for (E_Int n = 1; n <= nfld; n++)
        {
          E_Float* fcp = fc->begin(n);
          E_Float* f2p = f2->begin(n);
          #pragma omp for nowait
          for (E_Int ind = 0; ind < size; ind++) f2p[ind] = fcp[ind];
        }
      }
    }
    else
    {
      #pragma omp parallel
      {
        E_Int offset;
        E_Float val;
        for (E_Int n = 1; n <= nfld; n++)
        {
          E_Float* fcp = fc->begin(n);
          E_Float* f2p = f2->begin(n);
          for (E_Int noz = 0; noz < km-1; noz++)
          {
            offset = noz*imcjmc;
            #pragma omp for nowait
            for (E_Int ind = 0; ind < imcjmc; ind++)
            {
              val = fcp[ind];
              f2p[ind] = val;
              f2p[ind+offset] = val;
            }
          }
        }
      }
    }
  }
  else
  {
    E_Int npts = f->getSize();
    E_Int neltsC = fc->getSize();

    tpl = K_ARRAY::buildArray3(nfld, varStringc, npts,
                               *cn, eltType, true, api, true);
    K_ARRAY::getFromArray3(tpl, f2);

    #pragma omp parallel
    {
      E_Int ind;
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fcp = fc->begin(n);
        E_Float* f2p = f2->begin(n);
        #pragma omp for collapse(2) nowait
        for (E_Int noz = 0; noz < N; noz++)
        for (E_Int i = 0; i < neltsC; i++)
        {
          ind = i + noz*neltsC;
          f2p[ind] = fcp[i];
        }
      }
    }
  }

  RELEASESHAREDS(tpl, f2);
  RELEASESHAREDB(resc, arrayC, fc, cnc);
  RELEASESHAREDB(res, arrayK, f, cn);
  return tpl;
}
