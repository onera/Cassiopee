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
# include "converter.h"
#include "Nuga/include/ngon_t.hxx"
#include "Nuga/include/Triangulator.h"

using namespace K_FLD;
using ngon_type = ngon_t<K_FLD::IntArray>;

//=============================================================================
/* Convert a NGon NFace numpy to a NGON Parent Element numpy */
//=============================================================================
PyObject* K_CONVERTER::adaptNFace2PE(PyObject* self, PyObject* args)
{
  PyObject* arrayNF; PyObject* arrayNG;  
  E_Int nfaces; E_Int nelts;
  E_Int methodPE(0);
  PyObject* arrayX; PyObject* arrayY; PyObject* arrayZ; 

  if (!PYPARSETUPLE_(args, OOOO_ O_ III_, &arrayNF,  &arrayNG, 
                     &arrayX, &arrayY, &arrayZ, 
                     &nelts, &nfaces, &methodPE)) return NULL;

  // Check numpy (NFace)
  FldArrayI* cNFace;
  E_Int res = K_NUMPY::getFromNumpyArray(arrayNF, cNFace);
  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "adaptNFace2PE: cNFace numpy is invalid.");
    return NULL;
  }
  // Check numpy (NGon)
  FldArrayI* cNGon;
  res = K_NUMPY::getFromNumpyArray(arrayNG, cNGon);
  if (res == 0)
  {
    RELEASESHAREDN(arrayNF, cNFace);
    PyErr_SetString(PyExc_TypeError, 
                    "adaptNFace2PE: cNGon numpy is invalid.");
    return NULL;
  }
  // Check numpy (CoordinateX)
  FldArrayF* coordX;
  res = K_NUMPY::getFromNumpyArray(arrayX, coordX);
  if (res == 0)
  {
    RELEASESHAREDN(arrayNF, cNFace);
    RELEASESHAREDN(arrayNG, cNGon);
    PyErr_SetString(PyExc_TypeError, 
                    "adaptNFace2PE: coordX numpy is invalid.");
    return NULL;
  }
  // Check numpy (CoordinateY)
  FldArrayF* coordY;
  res = K_NUMPY::getFromNumpyArray(arrayY, coordY);
  if (res == 0)
  {
    RELEASESHAREDN(arrayNF, cNFace);
    RELEASESHAREDN(arrayNG, cNGon);
    RELEASESHAREDN(arrayX, coordX);
    PyErr_SetString(PyExc_TypeError, 
                    "adaptNFace2PE: coordY numpy is invalid.");
    return NULL;
  }
  // Check numpy (CoordinateZ)
  FldArrayF* coordZ;
  res = K_NUMPY::getFromNumpyArray(arrayZ, coordZ);
  if (res == 0)
  {
    RELEASESHAREDN(arrayNF, cNFace);
    RELEASESHAREDN(arrayNG, cNGon);
    RELEASESHAREDN(arrayX, coordX);
    RELEASESHAREDN(arrayY, coordY);
    PyErr_SetString(PyExc_TypeError, 
                    "adaptNFace2PE: coordZ numpy is invalid.");
    return NULL;
  }
  PyObject* tpl = K_NUMPY::buildNumpyArray(nfaces, 2, 1, 1);
  E_Int* cFE = K_NUMPY::getNumpyPtrI(tpl);

  bool is_3D = (cNGon->begin()[0] > 2);
  if (!is_3D) methodPE = 0;

  if (methodPE == 0) // GEOMETRIC APPROACH => assume all cells are centroid-star-shaped.
  {
    E_Int* facesp1 = cFE;
    E_Int* facesp2 = cFE + nfaces;

    E_Int* ptrNF = cNFace->begin();
    E_Int face, nf;

    #pragma omp simd
    for (E_Int i = 0; i < nfaces*2; i++) cFE[i] = 0;

    for (E_Int i = 0; i < nelts; i++)
    {
      nf = ptrNF[0]; // nb de faces pour l'elt
      for (E_Int j = 1; j <= nf; j++)
      {
        face = ptrNF[j]-1;
        if (facesp1[face] == 0) facesp1[face] = i+1;
        else facesp2[face] = i+1;
      }
      ptrNF += nf+1;
    }
  
    // post - reordonne les PEs pour que les normales soient exterieures a gauche
    E_Int* ptrNG = cNGon->begin();
    ptrNF = cNFace->begin();
    nf = ptrNF[0]; // nb de faces pour l'elt
    if (nf > 2) // NGON volumique
    { 
      FldArrayI posElts(nelts,1); FldArrayI posFaces(nfaces,1);
      K_CONNECT::getPosFacets(cNFace->begin(), 0, nelts, posElts);
      //E_Int* posEltsp = posElts.begin();
      K_CONNECT::getPosFacets(cNGon->begin(), 0, nfaces, posFaces);
      //E_Int* posFacesp = posFaces.begin();

      E_Int nvert, ind1, ind2, etg, etd;
      E_Float sx, sy, sz, ps;
      //E_Float ps1, ps2;
      E_Float l1x, l1y , l1z, l2x, l2y, l2z;
      E_Float surfnx, surfny, surfnz;
      E_Float xbf, ybf, zbf, xbg, ybg, zbg, xbd, ybd, zbd;
      E_Float dx, dy, dz;
      //E_Float dx1, dy1, dz1, dx2, dy2, dz2;
      E_Float* xt = coordX->begin();
      E_Float* yt = coordY->begin();
      E_Float* zt = coordZ->begin();
      E_Int* ptrNG2 = cNGon->begin();
      for (E_Int fa = 0; fa < nfaces; fa++)
      {     
        // 1. calcul de la normale a la face
        K_COMPGEOM::getNGONFaceBarycenter(0, ptrNG2, xt, yt, zt, xbf, ybf, zbf);

        nvert = ptrNG2[0];
        sx = 0.; sy = 0.; sz = 0.;
        for (E_Int nv = 1; nv < nvert; nv++)
        {
          ind1 = ptrNG2[nv]-1; ind2 = ptrNG2[nv+1]-1;
          // calcul de la normale au triangle (nv, nv+1, bf)
          l1x = xt[ind1]-xbf; l1y = yt[ind1]-ybf; l1z = zt[ind1]-zbf;
          l2x = xt[ind2]-xbf; l2y = yt[ind2]-ybf; l2z = zt[ind2]-zbf;
          surfnx = l1y*l2z-l1z*l2y;
          surfny = l1z*l2x-l1x*l2z;
          surfnz = l1x*l2y-l1y*l2x;
          sx += surfnx; sy += surfny; sz += surfnz;       
        }
        // le dernier et le premier pour boucler
        ind1 = ptrNG2[nvert]-1; ind2 = ptrNG2[1]-1;
        l1x = xt[ind1]-xbf; l1y = yt[ind1]-ybf; l1z = zt[ind1]-zbf;
        l2x = xt[ind2]-xbf; l2y = yt[ind2]-ybf; l2z = zt[ind2]-zbf;
        surfnx = l1y*l2z-l1z*l2y;
        surfny = l1z*l2x-l1x*l2z; 
        surfnz = l1x*l2y-l1y*l2x;
        sx += surfnx; sy += surfny; sz += surfnz;

        etg = facesp1[fa]-1; etd = facesp2[fa]-1;
        if (etg > -1 && etd > -1)
        {
          // 2. calcul du barycentre de l'elt gauche G
          K_COMPGEOM::getNGONEltBarycenter(etg, posElts.begin(), posFaces.begin(),
                                         ptrNF, ptrNG, xt, yt, zt, xbg, ybg, zbg);

          // 3. calcul du barycentre de l'elt droit D
          K_COMPGEOM::getNGONEltBarycenter(etd, posElts.begin(), posFaces.begin(),
                                         ptrNF, ptrNG, xt, yt, zt, xbd, ybd, zbd);

          // 4. produit scalaire de n.GD
          dx = xbd-xbg; dy = ybd-ybg; dz = zbd-zbg;
          ps = sx*dx+sy*dy+sz*dz;
          //dx1 = xbd-xbf; dy1 = ybd-ybf; dz1 = zbd-zbf;
          //ps1 = sx*dx1+sy*dy1+sz*dz1;
          //dx2 = xbf-xbg; dy2 = ybf-ybg; dz2 = zbf-zbg;
          //ps2 = sx*dx2+sy*dy2+sz*dz2;
          //if (ps*ps1 < 0 || ps*ps2 < 0) printf("face=%d, ps %f %f %f\n", fa, ps, ps1,ps2);
          if (ps < 0.)
          {
            facesp2[fa] = etg+1; facesp1[fa] = etd+1;
          }
        }
        else if (etg == -1 && etd > -1)
        {
          // force inversion (elsA)
          facesp1[fa] = etd+1; facesp1[fa] = etg+1;
          etg = etd; etd = -1;
          K_COMPGEOM::getNGONEltBarycenter(etg, posElts.begin(), posFaces.begin(),
                                         ptrNF, ptrNG, xt, yt, zt, xbg, ybg, zbg);
          dx = xbf-xbg; dy = ybf-ybg; dz = zbf-zbg;
          ps = sx*dx+sy*dy+sz*dz;
          if (ps < 0.)
          {
            // renumerote NGON Face
            nvert = ptrNG2[0];
            std::vector<E_Int> sav(nvert);
            for (E_Int nv = 0; nv < nvert; nv++) sav[nv] = ptrNG2[nv+1];
            for (E_Int nv = 0; nv < nvert; nv++) ptrNG2[nv+1] = sav[nvert-1-nv];
          }
        }
        else if (etg > -1 && etd == -1)
        {
          K_COMPGEOM::getNGONEltBarycenter(etg, posElts.begin(), posFaces.begin(),
                                           ptrNF, ptrNG, xt, yt, zt, xbg, ybg, zbg);
          dx = xbf-xbg; dy = ybf-ybg; dz = zbf-zbg;
          ps = sx*dx+sy*dy+sz*dz;
          if (ps < 0.)
          {
            // renumerote NGON Face
            nvert = ptrNG2[0];
            std::vector<E_Int> sav(nvert);
            for (E_Int nv = 0; nv < nvert; nv++) sav[nv] = ptrNG2[nv+1];         
            for (E_Int nv = 0; nv < nvert; nv++) ptrNG2[nv+1] = sav[nvert-1-nv];
          }
        }
        ptrNG2 += nvert+1;
      }
    }

    // Check
    bool check=false;
    if (check)
    {
      FldArrayF bilan(nelts,3); bilan.setAllValuesAtNull();
      E_Float* bilanx = bilan.begin(1);
      E_Float* bilany = bilan.begin(2);
      E_Float* bilanz = bilan.begin(3);
      E_Int nvert, ind1, ind2, etg, etd;
      E_Float sx, sy, sz;
      E_Float l1x, l1y , l1z, l2x, l2y, l2z;
      E_Float surfnx, surfny, surfnz;
      E_Float xbf, ybf, zbf;
      //E_Float dx, dy, dz;
      E_Float* xt = coordX->begin();
      E_Float* yt = coordY->begin();
      E_Float* zt = coordZ->begin();

      E_Int* ptrNG2 = cNGon->begin();
      for (E_Int fa = 0; fa < nfaces; fa++)
      {
        K_COMPGEOM::getNGONFaceBarycenter(0, ptrNG2, xt, yt, zt, xbf, ybf, zbf);
        nvert = ptrNG2[0];
        etg = facesp1[fa]-1;
        etd = facesp2[fa]-1;
        
        sx = 0.; sy = 0.; sz = 0.;
        for (E_Int nv = 1; nv < nvert; nv++)
        {
          ind1 = ptrNG2[nv]-1; ind2 = ptrNG2[nv+1]-1;
          // calcul de la normale au triangle (nv, nv+1, bf)
          l1x = xt[ind1]-xbf; l1y = yt[ind1]-ybf; l1z = zt[ind1]-zbf;
          l2x = xt[ind2]-xbf; l2y = yt[ind2]-ybf; l2z = zt[ind2]-zbf;
          surfnx = l1y*l2z-l1z*l2y;
          surfny = l1z*l2x-l1x*l2z;
          surfnz = l1x*l2y-l1y*l2x;
          sx += surfnx; sy += surfny; sz += surfnz;       
        }
        // le dernier et le premier pour boucler
        ind1 = ptrNG2[nvert]-1; ind2 = ptrNG2[1]-1;
        l1x = xt[ind1]-xbf; l1y = yt[ind1]-ybf; l1z = zt[ind1]-zbf;
        l2x = xt[ind2]-xbf; l2y = yt[ind2]-ybf; l2z = zt[ind2]-zbf;
        surfnx = l1y*l2z-l1z*l2y;
        surfny = l1z*l2x-l1x*l2z; 
        surfnz = l1x*l2y-l1y*l2x;
        sx += surfnx; sy += surfny; sz += surfnz;

        if (etg >= 0)
        {
          bilanx[etg] += sx;
          bilany[etg] += sy;
          bilanz[etg] += sz;
        }
        else
        {
          bilanx[etd] += 1000;
          bilany[etd] += 1000;
          bilanz[etd] += 1000;
        }
        if (etd >= 0)
        {
          bilanx[etd] -= sx;
          bilany[etd] -= sy;
          bilanz[etd] -= sz;
        }
        else
        {
          bilanx[etg] += 1000;
          bilany[etg] += 1000;
          bilanz[etg] += 1000;
        }
        ptrNG2 += nvert+1;
      }
      for (E_Int e = 0; e < nelts; e++)
      {
        printf(SF_D_ ": " SF_F3_ "\n", e, bilanx[e], bilany[e], bilanz[e]);
      }
    }
  }
  else if (methodPE == 1) // TOPOLOGICAL => GENERAL CASE
  {
    E_Float* xt = coordX->begin();
    E_Float* yt = coordY->begin();
    E_Float* zt = coordZ->begin();
    ngon_type ng(cNGon->begin(), cNGon->getSize(), nfaces, cNFace->begin(), cNFace->getSize(), nelts);
    K_FLD::FloatArray crd(3, coordX->getSize());
    #pragma omp simd
    for (E_Int i = 0; i < coordX->getSize(); ++i)
    {
      crd(0,i) = xt[i];
      crd(1,i) = yt[i];
      crd(2,i) = zt[i];
    }

    // reorient skin
    ng.flag_externals(1);
    DELAUNAY::Triangulator dt;
    bool has_been_reversed;
    E_Int err = ngon_type::reorient_skins(dt, crd, ng, has_been_reversed);
    if (err)
    {
      RELEASESHAREDN(arrayNF, cNFace);
      RELEASESHAREDN(arrayNG, cNGon);
      RELEASESHAREDN(arrayX, coordX);
      RELEASESHAREDN(arrayY, coordY);
      PyErr_SetString(PyExc_TypeError, 
                      "adaptNFace2PE: method 2: failed to reorient external faces.");
      return NULL;
    }

    // build oriented F2E
    K_FLD::IntArray F2E2;
    ngon_unit neighbors;
    ng.build_ph_neighborhood(neighbors);
    ng.build_F2E(neighbors, F2E2);

    // fill cFE
    #pragma omp simd
    for (E_Int i = 0; i < nfaces; ++i)
    {
      E_Int fli = F2E2(0,i);
      E_Int fri = F2E2(1,i);
      cFE[i] = (fli == E_IDX_NONE) ? 0 : fli+1;          //left elt
      cFE[i + nfaces] = (fri == E_IDX_NONE) ? 0 : fri+1; //right elt
    }

    // replace NGON/NFACE
    E_Int* ptrNG2 = cNGon->begin();
    std::copy (ng.PGs._NGON.begin()+2, ng.PGs._NGON.end(), ptrNG2);
    //E_Int* ptrNF = cNFace->begin();
    //std::copy ( ng.PHs._NGON.begin()+2, ng.PHs._NGON.end(), ptrNF );
  }
  else // Unoriented ParentElements
  {
    E_Int* facesp1 = cFE;
    E_Int* facesp2 = cFE + nfaces;
    E_Int* ptrNF = cNFace->begin();
    E_Int face, nf;

    #pragma omp simd
    for (E_Int i = 0; i < nfaces*2; i++) cFE[i] = 0;

    for (E_Int i = 0; i < nelts; i++)
    {
      nf = ptrNF[0]; // nb de faces pour l'elt
      for (E_Int j = 1; j <= nf; j++)
      {
        face = ptrNF[j]-1;
        if (facesp1[face] == 0) facesp1[face] = i+1;
        else facesp2[face] = i+1;
      }
      ptrNF += nf+1;
    }
  }
  
  RELEASESHAREDN(arrayNF, cNFace);
  RELEASESHAREDN(arrayNG, cNGon);
  RELEASESHAREDN(arrayX, coordX);
  RELEASESHAREDN(arrayY, coordY);
  RELEASESHAREDN(arrayZ, coordZ);

  // export du numpy de sortie
  return tpl;
}
