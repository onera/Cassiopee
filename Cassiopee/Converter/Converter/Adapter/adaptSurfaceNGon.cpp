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

# include <unordered_map>
# include "converter.h"

using namespace K_FLD;
using namespace K_FUNC;
using namespace std;

# define PRINTINFO 0

//=============================================================================
/* Convert a surface NGon */
/* if INPUT is a NGON=bars, NFACE=polygon (Type A) -> OUTPUT is NGON=polygon, NFACE=NULL (Type B)
   if INPUT is NGON=polygon, NFACE=NULL (Type B) -> OUTPUT a NGON=bars, NFACE=polygon (Type A) */
//=============================================================================
PyObject* K_CONVERTER::adaptSurfaceNGon(PyObject* self, PyObject* args)
{
  PyObject* o; 
  if (!PYPARSETUPLE_(args, O_, &o)) return NULL;

  // Check array
  E_Int ni, nj, nk;
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* c;
  char* varString; char* eltType;
  E_Int ret = K_ARRAY::getFromArray3(o, varString, f, ni, nj, nk, c, eltType);
  
  if (ret <= 0)
  { PyErr_SetString(PyExc_TypeError, "adaptSurfaceNGon: only for NGons."); return NULL; }

  if (ret == 1)
  { 
    PyErr_SetString(PyExc_TypeError, "adaptSurfaceNGon: only for NGons."); 
    RELEASESHAREDS(o, f);
    return NULL;
  }

  // Analyse input
  E_Int nelts = c->getNElts();
  E_Int nfaces = c->getNFaces();
  E_Int* ngon = c->getNGon();
  E_Int* nface = c->getNFace();
  E_Int* indPG = c->getIndPG(); 
  E_Int* indPH = c->getIndPH();

  E_Int sizeNGon = c->getSizeNGon();
  E_Int sizeNFace = c->getSizeNFace();

  E_Int sizeElt;
  E_Int sizeFace;

  E_Int nfld = f->getNfld();
  E_Int npts = f->getSize();

  PyObject* tpl = NULL;
  if (nelts == 0) // INPUT is type B
  {
    E_Int ngonType = c->getNGonType();
    E_Int nfacesA = 0;
    E_Int neltsA = nfaces;
    E_Int sizeNGonA = 0;
    E_Int sizeNFaceA = sizeNGon;

    // Initialize sizeNGonA & sizeNFaceA
    std::vector<E_Int> edge{0,0};
    Topology E;
    std::unordered_map<Topology, E_Int, BernsteinHash<Topology> > ngonMapA;
    std::vector<std::pair<E_Int,E_Int> > ngonListA;

    for (E_Int j = 0; j < nfaces; j++)
    {
      E_Int* face = c->getFace(j, sizeFace, ngon, indPG);
      for (E_Int i = 0; i < sizeFace; i++)
      {
        edge[0] = face[i];
        edge[1] = face[(i+1)%sizeFace];
        E.set(edge);
        // try to insert topological edge element (E) in ngonMapA
        auto res = ngonMapA.insert({E, nfacesA+1});
        if (res.second) // check insertion outcome
        {
          // if E was not already in ngonMapA, add the new bar face in ngonListA and increase nfacesA counter
          ngonListA.push_back({edge[0], edge[1]});
          nfacesA++;
        }
      }
    }

    sizeNGonA = 2*nfacesA;
    if (ngonType < 3) sizeNGonA = sizeNGonA+nfacesA;

    // TYPE B
    #if PRINTINFO
    printf("Info: adaptSurfaceNGon: universel NGON (type B): nbre de faces=" SF_D_ ", nbre d'elements=" SF_D_ "\n", nfaces, nelts);
    printf("Info: adaptSurfaceNGon: universel NGON (type B): sizeNGon=" SF_D_ ", sizeNFace=" SF_D_ "\n", sizeNGon, sizeNFace);
    #endif
    tpl = K_ARRAY::buildArray3(nfld, varString, npts, neltsA, nfacesA, 
                               "NGON", sizeNGonA, sizeNFaceA, ngonType, false, 3);
    K_FLD::FldArrayF* fo; K_FLD::FldArrayI* co;
    K_ARRAY::getFromArray3(tpl, fo, co);

    E_Int* ngonA = co->getNGon();
    E_Int* nfaceA = co->getNFace();
    E_Int* indPGA = co->getIndPG(); 
    E_Int* indPHA = co->getIndPH();

    // TYPE A
    #if PRINTINFO
    printf("Info: adaptSurfaceNGon: universel NGON (type A): nbre de faces=" SF_D_ ", nbre d'elements=" SF_D_ "\n", nfacesA, neltsA);
    printf("Info: adaptSurfaceNGon: universel NGON (type A): sizeNGon=" SF_D_ ", sizeNFace=" SF_D_ "\n", sizeNGonA, sizeNFaceA);
    #endif

    // fill ngonA
    E_Int shiftNGonA = 0;
    indPGA[0] = 0;
    for (E_Int i = 0; i < nfacesA; i++)
    {
      if (ngonType < 3) 
      {
        ngonA[shiftNGonA] = 2;
        if (i != nfacesA-1) indPGA[i+1] = indPGA[i]+3;
        shiftNGonA++;
      }
      else indPGA[i+1] = indPGA[i]+2;

      ngonA[shiftNGonA] = ngonListA[i].first;
      ngonA[shiftNGonA+1] = ngonListA[i].second;
      shiftNGonA += 2;
    }

    // fill nfaceA
    E_Int shiftNFaceA = 0;
    indPHA[0] = 0;
    for (E_Int j = 0; j < nfaces; j++)
    {
      E_Int* face = c->getFace(j, sizeFace, ngon, indPG);

      if (ngonType < 3) 
      {
        nfaceA[shiftNFaceA] = sizeFace;
        if (j != nfaces-1) indPHA[j+1] = indPHA[j]+sizeFace+1;
        shiftNFaceA++;
      }
      else indPHA[j+1] = indPHA[j]+sizeFace;

      for (E_Int i = 0; i < sizeFace; i++)
      {
        edge[0] = face[i];
        edge[1] = face[(i+1)%sizeFace];
        E.set(edge);
        nfaceA[shiftNFaceA] = ngonMapA[E]; // get the unique E_Int ID of the face
        shiftNFaceA++;
      }
    }

    // set fields
    #pragma omp parallel
    {
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fp = f->begin(n);
        E_Float* f2p = fo->begin(n);
        #pragma omp for
        for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
      }
    }

    RELEASESHAREDU(tpl, fo, co);
  }
  else // INPUT is type A
  {
    E_Int ngonType = c->getNGonType();
    E_Int nfacesB = nelts;
    E_Int neltsB = 0;
    E_Int sizeNGonB = sizeNFace;
    E_Int sizeNFaceB = 0;

    // TYPE A
    #if PRINTINFO
    printf("Info: adaptSurfaceNGon: convert type A to type B\n");
    printf("Info: adaptSurfaceNGon: universel NGON (type A): nbre de faces=" SF_D_ ", nbre d'elements=" SF_D_ "\n", nfaces, nelts);
    printf("Info: adaptSurfaceNGon: universel NGON (type A): sizeNGon=" SF_D_ ", sizeNFace=" SF_D_ "\n", sizeNGon, sizeNFace);
    #endif
    tpl = K_ARRAY::buildArray3(nfld, varString, npts, neltsB, nfacesB, 
                               "NGON", sizeNGonB, sizeNFaceB, ngonType, false, 3);
    K_FLD::FldArrayF* fo; K_FLD::FldArrayI* co;
    K_ARRAY::getFromArray3(tpl, fo, co);

    E_Int* ngonB = co->getNGon();
    E_Int* indPGB = co->getIndPG(); 

    // TYPE B
    #if PRINTINFO
    printf("Info: adaptSurfaceNGon: universel NGON (type B): nbre de faces=" SF_D_ ", nbre d'elements=" SF_D_ "\n", nfacesB, neltsB);
    printf("Info: adaptSurfaceNGon: universel NGON (type B): sizeNGon=" SF_D_ ", sizeNFace=" SF_D_ "\n", sizeNGonB, sizeNFaceB);
    #endif

    // fill ngonB
    E_Int shiftNGonB = 0;
    indPGB[0] = 0;
    for (E_Int k = 0; k < nelts; k++)
    {
      E_Int* elt = c->getElt(k, sizeElt, nface, indPH);
      if (ngonType < 3) 
      {
        ngonB[shiftNGonB] = sizeElt;
        if (k != nelts-1) indPGB[k+1] = indPGB[k]+sizeElt+1;
        shiftNGonB++;
      }
      else indPGB[k+1] = indPGB[k]+sizeElt;

      // First, reorder the first three vertices
      // Then, simply add new vertices based on the last entry
      E_Int* face1 = c->getFace(elt[0]-1, sizeFace, ngon, indPG);
      E_Int* face2 = c->getFace(elt[1]-1, sizeFace, ngon, indPG);

      if (face1[0] == face2[0])
      {
        ngonB[shiftNGonB+0] = face1[1];
        ngonB[shiftNGonB+1] = face1[0];
        ngonB[shiftNGonB+2] = face2[1];
        shiftNGonB += 3;
      }
      else if (face1[1] == face2[1])
      {
        ngonB[shiftNGonB+0] = face1[0];
        ngonB[shiftNGonB+1] = face1[1];
        ngonB[shiftNGonB+2] = face2[0];
        shiftNGonB += 3;
      }
      else if (face1[1] == face2[0])
      {
        ngonB[shiftNGonB+0] = face1[0];
        ngonB[shiftNGonB+1] = face1[1];
        ngonB[shiftNGonB+2] = face2[1];
        shiftNGonB += 3;
      }
      else
      {
        ngonB[shiftNGonB+0] = face1[1];
        ngonB[shiftNGonB+1] = face1[0];
        ngonB[shiftNGonB+2] = face2[0];
        shiftNGonB += 3;
      }

      for (E_Int i = 2; i < sizeElt-1; i++)
      {
          E_Int* face = c->getFace(elt[i]-1, sizeFace, ngon, indPG);
          if (face[0] == ngonB[shiftNGonB-1]) ngonB[shiftNGonB] = face[1];
          else ngonB[shiftNGonB] = face[0];
          shiftNGonB++;
      }
    }

    // set fields
    #pragma omp parallel
    {
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fp = f->begin(n);
        E_Float* f2p = fo->begin(n);
        #pragma omp for
        for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
      }
    }

    RELEASESHAREDU(tpl, fo, co);
  }

  RELEASESHAREDU(o, f, c);

  return tpl;
}
