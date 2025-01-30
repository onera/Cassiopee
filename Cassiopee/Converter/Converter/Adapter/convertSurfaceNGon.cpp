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

using namespace K_FLD;
using namespace K_FUNC;
using namespace std;

//define PRINTINFO

//=============================================================================
/* Convert a surface NGon */
/* if INPUT is a NGON=bars, NFACE=polygon (Type A) -> OUTPUT is NGON=polygon, NFACE=NULL (Type B)
   if INPUT is NGON=polygon, NFACE=NULL (Type B) -> OUTPUT a NGON=bars, NFACE=polygon (Type A) */
//=============================================================================
PyObject* K_CONVERTER::convertSurfaceNGon(PyObject* self, PyObject* args)
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
    E_Int isNGon = c->isNGon(); 
    E_Int nfacesA = 0;
    E_Int neltsA = nfaces;
    E_Int sizeNGonA = 0;
    E_Int sizeNFaceA = sizeNGon;
    #ifdef PRINTINFO
    printf("Info: convertSurfaceNGon: convert type B to type A\n");
    #endif
    // Initialize sizeNGonA & sizeNFaceA
    std::set<vector<E_Int>> ngonSetsA;
    std::vector<E_Int> edge{0,0};
    E_Int v1, v2;
    for (E_Int j = 0; j < nfaces; j++)
    {
      E_Int* face = c->getFace(j, sizeFace, ngon, indPG);
      for (E_Int i = 0; i < sizeFace-1; i++)
      {
        v1 = face[i];
        v2 = face[i+1];
        edge[0] = E_min(v1,v2);
        edge[1] = E_max(v1,v2);
        ngonSetsA.insert(edge);
      }
      v1 = face[0];
      v2 = face[sizeFace-1];
      edge[0] = E_min(v1,v2);
      edge[1] = E_max(v1,v2);
      ngonSetsA.insert(edge);
    }

    nfacesA = ngonSetsA.size();
    sizeNGonA = 2*nfacesA;
    if (isNGon < 3) sizeNGonA = sizeNGonA+nfacesA;

    // TYPE B
    #ifdef PRINTINFO
    printf("Info: convertSurfaceNGon: universel NGON (type B): nbre de faces=%d, nbre d'elements=%d\n", nfaces, nelts);
    printf("Info: convertSurfaceNGon: universel NGON (type B): sizeNGon=%d, sizeNFace=%d\n", sizeNGon, sizeNFace);
    #endif
    tpl = K_ARRAY::buildArray3(f->getNfld(), varString, f->getSize(), neltsA, nfacesA, 
                               "NGON", sizeNGonA, sizeNFaceA, isNGon, false, 3);
    K_FLD::FldArrayF* fo; K_FLD::FldArrayI* co;
    K_ARRAY::getFromArray3(tpl, fo, co);
    // co->setNGon(isNGon); // force NGonType

    E_Int* ngonA = co->getNGon();
    E_Int* nfaceA = co->getNFace();
    E_Int* indPGA = co->getIndPG(); 
    E_Int* indPHA = co->getIndPH();

    // TYPE A
    #ifdef PRINTINFO
    printf("Info: convertSurfaceNGon: universel NGON (type A): nbre de faces=%d, nbre d'elements=%d\n", nfacesA, neltsA);
    printf("Info: conertSurfaceNGon: universel NGON (type A): sizeNGon=%d, sizeNFace=%d\n", sizeNGonA, sizeNFaceA);
    #endif
    // fill ngonA
    E_Int shiftNGonA = 0;
    E_Int cptNGonA = 1;
    indPGA[0] = 0;
    for (auto it : ngonSetsA) 
    {
      if (isNGon < 3) 
      {
        ngonA[shiftNGonA] = 2;
        shiftNGonA++;
      }
      ngonA[shiftNGonA] = it[0];
      ngonA[shiftNGonA+1] = it[1];
      shiftNGonA = shiftNGonA+2;

      
      if ((isNGon < 3 && cptNGonA == nfacesA) == false) indPGA[cptNGonA] = shiftNGonA;
      cptNGonA++;
    }

    // fill nfaceA
    E_Int shiftNFaceA = 0;
    E_Int facePos = 0;
    indPHA[0] = 0;
    for (E_Int j = 0; j < nfaces; j++)
    {
      E_Int* face = c->getFace(j, sizeFace, ngon, indPG);
      std::vector<E_Int> faceA;
      for (E_Int i = 0; i < sizeFace-1; i++)
      {
        v1 = face[i];
        v2 = face[i+1];
        edge[0] = E_min(v1,v2);
        edge[1] = E_max(v1,v2);
        auto pos = ngonSetsA.find(edge);
        facePos = distance(ngonSetsA.begin(), pos);
        faceA.push_back(facePos+1);
      }
      v1 = face[0];
      v2 = face[sizeFace-1];
      edge[0] = E_min(v1,v2);
      edge[1] = E_max(v1,v2);
      auto pos = ngonSetsA.find(edge);
      facePos = distance(ngonSetsA.begin(), pos);
      faceA.push_back(facePos+1);

      if (isNGon < 3) 
      {
        nfaceA[shiftNFaceA] = faceA.size();
        shiftNFaceA++;
      }
      for (size_t i = 0; i < faceA.size(); i++) nfaceA[i+shiftNFaceA] = faceA[i];
      shiftNFaceA = shiftNFaceA+faceA.size();
      if ((isNGon < 3 && j == nfaces-1) == false) indPHA[j+1] = shiftNFaceA;
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
    E_Int isNGon = c->isNGon(); 
    E_Int nfacesB = nelts;
    E_Int neltsB = 0;
    E_Int sizeNGonB = sizeNFace;
    E_Int sizeNFaceB = 0;

    // TYPE A
    #ifdef PRINTINFO
    printf("Info: convertSurfaceNGon: convert type A to type B\n");
    printf("Info: convertSurfaceNGon: universel NGON (type A): nbre de faces=%d, nbre d'elements=%d\n", nfaces, nelts);
    printf("Info: convertSurfaceNGon: universel NGON (type A): sizeNGon=%d, sizeNFace=%d\n", sizeNGon, sizeNFace);
    #endif
    tpl = K_ARRAY::buildArray3(f->getNfld(), varString, f->getSize(), neltsB, nfacesB, 
                               "NGON", sizeNGonB, sizeNFaceB, isNGon, false, 3);
    K_FLD::FldArrayF* fo; K_FLD::FldArrayI* co;
    K_ARRAY::getFromArray3(tpl, fo, co);
    // co->setNGon(isNGon); // force NGonType

    E_Int* ngonB = co->getNGon();
    E_Int* indPGB = co->getIndPG(); 

    // TYPE B
    #ifdef PRINTINFO
    printf("Info: convertSurfaceNGon: universel NGON (type B): nbre de faces=%d, nbre d'elements=%d\n", nfacesB, neltsB);
    printf("Info: convertSurfaceNGon: universel NGON (type B): sizeNGon=%d, sizeNFace=%d\n", sizeNGonB, sizeNFaceB);
    #endif
    // fill indPGB & ngonB
    E_Int shift = 0;
    indPGB[0] = 0;
    for (E_Int k = 0; k < nelts; k++)
    {
      E_Int* elt = c->getElt(k, sizeElt, nface, indPH);
      std::vector<E_Int> faceB;
      for (E_Int i = 0; i < sizeElt-1; i++)
      {
        E_Int* face = c->getFace(elt[i]-1, sizeFace, ngon, indPG); // Faces are bar types (two vertices per face)
        if (faceB.size() < 1)
        {
          faceB.push_back(face[0]);
          faceB.push_back(face[1]);
        }
        else if (face[0]-1 == faceB[faceB.size()-1])
          faceB.push_back(face[1]);
        else
          faceB.insert(faceB.begin(), face[1]);
      }

      if (isNGon < 3) 
      {
        ngonB[shift] = faceB.size();
        shift++;
      }
      
      for (size_t i = 0; i < faceB.size(); i++) ngonB[i+shift] = faceB[i];

      shift = shift+faceB.size();
      if ((isNGon < 3 && k == nelts-1) == false) indPGB[k+1] = shift;
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
