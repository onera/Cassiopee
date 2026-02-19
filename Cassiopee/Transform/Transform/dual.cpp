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
#include <algorithm>

using namespace K_FLD;
using namespace std;

// ============================================================================
/* Cree le dual d'un maillage NGON conforme */
// ============================================================================
PyObject* K_TRANSFORM::dualNGon(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int extraPoints; // 0: pas de pts en plus des centres, 1: points en plus sur les faces externes
  if (!PYPARSETUPLE_(args, O_ I_, &array, &extraPoints)) return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km,
                                     cn, eltType);

  if (res != 2)
  {
    if (res == 1) RELEASESHAREDS(array,f);
    PyErr_SetString(PyExc_TypeError,
                    "dualNGon: invalid type of array.");
    return NULL;
  }
  if (strcmp(eltType, "NGON") != 0)
  {
    RELEASESHAREDU(array,f,cn);
    PyErr_SetString(PyExc_TypeError,
                    "dualNGon: array must be NGON.");
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "dualNGon: array must contain coordinates.");
    RELEASESHAREDU(array,f,cn); return NULL;
  }
  posx++; posy++; posz++;

  FldArrayF fd;
  FldArrayI cNGD;
  E_Int api = f->getApi();
  E_Int dim = cn->getDim();
  cNGD.setNGonType(1);
  if (dim == 1) dualNGON1D(*f, *cn, extraPoints, fd, cNGD);
  else if (dim == 2) dualNGON2D(*f, *cn, extraPoints, fd, cNGD);
  else if (dim == 3)
  {
    FldArrayI cn2; FldArrayF f2;
    if (extraPoints == 1) // ajoute les pts sur les faces externes
    {
      E_Int ext = K_TRANSFORM::createDegeneratedPrimalMesh3D(*f, *cn, f2, cn2);
      if (ext != 0) dualNGON3D(f2, cn2, fd, cNGD);
      else dualNGON3D(*f, *cn, fd, cNGD);
    }
    else dualNGON3D(*f, *cn, fd, cNGD);
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "dualNGon: array must be 1D, 2D or 3D.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }
  // build array
  RELEASESHAREDU(array, f, cn);
  PyObject* tpl = NULL;
  if (extraPoints == 1)
  {
    tpl = K_CONNECT::V_cleanConnectivityNGon(posx, posy, posz, varString,
                                             fd, cNGD, 1.e-10);
    return tpl;
  }
  else
  {
    tpl = K_ARRAY::buildArray3(fd, varString, cNGD, eltType, api);
  }
  return tpl;
}
//=============================================================================
/* dualNGON 3D */
//=============================================================================
void K_TRANSFORM::dualNGON3D(FldArrayF& f, FldArrayI& cn,
                             FldArrayF& fd, FldArrayI& cNGD)
{
  // Cree les centres du maillage primal = les noeuds du maillage dual
  E_Int nfld = f.getNfld();
  E_Int nptsp = f.getSize();  // nb de pts ds le maillage primal
  E_Int* ngon = cn.getNGon(); E_Int* indPG = cn.getIndPG();
  E_Int neltsp = cn.getNElts();
  E_Int sizeFNp = cn.getSizeNGon();
  E_Int sizeEFp = cn.getSizeNFace();

  FldArrayI cFNd(2*sizeFNp);
  FldArrayI cEFd(2*sizeEFp);
  E_Int* ptrFN = cFNd.begin();
  E_Int* ptrEF = cEFd.begin();

  fd.malloc(neltsp, nfld);  // nb de pts ds le dual = nb d'elts dans le primal
  K_LOC::node2centerNGon(f, cn, fd);

  // calcul de la connectivite Vertex/Faces primale
  vector<vector<E_Int> > cVFp(nptsp);
  K_CONNECT::connectNG2VF(cn, cVFp);

  // calcul de la connectivite faces/elts primale
  FldArrayI cFEp; K_CONNECT::connectNG2FE(cn, cFEp);
  E_Int* cFEp1 = cFEp.begin(1);
  E_Int* cFEp2 = cFEp.begin(2);

  E_Int currentFace = 0;  // nb de faces totales
  E_Int sizeEFd = 0; E_Int sizeFNd = 0;
  E_Int etg, etd, face0, face1, found, indface, nvert, vsize, sizeFacesPP0;
  E_Int next, start, cf;
  vector<E_Int> facesPP0, facesPP;
  FldArrayI dejaVu;
  E_Int compt = 0;
  E_Int neltsfin = 0;  // nb d elts ds le dual

  for (E_Int indv = 0; indv < nptsp; indv++)
  {
    E_Int noface = 0;
    vector<E_Int>& faces = cVFp[indv];
    E_Int nfacesv = faces.size();
    // Determination des sommets P' differents de P appartenant a des faces contenant P
    vector<E_Int> vertices; // sommets differents de indv appartenant aux faces
    for (E_Int nof = 0; nof < nfacesv; nof++)
    {
      indface = faces[nof]-1;
      E_Int* face = cn.getFace(indface, nvert, ngon, indPG);
      for (E_Int nov = 0; nov < nvert; nov++)
      {
        vertices.push_back(face[nov]-1);
      }
    }
    sort(vertices.begin(), vertices.end());
    vertices.erase(unique(vertices.begin(), vertices.end()), vertices.end());
    vsize = vertices.size(); // attention peut etre decremente

    for (size_t nov = 0; nov < vertices.size(); nov++)
    {
      E_Int indv2 = vertices[nov];
      facesPP0.clear();
      // trouver toutes les faces contenant l edge [indv-indv2]
      for (E_Int nof = 0; nof < nfacesv; nof++)
      {
        found = 0;
        indface = faces[nof]-1;
        E_Int* face = cn.getFace(indface, nvert, ngon, indPG);
        for (E_Int nov = 0; nov < nvert; nov++)
        {
          if (face[nov]-1 == indv || face[nov]-1 == indv2) found++;
        }
        if (found >= 2) facesPP0.push_back(indface);
      }

      // Ordonnnancement cyclique des faces de facesPP
      facesPP.clear();
      sort(facesPP0.begin(), facesPP0.end());
      facesPP0.erase(unique(facesPP0.begin(), facesPP0.end()), facesPP0.end());
      sizeFacesPP0 = facesPP0.size();

      if (sizeFacesPP0 < 2) { vsize--; goto skipFace; } // le point P' ne partage pas une arete avec P (sur une face quad = sommets opposes)
      for (E_Int j = 0; j < sizeFacesPP0; j++)
      {
        face0 = facesPP0[j];
        etg = cFEp1[face0];
        etd = cFEp2[face0];
        if (etg == 0 || etd == 0)
        {
          vsize--; goto skipFace;
        }
      }

      face0 = facesPP0[0]; //1ere face
      etg = cFEp1[face0];
      etd = cFEp2[face0];
      if (etd != 0) {start = etg; next = etd;}
      else {start = etd; next = etg;}
      facesPP.push_back(face0);
      cf = 1;

      dejaVu.malloc(sizeFacesPP0);
      dejaVu.setAllValuesAtNull(); dejaVu[0] = 1;

      compt = 1;
      while (cf < sizeFacesPP0)
      {
        for (E_Int j = 1; j < sizeFacesPP0; j++)
        {
          if (dejaVu[j] == 0)
          {
            face1 = facesPP0[j];
            etg = cFEp1[face1];
            etd = cFEp2[face1];
            if (etg == next && etd != start) {
              facesPP.push_back(face1); start = etg; next = etd; dejaVu[j] = 1; compt++; goto nextface;}
            else if (etd == next && etg != start) {
              facesPP.push_back(face1); start = etd; next = etg; dejaVu[j] = 1; compt++; goto nextface;}
          }
        }
        cf += sizeFacesPP0;// on n a pas trouve de face qui permet de continuer a boucler
        nextface:;
        if (compt == sizeFacesPP0) break;
      }
      // Parcours des faces
      ptrFN[0] = facesPP.size();
      face0 = facesPP[0]; etg = cFEp1[face0]; etd = cFEp2[face0];
      start = etg; next = etd; ptrFN[1] = etg; ptrFN[2] = etd;

      for (size_t nof = 2; nof <= facesPP.size(); nof++)
      {
        face1 = facesPP[nof-1]; etg = cFEp1[face1]; etd = cFEp2[face1];
        if (etg == next && etd != start) {ptrFN[nof+1] = etd; start = etg; next = etd;}
        else if (etd == next && etg != start) {ptrFN[nof+1] = etg; start = etd; next = etg;}
      }
      ptrFN += facesPP.size()+1; sizeFNd += facesPP.size()+1;
      ptrEF[noface+1] = currentFace+1; currentFace++; noface++;
      skipFace:;
    }  // fin boucle sur les sommets P' appartenant a une face contenant P
    if (vsize > 1)
    {
      ptrEF[0] = noface; ptrEF+= noface+1; sizeEFd += noface+1; neltsfin++;
    }
  }  // boucle sur les sommets P

  cFNd.reAlloc(sizeFNd); cEFd.reAlloc(sizeEFd);
  cNGD.malloc(4+sizeFNd+sizeEFd);
  E_Int* cNGDp = cNGD.begin();
  cNGDp[0] = currentFace; cNGDp[1] = sizeFNd; cNGDp += 2;
  for (E_Int v = 0; v < sizeFNd; v++) cNGDp[v] = cFNd[v];
  cNGDp += sizeFNd;
  cNGDp[0] = neltsfin; cNGDp[1] = sizeEFd; cNGDp += 2;
  for (E_Int v = 0; v < sizeEFd; v++) cNGDp[v] = cEFd[v];
  return;
}

//=============================================================================
/* dualNGON 2D */
//=============================================================================
void K_TRANSFORM::dualNGON2D(FldArrayF& f, FldArrayI& cn, E_Int extraPoints,
                             FldArrayF& fd, FldArrayI& cNGD)
{
  // Cree les centres du maillage primal = les noeuds du maillage dual
  E_Int nfld = f.getNfld();
  E_Int nptsp = f.getSize();  // nb de pts ds le maillage primal
  E_Int* ngon = cn.getNGon(); E_Int* indPG = cn.getIndPG();
  E_Int* nface = cn.getNFace(); E_Int* indPH = cn.getIndPH();
  E_Int nfacesp = cn.getNFaces();  // nb de faces ds le maillage primal
  E_Int neltsp = cn.getNElts();
  E_Int sizeFNp = cn.getSizeNGon();
  E_Int sizeEFp = cn.getSizeNFace();

  fd.malloc(neltsp, nfld); //nb de pts ds le dual=nb d'elts dans le primal
  K_LOC::node2centerNGon(f, cn, fd);

  // calcul de la connectivite Vertex/Faces primale
  vector<vector<E_Int> > cVFp(nptsp);
  K_CONNECT::connectNG2VF(cn, cVFp);

  // calcul de la connectivite faces/elts primale
  FldArrayI cFEp; K_CONNECT::connectNG2FE(cn, cFEp);
  E_Int* cFEp1 = cFEp.begin(1);
  E_Int* cFEp2 = cFEp.begin(2);
  E_Int nfacesd = 0;  // nb de faces du dual
  E_Int nvertp, nf;

  // parcours de chq noeud primal = elt dual
  for (E_Int nov = 0; nov < nptsp; nov++)
  {
    const vector<E_Int>& faces = cVFp[nov];
    nfacesd += faces.size();
  }

  // determination des faces externes du maillage primal
  vector<E_Int> facesExt;  // faces externes du maillage primal: demarre a 0
  vector<E_Int> pointsExt;  // sommets externes du maillage primal: demarre a 0
  FldArrayI dejaVu(nptsp); dejaVu.setAllValuesAtNull();  // 1 = point Ext
  E_Int* dejaVup = dejaVu.begin();
  for (E_Int nof = 0; nof < nfacesp; nof++)
  {
    E_Int* face = cn.getFace(nof, nvertp, ngon, indPG);
    E_Int etg = cFEp1[nof]; E_Int etd = cFEp2[nof];
    if (etd == 0 || etg == 0) // face exterieure
    {
      facesExt.push_back(nof);
      for (E_Int nov = 0; nov < nvertp; nov++)
      {
        E_Int indv = face[nov]-1;
        if (dejaVup[indv] == 0)
        {
          pointsExt.push_back(indv); dejaVup[indv] = 1;
        }
      }
    }
  }
  dejaVu.malloc(0);

  if (extraPoints == 0)
  {
    // Connectivite Elts/Faces duale
    FldArrayI cEFd(nfacesd+nptsp);
    // connectivite Faces/Noeuds duale
    FldArrayI cFNd(3*nfacesd);// (2+1) x nb de faces
    E_Int* cEFp = cEFd.begin();
    E_Int* cFNp = cFNd.begin();
    FldArrayI tag(nfacesd); tag.setAllValuesAtNull();
    E_Int* tagp = tag.begin();
    // parcours de chq noeud primal=elt dual
    E_Int currentFace = 0;
    E_Int offFaces = 0; E_Int offElts = 0; E_Int ne = 0;
    for (E_Int nov = 0; nov < nptsp; nov++)
    {
      const vector<E_Int>& faces = cVFp[nov];
      E_Int nfacesv = faces.size();

      // l'element est-il correct? Pas de pt exterieur
      for (E_Int nof = 0; nof < nfacesv; nof++)
      {
        E_Int indface = faces[nof]-1;
        E_Int etg = cFEp1[indface]; E_Int etd = cFEp2[indface];
        if (etg == 0 || etd == 0) goto bad;
      }

      // element ok
      cEFp[0] = nfacesv;
      // recuperation des elts G et D de chq face
      for (E_Int nof = 0; nof < nfacesv; nof++)
      {
        E_Int indface = faces[nof]-1;
        cEFp[nof+1] = currentFace+1;
        if (tagp[indface] == 0) tagp[indface] = 1;
        E_Int etg = cFEp1[indface]; E_Int etd = cFEp2[indface];

        // connectivite faces/vertex
        cFNp[0] = 2; cFNp[1] = etg; cFNp[2] = etd;
        cFNp += 3; offFaces += 3;
        currentFace++;
      }
      cEFp += nfacesv+1; offElts += nfacesv+1; ne++;
      bad: ;
    }

    tag.malloc(0);
    cEFd.reAlloc(offElts); cFNd.reAlloc(offFaces);
    cNGD.malloc(4+cEFd.getSize()+cFNd.getSize());
    E_Int* cNGDp = cNGD.begin();
    cNGDp[0] = currentFace;
    cNGDp[1] = cFNd.getSize(); cNGDp += 2;
    for (E_Int v = 0; v < cFNd.getSize(); v++) cNGDp[v] = cFNd[v];
    cNGDp += cFNd.getSize();
    cNGDp[0] = ne;
    cNGDp[1] = cEFd.getSize(); cNGDp += 2;
    for (E_Int v = 0; v < cEFd.getSize(); v++) cNGDp[v] = cEFd[v];
    return;
  }

  if (facesExt.size() == 0)
  {
    // Connectivite Elts/Faces duale
    FldArrayI cEFd(nfacesd+nptsp);
    // connectivite Faces/Noeuds duale
    FldArrayI cFNd(3*nfacesd);// (2+1) x nb de faces
    E_Int* cEFp = cEFd.begin();
    E_Int* cFNp = cFNd.begin();
    FldArrayI tag(nfacesd); tag.setAllValuesAtNull();
    E_Int* tagp = tag.begin();
    // parcours de chq noeud primal=elt dual
    E_Int currentFace = 0;
    for (E_Int nov = 0; nov < nptsp; nov++)
    {
      const vector<E_Int>& faces = cVFp[nov];
      E_Int nfacesv = faces.size();
      cEFp[0] = nfacesv;
      // recuperation des elts G et D de chq face
      for (E_Int nof = 0; nof < nfacesv; nof++)
      {
        E_Int indface = faces[nof]-1;
        cEFp[nof+1] = currentFace+1;
        if (tagp[indface] == 0) tagp[indface] = 1;
        E_Int etg = cFEp1[indface]; E_Int etd = cFEp2[indface];

        // connectivite faces/vertex
        cFNp[0] = 2; cFNp[1] = etg; cFNp[2] = etd;
        cFNp += 3;
        currentFace++;
      }
      cEFp += nfacesv+1;
    }
    tag.malloc(0);
    cNGD.malloc(4+cEFd.getSize()+cFNd.getSize());
    E_Int* cNGDp = cNGD.begin();
    cNGDp[0] = nfacesd;
    cNGDp[1] = cFNd.getSize(); cNGDp += 2;
    for (E_Int v = 0; v < cFNd.getSize(); v++) cNGDp[v] = cFNd[v];
    cNGDp += cFNd.getSize();
    cNGDp[0] = nptsp;
    cNGDp[1] = cEFd.getSize(); cNGDp += 2;
    for (E_Int v = 0; v < cEFd.getSize(); v++) cNGDp[v] = cEFd[v];
    return;
  }

  /* Frontieres Externes */
  // Construction des elts degeneres en BAR
  E_Int nfacesExt = facesExt.size();
  E_Int nptsExt = pointsExt.size();
  FldArrayI cFNp2(nfacesExt*3*(2+1)+nptsExt*1*(2+1));// on rajoute 3 faces par face exterieure et 1 faces par sommet exterieur
  FldArrayI cEFp2(nfacesExt*(4+1)+nptsExt*(3+1));// on rajoute 1 elt QUAD degenere par face ext et 1 elt TRI degenere par point ext
  E_Int* ptrFN = cFNp2.begin();
  E_Int* ptrEF = cEFp2.begin();
  FldArrayI facesopp(nptsp,2); facesopp.setAllValuesAt(-1);
  E_Int* foppg = facesopp.begin(1);// pour chq sommet externe on recupere les faces gauche et droite creees par les faces externes
  E_Int* foppd = facesopp.begin(2);
  E_Int cf = 0;//compteur sur les faces
  E_Int ce = 0; // compteur sur les elts

  for (E_Int nof = 0; nof < nfacesExt; nof++)
  {
    E_Int indface = facesExt[nof];
    /*1. construction de l'elt QUAD degenere en BAR a partir de la face externe */
    E_Int* face = cn.getFace(indface, nvertp, ngon, indPG);  // connectivite Faces/Noeuds primale

    // creation de la face opposee identique a la face indface
    ptrFN[0] = nvertp;
    for (E_Int nov = 1; nov <= nvertp; nov++) ptrFN[nov] = face[nov-1];
    ptrFN += nvertp+1;
    // creation des faces laterales
    for (E_Int nov = 1; nov <= nvertp; nov++)
    {
      ptrFN[0] = 2;
      E_Int indv = face[nov-1]; ptrFN[1] = indv; ptrFN[2] = indv;
      ptrFN += 3;
      if (foppg[indv-1] == -1) foppg[indv-1] = nfacesp+cf+nov+1;
      else if (foppd[indv-1] == -1) foppd[indv-1] = nfacesp+cf+nov+1;
    }
    // creation de l elt QUAD contenant ces faces degenerees
    ptrEF[0] = 4;
    ptrEF[1] = indface+1;
    ptrEF[2] = nfacesp+cf+1;//cf demarre a 0
    ptrEF[3] = nfacesp+cf+2;
    ptrEF[4] = nfacesp+cf+3;
    ptrEF += 5;
    cf += nvertp+1; // 3 faces ajoutees
    ce++;
    // passage a la facette externe suivante
  }

  // Creation des elts degeneres en NODE a partir des sommets externes
  for (E_Int nov = 0; nov < nptsExt; nov++)
  {
    E_Int indv = pointsExt[nov];
    // Creation de la troisieme face
    ptrFN[0] = 2;
    ptrFN[1] = indv+1;
    ptrFN[2] = indv+1;
    ptrFN += 3;

    // Creation de l elt associe
    ptrEF[0] = 3;
    ptrEF[1] = foppg[indv];
    ptrEF[2] = foppd[indv];
    ptrEF[3] = nfacesp+cf+1;
    ptrEF += 4;

    ce++; cf++;  // 1 face ajoutee
  }
  facesopp.malloc(0);

  // Creation de la nouvelle connectivite avec degenerescences
  sizeFNp += cFNp2.getSize(); sizeEFp += cEFp2.getSize();
  FldArrayI cNGon(4+sizeFNp+sizeEFp);
  cNGon.setNGonType(1);
  cNGon[0] = nfacesp + cf; cNGon[1] = sizeFNp;
  cNGon[2+sizeFNp] = neltsp + ce; cNGon[3+sizeFNp] = sizeEFp;
  E_Int* ptr1 = cNGon.begin()+2;
  for (E_Int i = 0; i < nfacesp; i++)
  {
    E_Int* face = cn.getFace(i, nvertp, ngon, indPG);
    ptr1[0] = nvertp; ptr1++;
    for (E_Int nov = 0; nov < nvertp; nov++) ptr1[nov] = face[nov];
    ptr1 += nvertp;
  }
  for (E_Int i = 0; i < cFNp2.getSize(); i++) { ptr1[0] = cFNp2[i]; ptr1++; }

  ptr1 += 2;  // cEF initial
  for (E_Int i = 0; i < neltsp; i++)
  {
    E_Int* elt = cn.getElt(i, nf, nface, indPH);
    ptr1[0] = nf; ptr1++;
    for (E_Int nof = 0; nof < nf; nof++) ptr1[nof] = elt[nof];
    ptr1 += nf;
  }
  for (E_Int i = 0; i < cEFp2.getSize(); i++) { ptr1[0] = cEFp2[i]; ptr1++; }

  fd.malloc(neltsp + ce, nfld);  // nb de pts ds le dual = nb d elts dans le primal
  K_LOC::node2centerNGon(f, cNGon, fd);

  /*------------------------------------------------*/
  /* On travaille sur la nouvelle connectivite NGON */
  /*------------------------------------------------*/
  // calcul de la connectivite Vertex/Faces primale
  for (size_t nov = 0; nov < cVFp.size(); nov++) cVFp[nov].clear();
  K_CONNECT::connectNG2VF(cNGon, cVFp);
  // calcul de la connectivite faces/elts primale
  K_CONNECT::connectNG2FE(cNGon, cFEp);
  cFEp1 = cFEp.begin(1);
  cFEp2 = cFEp.begin(2);
  nfacesd = 0; // nb de faces du dual
  // parcours de chq noeud primal = elt dual
  for (E_Int nov = 0; nov < nptsp; nov++) nfacesd +=  cVFp[nov].size();

  // Connectivite Elts/Faces duale
  FldArrayI cEFd(nfacesd+nptsp);
  E_Int neltsd = 0;
  // connectivite Faces/Noeuds duale
  FldArrayI cFNd(3*nfacesd); // (2+1) x nb de faces
  FldArrayI tag(nfacesd); tag.setAllValuesAtNull();
  E_Int* cEFp = cEFd.begin();
  E_Int* cFNp = cFNd.begin();
  E_Int sizeEF2 = 0; E_Int sizeFN2 = 0;
  // parcours de chq noeud primal = elt dual
  E_Int currentFace = 0;
  for (E_Int nov = 0; nov < nptsp; nov++)
  {
    vector<E_Int>& faces = cVFp[nov];
    E_Int nfacesv = faces.size();
    cEFp[0] = nfacesv;
    // recuperation des elts G et D de chq face
    for (E_Int nof = 0; nof < nfacesv; nof++)
    {
      E_Int indface = faces[nof]-1;
      if (tag[indface] == 0) tag[indface] = 1;
      E_Int etg = cFEp1[indface]; E_Int etd = cFEp2[indface];
      cEFp[nof+1] = currentFace+1;
      if (etd == 0) etd = etg;
      // connectivite faces/vertex
      cFNp[0] = 2; cFNp[1] = etg; cFNp[2] = etd;
      cFNp += 3;
      currentFace++;
      sizeFN2 += 3;
    }
    sizeEF2 += nfacesv+1;
    cEFp += nfacesv+1; neltsd++;
  }
  cEFd.resize(sizeEF2); cFNd.resize(sizeFN2);
  //cleanings
  tag.malloc(0); cFEp.malloc(0);
  for (size_t nov = 0; nov < cVFp.size(); nov++) cVFp[nov].clear();
  cVFp.clear();
  // sortie
  cNGD.malloc(4+sizeFN2+sizeEF2);
  E_Int* cNGDp = cNGD.begin();
  cNGDp[0] = currentFace;
  cNGDp[1] = cFNd.getSize(); cNGDp += 2;
  for (E_Int v = 0; v < cFNd.getSize(); v++) cNGDp[v] = cFNd[v];
  cNGDp += cFNd.getSize();
  cNGDp[0] = neltsd;
  cNGDp[1] = cEFd.getSize(); cNGDp += 2;
  for (E_Int v = 0; v < cEFd.getSize(); v++) cNGDp[v] = cEFd[v];
  return;
}

//=============================================================================
/* dualNGON 1D */
//=============================================================================
void K_TRANSFORM::dualNGON1D(FldArrayF& f, FldArrayI& cn, E_Int extraPoints,
                             FldArrayF& fd, FldArrayI& cNGD)
{
  E_Int nfld = f.getNfld();
  E_Int nptsp = f.getSize();
  E_Int neltsp = cn.getNElts();
  E_Int nptsd = neltsp; // nb de pts ds le dual = nb d'elts dans le primal
  // connectivite 1D fermee: neltsd = nptsd
  // connectivite 1D ouverte: neltsd = nptsd-1
  E_Int neltsd = nptsd - (nptsp - neltsp);
  fd.malloc(nptsd, nfld);
  K_LOC::node2centerNGon(f, cn, fd);

  // Connectivite duale
  E_Int shift = 1;
  E_Int sizeFNp = (1+shift)*nptsp;
  E_Int sizeFNd = (1+shift)*nptsd;
  E_Int sizeEFd = (2+shift)*neltsd;
  E_Int sizeConn = 4 + sizeFNd + sizeEFd;
  cNGD.malloc(sizeConn);

  cNGD[0] = nptsd; cNGD[1] = sizeFNd;
  cNGD[2+sizeFNd] = neltsd; cNGD[2+sizeFNd+1] = sizeEFd;

  #pragma omp parallel
  {
    E_Int ind, indp, offset, offsetp, fac;
    offset = 2; fac = (1+shift);
    #pragma omp for
    for(E_Int i = 0; i < nptsd; i++)
    {
      ind = offset + fac*i;
      cNGD[ind] = 1; cNGD[ind+shift] = cn[ind+shift];
    }

    offset = 4 + sizeFNd; offsetp = 4 + sizeFNp; fac = (2+shift);
    #pragma omp for
    for(E_Int i = 0; i < neltsd; i++)
    {
      ind = offset + fac*i; indp = offsetp + fac*i;
      cNGD[ind] = 2;
      cNGD[ind+shift] = cn[indp+shift];
      cNGD[ind+shift+1] = cn[indp+shift+1];
    }
  }

  if (nptsd == neltsd)
    cNGD[sizeConn-1] = cNGD[4+sizeFNd+shift]; // closed contour

  return;
}

//=============================================================================
/* Creation du maillage primal en ajoutant des elts degeneres a partir des
   faces exterieures et des aretes partageant 2 faces exterieures
   Le tableau de coordonnees reste identique, avec nptsp points
   cNGD est alloue ici, si pas de faces externes, de taille nulle */
//=============================================================================
E_Int K_TRANSFORM::createDegeneratedPrimalMesh3D(
  FldArrayF& fNG, FldArrayI& cNG,
  FldArrayF& fNGD, FldArrayI& cNGD)
{
  E_Int nptsp = fNG.getSize(); E_Int nfld = fNG.getNfld();
  cNGD.malloc(0);  fNGD.malloc(0);
  E_Int* ngon = cNG.getNGon(); E_Int* indPG = cNG.getIndPG();
  E_Int* nface = cNG.getNFace(); E_Int* indPH = cNG.getIndPH();
  E_Int nfacesp = cNG.getNFaces();  // nb de faces ds le maillage primal
  E_Int neltsp = cNG.getNElts();
  E_Int sizeFNp = cNG.getSizeNGon();
  E_Int sizeEFp = cNG.getSizeNFace();

  FldArrayI cFEp; K_CONNECT::connectNG2FE(cNG, cFEp);
  E_Int* cFEp1 = cFEp.begin(1);
  E_Int* cFEp2 = cFEp.begin(2);

  // determination des faces externes du maillage primal
  vector<E_Int> facesExt;  // demarre a 0
  E_Int sizeFN2 = 0; E_Int sizeEF2 = 0;
  // dimensionnements
  E_Int nptsd = 0;
  E_Int nvertp, nf;

  for (E_Int nof = 0; nof < nfacesp; nof++)
  {
    cNG.getFace(nof, nvertp, ngon, indPG);
    E_Int etg = cFEp1[nof]; E_Int etd = cFEp2[nof];
    if (etd == 0 || etg == 0) // face exterieure
    {
      facesExt.push_back(nof);
      sizeFN2 += nvertp+1; // face courante dupliquee avec nvertp sommets
      sizeFN2 += (4+1)*nvertp; // autant de faces laterales quad que d'aretes dans la face nof
      sizeEF2 += (nvertp+2+1); // face courante + face dupliquee + nvertp faces laterales + 1 pour la taille
      nptsd += nvertp;
    }
  }
  E_Int nfacesExt = facesExt.size();
  if (nfacesExt == 0) return 0;
  fNGD.malloc(nptsd+nptsp, nfld);
  for (E_Int eq = 1; eq <= nfld; eq++)
  {
    E_Float* fp = fNG.begin(eq);
    E_Float* fd = fNGD.begin(eq);
    for (E_Int ind = 0; ind < nptsp; ind++) fd[ind] = fp[ind];
  }
  FldArrayI cFN2(sizeFN2); FldArrayI cEF2(sizeEF2);
  E_Int* ptrFN2 = cFN2.begin();
  E_Int* ptrEF2 = cEF2.begin();
  E_Int sizeFN3 = 0; E_Int sizeEF3 = 0;
  E_Int ce = 0; // compteur sur les elts
  // Creation des elements degeneres a partir des frontieres externes
  E_Int compt = nptsp;
  FldArrayI indirExt(nptsp); indirExt.setAllValuesAt(-1);  // demarre a 0
  // Recherche des pts exterieurs
  E_Int nptsExt = 0;

  for (E_Int nof = 0; nof < nfacesExt; nof++)
  {
    E_Int indface = facesExt[nof];
    E_Int* face = cNG.getFace(indface, nvertp, ngon, indPG);
    for (E_Int nov = 0; nov < nvertp; nov++)
    {
      E_Int indv = face[nov]-1;
      if (indirExt[indv] == -1) { indirExt[indv] = nptsExt; nptsExt++;}
    }
  }

  vector<vector<E_Int> > facesLaterales(nptsExt);
  FldArrayI indir(nptsExt); indir.setAllValuesAtNull();
  E_Int nfacesTot = nfacesp;
  for (E_Int nof = 0; nof < nfacesExt; nof++)
  {
    E_Int indface = facesExt[nof];
    E_Int* face = cNG.getFace(indface, nvertp, ngon, indPG);
    ptrFN2[0] = nvertp;
    ptrEF2[0] = nvertp+2; // l elt degenere contient nvertp+2 faces
    ptrEF2[1] = indface+1;

    // 1.creation de la face opposee a indface
    for (E_Int nov = 1; nov <= nvertp; nov++)
    {
      E_Int indv = face[nov-1]-1;
      E_Int indExt = indirExt[indv]; // numero de noeud externe associe
      if (indir[indExt] == 0) // on cree le nouveau pt degenere
      {
        for (E_Int eq = 1; eq <= nfld; eq++) fNGD(compt,eq) = fNG(indv,eq);
        indir[indExt] = compt+1;  // demarre a 1
        compt++;
      }
      ptrFN2[nov] = indir[indExt];
    }// boucle sur les sommets
    ptrFN2 += nvertp+1; sizeFN3 += nvertp+1;
    ptrEF2[2] = nfacesTot+1; nfacesTot++;

    // 2.creation des faces laterales
    for (E_Int nov = 1; nov <= nvertp; nov++)
    {
      ptrFN2[0] = 4; // 4 sommets par face laterale
      E_Int indv1 = face[nov-1];
      E_Int indv2 = face[nov];
      if (nov == nvertp) indv2 = face[0];
      E_Int indExt1 = indirExt[indv1-1];  // demarre a 0
      E_Int indExt2 = indirExt[indv2-1];

      vector<E_Int>& facesVExt1 = facesLaterales[indExt1];
      vector<E_Int>& facesVExt2 = facesLaterales[indExt2];
      // recherche de la face laterale deja creee eventuellement associe a l'edge indExt1-indExt2
      for (size_t i1 = 0; i1 < facesVExt1.size(); i1++)
        for (size_t i2 = 0; i2 < facesVExt2.size(); i2++)
        {
          if (facesVExt1[i1] == facesVExt2[i2])
          {
            ptrEF2[2+nov] = facesVExt1[i1];  // numero de la face qui existe deja
            goto skipv;
          }
        }
      ptrFN2[1] = indv1;
      ptrFN2[2] = indv2;
      ptrFN2[3] = indir[indExt2];
      ptrFN2[4] = indir[indExt1];
      ptrFN2 += 5; sizeFN3 += 5;
      ptrEF2[2+nov] = nfacesTot+1;

      facesVExt1.push_back(nfacesTot+1); facesVExt2.push_back(nfacesTot+1);
      nfacesTot++;
      skipv:;
    }

    ptrEF2 += nvertp+3; ce++; sizeEF3 += nvertp+3;
  }
  fNGD.reAllocMat(compt,nfld); indir.malloc(0); indirExt.malloc(0);
  for (size_t nof = 0; nof < facesLaterales.size(); nof++)
    facesLaterales[nof].clear();
  facesLaterales.clear();
  sizeFN2 = sizeFN3; sizeEF2 = sizeEF3;

  // Creation de la nouvelle connectivite avec degenerescences
  sizeFNp += sizeFN2; sizeEFp += sizeEF2;
  cNGD.malloc(4+sizeFNp+sizeEFp);
  cNGD.setNGonType(1);

  cNGD[0] = nfacesTot; cNGD[1] = sizeFNp;
  cNGD[2+sizeFNp] = neltsp + ce; cNGD[3+sizeFNp] = sizeEFp;

  ptrFN2 = cFN2.begin(); ptrEF2 = cEF2.begin();

  E_Int* ptr1 = cNGD.begin()+2;
  for (E_Int i = 0; i < nfacesp; i++)
  {
    E_Int* face = cNG.getFace(i, nvertp, ngon, indPG);
    ptr1[0] = nvertp; ptr1++;
    for (E_Int nov = 0; nov < nvertp; nov++) ptr1[nov] = face[nov];
    ptr1 += nvertp;
  }
  for (E_Int i = 0; i < sizeFN2; i++) { ptr1[0] = ptrFN2[i]; ptr1++; }

  ptr1 += 2;  // cEF initial
  for (E_Int i = 0; i < neltsp; i++)
  {
    E_Int* elt = cNG.getElt(i, nf, nface, indPH);
    ptr1[0] = nf; ptr1++;
    for (E_Int nof = 0; nof < nf; nof++) ptr1[nof] = elt[nof];
    ptr1 += nf;
  }
  for (E_Int i = 0; i < sizeEF2; i++) { ptr1[0] = ptrEF2[i]; ptr1++; }

  return nfacesExt;
}
