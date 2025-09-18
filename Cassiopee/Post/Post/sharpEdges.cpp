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

#include "post.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Detect sharp edges on a surface defined by BAR, TRI or QUAD */
//=============================================================================
PyObject* K_POST::sharpEdges(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float dirVect[3];
  dirVect[0] = 0; dirVect[1] = 0; dirVect[2] = 1; 
  E_Float alpref;
  if (!PYPARSETUPLE_(args, O_ R_, &array, &alpref))
  {
    return NULL;
  }
  /*-----------------------------------------------*/
  /* Extraction des donnees du maillage surfacique */ 
  /*-----------------------------------------------*/
  char* varString0; char* eltType0;
  FldArrayF* f; FldArrayI* cn;
  E_Int nil, njl, nkl;
  E_Int res = K_ARRAY::getFromArray3(array, varString0, f, nil, njl, nkl, 
                                     cn, eltType0);
  
  if (res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "sharpEdges: array must be unstructured.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  E_Int type = 0;
  if (strcmp(eltType0, "BAR") == 0) type = 2;
  else if (strcmp(eltType0, "TRI") == 0) type = 3;
  else if (strcmp(eltType0, "QUAD") == 0) type = 4;
  else if (strcmp(eltType0, "NGON") == 0) type = 5;
  if (type == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "sharpEdges: array must be BAR, TRI or QUAD.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }  
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString0);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString0);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString0);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "sharpEdges: array must contain coordinates.");
    RELEASESHAREDU(array, f, cn); return NULL;  
  }
  posx++; posy++; posz++;
  
  E_Int nelts = cn->getSize();
  E_Int npts = f->getSize();
  E_Int nfld = f->getNfld();
  E_Int api = f->getApi();
  // pointers on field f
  vector<E_Float*> fp(nfld);
  for (E_Int p = 0; p < nfld; p++) fp[p] = f->begin(p+1);
  // coordonnees
  E_Float* x = f->begin(posx);
  E_Float* y = f->begin(posy);
  E_Float* z = f->begin(posz);

  // pointers on field fe (field associated with new connectivity cne)
  FldArrayF* fe = new FldArrayF(npts,nfld);
  vector<E_Float*> fep(nfld);
  for (E_Int p = 0; p < nfld; p++) fep[p] = fe->begin(p+1);

  E_Int indA1, indB1, indC1, indD1, indA2, indB2, indC2, indD2;
  vector<E_Float> ptA1(nfld); vector<E_Float> ptB1(nfld);
  vector<E_Float> ptC1(nfld); vector<E_Float> ptD1(nfld);
  vector<E_Float> ptA2(nfld); vector<E_Float> ptB2(nfld);
  vector<E_Float> ptC2(nfld); vector<E_Float> ptD2(nfld);
  vector< vector<E_Int> > cEEN(nelts);
  K_CONNECT::connectEV2EENbrs(eltType0, npts, *cn, cEEN); 

  E_Float alphamin = 180.-alpref; E_Float alphamax = 180.+alpref;
  E_Int nop = 0; E_Int noe = 0;
  E_Int found1 = 0; E_Int found2 = 1;
  if (type == 2) // BAR
  {
    E_Int* cn1 = cn->begin(1); E_Int* cn2 = cn->begin(2);
    FldArrayI* cne = new FldArrayI(0,1);// NODE

    for (E_Int et1 = 0; et1 < nelts; et1++)
    {
      indA1 = cn1[et1]-1; indB1 = cn2[et1]-1; 
      ptA1[0] = x[indA1];ptB1[0] = x[indB1];
      ptA1[1] = y[indA1];ptB1[1] = y[indB1];
      ptA1[2] = z[indA1];ptB1[2] = z[indB1];
      vector<E_Int>& eltsVoisins = cEEN[et1]; E_Int nvoisins = eltsVoisins.size();
      for (E_Int noet2 = 0; noet2 < nvoisins; noet2++)
      {
        E_Int et2 = eltsVoisins[noet2];
        indA2 = cn1[et2]-1; indB2 = cn2[et2]-1; 
        ptA2[0] = x[indA2]; ptB2[0] = x[indB2];
        ptA2[1] = y[indA2]; ptB2[1] = y[indB2];
        ptA2[2] = z[indA2]; ptB2[2] = z[indB2];
        E_Float alpha = K_COMPGEOM::getAlphaAngleBetweenBars(&ptA1[0], &ptB1[0], &ptA2[0], &ptB2[0], dirVect);
        if ((alpha < alphamin || alpha > alphamax) && alpha != -1000.) 
        {
          found1 = 0; found2 = 0;
          if (indA1 == indA2 || indA1 == indB2) 
          {for (E_Int p=0;p<nfld;p++) fep[p][nop] = fp[p][indA1];}
          else {for (E_Int p=0;p<nfld;p++) fep[p][nop] = fp[p][indB1];}
          nop++; noe++;
          if ( nop+2 > fe->getSize() ) 
          { fe->reAllocMat(2*fe->getSize(), nfld); for (E_Int p=0;p<nfld;p++) fep[p] = fe->begin(p+1);}
        }
      }
    }
    RELEASESHAREDU(array, f, cn);
    fe->reAllocMat(nop, nfld); 
    PyObject* tpl;
    if (fe->getSize() == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "sharpEdges: sharp edges set is empty.");
      tpl = NULL;
    }
    else tpl = K_ARRAY::buildArray3(*fe, varString0, *cne, "NODE", api);
    delete fe; delete cne;
    return tpl;
  }
  else if ( type == 3 ) // TRI 
  {
    E_Int* cn1 = cn->begin(1); E_Int* cn2 = cn->begin(2); E_Int* cn3 = cn->begin(3);
    FldArrayI* cne = new FldArrayI(nelts, 2);// BAR
    E_Int* cne1 = cne->begin(1); E_Int* cne2 = cne->begin(2);
    
    for (E_Int et1 = 0; et1 < nelts; et1++)
    {
      indA1 = cn1[et1]-1; indB1 = cn2[et1]-1; indC1 = cn3[et1]-1;
      ptA1[0] = x[indA1]; ptB1[0] = x[indB1]; ptC1[0] = x[indC1];
      ptA1[1] = y[indA1]; ptB1[1] = y[indB1]; ptC1[1] = y[indC1];
      ptA1[2] = z[indA1]; ptB1[2] = z[indB1]; ptC1[2] = z[indC1];
      for (E_Int p=0;p<nfld;p++){ptA1[p] = fp[p][indA1];ptB1[p] = fp[p][indB1];ptC1[p] = fp[p][indC1];}
      vector<E_Int>& eltsVoisins = cEEN[et1]; E_Int nvoisins = eltsVoisins.size();
      for (E_Int noet2 = 0; noet2 < nvoisins; noet2++)
      {
        E_Int et2 = eltsVoisins[noet2];
        indA2 = cn1[et2]-1; indB2 = cn2[et2]-1; indC2 = cn3[et2]-1;
        ptA2[0] = x[indA2];ptB2[0] = x[indB2]; ptC2[0] = x[indC2];
        ptA2[1] = y[indA2];ptB2[1] = y[indB2]; ptC2[1] = y[indC2];
        ptA2[2] = z[indA2];ptB2[2] = z[indB2]; ptC2[2] = z[indC2];
        E_Float alpha = K_COMPGEOM::getAlphaAngleBetweenTriangles(&ptA1[0], &ptB1[0], &ptC1[0], &ptA2[0], &ptB2[0], &ptC2[0]);

        if ((alpha < alphamin || alpha > alphamax ) && alpha != -1000.) 
        {
          found1 = 0; found2 = 0;
          if (indA1 == indA2 || indA1 == indB2 || indA1 == indC2) 
          {
            found1 = 1;
          }
          if (indB1 == indA2 || indB1 == indB2 || indB1 == indC2) 
          {
            if (found1 == 0) found1 = 2;
            else found2 = 2;
          }
          if (indC1 == indA2 || indC1 == indB2 || indC1 == indC2) 
          {
            found2 = 3;
          }
          if ( found1 == 1 ) {for (E_Int p=0;p<nfld;p++) fep[p][nop] = fp[p][indA1];}
          else {for (E_Int p=0;p<nfld;p++) fep[p][nop] = fp[p][indB1];}
          if ( found2 == 2 ) {for (E_Int p=0;p<nfld;p++) fep[p][nop+1] = fp[p][indB1];}
          else {for (E_Int p=0;p<nfld;p++) fep[p][nop+1] = fp[p][indC1];}
          cne1[noe] = nop+1; cne2[noe] = nop+2;
          nop = nop+2; noe++;
          if ( nop+2 >fe->getSize() ) 
          { fe->reAllocMat(2*fe->getSize(),nfld); for (E_Int p=0;p<nfld;p++) fep[p] = fe->begin(p+1);}
          if ( noe+2 > cne->getSize() ) {cne->reAllocMat(2*cne->getSize(),2); cne1 = cne->begin(1); cne2 = cne->begin(2);}
        }
      }
      
    }
    RELEASESHAREDU(array, f, cn);
    fe->reAllocMat(nop, nfld); cne->reAllocMat(noe, 2);
    PyObject* tpl;
    if (fe->getSize() == 0 || cne->getSize() == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "sharpEdges: sharp edges set is empty.");
      tpl = NULL;
    }
    else tpl = K_ARRAY::buildArray3(*fe, varString0, *cne, "BAR", api);
    delete fe; delete cne;
    return tpl;
  }
  else if (type == 4) //QUAD
  {
    E_Int* cn1 = cn->begin(1); E_Int* cn2 = cn->begin(2); E_Int* cn3 = cn->begin(3); E_Int* cn4 = cn->begin(4);
    FldArrayI* cne = new FldArrayI(nelts,2);// BAR
    E_Int* cne1 = cne->begin(1); E_Int* cne2 = cne->begin(2);

    for (E_Int et1 = 0; et1 < nelts; et1++)
    {
      indA1 = cn1[et1]-1; indB1 = cn2[et1]-1; indC1 = cn3[et1]-1; indD1 = cn4[et1]-1;
      ptA1[0] = x[indA1]; ptB1[0] = x[indB1]; ptC1[0] = x[indC1]; ptD1[0] = x[indD1];
      ptA1[1] = y[indA1]; ptB1[1] = y[indB1]; ptC1[1] = y[indC1]; ptD1[1] = y[indD1];
      ptA1[2] = z[indA1]; ptB1[2] = z[indB1]; ptC1[2] = z[indC1]; ptD1[2] = z[indD1];
      vector<E_Int>& eltsVoisins = cEEN[et1]; E_Int nvoisins = eltsVoisins.size();
      for (E_Int noet2 = 0; noet2 < nvoisins; noet2++)
      {
        E_Int et2 = eltsVoisins[noet2];
        indA2 = cn1[et2]-1; indB2 = cn2[et2]-1; indC2 = cn3[et2]-1; indD2 = cn4[et2]-1;
        ptA2[0] = x[indA2];ptB2[0] = x[indB2]; ptC2[0] = x[indC2]; ptD2[0] = x[indD2];
        ptA2[1] = y[indA2];ptB2[1] = y[indB2]; ptC2[1] = y[indC2]; ptD2[1] = y[indD2];
        ptA2[2] = z[indA2];ptB2[2] = z[indB2]; ptC2[2] = z[indC2]; ptD2[2] = z[indD2];
        for (E_Int p=0;p<nfld;p++){ptA2[p] = fp[p][indA2];ptB2[p] = fp[p][indB2];ptC2[p] = fp[p][indC2];ptD2[p] = fp[p][indD2];}
        E_Float alpha = 
          K_COMPGEOM::getAlphaAngleBetweenQuads(&ptA1[0], &ptB1[0], &ptC1[0], &ptD1[0], 
                                                &ptA2[0], &ptB2[0], &ptC2[0], &ptD2[0]);
        if ( (alpha < alphamin || alpha > alphamax ) && alpha != -1000.) 
        {
          found1 = 0; found2 = 0;
          if ( indA1 == indA2 || indA1 == indB2 || indA1 == indC2 || indA1 == indD2) 
          {found1 = 1;}
          if ( indB1 == indA2 || indB1 == indB2 || indB1 == indC2 || indB1 == indD2) 
          {
            if ( found1 == 0) found1 = 2;
            else found2 = 2;
          }
          if ( indC1 == indA2 || indC1 == indB2 || indC1 == indC2 || indC1 == indD2) 
          {
            if ( found1 == 0 ) found1 = 3; 
            else found2 = 3;
          }
          if ( indD1 == indA2 || indD1 == indB2 || indD1 == indC2 || indD1 == indD2) 
          {found2 = 4;}

          if ( found1 == 1 ) {for (E_Int p=0;p<nfld;p++) fep[p][nop] = fp[p][indA1];}
          else if ( found1 == 2 ) {for (E_Int p=0;p<nfld;p++) fep[p][nop] = fp[p][indB1];}
          else {for (E_Int p=0;p<nfld;p++) fep[p][nop] = fp[p][indC1];}

          if (found2 == 2) {for (E_Int p=0;p<nfld;p++) fep[p][nop+1] = fp[p][indB1];}
          else if ( found2 == 3 ) {for (E_Int p=0;p<nfld;p++) fep[p][nop+1] = fp[p][indC1];}
          else {for (E_Int p=0;p<nfld;p++) fep[p][nop+1] = fp[p][indD1];}

          cne1[noe] = nop+1; cne2[noe] = nop+2; nop = nop+2; noe++;
          if ( nop+2 > fe->getSize() ) {fe->reAllocMat(2*fe->getSize(),nfld);for (E_Int p=0;p<nfld;p++) fep[p] = fe->begin(p+1);}
          if ( noe+2 > cne->getSize() ) {cne->reAllocMat(2*cne->getSize(),2);  cne1 = cne->begin(1); cne2 = cne->begin(2);}       
        }
      }
    }
    RELEASESHAREDU(array, f, cn);
    fe->reAllocMat(nop, nfld); cne->reAllocMat(noe, 2);
    PyObject* tpl;
    if (fe->getSize() == 0 || cne->getSize() == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "sharpEdges: sharp edges set is empty.");
      tpl = NULL;
    }
    else tpl = K_ARRAY::buildArray3(*fe, varString0, *cne, "BAR", api);
    delete fe; delete cne;
    return tpl;  
  }
  else // NGON
  {
    E_Int* cnp = cn->begin();       // connectivite NGON
    E_Int sizeFN = cnp[1];         // taille de la connectivite face/noeuds
    nelts = cnp[sizeFN+2];         // nombre d elements
    // connectivite element/faces
    E_Int* ptrEF = cnp+sizeFN+4;
    // connectivite element/elements voisins cEEN
    FldArrayI cFE;
    K_CONNECT::connectNG2FE(*cn, cFE);
    vector< vector<E_Int> > cEEN(nelts);
    K_CONNECT::connectFE2EENbrs(cFE, cEEN);
    
    // tableau des positions des elts dans la connectivite
    FldArrayI pos; K_CONNECT::getPosElts(*cn, pos);
    // tableau des positions des faces dans la connectivite
    FldArrayI posFaces; K_CONNECT::getPosFaces(*cn, posFaces);
    FldArrayI dimElts; // tableau donnant pour chaque element sa dimension
    K_CONNECT::getDimElts(*cn, dimElts);
    E_Int* dimEltsp = dimElts.begin();
    // taille du champ f
    E_Int nfld = f->getNfld();
    E_Int api = f->getApi();
    // pointeur sur le champ f
    vector<E_Float*> fp(nfld);
    for (E_Int p = 0; p < nfld; p++) fp[p] = f->begin(p+1);

    // noeuds des elements elt1 et elt2 
    vector<E_Float*> pts1; vector<E_Float*> pts2;
    // nombre de faces de elt1 et elt2
    E_Int nbfaces1;
    // coordonnees de noeuds de fa1 et fa2
    E_Float x1, y1, z1, x2, y2, z2;
    // nof: indice d'une face/d'un noeud de la nouvelle connectivite
    E_Int nof;
    // vecteur d'indices des noeuds des element elt1 et elt2
    vector<E_Int> indices1; vector<E_Int> indices2;
    // vecteur des points coincidents entre la face fa1 et un element voisin
    vector<E_Int> matchingPoints;
    E_Int next1=0; // nbre de faces pour la nouvelle connectivite 
    E_Int next2=0; // nbre d'elements pour la nouvelle connectivite
    E_Int sizeco1 = 0; // taille de la nouvelle connectivite face/noeuds
    E_Int sizeco2 = 0; // taille de la nouvelle connectivite element/faces
    // autres variables locales de travail ...
    E_Int ind, indn, nvert, numface, posfa1, nbnodes1, indpt1;      
    // objet de sortie
    PyObject* tpl;

    // Verification de la dimensions des elements
    E_Int dim=0; E_Int error = 0;
    for (E_Int elt1 = 0; elt1 < nelts; elt1++)
    {
      dim = dimEltsp[elt1];
      if (dim == 3) {error = 1; break;}
    }
    if (error == 1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "sharpEdges: not implemented for 3D NGON elements.");
      return NULL;
    }

    // On suppose que la dimension des elements est homogene. 
    // La dimension est donc basee sur celle du dernier element
    if (dim == 1)
    {
      // estimation surevaluee des dimensions pour allouer les nouvelles connectivites newcnFN et newcnEF
      for (E_Int elt1 = 0; elt1 < nelts; elt1++)
      {
        for (size_t ev = 0; ev < cEEN[elt1].size(); ev++)
        {
          sizeco1 += 1;
        }
      }

      FldArrayI* cne = new FldArrayI(sizeco1);
      E_Int* cnep = cne->begin();

      // parcours des elements
      for (E_Int elt1 = 0; elt1 < nelts; elt1++)
      {
        K_CONNECT::getVertexIndices(cn->begin(), posFaces.begin(), pos[elt1], indices1);
        indA1 = indices1[0]-1; indB1 = indices1[1]-1;
        ptA1[0] = x[indA1];ptB1[0] = x[indB1];
        ptA1[1] = y[indA1];ptB1[1] = y[indB1];
        ptA1[2] = z[indA1];ptB1[2] = z[indB1];
        // parcours des elements voisins
        for (unsigned int ev = 0; ev < cEEN[elt1].size(); ev++)
        {
          // elt2 voisin
          E_Int elt2 = cEEN[elt1][ev];
          // recuperation des noeuds de elt2
          K_CONNECT::getVertexIndices(cn->begin(), posFaces.begin(), pos[elt2], indices2);
          indA2 = indices2[0]-1; indB2 = indices2[1]-1;
          ptA2[0] = x[indA2];ptB2[0] = x[indB2];
          ptA2[1] = y[indA2];ptB2[1] = y[indB2];
          ptA2[2] = z[indA2];ptB2[2] = z[indB2];
          E_Float alpha = K_COMPGEOM::getAlphaAngleBetweenBars(&ptA1[0], &ptB1[0], &ptA2[0], &ptB2[0], dirVect);
          if ((alpha < alphamin || alpha > alphamax ) && alpha != -1000.) 
         {
           found1 = 0; found2 = 0;
           if (indA1 == indA2 || indA1 == indB2) 
           {for (E_Int p=0;p<nfld;p++) fep[p][nop] = fp[p][indA1];}
           else {for (E_Int p=0;p<nfld;p++) fep[p][nop] = fp[p][indB1];}
           cnep[noe] = nop+1;  nop++; noe++;
           if (nop+2 > fe->getSize()) 
           { fe->reAllocMat(2*fe->getSize(), nfld); for (E_Int p=0;p<nfld;p++) fep[p] = fe->begin(p+1);}
         }
        }  // fin boucle elt2   
      } // fin boucle elt1

      RELEASESHAREDU(array, f, cn);
      fe->reAllocMat(nop, nfld); cne->reAllocMat(noe, 1);
      if (fe->getSize() == 0 || cne->getSize() == 0)
      {
        PyErr_SetString(PyExc_TypeError,
                        "sharpEdges: sharp edges set is empty.");
        tpl = NULL;
      }
      else tpl = K_ARRAY::buildArray3(*fe, varString0, *cne, "NODE", api);
      delete fe; delete cne;
    }
    else // dim == 3
    {
      // estimation surevaluee des dimensions pour allouer les nouvelles connectivites newcnFN et newcnEF
      E_Int nbVoisins = 0;
      for (E_Int elt1 = 0; elt1 < nelts; elt1++)
      {
        nbVoisins = cEEN[elt1].size();
        nbfaces1 = ptrEF[0];
        for (E_Int fa1 = 0; fa1 < nbfaces1; fa1++)
        {
          numface = ptrEF[fa1+1] -1;
          posfa1 = posFaces[numface]; 
          nbnodes1 = cnp[posfa1];
          sizeco2 += nbVoisins*(nbnodes1+1);
          sizeco1 += nbVoisins*(3*nbnodes1);
        }
        ptrEF += nbfaces1+1;
      }

      // nouvelle connectivite (representant les "sharp edges")
      FldArrayI newcnFN(sizeco1);
      FldArrayI newcnEF(sizeco2);
      E_Int* newptrFN = newcnFN.begin();
      E_Int* newptrEF = newcnEF.begin();
      sizeco1 = 0; sizeco2 = 0; ptrEF = cnp+sizeFN+4;
      nof =1;
      // parcours des elements
      for (E_Int elt1 = 0; elt1 < nelts; elt1++)
      {
        // dimension de l'element elt1
        dim = dimEltsp[elt1];
        // nombre de faces de elt1
        nbfaces1 = ptrEF[0];
        // recuperation des noeuds de elt1      
        K_CONNECT::getVertexIndices(cn->begin(), posFaces.begin(), pos[elt1], indices1);
        nvert = indices1.size();
        pts1.reserve(nvert);
        for (E_Int k = 0; k < nvert; k++)
        {
          E_Float* pt = new E_Float[3];
          ind = indices1[k]-1;
          pt[0] = x[ind]; pt[1] = y[ind]; pt[2] = z[ind];
          pts1.push_back(pt);
        }
        // parcours des elements voisins
        for (size_t ev = 0; ev < cEEN[elt1].size(); ev++)
        {
          // elt2 voisin
          E_Int elt2 = cEEN[elt1][ev];
          // nombre de faces de elt2
          //nbfaces2 = cnp[pos[elt2]];
          // recuperation des noeuds de elt2
          K_CONNECT::getVertexIndices(cn->begin(), posFaces.begin(), pos[elt2], indices2);
          nvert = indices2.size();
          pts2.reserve(nvert);
          for (E_Int k = 0; k < nvert; k++)
          {
            E_Float* pt = new E_Float[3];
            ind = indices2[k]-1;
            pt[0] = x[ind]; pt[1] = y[ind]; pt[2] = z[ind];
            pts2.push_back(pt);
          }
          // angle entre les polygones elt1 et elt2
          E_Float alpha;
          if (dim == 1)
          {
            E_Float dirVect[3];
            E_Float* ptA = new E_Float[3]; 
            E_Float* ptB = new E_Float[3]; 
            E_Float* ptC = new E_Float[3]; 
            E_Float* ptD = new E_Float[3]; 
            // determination de la face (donc d'un noeud en 1D) commune
            ptA[0] = pts1[0][0];  ptA[1] = pts1[0][1];  ptA[2] = pts1[0][2];
            ptB[0] = pts1[1][0];  ptB[1] = pts1[1][1];  ptB[2] = pts1[1][2];
            ptC[0] = pts2[0][0];  ptC[1] = pts2[0][1];  ptC[2] = pts2[0][2];
            ptD[0] = pts2[1][0];  ptD[1] = pts2[1][1];  ptD[2] = pts2[1][2];
            // produit vectoriel pour calculer dirVect
            E_Float dx1 = ptA[0]-ptB[0];
            E_Float dy1 = ptA[1]-ptB[1];
            E_Float dz1 = ptA[2]-ptB[2];
            E_Float dx2 = ptD[0]-ptC[0];
            E_Float dy2 = ptD[1]-ptC[1];
            E_Float dz2 = ptD[2]-ptC[2];
            E_Float p1 = dy1*dz2-dz1*dy2;
            E_Float p2 = dz1*dx2-dz2*dx1;
            E_Float p3 = dx1*dy2-dy1*dx2;
            E_Float n = sqrt(p1*p1+p2*p2+p3*p3);
            if (n > 1.e-12) 
            {
              dirVect[0] = p1/n; dirVect[1] = p2/n; dirVect[2] = p3/n; 
            }
            else
            {
              dirVect[0] = 0.; dirVect[1] = 0.; dirVect[2] = 1.;
            }
            alpha = K_COMPGEOM::getAlphaAngleBetweenBars(ptA,ptB,ptC,ptD,dirVect);
            if ((alpha < alphamin || alpha > alphamax) && alpha != -1000.) 
            {
              found1 = 0; found2 = 0;
              if ( indices1[0] == indices2[0] || indices1[0] == indices2[1] ) 
              {for (E_Int p=0;p<nfld;p++) fep[p][nop] = ptA[p];}
              else {for (E_Int p=0;p<nfld;p++) fep[p][nop] = ptC[p];}
              newptrFN[noe] = nop+1;  nop++; noe++;
            }
          }
          else
          {
            alpha = K_COMPGEOM::getAlphaAngleBetweenPolygons(pts1, pts2);
            if ( (alpha < alphamin || alpha > alphamax) && alpha != -1000.) 
            {
              // si une face fa1 coincide entierement avec des pts de elt2,
              // on ajoute cette face a la nouvelle connectivite
              for (E_Int fa1 = 0; fa1 < nbfaces1; fa1++)
              {
                numface = ptrEF[fa1+1] -1;
                posfa1 = posFaces[numface]; 
                nbnodes1 = cnp[posfa1];
                matchingPoints.clear();
                for (E_Int n = 0; n < nbnodes1; n++)
                {
                  indn = cnp[posfa1+1+n] -1;
                  x1 = x[indn]; y1 = y[indn]; z1 = z[indn];
                  for (unsigned int v2 = 0; v2 < pts2.size(); v2++)
                  {
                    x2 = pts2[v2][0]; y2 = pts2[v2][1]; z2 = pts2[v2][2];
                    if ((K_FUNC::fEqual(x1,x2,K_CONST::E_GEOM_CUTOFF) == true) &&
                        (K_FUNC::fEqual(y1,y2,K_CONST::E_GEOM_CUTOFF) == true) && 
                        (K_FUNC::fEqual(z1,z2,K_CONST::E_GEOM_CUTOFF) == true))
                    {
                      matchingPoints.push_back(n); break;
                    }
                  }
                } // fin boucle sur les noeuds n pour determiner le nb de points coincidents
                if (matchingPoints.size() == (size_t)nbnodes1)
                {
                  // connectivites Face/Noeuds et Elt/Faces : la face devient un element de la nouvelle connectivite
                  newptrEF[0] = nbnodes1; newptrEF += 1;
                  for (unsigned int v1 = 0; v1 < matchingPoints.size(); v1++)
                  {
                    E_Int pt = matchingPoints[v1];
                    indpt1 = cnp[posfa1+pt+1];
                    // nouvelle face : 
                    newptrEF[v1] = nof; nof++;
                    // noeud  de la nouvelle face
                    newptrFN[0] = 1;
                    newptrFN[1] = indpt1;
                    newptrFN += 2;
                  } 
                  newptrEF += matchingPoints.size();

                  sizeco1 += 2*matchingPoints.size(); next1 += matchingPoints.size();
                  sizeco2 += matchingPoints.size()+1; next2 += 1;

                } // fin test sur matchingPoints
              } // fin boucle sur les faces fa1
            } // fin test angle alpha
          }      
          for (size_t k = 0; k < pts2.size(); k++) delete [] pts2[k];
          pts2.clear();
        } // fin boucle elt2
        ptrEF += nbfaces1+1;
        for (size_t k = 0; k < pts1.size(); k++) delete [] pts1[k];
        pts1.clear();
      } // fin boucle elt1

      // reallocation pour la nouvelle connectivite
      newcnFN.resize(sizeco1); newcnEF.resize(sizeco2);

      // Construction des arrays de sortie + cleanConnectivity
      newptrFN = newcnFN.begin();
      newptrEF = newcnEF.begin();
    
      // aucune face exterieur : renvoi un objet nul
      if (sizeco1+sizeco2 == 0)
      {
        PyErr_SetString(PyExc_TypeError,
                        "sharpEdges: sharp Edges set is empty.");
        return NULL;
      }

      FldArrayI* cne = new FldArrayI(sizeco1+sizeco2+4);
      E_Int* cnep = cne->begin();
      cnep[0] = next1; // nombre de faces
      cnep[1] = sizeco1; // taille du tableau de faces
      cnep += 2;
      for (E_Int i = 0; i < sizeco1; i++) cnep[i] = newptrFN[i];
      cnep += sizeco1;
      cnep[0] = next2; // nombre d'elements
      cnep[1] = sizeco2; // taille du tableau d'elements
      cnep += 2;
      for (E_Int i = 0; i < sizeco2; i++) cnep[i] = newptrEF[i];
      K_CONNECT::cleanConnectivityNGon(posx, posy, posz, 1.e-10,
                                       *f, *cne);
      cne->setNGon(1);
      tpl = K_ARRAY::buildArray3(*f, varString0, *cne, "NGON", api);
      delete fe; delete cne;
      RELEASESHAREDU(array, f, cn);
    }
    return tpl;
  } // fin test NGON
}

