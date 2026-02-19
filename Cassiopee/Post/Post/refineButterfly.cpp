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

// refine elements of a surface triangular mesh using butterfly interpolation

# include "post.h"
using namespace std;
using namespace K_FLD;
E_Int getOtherVertex(E_Int* cn1, E_Int* cn2, E_Int* cn3, E_Int e, E_Int v1, E_Int v2);
void getCoeffs(E_Float w, E_Int p1, E_Int p2, E_Int p3, E_Int nume, FldArrayI cEV, E_Int np,
               vector< vector<E_Int> >& cEEN, vector< vector<E_Int> >& cVE,
               vector<E_Int>& neigh1,
               E_Int* cn1, E_Int* cn2, E_Int* cn3,
               E_Int& p4, E_Int& p5, E_Int& p6, E_Int& p7, E_Int& p8, E_Int& p9, E_Int& p10,
               E_Float& w1, E_Float& w2, E_Float& w3, E_Float& w4,
               E_Float& w5, E_Float& w6, E_Float& w7, E_Float& w8, E_Float& w9, E_Float& w10);

//=============================================================================
/* Raffine un maillage surfacique TRI partour en utilisant l'algorithme
   butterfly.
   Retourne un maillage surfacique TRI modifie */
//=============================================================================
PyObject* K_POST::refineButterfly(PyObject* self, PyObject* args)
{
  // surf: maillage a raffiner (x,y,z+sol)
  // w: parametre pour le butterfly
  PyObject* surf; E_Float w;
  if (!PYPARSETUPLE_(args, O_ R_, &surf, &w))
  {
    return NULL;
  }

  /*-----------------------------------------------*/
  /* Extraction des donnees du maillage surfacique */
  /*-----------------------------------------------*/
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int nil, njl, nkl;
  E_Int res =
    K_ARRAY::getFromArray3(surf, varString, f, nil, njl, nkl, cn, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "refine: input array is invalid.");
    return NULL;
  }
  if (res != 2 || strcmp(eltType, "TRI") != 0)
  {
    RELEASESHAREDB(res, surf, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "refine: input array must be TRI.");
    return NULL;
  }

  // Check coordinates
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "refine: coordinates not found in array.");
    RELEASESHAREDB(res, surf, f, cn);
    return NULL;
  }
  posx++; posy++; posz++;

  E_Int api = f->getApi();
  FldArrayF* fo; FldArrayI* cno;
  refineButterfly(*f, *cn, w, fo, cno);

  E_Float tolc = 1.e-12;
  K_CONNECT::cleanConnectivity(posx, posy, posz, tolc, "TRI", *fo, *cno);

  PyObject* t = K_ARRAY::buildArray3(*fo, varString, *cno, "TRI", api);
  delete fo; delete cno;
  return t;
}

//=============================================================================
void K_POST::refineButterfly(FldArrayF& f, FldArrayI& cn, E_Float w,
  FldArrayF*& fo, FldArrayI*& co)
{
  E_Int i, v, n1, n2, n3;
  // Nombre de champs
  E_Int nf = f.getNfld();
  // Nombre d'elements sur le maillage d'origine
  E_Int ne = cn.getSize();
  // Nombre de noeuds sur le maillage d'origine
  E_Int np = f.getSize();

  // Connectivite des elements voisins
  vector< vector<E_Int> > cEEN(ne);
  K_CONNECT::connectEV2EENbrs("TRI", np, cn, cEEN);
  // Connectivite des elements voisins
  vector< vector<E_Int> > cVE(np);
  K_CONNECT::connectEV2VE(cn, cVE);

  // Nouvelle connectivite
  co = new FldArrayI(4*ne, cn.getNfld());

  // Nouveaux coord des points
  fo = new FldArrayF(np+3*ne, nf);
  for (v = 1; v <= nf; v++)
  {
    E_Float* fop = fo->begin(v);
    E_Float* fp = f.begin(v);
    for (i = 0; i < np; i++) fop[i] = fp[i];
  }

  // no de l'element courant dans le nouveau maillage
  E_Int nelt = 0;
  // no du point courant dans le nouveau maillage
  E_Int npt = np;
  E_Int* cn1 = cn.begin(1); E_Int* cn2 = cn.begin(2); E_Int* cn3 = cn.begin(3);
  E_Int* co1 = co->begin(1); E_Int* co2 = co->begin(2); E_Int* co3 = co->begin(3);
  E_Float* fop; E_Float* fp;
  E_Int p1,p2,p3,p4,p5,p6,p7,p8,p9,p10;
  E_Float w1,w2,w3,w4,w5,w6,w7,w8,w9,w10;

  for (i = 0; i < ne; i++)
  {
    // printf("element %d\n",i );
    // vertices
    n1 = cn1[i]-1; n2 = cn2[i]-1; n3 = cn3[i]-1;
    // neighbours (level1)
    vector<E_Int>& neigh1 = cEEN[i];

    // Point X1
    p1 = n1; p2 = n2; p3 = n3;
    getCoeffs(w,p1,p2,p3,i,cn,np,cEEN,cVE,neigh1,cn1,cn2,cn3,p4,p5,p6,p7,p8,p9,p10,
              w1,w2,w3,w4,w5,w6,w7,w8,w9,w10);

    for (v = 1; v <= nf; v++)
    {
      fop = fo->begin(v); fp = f.begin(v);
      fop[npt] = w1*fp[p1]+w2*fp[p2]+w3*fp[p3]+w4*fp[p4]+w5*fp[p5]+
      w6*fp[p6]+w7*fp[p7]+w8*fp[p8]+w9*fp[p9]+w10*fp[p10];
    }

    // Point X2
    p1 = n1; p2 = n3; p3 = n2;
    getCoeffs(w,p1,p2,p3,i,cn,np,cEEN,cVE,neigh1,cn1,cn2,cn3,p4,p5,p6,p7,p8,p9,p10,
              w1,w2,w3,w4,w5,w6,w7,w8,w9,w10);

    for (v = 1; v <= nf; v++)
    {
      fop = fo->begin(v); fp = f.begin(v);
      fop[npt+1] = w1*fp[p1]+w2*fp[p2]+w3*fp[p3]+w4*fp[p4]+w5*fp[p5]+
      w6*fp[p6]+w7*fp[p7]+w8*fp[p8]+w9*fp[p9]+w10*fp[p10];
    }

    // Point X3
    p1 = n2; p2 = n3; p3 = n1;
    getCoeffs(w,p1,p2,p3,i,cn,np,cEEN,cVE,neigh1,cn1,cn2,cn3,p4,p5,p6,p7,p8,p9,p10,
              w1,w2,w3,w4,w5,w6,w7,w8,w9,w10);    

    for (v = 1; v <= nf; v++)
    {
      fop = fo->begin(v); fp = f.begin(v);
      fop[npt+2] = w1*fp[p1]+w2*fp[p2]+w3*fp[p3]+w4*fp[p4]+w5*fp[p5]+
      w6*fp[p6]+w7*fp[p7]+w8*fp[p8]+w9*fp[p9]+w10*fp[p10];
    }

    // Connectivite
    co1[nelt] = n1+1;
    co2[nelt] = npt+1;
    co3[nelt] = npt+2;
    co1[nelt+1] = npt+2;
    co2[nelt+1] = npt+3;
    co3[nelt+1] = n3+1;
    co1[nelt+2] = npt+1;
    co2[nelt+2] = npt+3;
    co3[nelt+2] = npt+2;
    co1[nelt+3] = npt+1;
    co2[nelt+3] = n2+1;
    co3[nelt+3] = npt+3;
    npt += 3;
    nelt += 4;
  }
}

//=============================================================================
// Soit deux no de vertex d'un element, retourne le dernier vertex de l'element
//=============================================================================
E_Int getOtherVertex(E_Int* cn1, E_Int* cn2, E_Int* cn3, E_Int e,
  E_Int v1, E_Int v2)
{
  E_Int n1, n2, n3;
  n1 = cn1[e]-1; n2 = cn2[e]-1; n3 = cn3[e]-1;
  if (v1 == n1)
  {
    if (v2 == n2) return n3;
    else if (v2 == n3) return n2;
  }
  else if (v1 == n2)
  {
    if (v2 == n1) return n3;
    else if (v2 == n3) return n1;
  }
  else if (v1 == n3)
  {
    if (v2 == n1) return n2;
    else if (v2 == n2) return n1;
  }
  return -1;
}


//==============================================================================
void getCoeffs(E_Float w, E_Int p1, E_Int p2, E_Int p3, E_Int nume, FldArrayI cEV, E_Int np,
               vector< vector<E_Int> >& cEEN, vector< vector<E_Int> >& cVE,
               vector<E_Int>& neigh1,
               E_Int* cn1, E_Int* cn2, E_Int* cn3,
               E_Int& p4, E_Int& p5, E_Int& p6, E_Int& p7, E_Int& p8, E_Int& p9, E_Int& p10,
               E_Float& w1, E_Float& w2, E_Float& w3, E_Float& w4,
               E_Float& w5, E_Float& w6, E_Float& w7, E_Float& w8, E_Float& w9, E_Float& w10)
{
  E_Int sizeNeigh1,e,ret,sizeNeigh2,e2;
  p4 = -1; p5 = -1; p6 = -1; p7 = -1; p8 = -1; p9 = -1; p10 = -1;
  sizeNeigh1 = neigh1.size();
  for (E_Int i = 0; i < sizeNeigh1; i++)
  {
    // check for P1, P2
    e = neigh1[i];
    ret = getOtherVertex(cn1, cn2, cn3, e, p1, p2);
    if (ret != -1)
    {
      p4 = ret;
      vector<E_Int>& neigh2 = cEEN[e];
      sizeNeigh2 = neigh2.size();
      for (E_Int j = 0; j < sizeNeigh2; j++)
      {
        e2 = neigh2[j];
        ret = getOtherVertex(cn1, cn2, cn3, e2, p1, p4);
        if (ret != -1) p6 = ret;
        ret = getOtherVertex(cn1, cn2, cn3, e2, p2, p4);
        if (ret != -1) p5 = ret;
      }
    }
    else
    {
      ret = getOtherVertex(cn1, cn2, cn3, e, p1, p3);
      if (ret != -1) p7 = ret;
      else
      {
        p8 = getOtherVertex(cn1, cn2, cn3, e, p2, p3);
      }
    }
  }

  w1 = 0.5; w2 = 0.5; w3 = 2*w; w4 = 2*w;
  w5 = -w; w6 = -w; w7 = -w; w8 = -w; w9 = 0.; w10 = 0.;
  //w3 = 0; w4 = 0; w5 = 0; w6 = 0; w7 = 0; w8 = 0;
  //printf("%d %d: %f %f %f %f %f %f %f %f \n",p1,p2,w1,w2,w3,w4,w5,w6,w7,w8);
  //printf("A - %d %d: %d %d %d %d %d %d\n",p1,p2,p3,p4,p5,p6,p7,p8);

  if (p7 == -1) p7 = p1;
  if (p8 == -1) p8 = p2;
  if (p5 == -1) p5 = p2;
  if (p6 == -1) p6 = p1;
  if (p4 == -1) // boundary
  {
    E_Int sizep1 = cVE[p1].size(); // nb elements avec noeud p1
    E_Int sizep2 = cVE[p2].size(); // nb elements avec noeud p2
    // printf("%d %d %d %d\n", p1,p2,sizep1,sizep2);
    std::list<int> list_ind_orig; // liste indices noeuds element courant
    list_ind_orig.push_back(p1);list_ind_orig.push_back(p2);list_ind_orig.push_back(p3);
    std::list<int> list_ind_mod; // liste indices noeuds elements voisins
    list_ind_mod.push_back(p1);list_ind_mod.push_back(p2);list_ind_mod.push_back(p3);
    // printf("element %d avec noeuds %d %d %d\n", nume,p1,p2,p3); // nume = indice element courant
    // traitement a gauche de p1
    for (E_Int i = 0; i < sizep1; i++)
    {
      E_Int numev = cVE[p1][i]; // numev = indice element voisin p1
      if (numev != nume)
      {
        E_Int nv1 = cEV.begin(1)[numev]-1;
        E_Int nv2 = cEV.begin(2)[numev]-1;
        E_Int nv3 = cEV.begin(3)[numev]-1;
        // printf("noeuds de l element voisin %d %d %d\n", nv1,nv2,nv3);
        if ( std::find(list_ind_mod.begin(), list_ind_mod.end(), nv1) != list_ind_mod.end() )
        {
          list_ind_mod.remove(nv1);
        }
        else
        {
          list_ind_mod.push_back(nv1);
        }
        if ( std::find(list_ind_mod.begin(), list_ind_mod.end(), nv2) != list_ind_mod.end() )
        {
          list_ind_mod.remove(nv2);
        }
        else
        {
          list_ind_mod.push_back(nv2);
        }
        if ( std::find(list_ind_mod.begin(), list_ind_mod.end(), nv3) != list_ind_mod.end() )
        {
          list_ind_mod.remove(nv3);
        }
        else
        {
          list_ind_mod.push_back(nv3);
        }
      }
    }
    // comparaison liste originale et liste modifiee
    for (std::list<int>::iterator it=list_ind_mod.begin(); it != list_ind_mod.end(); ++it)
    {
        if ( std::find(list_ind_orig.begin(), list_ind_orig.end(), *it) != list_ind_orig.end() )
        {
        }
        else
        {
          p9 = *it;
          // printf("noeud voisin gauche de %d est %d\n", p1,p9);
        }
    }
    if (sizep1 <= 2) p9 = p1;

    // traitement a droite de p2
    list_ind_mod.clear();list_ind_mod.push_back(p1);list_ind_mod.push_back(p2);list_ind_mod.push_back(p3);
    for (E_Int i = 0; i < sizep2; i++)
    {
      E_Int numev = cVE[p2][i]; // numev = indice element voisin p1
      if (numev != nume)
      {
        E_Int nv1 = cEV.begin(1)[numev]-1;
        E_Int nv2 = cEV.begin(2)[numev]-1;
        E_Int nv3 = cEV.begin(3)[numev]-1;
        // printf("noeuds de l element voisin %d %d %d\n", nv1,nv2,nv3);
        if ( std::find(list_ind_mod.begin(), list_ind_mod.end(), nv1) != list_ind_mod.end() )
        {
          list_ind_mod.remove(nv1);
        }
        else
        {
          list_ind_mod.push_back(nv1);
        }
        if ( std::find(list_ind_mod.begin(), list_ind_mod.end(), nv2) != list_ind_mod.end() )
        {
          list_ind_mod.remove(nv2);
        }
        else
        {
          list_ind_mod.push_back(nv2);
        }
        if ( std::find(list_ind_mod.begin(), list_ind_mod.end(), nv3) != list_ind_mod.end() )
        {
          list_ind_mod.remove(nv3);
        }
        else
        {
          list_ind_mod.push_back(nv3);
        }
      }
    }
    // comparaison liste originale et liste modifiee
    for (std::list<int>::iterator it=list_ind_mod.begin(); it != list_ind_mod.end(); ++it)
    {
        if ( std::find(list_ind_orig.begin(), list_ind_orig.end(), *it) != list_ind_orig.end() )
        {
        }
        else
        {
          p10=*it;
          // printf("noeud voisin droite de %d est %d\n", p2,p10);
        }
    }
    if (sizep2 <= 2) p10 = p2;

    w1 = 9./16.; w2 = 9./16.; w3 = 0; w4 = 0; w5 = 0; w6 = 0; w7 = 0; w8 = 0; w9 = -1./16.; w10 = -1./16;
    // p4 = p3; w3 = 0; w4 = 0; w5 = 0; w6 = 0; w7 = 0; w8 = 0; // lineaire
  }
}
