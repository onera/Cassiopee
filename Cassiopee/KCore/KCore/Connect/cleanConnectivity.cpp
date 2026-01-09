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
#include <string.h>
#include "String/kstring.h"
#include "Connect/connect.h"
# include "Nuga/include/ArrayAccessor.h"
#include "Nuga/include/merge.h"
# include "Nuga/include/BbTree.h"
#include <iostream>
#include "Nuga/include/ngon_t.hxx"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* 
   Nettoyage de la connectivite et du tableau des sommets.
*/
//=============================================================================
void K_CONNECT::cleanConnectivity(E_Int posx, E_Int posy, E_Int posz, 
                                  E_Float eps,  const char* eltType, 
                                  FldArrayF& f, FldArrayI& cn, bool remove_degen, bool ordered_merge)
{
  if (K_STRING::cmp(eltType, "NGON") == 0 || K_STRING::cmp(eltType, "NGON*") == 0) 
    cleanConnectivityNGon(posx, posy, posz, eps, f, cn, remove_degen, ordered_merge);
  else
    cleanConnectivityBasic(posx, posy, posz, eps, eltType, f, cn, ordered_merge);
}

//=============================================================================
/* 
   Nettoyage de la connectivite et du tableau des sommets: optim openmp 
   corse grain
*/
//=============================================================================
void K_CONNECT::cleanConnectivity_opt(E_Int posx, E_Int posy, E_Int posz, 
                                  E_Float eps,  const char* eltType, 
                                  FldArrayF& f, FldArrayI& cn, bool remove_degen, bool ordered_merge)
{
  if (K_STRING::cmp(eltType, "NGON") == 0 || K_STRING::cmp(eltType, "NGON*") == 0) 
    cleanConnectivityNGon(posx, posy, posz, eps, f, cn, remove_degen, ordered_merge);
  else
    cleanConnectivityBasic_opt(posx, posy, posz, eps, eltType, f, cn, ordered_merge);
}


//=============================================================================
/* 
   Nettoyage de la connectivite et du tableau des sommets associe 
   pour les elements basiques (TRI, TETRA, ...): 
   - suppression des noeuds doubles dans coord a eps pres.
   - suppression des noeuds non utilises
   - suppression des elements degeneres.
   - suppression des elements doubles pour les BARs.
   - Manque: suppression des elements doubles pour tous les autres elts.
   Optimise par KdTree.
*/
//=============================================================================
void K_CONNECT::cleanConnectivityBasic(E_Int posx, E_Int posy, E_Int posz, 
                                       E_Float eps,  const char* eltType, 
                                       FldArrayF& f, FldArrayI& cn,
                                       bool ordered_merge)
{
  E_Int nelts = cn.getSize(); E_Int nvert = cn.getNfld();
  E_Int npts = f.getSize(); E_Int nfld = f.getNfld();
  E_Int et1, ind1, np, indp;
  FldArrayI indir(npts); //indir.setAllValuesAt(-1);
  E_Int* indirp = indir.begin();
  FldArrayI* indir2 = new FldArrayI(npts);
  E_Int* indir2p = indir2->begin();

  // Recherche des vertex non utilises dans la connectivite
  // Construit f2 a partir de f, cn est modifie en cn2
  FldArrayI* used = new FldArrayI(npts);
  if (nelts == 0) used->setAllValuesAt(1); // cas NODE
  else used->setAllValuesAtNull(); // autres cas
  E_Int* usedp = used->begin();
  for (E_Int n1 = 1; n1 <= nvert; n1++)
  {
    E_Int* cnp = cn.begin(n1);
    for (et1 = 0; et1 < nelts; et1++)
    {
      ind1 = cnp[et1]-1; usedp[ind1]++;
    }
  }

  np = 0;
  for (E_Int i = 0; i < npts; i++)
  {
    if (usedp[i] > 0) // used
    {
      indirp[i] = np; indir2p[np] = i; np++;
    }
  }
  delete used;

  // Copie de la solution
  FldArrayF f2(np, nfld);
#pragma omp parallel default(shared)
  {
    for (E_Int n1 = 1; n1 <= nfld; n1++)
    {
      E_Float* fp = f.begin(n1);
      E_Float* f2p = f2.begin(n1);
#pragma omp for nowait
      for (E_Int i = 0; i < np; i++) { f2p[i] = fp[indir2p[i]]; }
    }
  }

  FldArrayI cn2(nelts, nvert);
#pragma omp parallel default(shared) private(et1, indp)
  {
    for (E_Int n1 = 1; n1 <= nvert; n1++)
    {
      E_Int* cnp = cn.begin(n1);
      E_Int* cn2p = cn2.begin(n1);
#pragma omp for nowait
      for (et1 = 0; et1 < nelts; et1++)
      {
        indp = cnp[et1]-1;
        cn2p[et1] = indirp[indp]+1;
      }
    }
  }

  // Recherche des doublons par kdtree
  // f2 devient f et cn2 -> cn
  np = 0;
  npts = f2.getSize();
  nelts = cn2.getSize();

  vector<E_Int> newId;
  {
    ArrayAccessor<FldArrayF> coordAcc(f2, posx, posy, posz);

    if (ordered_merge)
      ::merge(coordAcc, eps, newId, true /*do omp*/);
    else
      ::merge_no_order_omp(coordAcc, eps, newId);
    
    for (E_Int i = 0; i < npts; i++)
    {
      ind1 = newId[i];
      if (ind1 == i) // not merged
      {
        indirp[ind1] = np; indir2p[np] = ind1; np++;
      }
    }
  }

  // Copie de la solution
  f.malloc(np, nfld);
#pragma omp parallel default(shared)
  {
    for (E_Int n1 = 1; n1 <= nfld; n1++)
    {
      E_Float* fp = f.begin(n1);
      E_Float* f2p = f2.begin(n1);
#pragma omp for
      for (E_Int i = 0; i < np; i++) fp[i] = f2p[indir2p[i]];
    }
  }
  delete indir2; 

  // Creation de la nouvelle connectivite  Elements -> Noeuds
  K_CONNECT::createConnectEV(cn2, newId, indirp, cn);
    
  // Elimination des elements degeneres
  // nmax+1 = nbre de pts qui match pour enlever l'element
  // on tolere toute degenerescence sauf la BAR (pour les elts 2D et 3D)
  E_Int nmax; 
  if (K_STRING::cmp(eltType, "QUAD") == 0 || 
      K_STRING::cmp(eltType, "TETRA") == 0)
    nmax = 1;
  else if (K_STRING::cmp(eltType, "PYRA") == 0 || 
           K_STRING::cmp(eltType, "PENTA") == 0)
    nmax = 3;
  else if (K_STRING::cmp(eltType, "HEXA") == 0) nmax = 5;
  else nmax = 0;
  E_Int ec = 0;
  E_Int c, n2, ind, ind2;
    
  for (et1 = 0; et1 < nelts; et1++)
  {
    c = 0; // compte les elements degeneres
    for (E_Int n1 = 1; n1 <= nvert-1; n1++)
    {
      ind = cn(et1, n1);
      for (n2 = n1+1; n2 <= nvert; n2++)
      {
        ind2 = cn(et1, n2);
        if (ind2 == ind) c++;
      }
    }
    if (c <= nmax)
    {
      for (E_Int n1 = 1; n1 <= nvert; n1++) cn2(ec, n1) = cn(et1, n1);
      ec++;
    }
  }
  cn2.reAllocMat(ec, nvert); cn = cn2;
  if (K_STRING::cmp(eltType, "BAR") == 0)
  {
    FldArrayI cnout = cn;
    removeDoubleBARElts(f.begin(posx), f.begin(posy), f.begin(posz), 
                        cn, cnout);
    cn = cnout;
  }
}


//=============================================================================
/* 
   Nettoyage de la connectivite et du tableau des sommets associe 
   pour les elements basiques (TRI, TETRA, ...): 
   - suppression des noeuds doubles dans coord a eps pres.
   - suppression des noeuds non utilises
   - suppression des elements degeneres.
   - suppression des elements doubles pour les BARs.
   - Manque: suppression des elements doubles pour tous les autres elts.
   Optimise par KdTree.
*/
//=============================================================================
void K_CONNECT::cleanConnectivityBasic_opt(E_Int posx, E_Int posy, E_Int posz, 
                                           E_Float eps,  const char* eltType, 
                                           FldArrayF& f, FldArrayI& cn,
                                           bool ordered_merge)
{
  E_Int nelts = cn.getSize(); 
  E_Int nvert = cn.getNfld();
  E_Int npts  = f.getSize(); 
  E_Int nfld  = f.getNfld();

  E_Int et1, ind1, np, indp;

  FldArrayI indir(npts); E_Int* indirp = indir.begin();     //indir.setAllValuesAt(-1);
  FldArrayI* indir2 = new FldArrayI(npts);
  E_Int* indir2p = indir2->begin();

  // Recherche des vertex non utilises dans la connectivite
  // Construit f2 a partir de f, cn est modifie en cn2
  FldArrayI* used = new FldArrayI(npts);
  if (nelts == 0) used->setAllValuesAt(1); // cas NODE
  else used->setAllValuesAtNull(); // autres cas
  E_Int* usedp = used->begin();
  for (E_Int n1 = 1; n1 <= nvert; n1++)
  {
    E_Int* cnp = cn.begin(n1);
    for (et1 = 0; et1 < nelts; et1++)
    {
      ind1 = cnp[et1]-1; usedp[ind1]++;
    }
  }

  np = 0;
  for (E_Int i = 0; i < npts; i++)
  {
    if (usedp[i] > 0) // used
    {
      indirp[i] = np; indir2p[np] = i; np++;
    }
  }
  delete used;

  // Copie de la solution
  FldArrayF f2(np, nfld);
    for (E_Int n1 = 1; n1 <= nfld; n1++)
    {
      E_Float* fp = f.begin(n1);
      E_Float* f2p = f2.begin(n1);
      for (E_Int i = 0; i < np; i++) { f2p[i] = fp[indir2p[i]]; }
    }

  FldArrayI cn2(nelts, nvert);
    for (E_Int n1 = 1; n1 <= nvert; n1++)
    {
      E_Int* cnp = cn.begin(n1);
      E_Int* cn2p = cn2.begin(n1);
      for (et1 = 0; et1 < nelts; et1++)
      {
        indp = cnp[et1]-1;
        cn2p[et1] = indirp[indp]+1;
      }
    }

  // Recherche des doublons par kdtree
  // f2 devient f et cn2 -> cn
  np = 0;
  npts = f2.getSize();
  nelts = cn2.getSize();

  vector<E_Int> newId;
  {
    ArrayAccessor<FldArrayF> coordAcc(f2, posx, posy, posz);

    if (ordered_merge)
      ::merge(coordAcc, eps, newId, true/*do omp*/);
    else
      ::merge_no_order_omp(coordAcc, eps, newId);

    
    for (E_Int i = 0; i < npts; i++)
    {
      ind1 = newId[i];
      if (ind1 == i) // not merged
      {
        indirp[ind1] = np; indir2p[np] = ind1; np++;
      }
    }
  }

  // Copie de la solution
  f.malloc(np, nfld);

  for (E_Int n1 = 1; n1 <= nfld; n1++)
  {
    E_Float* fp = f.begin(n1);
    E_Float* f2p = f2.begin(n1);
    for (E_Int i = 0; i < np; i++) fp[i] = f2p[indir2p[i]];
  }
  delete indir2; 

  // Creation de la nouvelle connectivite  Elements -> Noeuds
  K_CONNECT::createConnectEV_opt(cn2, newId, indirp, cn);
    
  // Elimination des elements degeneres
  // nmax+1 = nbre de pts qui match pour enlever l'element
  // on tolere toute degenerescence sauf la BAR (pour les elts 2D et 3D)
  E_Int nmax; 
  if (K_STRING::cmp(eltType, "QUAD") == 0 || 
      K_STRING::cmp(eltType, "TETRA") == 0)
    nmax = 1;
  else if (K_STRING::cmp(eltType, "PYRA") == 0 || 
           K_STRING::cmp(eltType, "PENTA") == 0)
    nmax = 3;
  else if (K_STRING::cmp(eltType, "HEXA") == 0) nmax = 5;
  else nmax = 0;
  E_Int ec = 0;
  E_Int c, n2, ind, ind2;
    
  for (et1 = 0; et1 < nelts; et1++)
  {
    c = 0; // compte les elements degeneres
    for (E_Int n1 = 1; n1 <= nvert-1; n1++)
    {
      ind = cn(et1, n1);
      for (n2 = n1+1; n2 <= nvert; n2++)
      {
        ind2 = cn(et1, n2);
        if (ind2 == ind) c++;
      }
    }
    if (c <= nmax)
    {
      for (E_Int n1 = 1; n1 <= nvert; n1++) cn2(ec, n1) = cn(et1, n1);
      ec++;
    }
  }
  cn2.reAllocMat(ec, nvert); cn = cn2;
  if (K_STRING::cmp(eltType, "BAR") == 0)
  {
    FldArrayI cnout = cn;
    removeDoubleBARElts(f.begin(posx), f.begin(posy), f.begin(posz), 
                        cn, cnout);
    cn = cnout;
  }
}

//=======================================================================
// Elimination des vertex non utilise par une connectivite basique
// Retourne un nouvel fout, un nouvel cnout
//=======================================================================
void K_CONNECT::cleanUnreferencedVertexBasic(FldArrayF& f, FldArrayI& cn,
                                             FldArrayF& fout, FldArrayI& cnout)
{
  E_Int nelts = cn.getSize(); E_Int nvert = cn.getNfld();
  E_Int npts = f.getSize(); E_Int nfld = f.getNfld();
  E_Int et1, ind1, np, indp;
  FldArrayI indir(npts); //indir.setAllValuesAt(-1);
  E_Int* indirp = indir.begin();
  FldArrayI* indir2 = new FldArrayI(npts);
  E_Int* indir2p = indir2->begin();

  // Recherche des vertex non utilises dans la connectivite
  // Construit fout a partir de f, cn est modifie en cnout
  FldArrayI* used = new FldArrayI(npts);
  if (nelts == 0) used->setAllValuesAt(1); // cas NODE
  else used->setAllValuesAtNull(); // autres cas
  E_Int* usedp = used->begin();
  for (E_Int n1 = 1; n1 <= nvert; n1++)
  {
    E_Int* cnp = cn.begin(n1);
    for (et1 = 0; et1 < nelts; et1++)
    {
      ind1 = cnp[et1]-1; usedp[ind1]++;
    }
  }

  np = 0;
  for (E_Int i = 0; i < npts; i++)
  {
    if (usedp[i] > 0) // used
    {
      indirp[i] = np; indir2p[np] = i; np++;
    }
  }
  delete used;

  // Copie de la solution
  fout.malloc(np, nfld);
#pragma omp parallel default(shared)
  {
    for (E_Int n1 = 1; n1 <= nfld; n1++)
    {
      E_Float* fp = f.begin(n1);
      E_Float* f2p = fout.begin(n1);
#pragma omp for nowait
      for (E_Int i = 0; i < np; i++) { f2p[i] = fp[indir2p[i]]; }
    }
  }

  cnout.malloc(nelts, nvert);
#pragma omp parallel default(shared) private(et1, indp)
  {
    for (E_Int n1 = 1; n1 <= nvert; n1++)
    {
      E_Int* cnp = cn.begin(n1);
      E_Int* cn2p = cnout.begin(n1);
#pragma omp for nowait
      for (et1 = 0; et1 < nelts; et1++)
      {
        indp = cnp[et1]-1;
        cn2p[et1] = indirp[indp]+1;
      }
    }
  }
  delete indir2;
}

//=============================================================================
// create de la connection Elements -> Noeuds
// ------------------------------------------
// IN: cEV: tableau de connectivite regroupant la connection Elements->Vertex
// newId
// indirp
// --------------
// OUT: cEVout: tableau de connectivite nettoye, regroupant la connection 
// Elements -> Vertex
// cnout alloue en externe
//=============================================================================
void K_CONNECT::createConnectEV(FldArrayI& cn, vector<E_Int> newId, 
                                E_Int* indirp, FldArrayI& cnout)
{
  E_Int nelts = cn.getSize(); E_Int nvert = cn.getNfld();

#pragma omp parallel default(shared)
  {
    E_Int indp, indn;
    E_Int n1, et1;
    
    for (n1 = 1; n1 <= nvert; n1++)
    {
      E_Int* cnp = cn.begin(n1);
      E_Int* cnoutp = cnout.begin(n1);
#pragma omp for nowait   
      for (et1 = 0; et1 < nelts; et1++)
      {
        indp = cnp[et1]-1; indn = newId[indp];
        cnoutp[et1] = indirp[indn]+1;
      }
    }
  }
}

//=============================================================================
// create de la connection Elements -> Noeuds: optim openmp coarsegrain
// ------------------------------------------
// IN: cEV: tableau de connectivite regroupant la connection Elements->Vertex
// newId
// indirp
// --------------
// OUT: cEVout: tableau de connectivite nettoye, regroupant la connection 
// Elements -> Vertex
// cnout alloue en externe
//=============================================================================
void K_CONNECT::createConnectEV_opt(FldArrayI& cn, vector<E_Int> newId, 
                                    E_Int* indirp, FldArrayI& cnout)
{
  E_Int nelts = cn.getSize(); E_Int nvert = cn.getNfld();

  E_Int indp, indn;
  E_Int n1, et1;
    
  for (n1 = 1; n1 <= nvert; n1++)
  {
    E_Int* cnp = cn.begin(n1);
    E_Int* cnoutp = cnout.begin(n1);
    for (et1 = 0; et1 < nelts; et1++)
    {
      indp = cnp[et1]-1; indn = newId[indp];
      cnoutp[et1] = indirp[indn]+1;
    }
  }
}
//=============================================================================
/* 
   CleanConnectivity specifique aux NGons
   0. elimination des noeuds non references par les faces
   1. elimination des noeuds doubles dans f a eps pres et maj dans la 
   connectivite Faces/Noeuds cFN
   2. supp faces degenerees ->modif de cFN 
   (ex une face quad degeneree passe en tri)
   3. elimination des faces doubles et maj dans cEF
   4. elimination des faces non utilisees dans les elements
   [5. supp des elemts doubles] : INACTIF
   [6. elimination, des elemts degeneres]

*/
//=============================================================================
void K_CONNECT::cleanConnectivityNGon(E_Int posx, E_Int posy, E_Int posz, 
                                      E_Float eps, FldArrayF& f, FldArrayI& cn, bool remove_degen, bool ordered_merge)
{
  using acrd_t = K_FLD::ArrayAccessor<FldArrayF> ;
  acrd_t fcA(f, posx, posy, posz);
  typedef ngon_t<K_FLD::FldArrayI> ngon_type;
  ngon_type NG(cn);
  
  NG.PGs.updateFacets();
  NG.PHs.updateFacets();
  
  // type of NGON: volumic , surfacic or lineic?
  E_Int ngon_dim = 3; // volumic
  E_Int max_stride=0;
  for (E_Int i=0; (i<NG.PGs.size()); ++i)
    max_stride = std::max(max_stride, NG.PGs.stride(i));
  
  if (max_stride == 1) // lineic ngon
    ngon_dim=1;
  else if (max_stride == 2) //surfacic ngon
    ngon_dim=2;

  // 1- Referencement unique pour les noeuds confondus par kdtree (les doublons seront supprimes a la fin)
  // si la tolerance est :
  // * positive : absolue
  // * nulle => pas de merge
  // * negative => relative
  if (eps >= 0.)
  {
    Vector_t<E_Int> nids;

    if (ordered_merge)
      ::merge(fcA, eps, nids, true/*do omp*/);
    else
      ::merge_no_order_omp(fcA, eps, nids);

    NG.PGs.change_indices(nids);
  }
  else // if (eps < 0.) 
  {
    /*std::cout << "RELATIVE TOL !!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    NG.flag_external_pgs(INITIAL_SKIN);
    std::vector<double> nodal_metric2;
    NUGA::MeshTool::computeNodalDistance2<acrd_t, ngon_unit>(fcA, NG.PGs, nodal_metric2);
    
    double RTOL = -eps;
    NG.join_phs2(f, posx, posy, posz, nodal_metric2, RTOL);*/

    std::vector<double> nodal_metric2;
    NUGA::MeshTool::computeNodalDistance2<acrd_t, ngon_unit>(fcA, NG.PGs, nodal_metric2);
    
    double RTOL = -eps;

    Vector_t<E_Int> nids;

    ::merge(fcA, nodal_metric2, RTOL, nids, true /*do_omp*/);
    // merge_no_order_omp not implemented for relative tol

    NG.PGs.change_indices(nids);
  }
  
  // 2- Elimination des faces degenerees
  std::vector<E_Int> pgids, phnids;
  if (ngon_dim != 1)
  {
    NG.remove_degenerated_pgs(ngon_dim, pgids, phnids);
    NG.PGs.remove_consecutive_duplicated(); //compact representation
  }
   
  // 3- Elimination des faces confondues
  if (ngon_dim == 3) //volumic
  {
    Vector_t<E_Int> pgnidstmp;
    NG.remove_duplicated_pgs(fcA,pgnidstmp, true/*do_omp*/);
  }
  else if (ngon_dim == 2) //surfacic
    NG.remove_duplicated_edges();
  else // lineic
    NG.remove_duplicated_nodes();

  // remove duplicated references to PGs within each elements
  NG.PHs.remove_duplicated();

  NG.PHs.updateFacets();

  // 4- Elimination des elts degeneres
  Vector_t<E_Int> toremove;
  if (remove_degen)
  {
    //std::cout << "REMOVE DEGN ENABLED" << std::endl;
    E_Int min_nb_facets = ngon_dim + 1;
    NG.PHs.get_degenerated(min_nb_facets, toremove);
  }

  // 5- Elimination des elts doubles : do not care of multiple occ in toremove as remove_entities handles it.
  /*if (remove_dups)
  {
    std::vector<E_Int> duphnids;
    NG.detect_phs_with_same_centroid(f, duphnids);
    for (size_t k = 0; k < duphnids.size(); ++k)
    {
      if (duphnids[k] != (E_Int)k)
        toremove.push_back(k);
    }
  }*/

  NG.PHs.remove_entities(toremove, phnids);
  
  // 6- Suppression des faces non referencees
  NG.remove_unreferenced_pgs(pgids, phnids);
  
  // 7- Compression du ngon aux seuls noeuds utilises
  FldArrayF fcopy = f;
  ngon_type::compact_to_used_nodes(NG.PGs, fcopy);
  f = fcopy; //fixme: FldArrayF not ready for dynamic
  
  // 8- Mise a disposition de la cn de sortie au format FldArrayI
  NG.export_to_array(cn);
}
//=============================================================================
 /*
    Supprime les elements definis 2 fois dans une BAR
    IN: xt, yt,zt: coords
    IN: cEV: connectivite element/vertex BAR
    IN: cEVout: connectivite element/vertex sans elements doublons.
  */
//=============================================================================
void K_CONNECT::removeDoubleBARElts(E_Float* xt, E_Float*yt, E_Float* zt, 
                                    FldArrayI& cn, FldArrayI& cnout)
{
  E_Int nelts = cn.getSize();
  E_Int* cn1 = cn.begin(1); E_Int* cn2 = cn.begin(2);
  short found = 0;
  E_Int noet = 0;
  FldArrayI dejaVu(nelts); dejaVu.setAllValuesAtNull();
  for (E_Int et1 = 0; et1 < nelts; et1++)
  {
    for (E_Int et2 = et1+1; et2 < nelts; et2++)
    {
      found = 0;
      if (dejaVu[et1] == 0 && dejaVu[et2] == 0) 
      {
        if (cn1[et1] == cn1[et2]) found += 1;
        else if (cn1[et1] == cn2[et2]) found += 1;
        if (cn2[et1] == cn1[et2]) found += 2;
        else if (cn2[et1] == cn2[et2]) found += 2;
        if (found == 3)
        {
          cnout(noet,1) = cn1[et1]; cnout(noet,2) = cn2[et1]; noet++; 
          dejaVu[et1] = 1; dejaVu[et2] = 1; goto end;
        }
      }
    }

    if (dejaVu[et1] == 0) 
    {  
      cnout(noet,1) = cn1[et1]; cnout(noet,2) = cn2[et1]; noet++;
      dejaVu[et1] = 1;
    }
    end:;
  }
  cnout.reAllocMat(noet, 2);
}

//==============================================================================
// Supprime les faces non referencees par les elements dans une connectivite
// NGON
//==============================================================================
void K_CONNECT::cleanUnreferencedFacesNGon(FldArrayI& cn)
{
  typedef ngon_t<K_FLD::FldArrayI> ngon_type;
  ngon_type NG(cn);
  // Suppression des faces non referencees
  std::vector<E_Int> pgnids, phnids;
  NG.remove_unreferenced_pgs(pgnids, phnids);
  // Mise a disposition de la cn de sortie au format FldArrayI
  NG.export_to_array(cn);
}
