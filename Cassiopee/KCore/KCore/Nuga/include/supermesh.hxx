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
//Authors : Sam Landier (sam.landier@onera.fr)
#ifndef NUGA_SUPERMESH_HXX
#define NUGA_SUPERMESH_HXX

//#define XMATCH_DBG

#include "Nuga/include/clipper.hxx"


#ifdef XMATCH_DBG
#include "Nuga/include/medit.hxx"
#endif


namespace NUGA
{

struct field
{
  E_Float* f;
  E_Float* gradf[3];
  field() :f(nullptr) { gradf[0] = gradf[1] = gradf[2] = nullptr; }
};

inline void __get_local_edge_ids(long szudzik_val, E_Int n0, E_Int& le0, E_Int& le1)
{
  // on decode les deux valeurs qui sont les indices locaux d'arete en intersection
  
  NUGA::szudzik_unpairing(szudzik_val, le0, le1);

  if (n0 != IDX_NONE) // ae1 has been reversed so need to convert local edge numbering
  {
    // conversion de l'indice local d'arete de l'element inverse vers l'element initial
    if (le1 < (n0 - 1))
      le1 = (n0 - 2) - le1;
  }
}

inline void __manage_x_points(aPolygon& bit, E_Int n0,
                              const E_Int* ge0s, const E_Int* ge1s,
                              K_FLD::FloatArray& xmcrd,
                              std::map<K_MESH::Edge, E_Int>& key_to_ptxid,
                              std::set<std::pair<E_Int, E_Int>>& edge_to_xnodes)
{
  for (size_t p = 0; p < bit.m_poids.size(); ++p)
  {
    bool is_an_x_point = (bit.m_poids[p] < 0);
    // valuation bits[k].m_poids[p] : soit avec un id dans xm.crd (le point est cree si inexistant)
    if (is_an_x_point)
    {
      // algo de szudzik : permet d'encoder 2 entiers en 1 seul
      // ici on recupere la valeur encodee
      E_Int szudzik_val = -bit.m_poids[p] - 1;

      E_Int le0, le1; //local edge ids
      __get_local_edge_ids(szudzik_val, n0, le0, le1);

      K_MESH::Edge key(ge0s[le0], ge1s[le1]);
      auto it = key_to_ptxid.find(key);

      if (it == key_to_ptxid.end())
      {
        // nouveau point a creer et a stocker : xpt aura new_ptid comme id dans xm.crd
        E_Int new_ptid = xmcrd.cols();
        E_Float* xpt = bit.m_crd.col(p);

        xmcrd.pushBack(xpt, xpt + 3); // stockage du point
        key_to_ptxid[key] = new_ptid; // association de la paire d'arete a ce point
        bit.m_poids[p] = new_ptid; // on fait reference ce point dans l'hsito => xm.add va utiliser ce point
      }
      else
        bit.m_poids[p] = it->second;

      edge_to_xnodes.insert(std::make_pair(le0, bit.m_poids[p])); // on recense le point pour construire a posteriori le ae0 impact�
    }
  }
}

inline void __impact_ae0(aPolygon& ae0, const std::set<std::pair<E_Int, E_Int>>& edge_to_xnodes, const K_FLD::FloatArray& xmcrd)
{
  std::vector<std::pair<E_Float, E_Int>> lambda_to_node;

  E_Int nnodes = ae0.nb_nodes();

  for (E_Int n = 0; n < nnodes; ++n)
    lambda_to_node.push_back(std::make_pair(E_Float(n), ae0.m_poids[n]));

  E_Float P0PX[3], P0P1[3];
  for (auto& e2n : edge_to_xnodes)
  {
    int n0 = e2n.first;

    const E_Float* P0 = xmcrd.col(ae0.m_poids[n0]);
    const E_Float* PX = xmcrd.col(e2n.second);
    const E_Float* P1 = xmcrd.col(ae0.m_poids[(n0 + 1) % nnodes]);

    NUGA::diff<3>(P1, P0, P0P1);
    NUGA::diff<3>(PX, P0, P0PX);

    E_Float Lx = sqrt(NUGA::sqrNorm<3>(P0PX));
    E_Float L = sqrt(NUGA::sqrNorm<3>(P0P1));

    assert(Lx <= L);

    E_Float lambda = (Lx / L) + E_Float(n0);

    lambda_to_node.push_back(std::make_pair(lambda, e2n.second));
  }

  std::sort(ALL(lambda_to_node));

  // from lambda_to_node to ae0
  K_FLD::FloatArray impact_crd;;
  std::vector<long> poids;

  E_Float lambda_prev = -1.;
  for (auto& l2n : lambda_to_node)
  {
    E_Float& lambda = l2n.first;
    if (lambda - lambda_prev > ZERO_M)
    {
      impact_crd.pushBack(xmcrd.col(l2n.second), xmcrd.col(l2n.second) + 3);
      poids.push_back(l2n.second);
    }
    lambda_prev = lambda;
  }

  ae0 = aPolygon(std::move(impact_crd));
  ae0.m_poids = poids;

#ifdef SUPERMESH_DBG
  std::ostringstream o;
  o << "modified_ae0";
  medith::write(o.str().c_str(), ae0.m_crd, &ae0.m_nodes[0], ae0.m_nodes.size(), ae0.shift());
  medith::write<DELAUNAY::Triangulator>("bits", xm.crd, xm.cnt, &xbit_ids, 0);
#endif
}

template <typename zmesh_t> inline
void xmatch(const zmesh_t& m0, const zmesh_t& m1, E_Float ARTOL, std::vector<E_Int>& anc0, std::vector<E_Int>& anc1, zmesh_t& xm, bool proj_on_first)
{
  using aelt_t = typename zmesh_t::aelt_t;

  xm.clear();
  anc0.clear();
  anc1.clear();

  auto loc1 = m1.get_localizer();

  // chaque pair d'arete en intersection est consideree jusqu'a 4 fois (chaque arete pouvant etre partagee par 2 faces)
  // les faces ayant des normales differentes, la projection pour se ramener au 2D peut etre differente, ces points d'intersection peuvent donc etre differents 
  // on decide donc de le calculer une seule fois et de l'utiliser pour les 4 decoupages => on force ainsi la conformite
  // on a donc besoin de connaitre l'indice global de chaque arete
  ngon_unit glob_edge_ids0;
  m0.build_global_edge_ids(glob_edge_ids0);
  ngon_unit glob_edge_ids1;
  m1.build_global_edge_ids(glob_edge_ids1);
  // on le stocke donc une map "cle vers id du pt d'intersection dans xm.crd". La cle etant formee a partir des indices globaux d'arete
  std::map<K_MESH::Edge, E_Int> key_to_ptxid;
  // Pour pouvoir reutiliser un maximum de points existants (voir les appels aux fonctions add plus bas), on concatene au depart tous les points de m0 et m1 dans xm.crd
  // cela minimise le besoin de fusion en aval
  xm.crd = m0.crd;
  xm.crd.pushBack(m1.crd); // => en extrayant un element de m1, on peut convertir son historique (ae1.m_poids) pour faire reference a des points de xm en decalant les indices (voir appel a  IdTool::shift dans la boucle)

  std::vector<E_Int> cands; // ids des candidats dans m2 en potention intersection avec ae0
  std::vector<aelt_t> bits;    // moreceaaux de l'elt courant ae0

  // donnees pour CLASSIFIER ae0 : classifier c'est dire si ae0 est completement dedans ou hors de ae1 s'il n'y a pas d'intersections retournee par :CLIP::isolated_clip. Dans ce cas true_clip vaut false.
  // pour classifier on pass en 2D : dans le referentiel de ae0 (2 axes dans le plan moyen de ae0 + sa normale)
  aPolygon ae0_2D, ae1_2D;
  bool true_clip(false); // pour distinguer les vrais intersections (celles avec un recouvrement non nul) des fausses (pas d'intersection du tout ou deux polygones en butte)
  bool ae0_crd_2D_computed = false; // booleen pour eviter de calculer plusieurs fois la transfo de ae0 dans son plan moyen
  K_FLD::FloatArray P(3, 3), iP(3, 3); // matrice de passage (et son inverse) pour passer du repere 2D au 3D (et vice versa)

  // donnees pour calculer le COMPLEMENTAIRE de ae0 : ae0 - bits
  std::vector<E_Int> xbit_ids; // ids dans xm des morceaux de ae0 : stockes pour calculer le complementaire de ae0
  // IMPACTED ae0  =  ae0 + points d'intersection sur ses aretes : on ajoute sur ae0 tous les points du contour de tous les morceaux : plus robuste et precis pour calculer le complementaite
  std::set<std::pair<E_Int, E_Int>> edge_to_xnodes;      // pour calculer impacted ae0 : indice local d'arete vers id de point dans xm.crd
  std::vector<std::vector<E_Int>> PGbs; // frontier des paquets de morceaux : il peut y en avoir plusieurs
  

  // structure de donnees declarees ici pour eviter de les instancier dans les boucles : faire un clear est bcp moins couteux
  // pour la fonction get_boundary
  std::vector<E_Int> orient;
  std::set<K_MESH::Edge> w_oe_set;
  std::map<E_Int, E_Int> w_n_map;
 

  E_Int nbcells0 = m0.ncells();
  E_Int nbpts0 = m0.crd.cols();

  ///
  for (E_Int i = 0; i < nbcells0; ++i)
  {
    auto ae0 = m0.aelement(i);
    //int nnodes = ae0.nb_nodes();

    E_Float surf0 = ae0.metrics();
    const E_Float* norm0 = ae0.get_normal();

    // on recupere tous les candidats : ceux dans la bounding box colisionne celle de ae0 
    cands.clear();
    loc1->get_candidates(ae0, ae0.m_crd, cands, 1/*return as 1-based*/, 1.e-2/*RTOL*/);// on etend un peu la boite de ae0 avec RTOL
    if (cands.empty()) continue;

#ifdef SUPERMESH_DBG
    medith::write("cands", m1.crd, m1.cnt, &cands, 1);
    medith::write("ae0", ae0.m_crd, &ae0.m_nodes[0], ae0.m_nodes.size(), ae0.shift());
#endif

    ae0_crd_2D_computed = false;
    edge_to_xnodes.clear();
    xbit_ids.clear();

    // polyclip de chaque candidat
    for (size_t n = 0; n < cands.size(); ++n)
    {
      //
      E_Int i2 = cands[n] - 1;
      auto ae1 = m1.aelement(i2);

      K_CONNECT::IdTool::shift(ae1.m_poids, nbpts0); //pour faire reference aux points de xm plutot que de m1

      E_Float surfc = ae1.metrics();
      const E_Float* normc = ae1.get_normal();

      E_Float ps = NUGA::dot<3>(norm0, normc);
      if (fabs(ps) < 0.9) continue; // doivent �tre colineaire

      bool ae1_reversed = false;
      if (ps < 0.) // CLIP::isolated_clip needs input polygons to have same orientation
      {
        ae1.reverse_orient();
        ps = -ps;
        ae1_reversed = true; // on garde l'info pour 1) retrouver les bons indices locaux d'arete 2) pour faire aussi le reverse sur la version 2D de ae1
      }

      bool check_for_equality = (fabs(surf0 - surfc) < EPSILON) && (ps > 1. - EPSILON); // filtre pour eviter de faire le test complet systematiquement

#ifdef SUPERMESH_DBG
      std::ostringstream o; o << "ae1_" << n;
      medith::write(o.str().c_str(), ae1.m_crd, &ae1.m_nodes[0], ae1.m_nodes.size(), ae1.shift());
#endif

      if (check_for_equality)
      {
        if (ae0 == ae1) // test complet (aPolygon::operator==)
        {
          xm.add(ae0, false/*do capitalize crds*/); // => aucun points nouveaux ajoutes a xm.xrd (puisque pas de points d'intersection)
          anc0.push_back(i);
          anc1.push_back(i2);
          break; // pure match found on a conformal mesh => no more candidate to check
        }
      }

      bits.clear();
     
      aPolygon subj(ae0); // restart from a clean face as isolated_clip might erase subj and we are not in an incremental-boolean loop

      if (proj_on_first)
        // reset cutter (ae1) normal to tell aPolygon specialization of isolated_clip to project on subj instead of cutter
        // "finte" : facon de passer une info a isolated_clip sans passer un argument a isolated_clip et donc sans modifier sa signature
        // car cette info n'est pertinenente que pour le surfacique, et pas le volumique.
        ae1.m_normal[0] = NUGA::FLOAT_MAX;

      int err = NUGA::CLIP::isolated_clip<aelt_t, aelt_t>(subj, ae1, NUGA::INTERSECT::INTERSECTION, ARTOL, bits, true_clip);

      if (err) std::cout << "clipping error : " << i << "-th cell with " << n << "-th candidates" << std::endl;
      assert(!err);

      // check that bits are nearly colinear with subj : fixme : not found yet a better way to detect parasite bits

      for (size_t b = 0; b < bits.size(); ++b)
      {
        const E_Float* normb = bits[b].get_normal();
        E_Float ps = NUGA::dot<3>(norm0, normb);
        if (fabs(ps) < 0.9) // must be nearly colinear
        {
          bits.clear(); 
          true_clip = false;
          break;
        }
      }

      if (true_clip)
      {
        // ajout des morceaux et mise a jour de leur historique
      
        auto ge0s = glob_edge_ids0.get_facets_ptr(i);
        auto ge1s = glob_edge_ids1.get_facets_ptr(i2);

        for (size_t k = 0; k < bits.size(); ++k)
        {
          // on gere au prealable les points d'intersection, ceux avec un poid negatif
          // si le point existe dans key_to_ptxid, on recupere son id
          // sinon on le stocke dans xm.crd et on l'ajoute dans la map
          // on l'ajoute aussi dans edge_to_xnodes pour calculer ae0 impacte

          // la cle est est definie par la paire des indices globaux des aretes
          // on retrouve les indices globaux  avec les indices locaux (le0 de ae0 et le1 de ae1)
          // si ae1 a ete reoriiente, l'indice local doit etre modifi�
         
          // n0 pour conversion de l'indice local d'arete le1
          E_Int n0 = ae1_reversed ? ae1.m_crd.cols() : IDX_NONE;

          __manage_x_points(bits[k], n0, ge0s, ge1s, xm.crd, key_to_ptxid, edge_to_xnodes);
          
          xbit_ids.push_back(anc0.size());
          xm.add(bits[k], false/*do capitalize crds*/);
          anc0.push_back(i);
          anc1.push_back(i2);
	  // std::cout << "anc0 : " << i << " - anc1 : " << i2 << std::endl;

#ifdef SUPERMESH_DBG
          std::ostringstream o;
          o << "bit_" << n << "_" << k;
          medith::write(o.str().c_str(), bits[k].m_crd, &bits[k].m_nodes[0], bits[k].m_nodes.size(), bits[k].shift());
#endif
        }

        continue;
      }

      // CLASSIFICATION : check for fully-in case : is ae0 fully inside ae1 ? or vice & versa ? => go 2D

      if (!ae0_crd_2D_computed)
      {
        ae0_crd_2D_computed = true;
        NUGA::computeAFrame(norm0, P);
        iP = P;
        K_FLD::FloatArray::inverse3(iP);

        ae0_2D = ae0;
        NUGA::transform(ae0_2D.m_crd, iP); // on est maintenant effectivement 2D i.e dans le plan moyen de ae0
      }

      ae1_2D = ae1;
      NUGA::transform(ae1_2D.m_crd, iP); // on passe dans le repere de ae0      

      if (ae1_reversed) // on fait subir le meme sort a la version 2D
        ae1_2D.reverse_orient();

      E_Float Lref2 = ae0.Lref2(); //ae0_2D ???
      E_Float ABSTOL = sqrt(Lref2) * 1.e-2;

      NUGA::eClassify c = NUGA::CLASSIFY::classify2D(ae0_2D, ae1_2D, ABSTOL);
      assert(c != AMBIGUOUS);

      //
      if (c == IN) // ae0 fully inside ae1
      {
        xm.add(ae0, false/*do capitalize crds*/);
        xbit_ids.push_back(anc0.size());
        anc0.push_back(i);
        anc1.push_back(i2);

#ifdef SUPERMESH_DBG
        std::ostringstream o;
        o << "bit__ae0_in_" << n;
        medith::write(o.str().c_str(), ae0.m_crd, &ae0.m_nodes[0], ae0.m_nodes.size(), ae0.shift());
#endif

        break; // ae0 cannot be in several PG
      }
      else if (c == OUT)
      {
        // empty
      }
      else if (c == IN_1) // ae1 fully inside ae0
      {
        if (proj_on_first)
        {
          // project ae1 on ae0
          E_Float zmean = 0;
          for (E_Int u = 0; u < ae0_2D.m_crd.cols(); ++u) zmean += ae0_2D.m_crd(2, u);
          zmean /= ae0_2D.m_crd.cols();
          for (E_Int u = 0; u < ae1_2D.m_crd.cols(); ++u) ae1_2D.m_crd(2, u) = zmean;
          NUGA::transform(ae1_2D.m_crd, P); // back to original ref frame  
          ae1_2D.m_poids = ae1.m_poids;
        }

        xm.add(ae1_2D, false/*do capitalize crds*/);
        xbit_ids.push_back(anc0.size());
        anc0.push_back(i);
        anc1.push_back(i2);

#ifdef SUPERMESH_DBG
        std::ostringstream o;
        o << "bit__" << n << "_in_ae0";
        medith::write(o.str().c_str(), ae1.m_crd, &ae1.m_nodes[0], ae1.m_nodes.size(), ae1.shift());
#endif
      }
    }

    
    // COMPUTE THE RESIDUAL BIT : impacted a0 minus {bits}
    
    // impacted ae0 : update ae0 with refining points
    if (!edge_to_xnodes.empty())
    {
      __impact_ae0(ae0, edge_to_xnodes, xm.crd);
#ifdef SUPERMESH_DBG
      std::ostringstream o;
      o << "impacted_ae0";
      medith::write(o.str().c_str(), ae0.m_crd, &ae0.m_nodes[0], ae0.m_nodes.size(), ae0.shift());
#endif
    }

    // boolean diff ae0 \ bound{bits} : add missing bits if any
 
    orient.clear();
    orient.resize(xbit_ids.size(), 1);
    PGbs.clear();
    int err = K_MESH::Polygon::get_boundary(xm.crd, xm.cnt, xbit_ids/*0 based*/, PGbs, orient, w_oe_set, w_n_map); // can have more than one element in PGbs !
    if (err) continue;

    assert(err == 0);
    assert(!PGbs.empty());

    bits.clear();
    bits.push_back(std::move(ae0));

    std::vector<aelt_t> tmpbits;
    
    // on part de ae0, on fait un diff avec chaque contour de PGbs
    // chacune des ces operation peut donner plusieurs morceaux
    std::vector<E_Int> tmp;
    for (size_t p = 0; (p < PGbs.size()) && !bits.empty(); ++p)
    {
      const auto& PGb = PGbs[p]; /*1-based*/
     
      // construction d'un aPolygon a prtir de la frontiere
      NUGA::aPolygon ae1f(&PGb[0], PGb.size(), 1/*idx_start*/, xm.crd);

      int nbits = bits.size();
      for (int b = 0; b < nbits; ++b)
      {
        if (bits[b].empty()) continue;

        if (bits[b] == ae1f)
        {
          // exact match for a diff => empty answer
          bits.clear();
          break;
        }

#ifdef SUPERMESH_DBG
        medith::write("ae1f", ae1f.m_crd, &ae1f.m_nodes[0], ae1f.m_nodes.size(), 0);
        {
          std::ostringstream o;
          o << "bitcur_" << i << "_" << b;
          medith::write(o.str().c_str(), bits[b].m_crd, &bits[b].m_nodes[0], bits[b].m_nodes.size(), bits[b].shift());
        }

#endif

        tmpbits.clear();
        bool true_clip(false);
        NUGA::CLIP::isolated_clip<aelt_t, aelt_t>(bits[b], ae1f, NUGA::INTERSECT::DIFFERENCE, ARTOL, tmpbits, true_clip);

        if (!true_clip) continue;

        // COMPRESS STRATEGY (with MOVE SEMANTICS) for bits => b can be decremented to treat the replaced at next iter 
        // if tmpbits is empty (IN) => compress (erase bits if single, put the last instead of current otherwise)
        // else replace the first bit, append the other bits  . 
        NUGA::CLIP::__replace_append_and_next_iter(bits, b, tmpbits);
        if ((E_Int)bits.size() != nbits) nbits = bits.size(); //update iff the last have replaced the current 
      }
    }

    if (!bits.empty())
    {
      for (size_t b = 0; b < bits.size(); ++b)
      {
        if (bits[b].empty()) continue;
        //fixme : m_poids might be inconsitent or inexistent for the residual bit => so force to be in 'coordinate appending' mode in add
        xm.add(bits[b], true/*append vertices*/);
        anc0.push_back(i);
	anc1.push_back(-1);
	// std::cout << "elt orphelin : " << i << std::endl;

#ifdef SUPERMESH_DBG
        std::ostringstream o;
        o << "remain_bit_" << i << "_" << b;
        medith::write(o.str().c_str(), bits[b].m_crd, &bits[b].m_nodes[0], bits[b].m_nodes.size(), bits[b].shift());
#endif
      }
      //std::cout << i << " has a remaining bit" << std::endl;
    }
  }

  xm.cnt.updateFacets();

#ifdef SUPERMESH_DBG
  K_FLD::IntArray cnto;
  ngon_type ngo(xm.cnt, false);
  ngo.export_to_array(cnto);
  tp::write("D:\\slandier\\DATA\\tmp\\interpol\\SURF\\tmp\\xm.tp", xm.crd, cnto, "NGON");
#endif

}

template <typename zmesh_t>
void interpol_coeffs_for_first(
  const zmesh_t& m0,
  const zmesh_t& m1, 
  E_Float RTOL,
  std::vector<int>& dindices,
  std::vector<E_Float>& dcoeffs,
  std::vector<int>& xr, bool do_omp=false)
{
  dindices.clear();
  dcoeffs.clear();
  xr.clear();

  using aelt_t = typename zmesh_t::aelt_t;

  auto loc1 = m1.get_localizer();

  m0.get_nodal_metric2();
  m1.get_nodal_metric2();

  E_Int nbcells0 = m0.ncells();
  std::vector<E_Int> cands;

  std::vector<aelt_t> bits;

  std::vector<std::vector<int>>    dindices_per_cell(nbcells0);
  std::vector<std::vector<E_Float>> dcoeffs_per_cell(nbcells0);
  
  size_t n, k;
  E_Int i, i2;

#pragma omp parallel for private(cands, bits, n, k, i, i2) if(do_omp)
  for (i = 0; i < nbcells0; ++i)
  {
    auto ae0 = m0.aelement(i);
    
    cands.clear();
    loc1->get_candidates(ae0, ae0.m_crd, cands, 1, RTOL); //return as 1-based
    if (cands.empty()) continue;

    for (n = 0; n < cands.size(); ++n)
    {
      i2 = cands[n] - 1;
      auto ae1 = m1.aelement(i2);
      bits.clear();
      bool true_clip = false;
      NUGA::CLIP::isolated_clip<aelt_t, aelt_t>(ae0, ae1, NUGA::INTERSECT::INTERSECTION, RTOL, bits, true_clip);

      if (bits.empty()) continue;

      E_Float si2 = 0.;
      for (k = 0; k < bits.size(); ++k)
        si2 += bits[k].extent();

      dindices_per_cell[i].push_back(i2);
      dcoeffs_per_cell[i].push_back(si2);
    }
  }

  // normalize coeffs with total covered surface
  for (size_t i = 0; i < nbcells0; ++i)
  {
    E_Float stot = 0.;
    for (size_t k = 0; k < dcoeffs_per_cell[i].size(); ++k)
      stot += dcoeffs_per_cell[i][k];

    //E_Float s0 = m0.aelement(i).extent();
    //if (stot / s0 < 0.99) std::cout << "cell " << i << " has " << stot / s0 << std::endl;

    for (size_t k = 0; k < dcoeffs_per_cell[i].size(); ++k)
      dcoeffs_per_cell[i][k] /= stot;
  }

  // concatenate info for exit
  xr.push_back(0);
  for (size_t i = 0; i < nbcells0; ++i)
  {
    dindices.insert(dindices.end(), ALL(dindices_per_cell[i]));
    dcoeffs.insert(dcoeffs.end(), ALL(dcoeffs_per_cell[i]));
    xr.push_back(dindices.size());
  }
}

inline E_Float transfer_mass(aPolyhedron<0>& pbit, aPolyhedron<0>& pdon, E_Float Vbit, E_Float rho, 
                             E_Float* gradx=nullptr, E_Float* grady=nullptr, E_Float* gradz = nullptr)
{
  E_Float m = Vbit * rho;

  if (gradx == nullptr) // order1
    return m;

  E_Float gradf[] = { *gradx, *grady, *gradz};
  const E_Float* Gbit = pbit.get_centroid();
  const E_Float* Gdon = pdon.get_centroid();
  E_Float GdonGbit[3];
  NUGA::diff<3>(Gbit, Gdon, GdonGbit);
  m += NUGA::dot<3>(GdonGbit, gradf) * Vbit;
  return m;
}

template <typename zmesh_t>
int interpolate(
  const zmesh_t& mrec,
  const zmesh_t& mdon,
  E_Float RTOL,
  const std::vector<field> & don_fields,
  std::vector<std::vector<E_Float>>& rec_fields,
  bool do_omp = false)
{
  rec_fields.clear();

  int nfields = don_fields.size();

  rec_fields.resize(nfields);
  for (size_t i = 0; i < don_fields.size(); ++i)
    rec_fields[i].resize(mrec.ncells(), 0.);

  std::vector<E_Float> covered_vol(mrec.ncells(), 0.);

  //if (mrec.oriented == 0) return;
  //if (mdon.oriented == 0) return;

  using aelt_t = typename zmesh_t::aelt_t;

  auto loc1 = mdon.get_localizer();

  mrec.get_nodal_metric2();
  mdon.get_nodal_metric2();

  E_Int nbcells0 = mrec.ncells();
  std::vector<E_Int> cands;
  std::vector<aelt_t> bits;

  size_t n, b;
  E_Int i, i2;
  E_Int f;

#pragma omp parallel for private(cands, bits, n, i, i2, f, b) if(do_omp)
  for (i = 0; i < nbcells0; ++i)
  {
    auto ae0 = mrec.aelement(i);

    cands.clear();
    loc1->get_candidates(ae0, ae0.m_crd, cands, 1, RTOL); //return as 1-based
    if (cands.empty()) continue;

#ifdef SUPERMESH_DBG
    medith::write<>("ae0", ae0.m_crd, ae0.m_pgs);
#endif

    for (n = 0; n < cands.size(); ++n)
    {
      i2 = cands[n] - 1;
      auto ae1 = mdon.aelement(i2);


#ifdef SUPERMESH_DBG
      medith::write<>("ae1", ae1.m_crd, ae1.m_pgs);
#endif

      bits.clear();
      bool just_io = !NUGA::CLIP::compute(ae0, ae1, NUGA::INTERSECT::INTERSECTION, bits); //robust_clip returns true if true clip
                                                                                         // IO : current bit does not intersect front (or only along its boundaries)
      if (just_io)
      {
        assert(bits.empty());
        NUGA::eClassify wher = NUGA::CLASSIFY::classify(ae0, ae1, true);

        if (wher == OUT) continue;
        
        if (wher == IN) // ae0 is fully inside ae1
        {
          E_Float Vbit = ae0.extent();
          // transfer fields
          for (f = 0; f < nfields; ++f)
          {
            E_Float* gradf[] = { nullptr, nullptr, nullptr };
            if (don_fields[f].gradf[0] != nullptr)//this field has gradients
            {
              gradf[0] = &don_fields[f].gradf[0][i2];
              gradf[1] = &don_fields[f].gradf[1][i2];
              gradf[2] = &don_fields[f].gradf[2][i2];
            }

            rec_fields[f][i] += transfer_mass(ae0, ae1, Vbit, don_fields[f].f[i2], gradf[0], gradf[1], gradf[2]);
              covered_vol[i] += Vbit;
          }
          break; // fully IN, so stop checking other candidates
        }
        else // IN_1 : ae1 is fully inside ae0
        {
          //gradients doesnt matter here since Vbit = Vdonnor
          E_Float Vbit = ae1.extent();
          // transfer fields
          for (f = 0; f < nfields; ++f)
          {
            rec_fields[f][i] += transfer_mass(ae0, ae1, Vbit, don_fields[f].f[i2]);
            covered_vol[i] += Vbit;
          }
        }
      }
      else // true clip
      {
        for (b = 0; b < bits.size(); ++b)
        {
          E_Float Vbit = bits[b].extent();
          // transfer fields
          for (f = 0; f < nfields; ++f)
          {
            E_Float* gradf[] = { nullptr, nullptr, nullptr };
            if (don_fields[f].gradf[0] != nullptr)//this field has gradients
            {
              gradf[0] = &don_fields[f].gradf[0][i2];
              gradf[1] = &don_fields[f].gradf[1][i2];
              gradf[2] = &don_fields[f].gradf[2][i2];
            }

            rec_fields[f][i] += transfer_mass(bits[b], ae1, Vbit, don_fields[f].f[i2], gradf[0], gradf[1], gradf[2]);
            covered_vol[i] += Vbit;
          }
        }
      }
    }
  }
  
#pragma omp parallel for private(i) if (do_omp)
  for (E_Int i = 0; i < nbcells0; ++i)
  {
    for (E_Int f = 0; f < nfields; ++f)
    {
      if (covered_vol[i] != 0.) // discard non interpolated cells
        rec_fields[f][i] /= covered_vol[i];
    }
  }
  return 0;
}


}

#endif // NUGA_SUPERMESHHXX
