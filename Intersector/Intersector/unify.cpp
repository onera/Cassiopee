
/*    
    Copyright 2013-2019 Onera.

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

#include "intersector.h"
#include <vector>
#include "Fld/DynArray.h"
#include "Fld/ngon_t.hxx"
#include "MeshElement/Hexahedron.h"
#include "MeshElement/Polyhedron.h"
#include "Search/BbTree.h"
#include "Nuga/include/localizer.hxx"
#include "Nuga/Delaunay/Triangulator.h"
#include "Nuga/include/ph_clipper.hxx"
#include "Nuga/include/polyhedron.hxx"

//#define FLAG_STEP
//#define OUTPUT_XCELLN

#ifdef FLAG_STEP
#include "chrono.h"
int chrono::verbose = 2;
#endif

#ifdef NETBEANSZ
//#define DEBUG_UNIFY
#endif

#ifdef OUTPUT_XCELLN
#include <sstream>
#include "Nuga/include/medit.hxx"
#endif

#ifdef DEBUG_UNIFY
#include "Nuga/include/medit.hxx"
#include "Nuga/Boolean/NGON_debug.h"
using NGDBG  = NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>;
#endif

enum eCellType { UNSET=-9, COLLIDE=-7, AMBIGUOUS=-1, IN_BODY=-5, IN=0, OUT=1};

using ngon_type = ngon_t<K_FLD::IntArray>;

using ELT1 = K_MESH::Polyhedron<UNKNOWN>;
using aELT1 = NUGA::haPolyhedron<UNKNOWN>;
using PG1_t = K_MESH::Polygon;//typename ELT1::boundary_type;
using ELT2 = K_MESH::Polyhedron<UNKNOWN>;
using aELT2 = NUGA::haPolyhedron<UNKNOWN>;
using PG2_t = K_MESH::Polygon;//typename ELT2::boundary_type;
using tree_t = K_SEARCH::BbTree3D;
using bbox_t = K_SEARCH::BBox3D;
using acrd_t = K_FLD::ArrayAccessor<K_FLD::FloatArray>;
using acnt_t = K_FLD::ArrayAccessor<K_FLD::IntArray>;
using loc_t = NUGA::localizer<tree_t, acrd_t, acnt_t>;

using prior_t = std::map< E_Int, std::vector<E_Int>>;

bool has_flag(ngon_unit& pgs, E_Int flag)
{
  bool ret = false;
  E_Int nbf = pgs.size();
  for (E_Int f=0; f < nbf; ++f)
  {
    if (pgs._type[f] == flag)
      return true;
  }
  return false;   
}

bool PH_has_flag(ngon_type& ng, E_Int PHi, E_Int flag)
{
  E_Int nbf = ng.PHs.stride(PHi);
  const E_Int* faces = ng.PHs.get_facets_ptr(PHi);
  for (E_Int f=0; f < nbf; ++f)
  {
    if (ng.PGs._type[faces[f]-1] == flag)
      return true;
  }
  return false;   
  
}

struct aIsLessThanb : std::binary_function<int, int, bool>
{
  
  aIsLessThanb(const prior_t p):_p(p){}
  
  bool operator()(int a, int b) const { 
    
    auto it = _p.find(a);
    if (it == _p.end()) return false;
    
    for (size_t i=0; i < it->second.size(); ++i)
      if (it->second[i] == b) return true;
    
    return false;
  }
  
  prior_t _p;
};
   
void comp_priorities(std::vector<std::pair<E_Int, E_Int>> & priority,  std::map<E_Int, std::vector<E_Int>> & prior_per_comp)
{
  prior_per_comp.clear();
  
  std::set<std::pair<E_Int, E_Int>> upairs(priority.begin(), priority.end());
  priority.clear();
  priority.insert(priority.end(), upairs.begin(), upairs.end());
    
  for (size_t i = 0; i < priority.size(); ++i)
  {
    prior_per_comp[priority[i].first].push_back(priority[i].second);
  }
  
  aIsLessThanb pred(prior_per_comp);
  for (auto it = prior_per_comp.begin(); it != prior_per_comp.end(); ++it)
    std::sort(ALL(it->second), pred);
}

eCellType is_inside(const E_Float* pt, const K_FLD::FloatArray& crdS, const E_Float* ptS, const tree_t* treeS, const ngon_unit& PGsS)
{
  E_Int Ti[3], tx;
  E_Float tol(E_EPSILON), tol2(E_EPSILON*E_EPSILON), d2;

  std::vector<E_Int> candidates;
  std::set<E_Int> valid_cands;
  treeS->getIntersectingBoxes(pt, ptS, candidates, E_EPSILON, true);//segment intersection

  size_t bx_sz=candidates.size();

#ifdef DEBUG_MASK
  assert (!candidates.empty());
  K_FLD::FloatArray crd;
  crd.pushBack(pt,pt+3);
  crd.pushBack(ptS, ptS+3);
  K_FLD::IntArray cnt(2,1,0); cnt(1,0)=1;
  
  std::ostringstream o;
  o << "/home/slandier/projects/UNIFY/2cubs/tmp/ray.mesh";
  medith::write(o.str().c_str(), crd, cnt, "BAR");
  {
    E_Int Tii[3];
    K_FLD::IntArray cnt2;
    for (size_t i=0; i < bx_sz; ++i)
    {
      const E_Int* nodes = PGsS.get_facets_ptr(candidates[i]);
      E_Int nb_nodes = PGsS.stride(candidates[i]);

      E_Int Tii[] = { nodes[0]-1, nodes[1]-1, nodes[2]-1};

      cnt2.pushBack(Tii, Tii+3);
    }
    //const K_FLD::FloatArray crd2(crdS);
    o.str("");
    o << "/home/slandier/projects/UNIFY/2cubs/tmp/cands.mesh";
    medith::write(o.str().c_str(), crdS, cnt2, "TRI");
  }
#endif

  E_Float /*Q0[3], Q1[3], Q2[3],*/ u0(K_CONST::E_MAX_FLOAT), u1, mindp(K_CONST::E_MAX_FLOAT), ptA[3];
  E_Bool overlap, tol_is_abs(1), coincident;
  size_t nb_visibles(0);

  K_FUNC::diff<3>(pt, ptS, ptA);
  K_FUNC::normalize<3>(ptA);
  E_Int T[3];

  acrd_t acrdS(crdS);
  DELAUNAY::Triangulator dt;
  
  for (size_t b=0; b < bx_sz; ++b)
  {     
    E_Int PGi = candidates[b];

    const E_Int* nodes = PGsS.get_facets_ptr(PGi);
    E_Int nb_nodes = PGsS.stride(PGi);
    
    K_MESH::Polygon PG(nodes, nb_nodes, -1/*to make it zero based*/);
    
    PG.triangulate(dt, acrdS);
    
    E_Int nb_tris = PG.nb_tris();
    for (E_Int t = 0; t < nb_tris; ++t)
    {
      PG.triangle(t, T); //watchme : base ?
      const E_Float* Q0 = crdS.col(T[0]);
      const E_Float* Q1 = crdS.col(T[1]);
      const E_Float* Q2 = crdS.col(T[2]);

      { 
        if (!K_MESH::Triangle::intersectv2<3>(Q0, Q1, Q2, pt, ptS, tol, tol_is_abs,  u0, u1, tx, overlap, coincident))
          continue;

        if (IS_ZERO(u0, E_EPSILON)) // lying on the mask surface
          return AMBIGUOUS;

        if (overlap || coincident) //ambiguous : overlap (but pt is not inside - otherwise u0=0. -) or coincident : touching a triangle's border
          return AMBIGUOUS;
      }
      // u0 has necessarily a value here.
      if (u0 < mindp)
      {
        //nb_visibles=1;
        //candidates[0]=candidates[b]; //overwrite candidates to have real candidates between rank 0 and (nb_visibles-1)
        valid_cands.clear();
        valid_cands.insert(candidates[b]);
        mindp=u0;
      }
      else if (IS_ZERO(u0-mindp, E_EPSILON))
      {
        valid_cands.insert(candidates[b]);
        //candidates[nb_visibles++]=candidates[b];
      }
    }
  }

  nb_visibles = valid_cands.size();

  if (nb_visibles == 0)
    return AMBIGUOUS;

  auto itPGi = valid_cands.begin();//[0];
  E_Int PGi = *itPGi;
  const E_Int* nodes = PGsS.get_facets_ptr(PGi);
  E_Int nb_nodes = PGsS.stride(PGi);
  const E_Float* Q0 = crdS.col(nodes[0]-1);
  const E_Float* Q1 = crdS.col(nodes[1]-1);
  const E_Float* Q2 = crdS.col(nodes[2]-1);
  E_Float normal[3];
  K_MESH::Triangle::normal(Q0,Q1,Q2, normal);
  E_Float ps = K_FUNC::dot<3>(normal, ptA);
  bool is_out= (ps > 0.); // true iff OUT

  if (nb_visibles > 1)
  {
    //check if all canddidates give the same orientation
    itPGi++;
    for (auto it = itPGi; it != valid_cands.end(); ++it)//size_t i=1; i < nb_visibles; ++i)
    {
      E_Int PGi = *(it);//valid_cands[i];
      const E_Int* nodes = PGsS.get_facets_ptr(PGi);
      E_Int nb_nodes = PGsS.stride(PGi);
      const E_Float* Q0 = crdS.col(nodes[0]-1);
      const E_Float* Q1 = crdS.col(nodes[1]-1);
      const E_Float* Q2 = crdS.col(nodes[2]-1);
      E_Float normal[3];
      K_MESH::Triangle::normal(Q0,Q1,Q2, normal);

      if ((K_FUNC::dot<3>(normal, ptA) > 0.) != is_out)
        return AMBIGUOUS;
    }
  }

  eCellType ret = is_out ? OUT : IN; // ASSUME WALLS ARE INWARD BODIES (i.e each fluid component is point outward)
  return ret;
}

eCellType cell_is_inside_mask(const K_FLD::FloatArray& crd1,const E_Int* unodes1, E_Int nb_nods1, E_Int idx_start, const K_FLD::FloatArray& mask_crd, const ngon_unit& mask_pgs, const tree_t* tree)
{
  eCellType ret = AMBIGUOUS; //-1 : cannot tell ; +1 is fully outside; 0 : is fully in

  E_Int nb_pgs_mask = mask_pgs.size();
  E_Float Gmask[3];
  ///
  for (E_Int n=0; n < nb_pgs_mask && (ret == AMBIGUOUS); ++n)
  {
    K_MESH::Polygon::centroid<3>(mask_crd, mask_pgs.get_facets_ptr(n), 3, 1, Gmask);
    const E_Float* pt = crd1.col(unodes1[n%nb_nods1] - idx_start); // try different nodes on the input cell
    ret = is_inside(pt, mask_crd, Gmask, tree, mask_pgs);
  }
  return ret;
}

eCellType cell_is_inside_mask(const aELT1& ph, const K_FLD::FloatArray& mask_crd, const ngon_unit& mask_pgs, const tree_t* tree)
{
  eCellType ret = AMBIGUOUS; //-1 : cannot tell ; +1 is fully outside; 0 : is fully in

  E_Int nb_pgs_mask = mask_pgs.size();
  E_Float Gmask[3];
  ///
  for (E_Int n=0; n < nb_pgs_mask && (ret == AMBIGUOUS); ++n)
  {
    K_MESH::Polygon::centroid<3>(mask_crd, mask_pgs.get_facets_ptr(n), 3, 1, Gmask);
    const E_Float* pt = ph.m_crd.col(n%ph.m_crd.cols()); // try different nodes on the input cell
    ret = is_inside(pt, mask_crd, Gmask, tree, mask_pgs);
  }
  return ret;
}

void init(const std::vector<K_FLD::FloatArray*> &crds, const std::vector<K_FLD::IntArray*>& cnts,
          const std::vector<K_FLD::FloatArray*> &mask_crds, const std::vector<K_FLD::IntArray*>& mask_cnts,
          std::vector<ngon_type*>& ngs, std::vector<acrd_t*>& acrds, std::vector<ngon_unit*>& pgfs, 
          const std::vector<std::vector<E_Int>> & base_wall_ids,
          std::vector<loc_t*>& locs)
{
  E_Int nb_zones = crds.size();
  assert (nb_zones == cnts.size());
  E_Int nb_comps = mask_crds.size();
  assert (nb_comps == mask_cnts.size());
    
  ngs.clear();
  ngs.resize(nb_zones);
  for (size_t i=0; i < nb_zones; ++i)ngs[i] = new ngon_type(*cnts[i]);

  acrds.clear();
  acrds.resize(nb_zones);
  for (size_t i=0; i < nb_zones; ++i)acrds[i] = new acrd_t(*crds[i]);
  
  pgfs.clear();
  pgfs.resize(nb_comps);

  for (size_t i=0; i < nb_comps; ++i)
  {
    pgfs[i] = new ngon_unit(mask_cnts[i]->begin());

    if (!base_wall_ids[i].empty()) pgfs[i]->_type.resize(pgfs[i]->size(), UNSET);

    for (size_t j=0; j < base_wall_ids[i].size(); ++j)
    {
      E_Int wi = base_wall_ids[i][j] - 1;
      //std::cout << "wi : " << wi << std::endl;
      assert(wi < pgfs[i]->_type.size());
      pgfs[i]->_type[wi] = BCWALL;
    }
  }

  locs.clear();
  locs.resize(nb_comps);
  for (size_t i=0; i < nb_comps; ++i)
  {
    tree_t* t = new tree_t(*mask_crds[i], *pgfs[i]);  //polygons !
    locs[i]   = new loc_t(t, E_EPSILON/*fixme*/); //not a leak : locs[i] owes t
  }
}


bool compute_zone_blanking(const K_FLD::FloatArray& crd, ngon_type& ng,
                           loc_t** mask_locs, ngon_unit** mask_skins, K_FLD::FloatArray** mask_crds, E_Int nb_masks, 
                           std::vector<E_Float>& xcelln)
{
  // per-mask loop / per-cell loop
  E_Int nb_phs(ng.PHs.size());

#ifdef DEBUG_UNIFY
  E_Int badPH = 3;
#endif
  
  xcelln.clear();
  xcelln.resize(nb_phs, UNSET);
  
  std::vector<E_Int> candsp;  
  acrd_t acrd1(crd);
  
  DELAUNAY::Triangulator dt;
  
  std::vector<E_Float> xcelln_cur;
  
  for (E_Int m=0; (m < nb_masks) ; ++m)
  {
#ifdef FLAG_STEP
    std::cout << "zone blanking for mask : " << m << " over " << nb_masks << std::endl;
#endif
    
    bool collides = false;
    const loc_t* locp = mask_locs[m];
    ngon_unit* skinp = mask_skins[m];
    const K_FLD::FloatArray& crdp = *mask_crds[m];
    acrd_t acrd2(crdp);
    
    xcelln_cur.clear();
    xcelln_cur.resize(nb_phs, UNSET);
          
    for (E_Int i = 0; i < nb_phs; ++i)
    {     
//      if (xcelln[i] == 0.) {
//        //xcelln_cur[i] = 0.;
//        continue; //fully inside a previous mask
//      }
      
      E_Int* faces1 = ng.PHs.get_facets_ptr(i);
      E_Int nb_faces1 = ng.PHs.stride(i);
      
#ifdef DEBUG_UNIFY
     //NGDBG::draw_PH("ecur.plt", crd, ng, i);
#endif

      ELT1 e1(&ng.PGs, faces1, nb_faces1); // global PH1
      
      candsp.clear();
      locp->get_candidates(e1, acrd1, candsp); //0-BASED
      if (candsp.empty()) continue;
      
#ifdef DEBUG_UNIFY
     //NGDBG::draw_PGs("shape", crdp, *skinp, candsp);
#endif
      
      E_Int nb_faces2 = candsp.size();
      K_CONNECT::IdTool::shift(candsp, 1);
      ELT2 e2(skinp, &candsp[0], nb_faces2);
      E_Int it1, it2;
      bool is_x = NUGA::COLLIDE::simplicial_colliding<acrd_t, 3>(acrd1, e1, acrd2, e2, K_MESH::Triangle::fast_intersectT3<3>, -1., it1, it2); //do not pass the tol as we use here a predicate
      
      if (is_x)
      {
        xcelln_cur[i] = COLLIDE;
        collides = true;
      }
    }
    
    if (!collides) continue;
  
    //coloring by neighborhood
    ngon_unit neighbors;
    ng.build_ph_neighborhood(neighbors);
    E_Float FIRST_COL = 1.;
    K_CONNECT::EltAlgo<K_MESH::Polyhedron<0>/*dummy*/>::coloring(neighbors, xcelln_cur, (E_Float)UNSET, (E_Float)FIRST_COL);
    
    E_Int nb_bits;
    {
      std::set<E_Int> cols(xcelln_cur.begin(), xcelln_cur.end());//fixme : hpc
      nb_bits = cols.size();
    }
    
#ifdef DEBUG_UNIFY
    {
      std::vector<ngon_type> ngbits(nb_bits);
      for (size_t b=0; b < nb_bits;++b)
      {
        ngbits[b].PGs = ng.PGs;
      }

      for (E_Int i =0; i <  ng.PHs.size(); ++i)
      {
        if (xcelln_cur[i] == COLLIDE) ngbits[0].PHs.add(ng.PHs.stride(i), ng.PHs.get_facets_ptr(i));
        else /*if (xcelln_cur[i] == 0.)*/ ngbits[E_Int(xcelln_cur[i])].PHs.add(ng.PHs.stride(i), ng.PHs.get_facets_ptr(i));
      }
      
      for (size_t b = 0; b < nb_bits; ++b)
      {
        K_FLD::IntArray cnt;
        ngbits[b].export_to_array(cnt);
        std::ostringstream o;
        o << "bit_" << b << ".mesh"; 
        medith::write(o.str().c_str(), crd, cnt, "NGON");
      }
    }
#endif
    
    //std::cout << "which-is-zero phase (nb bits: " <<  nb_bits << ")" << std::endl;

    // detect colors among 1 to n that are inside masks (and therefore are zero color)
    STACK_ARRAY(E_Int, nb_bits, is_zero);
    for(size_t u=0; u< nb_bits; ++u)is_zero[u]=E_IDX_NONE;
    E_Int remains_zones(nb_bits-1);//-1 : discard collide color of the test
    //
    for (E_Int i = 0; (i < nb_phs) && (remains_zones> 0); ++i)
    {
      E_Int nb_pgs = ng.PHs.stride(i);
      const E_Int* faces = ng.PHs.get_facets_ptr(i);
      
      E_Int col = xcelln_cur[i];
      
      if (col == COLLIDE || col == IN) continue;
      if (is_zero[col] != E_IDX_NONE) continue;
      
      //debugg=0;
            
      const E_Int* neighs = neighbors.get_facets_ptr(i);
      
      for (size_t v=0; (v<nb_pgs); ++v)
      {
        E_Int Fi = faces[v] - 1;
        E_Int PHn = neighs[v]; //0-based
        if (PHn == E_IDX_NONE) continue;
        if (xcelln_cur[PHn] != COLLIDE) continue;
        
        //do the classification test of each bit regarding mask m          
        E_Int ret = cell_is_inside_mask(crd, ng.PGs.get_facets_ptr(Fi), ng.PGs.stride(Fi), 1, crdp, *skinp, locp->get_tree());
        if (ret == AMBIGUOUS) continue;

        is_zero[col] = (ret == IN);
        --remains_zones;
        break;
      }
    }

    // reset non-colliding and non-zero colors to UNSET
    for (size_t i=0; i < xcelln_cur.size(); ++i)
    {
      E_Float& xc = xcelln_cur[i];
      if (xc == COLLIDE || xc == IN) continue;
      xc = (is_zero[xc] == 1) ? IN : UNSET;
    }
    
    // now update the output mask : IN (0) > COLLIDE(-7) > UNSET(-9)
    for (size_t i=0; i < xcelln.size(); ++i)
      xcelln[i] = std::max(xcelln[i], xcelln_cur[i]);
    
#ifdef DEBUG_UNIFY
    {      
      ngon_type ng1,ng2;
      ng1.PGs = ng2.PGs = ng.PGs;
      for (E_Int i =0; i <  ng.PHs.size(); ++i)
      {
        if (xcelln[i] == COLLIDE) ng1.PHs.add(ng.PHs.stride(i), ng.PHs.get_facets_ptr(i));
        else if (xcelln[i] != UNSET) ng2.PHs.add(ng.PHs.stride(i), ng.PHs.get_facets_ptr(i));
      }
      
      K_FLD::IntArray cnt;
      ng1.export_to_array(cnt);
      std::ostringstream o;
      o << "update_0.mesh"; 
      medith::write(o.str().c_str(), crd, cnt, "NGON");
      ng2.export_to_array(cnt);
      o.str("");
      o << "update_x.mesh"; 
      medith::write(o.str().c_str(), crd, cnt, "NGON");
    }
#endif

  } // masks loop

  return true;
}

void build_current(ngon_type& ng, E_Int i, const K_FLD::FloatArray& crd, aELT1& aPHcur, E_Float& V0, bool is_inward1, bool make_inward1)
{
  aPHcur.set(ng,i, crd); // self-contained with pt histo
  // make orientation consistent
  aPHcur.reorient();

  // Compute initial volume
  E_Int err = aPHcur.volume<DELAUNAY::Triangulator>(V0, false/*need reorient*/);

  is_inward1=false;
  if (V0 < 0.) //point inward
  {
    is_inward1 = true;
    V0 = -V0;
  }

  //classification based on operation
  if (make_inward1 && !is_inward1)
  {
    aPHcur.reverse();
    is_inward1 = true;
  }
  
#ifdef DEBUG_UNIFY
     //NGDBG::draw_PGs("ecur", aPHcur.m_crd, aPHcur.m_pgs);
#endif
}

void build_cutter(ngon_unit* skinp, std::vector<E_Int>& candsp, const K_FLD::FloatArray& crdp, bool make_inward, aELT1& aPHcutter)
{
#ifdef DEBUG_UNIFY
  //NGDBG::draw_PGs("shapeG", crdp,  *skinp, candsp);
#endif

  // gather the polygons as a single shape
  K_CONNECT::IdTool::shift(candsp, 1);//because it was 0-based
  ELT2 gPH2(skinp, &candsp[0], candsp.size());

  aPHcutter.set(gPH2, crdp); // self-contained with pt histo
  
  //WARNING : supposed outward at init (reorient skins) => nothing to do if DIFF or UNION (duplication)
   if (make_inward) aPHcutter.reverse();
      
#ifdef DEBUG_UNIFY
  NGDBG::draw_PGs("cutter", aPHcutter.m_crd, aPHcutter.m_pgs);    
#endif
}

E_Int robust_clip(acrd_t& acrdcur, aELT1& aPHcur, bool is_inward1, acrd_t& acrdcut, aELT2& aPHcutter, aELT1& aPHclip, bool dbg)
{
  E_Float PS_MIN0 = ::cos(0.1*K_CONST::E_1_PI/180.); // 0.1°
  
  E_Float RTOL=1.e-6;
  E_Float PS_MIN = PS_MIN0;
  E_Int err(0), contact; 

  do
  {
    if (dbg) std::cout << "isolated_clip : PS_MIN/RTOL : " << PS_MIN << "/" << RTOL << std::endl;

    contact = 0;
    err = NUGA::INTERSECT::isolated_clip(acrdcur, aPHcur, is_inward1, acrdcut, aPHcutter, PS_MIN, RTOL, aPHclip, contact, dbg);

    RTOL *=10;
    PS_MIN -=0.01;   
  }
  while (err && RTOL <= 1.e-1);
  
  return err;
  
}

bool update_volume(aELT1& aPHcur, const K_FLD::FloatArray&crdp, ngon_unit* skinp, const loc_t* locp, aELT1& aPHclip, E_Float& vcur, bool classify_geom)
{
  bool true_clip(true);
  E_Float v;
  if (!classify_geom) // non_empy polyclip answer
  {
    E_Int er = aPHclip.volume<DELAUNAY::Triangulator>(v, false/*need reorient*/);
    v = (v < 0.) ? -v : v;
    classify_geom = (v/vcur >= 1. - E_EPSILON); // if the volume is the same => do the test in/out
  }

  if (classify_geom)
  {        
    eCellType ret = cell_is_inside_mask(aPHcur, crdp, *skinp, locp->get_tree());
    assert (ret != AMBIGUOUS);
    v = (ret == IN) ? 0. : vcur;
  }
  
  true_clip = !classify_geom;
  if (true_clip) assert (v < vcur && 0. < v);
  
  vcur=v;
  
  return true_clip;
}

bool incremental_clip
(ngon_unit* skinp, const K_FLD::FloatArray& crdp, const loc_t* locp, bool is_inward1, bool make_inward2, aELT1& aPHclip, aELT2& aPHcutter, aELT1& aPHcur, E_Float& vcur, bool dbg)
{     
  acrd_t acrdcur(aPHcur.m_crd); //necessary to be in the loop as crd change
  
  std::vector<E_Int> candsp, poids;
  locp->get_candidates(aPHcur, acrdcur, candsp); //0-BASED

  bool classify_geom = candsp.empty(); // no collision => fully in or out ?

  E_Int err(0);
  if (!classify_geom)
  {
    build_cutter(skinp, candsp, crdp, make_inward2, aPHcutter);
    acrd_t acrdcut(aPHcutter.m_crd);
    err = robust_clip(acrdcur, aPHcur, is_inward1, acrdcut, aPHcutter, aPHclip, dbg);
    if (err)
    {
      vcur=0.;
      return true;
    }
    classify_geom = (aPHclip.m_crd.cols() == 0); //empty answer
  }

  //update the current volume and tells if it was a true clipping
  bool true_clip = update_volume(aPHcur, crdp, skinp, locp, aPHclip, vcur, classify_geom);

  if (true_clip)
  {
    if (!aPHcur.poids.empty())
      poids = aPHcur.poids; //sauvegrade avant excrasement (todo : spécialiser operator==

    aPHcur = aPHclip;// next patch. checkme : required ? 

    if (!poids.empty())
    {
      poids.resize(aPHclip.m_crd.cols(), E_IDX_NONE);
      for (size_t u=0; u < aPHclip.m_crd.cols(); ++u)
        if (aPHclip.poids[u] != E_IDX_NONE) poids[u] = poids[aPHclip.poids[u]];
      aPHcur.poids = poids;
    }

#ifdef DEBUG_UNIFY
  NGDBG::draw_PGs("lastcut", aPHcur.m_crd, aPHcur.m_pgs);
#endif
  }
  
  return false;
}



E_Float compute_coeff(aELT1& aPHcur, E_Int nb_pts0, E_Float vcur, E_Float V0, std::vector<bool>& in_points)
{
  E_Float coeff = vcur/V0;
  E_Float TOL = 1.e-6;//1.e-4; //fixme
  coeff = (coeff > 1. - TOL) ? OUT : coeff;
  coeff = (coeff < TOL) ? IN : coeff;

  // check if the PH after the complete polyclip doesn't contain any mask wall (other wall are not taken into account )
  // if it does, flag any point from initial cell as "in"
  if (aPHcur.m_pgs._type.empty())
    return coeff;
  
  bool has_in_body_bits = false;

  for (E_Int f=0; f < aPHcur.nb_faces(); ++f)
  {
    if (aPHcur.m_pgs._type[f] == BCWALL)
    {
      has_in_body_bits = true;
      break;
    }
  }
  
  if (!has_in_body_bits)
    return coeff;

  // multi-bit, check which one is in-body. if some, modify the coeff, if all, coeff := IN_BODY

  ngon_unit neighbors;
  std::vector<E_Int> colors;
  K_MESH::Polygon::build_pg_neighborhood(aPHcur.m_pgs, neighbors);
  K_CONNECT::EltAlgo<K_MESH::Polygon>::coloring (neighbors, colors);
  E_Int nb_connex = 1+*std::max_element(colors.begin(), colors.end());

  std::vector<bool> is_wall_free(nb_connex, true);

  if (nb_connex == 1)
  {
    coeff = IN_BODY;
    is_wall_free[0]=false;
  }
  else //detect any IN_BODY bit
  {
    for (E_Int f=0; f < aPHcur.nb_faces(); ++f)
    {
      if (aPHcur.m_pgs._type[f] == BCWALL)
       is_wall_free[colors[f]] = false;
    }

    bool has_good=false;
    for (E_Int c=0; c < nb_connex; ++c)
      has_good |= (is_wall_free[c] == true);

    if (!has_good)
      coeff = IN_BODY;
    else // compute the coef based on OUT pieces
    {
      std::vector<E_Int> facs;
      for (size_t u=0; u < aPHcur._nb_faces; ++u)
        if (is_wall_free[colors[u]])facs.push_back(aPHcur._faces[u]);

      E_Float v=0.;
      aPHcur.volume<DELAUNAY::Triangulator>(v, false/*need reorient*/);
      v = (v < 0.) ? -v : v;
      coeff = v/V0;
      coeff = (coeff > 1.) ? 1. : coeff;
      coeff = (coeff < 0.) ? 0. : coeff;
    }
  }

  // for NON AMBIGUOUS inside bits, flag orignal node to detect false OUT
  const ngon_unit* pgs1 = &aPHcur.m_pgs;
  bool ambiguous = false;
  for (E_Int u=0; (u < aPHcur._nb_faces) && !ambiguous; ++u)
  {
    if (is_wall_free[colors[u]]) continue;
    if (aPHcur.m_pgs._type[u] == WALL1)
    {
      ambiguous = true; break;
    }
  }

  //if (!ambiguous) std::cout << "good catch!!" << std::endl;

  for (E_Int u=0; (u < aPHcur._nb_faces) && !ambiguous; ++u)
  {
    if (is_wall_free[colors[u]]) continue;

    E_Int Fi = aPHcur._faces[u]-1;
    E_Int nbn = pgs1->stride(Fi);
    const E_Int * nodes = pgs1->get_facets_ptr(Fi);
    for (E_Int n=0; n < nbn; ++n)
    {
      E_Int Ni = nodes[n] - 1;
      E_Int lid = aPHcur.poids[Ni];
      if (lid >= nb_pts0) 
      {
        //std::cout << lid << " is greater than " << nb_pts1 << std::endl;
        continue;
      }

      //std::cout << "CAUGHT POINT : " << lid << std::endl;
      E_Int gpid = aPHcur.poids[lid];
      in_points[gpid] = true;
    }
  }
  
  return coeff;
}

void compute_zone_overlaps(const K_FLD::FloatArray& crd, ngon_type& ng, loc_t** mask_locs, ngon_unit** mask_skins, K_FLD::FloatArray** mask_crds, E_Int nb_masks, 
                           std::vector<E_Float>& xcelln, K_FLD::FloatArray& crdo, ngon_type& ngo, std::vector<bool>& in_points, bool make_inward1, bool make_inward2)
{
  // per-cell loop / per-mask loop
  
#ifdef DEBUG_UNIFY
    E_Int badPH = 102212;//102212;//618798;630963
    E_Int badZ = 1;
#endif
  
  E_Int nb_phs = ng.PHs.size();
    
  in_points.clear();
  in_points.resize(crd.cols(), false);
  
  aELT1 aPHcur, aPHclip, aPHcutter; // self_contained : current, clipped, open-cutter-cell
  
  ///  
  for (E_Int i = 0; i < nb_phs; ++i)
  {

#ifdef DEBUG_UNIFY
    //if (i != badPH) continue;
    //std::cout << "PH : " << i << " over " << nb_phs << std::endl;
    //std::cout << "xcelln[PH] : " << xcelln[i] << std::endl;
#endif

    if (xcelln[i] != COLLIDE) continue;
    
    E_Float V0;
    bool is_inward1;
    build_current(ng, i, crd, aPHcur, V0, is_inward1, make_inward1);
        
    E_Int nb_pts1 = aPHcur.m_crd.cols();

    bool dbg = false;
#ifdef DEBUG_UNIFY
    dbg = (i == badPH /*&& Z == badZ*/);
#endif

    E_Float vcur = V0;
    /// incremental clipping
    for (E_Int m=0; (m < nb_masks) && (vcur > E_EPSILON); ++m)  // multi-overlap loop : by decsreasing priority (sorted by the caller)
    {
      const loc_t* locp = mask_locs[m];
      ngon_unit* skinp = mask_skins[m];
      const K_FLD::FloatArray& crdp = *mask_crds[m];

      bool err = incremental_clip(skinp, crdp, locp, is_inward1, make_inward2, aPHclip, aPHcutter, aPHcur, vcur, dbg);

#ifdef DEBUG_UNIFY
      if (err) std::cout << "ERROR for PHS : " << i /*<< " in zone : " << z */<< std::endl;
#endif
    }
    
    // compute the coeff
    E_Float coeff = compute_coeff(aPHcur, nb_pts1, vcur, V0, in_points);
    
    if (coeff == OUT) coeff = UNSET; // reset any supposedly VISIBLE at this stage. might be a false and preturbate the logic for coloring IN_BODY zones.

    xcelln[i] = coeff;

#ifdef OUTPUT_XCELLN

    if (coeff <= 0.) continue; // empty intersection or error

    E_Int nb_pts = crdo.cols();

    aPHcur.m_pgs.shift(nb_pts);
    E_Int nb_pgs = ngo.PGs.size(); 
    ngo.PGs.append(aPHcur.m_pgs);
    assert (aPHcur.m_pgs.size() == aPHcur._nb_faces);

    std::vector<E_Int> fs(aPHcur._nb_faces);
    for (size_t k=0; k < aPHcur._nb_faces; ++k)
      fs[k] = aPHcur._faces[k] + nb_pgs;

    crdo.pushBack(aPHcur.m_crd);
    ngo.PHs.add(aPHcur._nb_faces, &fs[0]);
    ngo.PHs.updateFacets();
    
    //std::cout << "coeff pour PH : " << i << " : " << coeff << std::endl;

#endif

  } // xcelln loop  
}

///
void MOVLP_unify(const std::vector<K_FLD::FloatArray*> &crds, const std::vector<K_FLD::IntArray*>& cnts, const std::vector<E_Int>& comp_id, 
                 std::vector<std::pair<E_Int, E_Int>> & priority, 
                 const std::vector<K_FLD::FloatArray*> &mask_crds, const std::vector<K_FLD::IntArray*>& mask_cnts,
                 const std::vector< std::vector<E_Int> > &mask_wall_ids, 
                 std::vector< std::vector<E_Float>>& xcelln)
{
  NUGA::INTERSECT::eOPER oper = NUGA::INTERSECT::eOPER::DIFFERENCE; // for XcellN

  E_Int nb_zones = crds.size();
  E_Int nb_comps = mask_crds.size();

  std::vector<ngon_type*> ngs;
  std::vector<ngon_unit*> pgfs;
  std::vector<acrd_t*> acrds;
  std::vector<loc_t*> locs/*comp boundaries*/;
  
#ifdef FLAG_STEP
  chrono c;
  c.start();
  std::cout << "unify : init (nb zones and masks : " <<  nb_zones << "/" << nb_comps << ") ..." << std::endl;
#endif
  
  init(crds, cnts, mask_crds, mask_cnts, ngs, acrds, pgfs, mask_wall_ids, locs);

  // all start in a visible state //fixme : cuirrently unset , but we should be able to replace UNSET concept by OUT
  xcelln.resize(nb_zones);
  for (size_t i=0; i < nb_zones; ++i)xcelln[i].resize(ngs[i]->PHs.size(), UNSET);

  std::map<E_Int, std::vector<E_Int>> sorted_comps_per_comp;
  comp_priorities(priority, sorted_comps_per_comp);
  
//  for (auto it = sorted_comps_per_comp.begin(); it != sorted_comps_per_comp.end(); ++it)
//  {
//    std::cout << "priorized compos for comps : " << it->first << " => ";
//    for (size_t i=0; i < it->second.size(); ++i) std::cout << it->second[i] << "--";
//    std::cout << std::endl;
//  }
  
  
#ifdef FLAG_STEP
  std::cout << "unify : init : CPU : " << c.elapsed() << std::endl;
#endif
  
  std::vector<bool> process_overlap_pass(nb_zones, false);
  
#ifdef FLAG_STEP
  c.start();
  std::cout << "unify : first pass (blanking) ..." << std::endl;
#endif
  
  bool is_prior = true; // i.e. no mask is hiding it

  //
#pragma omp parallel for shared(process_overlap_pass) //private()
  for (size_t z=0; z < nb_zones; ++z)
  {
    E_Int compid = comp_id[z];
    auto it = sorted_comps_per_comp.find(compid);
    if (it == sorted_comps_per_comp.end()) continue;
    
    is_prior = false;

    std::vector<E_Int>& hiding_comps = it->second;
    E_Int sz(hiding_comps.size());

    STACK_ARRAY(loc_t*, sz, hiding_locs);
    STACK_ARRAY(ngon_unit*, sz, hiding_skins);
    STACK_ARRAY(K_FLD::FloatArray*, sz, hiding_crds);
    
    for (size_t i=0; i < sz; ++i)
    {
      E_Int Cp = hiding_comps[sz-1-i]; //store them in decreasing priority order
      hiding_locs[i] = locs[Cp];
      hiding_skins[i] = pgfs[Cp];
      hiding_crds[i] = mask_crds[Cp];
    }
    
    const K_FLD::FloatArray& crd = *crds[z];
    ngon_type& ng = *ngs[z];
    
#ifdef FLAG_STEP
    std::cout << "unify :  -- zone : " << z << std::endl;
#endif

    // BLANKING => IN, COLLIDE
    process_overlap_pass[z] = compute_zone_blanking(crd, ng, hiding_locs.get(), hiding_skins.get(), hiding_crds.get(), sz, xcelln[z]);
    
#ifdef OUTPUT_XCELLN
    {
      const std::vector<E_Float> & xc = xcelln[z];
      ngon_type ngin, ngx, ngunset;
      ngunset.PGs = ngin.PGs = ngx.PGs = ng.PGs;
      for (E_Int i =0; i <  ng.PHs.size(); ++i)
      {
        if (xc[i] == COLLIDE) ngx.PHs.add(ng.PHs.stride(i), ng.PHs.get_facets_ptr(i));
        else if (xc[i] ==IN) ngin.PHs.add(ng.PHs.stride(i), ng.PHs.get_facets_ptr(i));
        else /*if (xc[i] == 0.)*/ ngunset.PHs.add(ng.PHs.stride(i), ng.PHs.get_facets_ptr(i));
      }
      
      std::vector<E_Int> pgnids, phnids;
      ngx.remove_unreferenced_pgs(pgnids, phnids);
      ngin.remove_unreferenced_pgs(pgnids, phnids);
      ngunset.remove_unreferenced_pgs(pgnids, phnids);

      K_FLD::IntArray cntin, cntx, cntunset;
      ngin.export_to_array(cntin);
      ngx.export_to_array(cntx);
      ngunset.export_to_array(cntunset);
     
      std::ostringstream o;
      o << "/home/slandier/tmp/blanking_zone_" << z << "_IN.mesh";
      //MIO::write(o.str().c_str(), crd, cntin, "NGON");
      medith::write(o.str().c_str(), crd, ngin.PGs);
     
      o.str("");
      o << "/home/slandier/tmp/blanking_zone_" << z << "_X.mesh";
      //MIO::write(o.str().c_str(), crd, cntx, "NGON");
      medith::write(o.str().c_str(), crd, ngx.PGs);
      o.str("");
      o << "/home/slandier/tmp/blanking_zone_" << z << "_UNSET.mesh";
      //MIO::write(o.str().c_str(), crd, cntunset, "NGON");
      medith::write(o.str().c_str(), crd, ngunset.PGs);
    }
  
#endif
    
  } // coarse blanking loop
  
#ifdef FLAG_STEP

  if (is_prior)
    std::cout << "unify : this component is fully visible" << std::endl;
  else
  {
    std::cout << "unify : first pass blanking : CPU : " << c.elapsed() << std::endl;
    std::cout << "unify : second pass blanking (polyclip) ..." << std::endl;
    c.start();
  }

#endif
  
  if (!is_prior){
    
  // compute fine blanking
#pragma omp parallel for shared(process_overlap_pass) //private()
  for (size_t z=0; (z < nb_zones); ++z)
  {
    if (!process_overlap_pass[z]) continue; // zone is fully in or out => need a component wie coloration a posteriori
    
    E_Int compid = comp_id[z];
    auto it = sorted_comps_per_comp.find(compid);
    if (it == sorted_comps_per_comp.end()) continue;

    std::vector<E_Int>& hiding_comps = it->second;
    E_Int sz(hiding_comps.size());

    STACK_ARRAY(loc_t*, sz, hiding_locs);
    STACK_ARRAY(ngon_unit*, sz, hiding_skins);
    STACK_ARRAY(K_FLD::FloatArray*, sz, hiding_crds);
    
    for (size_t i=0; i < sz; ++i)
    {
      E_Int Cp = hiding_comps[sz-1-i]; //stor them in decreasing priority order
      hiding_locs[i] = locs[Cp];
      hiding_skins[i] = pgfs[Cp];
      hiding_crds[i] = mask_crds[Cp];
    }
    
    // optionnal output
    K_FLD::FloatArray crdo;
    ngon_type ngo;
    
    const K_FLD::FloatArray& crd = *crds[z];
    ngon_type& ng = *ngs[z];
    
    ng.flag_externals(WALL1);

#ifdef FLAG_STEP
    std::cout << " compute_zone_overlaps zone : " << z << std::endl;
    //std::cout << "check non-OUT before " << std::endl;
    //for (E_Int i=0; i < ng.PHs.size(); ++i) assert (xcelln[z][i] != OUT);
#endif

    // MULTI POLYCLIP": IN, IN_BODY, true partial coeff 
    std::vector<bool> in_points;
    bool make_inward1 = (oper == NUGA::INTERSECT::eOPER::INTERSECTION || oper == NUGA::INTERSECT::eOPER::DIFFERENCE);
    bool make_inward2 = (oper == NUGA::INTERSECT::eOPER::INTERSECTION);//outward at init (reorient skins) => nothing to do if DIFF or UNION (duplication)
    compute_zone_overlaps(crd, ng, hiding_locs.get(), hiding_skins.get(), hiding_crds.get(), sz, xcelln[z], crdo, ngo, in_points, make_inward1, make_inward2);

#ifdef OUTPUT_XCELLN
    std::ostringstream o;
    o << "/home/slandier/tmp/overlap_z_" << z << "_polyclip.plt";
//    K_FLD::IntArray cnto;
//    ngo.export_to_array(cnto);
    medith::write(o.str().c_str(), crdo, ngo);
    
    //std::cout << crdo << std::endl;
    //std::cout << cnto << std::endl;
    
    //MIO::write(o.str().c_str(), crdo, cnto, "NGON");
    
    // IN_BODY, IN and TRUE PARTIAL
    {
      std::vector<E_Int> in_body, in, true_x, out, unset, border_x_0, border_x_1;
      for (E_Int i=0; i < ng.PHs.size(); ++i)
      {
        if (xcelln[z][i] == IN_BODY) in_body.push_back(i);
        else if (xcelln[z][i] == IN) in.push_back(i);
        else if (xcelln[z][i] > 0. && xcelln[z][i] < 1.e-9) border_x_0.push_back(i);
        else if (xcelln[z][i] < 1. && xcelln[z][i] > 1. - 1.e-9) border_x_1.push_back(i);
        else if (xcelln[z][i] > 0. && xcelln[z][i] < 1.) true_x.push_back(i);
        else if (xcelln[z][i] >= OUT) out.push_back(i);
        else if (xcelln[z][i] == UNSET) unset.push_back(i);
      }
      std::ostringstream o;
      o << "/home/slandier/tmp/fine_blanking_for_zone_" << z << "_IN_BODY.mesh";
      medith::write(o.str().c_str(), crd, ng, in_body);

      o.str("");
      o << "/home/slandier/tmp/fine_blanking_for_zone_" << z << "_IN.mesh";
      medith::write(o.str().c_str(), crd, ng, in);

      o.str("");
      o << "/home/slandier/tmp/fine_blanking_for_zone_" << z << "_X.mesh";
      medith::write(o.str().c_str(), crd, ng, true_x);
      
      o.str("");
      o << "/home/slandier/tmp/fine_blanking_for_zone_" << z << "_BORDER_X_0.mesh";
      medith::write(o.str().c_str(), crd, ng, border_x_0);
      
      o.str("");
      o << "/home/slandier/tmp/fine_blanking_for_zone_" << z << "_BORDER_X_1.mesh";
      medith::write(o.str().c_str(), crd, ng, border_x_1);
      
      o.str("");
      o << "/home/slandier/tmp/fine_blanking_for_zone_" << z << "_UNSET.mesh";
      medith::write(o.str().c_str(), crd, ng, unset);
      
      o.str("");
      o << "/home/slandier/tmp/fine_blanking_for_zone_" << z << "_OUT.mesh";
      medith::write(o.str().c_str(), crd, ng, out);
      
      o.str("");
      o << "/home/slandier/tmp/fine_blanking_for_zone_" << z << "_UNSET.mesh";
      medith::write(o.str().c_str(), crd, ng, unset);
    }
#endif
    
    // COLOR any UNSET zone attached to IN_BODY into IN
    {
//      //coloring by neighborhood
      ngon_unit neighbors;
      ng.build_ph_neighborhood(neighbors);
//      
      std::vector<E_Int> toprocess; // elements attached to IN_BODY      
      for (E_Int i=0; i < ng.PHs.size(); ++i)
      {
        if (xcelln[z][i] != UNSET) continue;
        
        const E_Int* faces = ng.PHs.get_facets_ptr(i);
        E_Int nb_faces = ng.PHs.stride(i);
        
        bool found = false;
        
        for (E_Int j=0; (j < nb_faces) && !found; ++j)
        {
          const E_Int* nodes = ng.PGs.get_facets_ptr(faces[j]-1);
          E_Int nb_nodes = ng.PGs.stride(faces[j]-1);
          
          for (E_Int n=0; (n < nb_nodes) && !found; ++n)
          {
            if (in_points[nodes[n] - 1])
            {
              toprocess.push_back(i);
              found = true;
              break;
            }
          }
        }
      }
  
#ifdef OUTPUT_XCELLN
      {
        K_FLD::FloatArray toto;
        for (size_t u=0; u< in_points.size(); ++u)
        {
          if (in_points[u])toto.pushBack(crd.col(u), crd.col(u)+3);
        }

        std::cout << "nb in points for zone " << z << " : " << toto.cols() << std::endl;
        if (toto.cols())
        {
          K_FLD::IntArray cnt(2,1,0);cnt(1,0)=1;
          std::ostringstream o;
          o << "/home/slandier/tmp/in_point_for_" << z << ".mesh";

          medith::write(o.str().c_str(), toto, cnt, "BAR");
        }
      }
    // attached to IN_BODY
      {
        ngon_type titi;
        for (size_t k=0; k < toprocess.size(); ++k)
        {
          titi.addPH(ng, toprocess[k], true);
        }
        K_FLD::IntArray cT3;
        DELAUNAY::Triangulator dt;
        for (size_t k=0; k < titi.PGs.size();++k)
          K_MESH::Polygon::triangulate<DELAUNAY::Triangulator>(dt, crd, titi.PGs.get_facets_ptr(k), titi.PGs.stride(k), 1, cT3);
        std::ostringstream o;
        o << "/home/slandier/tmp/attached_to_IN_BODY_for_zone_" << z << ".mesh";
        medith::write(o.str().c_str(), crd, cT3, "TRI");
      }
      
      
#endif
      
      // now do the coloring
      for (size_t i=0; i < toprocess.size(); ++i)
      {
        E_Int Kseed = toprocess[i];
        E_Float bad_col;
        bool good_dom = K_CONNECT::EltAlgo<K_MESH::Triangle/*dummy*/>::coloring_one_connex_homogeneous (neighbors, xcelln[z], Kseed, (E_Float)UNSET, (E_Float)IN, (E_Float)IN_BODY, bad_col);
        
      }
    }
    
    // CONVERT any IN_BODY into IN and any remaining UNSET to OUT
    for (E_Int i=0; i < ng.PHs.size(); ++i)
    {
      if (xcelln[z][i] == UNSET) xcelln[z][i] = OUT;
      else if (xcelln[z][i] == IN_BODY) xcelln[z][i] = IN;
    }
    
    
#ifdef FLAG_STEP
    E_Int nb_in(0), nb_out(0), nb_x(0);
    for (E_Int i=0; i < ng.PHs.size(); ++i)
    {
      if (xcelln[z][i] == OUT) ++nb_out;
      else if (xcelln[z][i] == IN) ++nb_in;
      else ++nb_x;
    }
    std::cout << "unify : stats for zone : " << z << " : nb_in/nb_x/nb_out : " << nb_in << "/" << nb_x << "/" << nb_out << " over " << ng.PHs.size() << " cells." << std::endl;
#endif
    

    
#ifdef OUTPUT_XCELLN
    //visibles
    
    {
      ngon_type ngout, ngin, ngx;
      ngx.PGs = ngout.PGs = ngin.PGs = ng.PGs;
      
      for (E_Int i=0; i < ng.PHs.size(); ++i)
      {
        if (xcelln[z][i] == OUT)
          ngout.PHs.add(ng.PHs.stride(i), ng.PHs.get_facets_ptr(i));
        else if (xcelln[z][i] == IN)
          ngin.PHs.add(ng.PHs.stride(i), ng.PHs.get_facets_ptr(i));
        else
          ngx.PHs.add(ng.PHs.stride(i), ng.PHs.get_facets_ptr(i));
      }

      std::vector<E_Int> pgnids, phnids;
      ngout.remove_unreferenced_pgs(pgnids, phnids);
      K_FLD::IntArray cnto;
      ngout.export_to_array(cnto);
      std::ostringstream o;
      o << "overlap_z_" << z << "_out.plt";
      //MIO::write(o.str().c_str(), crd, cnto, "NGON");
      medith::write(o.str().c_str(), crd, ngout);
      
      ngin.export_to_array(cnto);
      o.str("");
      o << "overlap_z_" << z << "_in.plt";
      //MIO::write(o.str().c_str(), crd, cnto, "NGON");
      medith::write(o.str().c_str(), crd, ngin);
      
      ngx.export_to_array(cnto);
      o.str("");
      o << "overlap_z_" << z << "_x.plt";
      //MIO::write(o.str().c_str(), crd, cnto, "NGON");
      medith::write(o.str().c_str(), crd, ngx);
    }
    
#endif
    
  } // zones loop
  }
  
#ifdef FLAG_STEP
  if (!is_prior) std::cout << "unify : second pass blanking (polyclip) : CPU : " << c.elapsed() << std::endl;
#endif
  
  if (is_prior)
  {
    for (size_t i=0; i < nb_zones; ++i){
      xcelln[i].clear();
      xcelln[i].resize(ngs[i]->PHs.size(), 1.);
    }
  }
  
  for (size_t i=0; i < nb_zones; ++i)
    delete ngs[i];
  
  for (size_t i=0; i < nb_comps; ++i)
  {
    delete pgfs[i];
    delete locs[i];
  } 
}

#ifndef NETBEANSZ
PyObject* K_INTERSECTOR::unify(PyObject* self, PyObject* args)
{
  PyObject* zones, *basenum, *masks, *priorities, *basewallfaces;
  if (!PyArg_ParseTuple(args, "OOOOO", &zones, &basenum, &masks, &priorities, &basewallfaces)) return NULL;

  E_Int nb_zones = PyList_Size(zones);
  E_Int nb_basenum = PyList_Size(basenum);
  E_Int nb_masks = PyList_Size(masks);
  E_Int nb_priority_pairs = PyList_Size(priorities);

  if (nb_zones != nb_basenum)
  {
    //std::cout << "nb zones vs nb_basenum : " << nb_zones << "/" << nb_basenum << std::endl;
    PyErr_SetString(PyExc_ValueError,
       "unify: must have as many base ids as zones.");
      return NULL;
  }
  
  //std::cout << "zones/masks/prior_pairs/basenum : " << nb_zones << "/" << nb_masks << "/" << nb_priority_pairs << "/" << nb_basenum << std::endl;

  std::vector<K_FLD::FloatArray*> crds(nb_zones, nullptr), mask_crds(nb_masks, nullptr);
  std::vector<K_FLD::IntArray*>   cnts(nb_zones, nullptr), mask_cnts(nb_masks, nullptr);

  std::vector<std::pair<E_Int, E_Int>> priority;
  std::vector<E_Int> comp_id;
  std::vector< std::vector<E_Float>> xcelln(nb_zones);
 
  // get the zones
  for (E_Int i=0; i < nb_zones; ++i)
  {
    //std::cout << "getting zone in list : " << i << std::endl;
    PyObject* py_zone = PyList_GetItem(zones, i);

    char* varString, *eltType;
    
    E_Int err = check_is_NGON(py_zone, crds[i], cnts[i], varString, eltType);
    if (err)
    {
      for (E_Int i=0; i < nb_zones; ++i)
      {
        delete crds[i];
        delete cnts[i];
      }
      return NULL;
    }

    //std::cout << "zone sizes : " << crds[i]->cols() << " points" << std::endl;
    //std::cout << "zone sizes : " << cnts[i]->cols() << " cells" << std::endl;
  }

  // get the masks
  char* varString, *eltType;
  for (E_Int i=0; i < nb_masks; ++i)
  {
    PyObject* py_mask = PyList_GetItem(masks, i);

    E_Int err = check_is_NGON(py_mask, mask_crds[i], mask_cnts[i], varString, eltType);
    if (err)
    {
      for (E_Int i=0; i < nb_zones; ++i)
      {
        delete crds[i];
        delete cnts[i];
      }

      return NULL;
    }

    //std::cout << "mask sizes : " << mask_crds[i]->cols() << " points" << std::endl;
    //std::cout << "mask sizes : " << mask_cnts[i]->cols() << " cells" << std::endl;
  }

  // get the priority pairs
  for (E_Int i=0; i < nb_priority_pairs; ++i)
  {
    //std::cout << "getting priority pairs in list : " << i << std::endl;
    PyObject* py_priority = PyList_GetItem(priorities, i);
    if (!PyTuple_Check(py_priority)) continue;

    PyObject* s = PyTuple_GetItem(py_priority, 0);
    E_Int p0 = PyLong_AsLong(s);
    s = PyTuple_GetItem(py_priority, 1);
    E_Int p1 = PyLong_AsLong(s);

    priority.push_back(std::make_pair(p0,p1));
  }

  // get the basenum (zone 's component id')
  for (E_Int i=0; i < nb_zones; ++i)
  {
    //std::cout << "getting priority pairs in list : " << i << std::endl;
    PyObject* py_comp_id = PyList_GetItem(basenum, i);
    E_Int c_id = PyLong_AsLong(py_comp_id);
    comp_id.push_back(c_id);
    //std::cout << "zone/comp : " << i << c_id << std::endl;
  }

  std::vector <std::vector<E_Int>>  mask_wall_ids(nb_masks);

  if (basewallfaces != Py_None)
  {
    E_Int ok(1);
    for (E_Int i=0; (i < nb_masks) /*&& ok*/; ++i)
    {
      PyObject* py_base_wall_ids = PyList_GetItem(basewallfaces, i);
      
      E_Int sz, nfld, *pw(nullptr);
      
      ok =  K_NUMPY::getFromNumpyArray(py_base_wall_ids, pw, sz, nfld, 1/*shared*/, 0 /*inverse*/);
      if (!ok || pw == nullptr) continue;

      mask_wall_ids[i].insert(mask_wall_ids[i].end(), pw, pw+sz);
      //std::cout << "some walls has been passed for comp : " << i << std::endl;
    }
  }

  MOVLP_unify (crds, cnts, comp_id, priority, mask_crds, mask_cnts, mask_wall_ids, xcelln);

  PyObject *l(PyList_New(0)), *tpl;

  for (E_Int i = 0; i < nb_zones; ++i)
  {
    E_Int sz = xcelln[i].size();
    K_FLD::FloatArray xcellno(1,sz);
    for (size_t j = 0; j < sz; ++j)xcellno(0,j)=xcelln[i][j];
    tpl = K_ARRAY::buildArray(xcellno, "xcelln", *cnts[i], -1, eltType, false);
    //tpl = K_NUMPY::buildNumpyArray(&xcelln[i][0], xcelln[i].size(), 1, 0);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  for (E_Int i=0; i < nb_zones; ++i)
  {
    delete crds[i];
    delete cnts[i];

  }

  return l;

}
#endif

