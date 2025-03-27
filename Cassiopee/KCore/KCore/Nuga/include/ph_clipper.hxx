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
//Authors : Sam Landier (sam.landier@onera.fr)

#ifndef NUGA_INTERSECT_HXX
#define NUGA_INTERSECT_HXX

#include <memory>

#include "Nuga/include/macros.h"
#include "Nuga/include/polyhedron.hxx"
#include "Nuga/include/collider.hxx"
#include "Nuga/include/TRI_Conformizer.h"

//#define DEBUG_CLIPPER

#ifdef DEBUG_CLIPPER
//#include "Nuga/include/NGON_debug.h"
//using NGDBG = NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>;
#include "medit.hxx"
#endif


enum ePGType { ANY=-8, BCWALL=1, WALL1 = 7};

#define THIRD_NODE_POS(pS, E0, E1) ( (*pS != E0 && *pS != E1) ? 0 : (*(pS+1) != E0 && *(pS+1) != E1) ? 1 : 2) 

namespace NUGA
{
  namespace INTERSECT
  {
    
    enum eOPER { INTERSECTION, DIFFERENCE, UNION };
        
    template <typename acrd_t, typename ELT1, typename ELT2>
    void discard_overlaps(const acrd_t& acrd1, ELT1& oe1, bool inward1, const acrd_t& acrd2, ELT2& oe2, E_Float RTOL, const K_FLD::FloatArray& L, E_Float ps_min, E_Int& contact, bool dbg)
    {
      using PG1_t = typename ELT1::boundary_type;
      using PG2_t = typename ELT2::boundary_type;
      
      E_Int nb_faces1 = oe1.nb_faces();
      E_Int nb_faces2 = oe2.nb_faces();
      
      contact = 0;
      
      STACK_ARRAY(bool, nb_faces2, keep);
        
      for (E_Int i=0; i < nb_faces2; ++i)keep[i]=true;
        
      for (E_Int i1=0; i1 < nb_faces1; ++i1)
      {
        const E_Int* nodes1 = oe1.nodes(i1);
        E_Int nb_nodes1 = oe1.nb_nodes(i1);

        PG1_t PG1(nodes1, nb_nodes1, -1);

        E_Float n1[3];
        PG1.template normal<acrd_t, 3>(acrd1, nodes1, nb_nodes1, 1, n1); //watchme : base

        // Tolerance hypothese : constant per PG : take the min over the polygon nodes
        E_Float Lref = NUGA::FLOAT_MAX;
        for (E_Int n=0; n < nb_nodes1;++n)
          Lref = MIN(Lref, L(0,nodes1[n]-1));

        E_Float abstol = MAX(EPSILON, RTOL*::sqrt(Lref));

        for (E_Int i2=0; i2 < nb_faces2; ++i2)
        {
          const E_Int* nodes2 = oe2.nodes(i2);
          E_Int nb_nodes2 = oe2.nb_nodes(i2);

          PG2_t PG2(nodes2, nb_nodes2, -1);

          E_Float n2[3];
          PG2.template normal<acrd_t, 3>(acrd2, nodes2, nb_nodes2, 1, n2); //watchme : base

          double ps = NUGA::dot<3>(n1,n2);
          if (::fabs(ps) < ps_min) continue;

          //std::cout << "ps : " << ps << std::endl;

#ifdef DEBUG_CLIPPER
          /*if (dbg)
          {
            {
            ngon_unit pg;
            pg.add(nb_nodes1, nodes1);
            ngon_type ng(pg);
            K_FLD::IntArray cnt;
            ng.export_to_array(cnt);
            medith::write("pg1.plt", acrd1.array(), cnt, "NGON");
            }
            {
            ngon_unit pg;
            pg.add(nb_nodes2, nodes2);
            ngon_type ng(pg);
            K_FLD::IntArray cnt;
            ng.export_to_array(cnt);
            medith::write("pg2.plt", acrd2.array(), cnt, "NGON");
            }
          }*/
#endif
         E_Int it1, it2;
         bool ov = NUGA::COLLIDE::simplicial_colliding<K_FLD::FloatArray, 3>(acrd1.array(), PG1, acrd2.array(), PG2, K_MESH::Triangle::overlap2, abstol, it1, it2);

          if (ov)
            contact = ( (ps > 0. && inward1) || (ps < 0. && !inward1)) ? 1 : -1;

          keep[i2] &= !ov;
        }
      }
      
      K_CONNECT::IdTool::compact (oe2.faces(), oe2.nb_faces()/*IO!*/, keep.get()/*, Vector2& new_Ids*/);

#ifdef DEBUG_CLIPPER
//        {
//          std::vector<E_Int> ids;
//          const E_Int* faces = oe2.faces();
//          for (size_t i=0; i <oe2.nb_faces(); ++i)ids.push_back(faces[i]-1);
//          NGDBG::draw_PGs("shape_wo_ovlp", acrd2.array(), *oe2.pgs(), ids);
//        }
#endif
    }
    
    ///
    template <typename acrd_t, typename ELT1, typename ELT2>
    E_Int isolated_clip(const acrd_t& acrd1, const ELT1& subj, bool inward1, const acrd_t& acrd2, ELT2& cutter, E_Float ps_min, E_Float RTOL, NUGA::haPolyhedron<UNKNOWN>& result, E_Int& contact, bool& true_clip, bool dbg)
    {
      E_Int err(0);

      true_clip = false;
      
      result.clear();

      // for nodal tolerance
      K_FLD::FloatArray L;
      NUGA::MeshTool::computeIncidentEdgesSqrLengths(acrd1.array(), *subj.pgs(), L);
      if (L.cols() == 0) return 1;

      //E_Int nb_faces1 = subj.nb_faces();
      E_Int &nb_faces2 = cutter.nb_faces();

      if (ps_min > 0.) //discard any "overlapping" PG in e2 (the more ps_min tends to 1, the more 'overlapping' has a meaning)
        discard_overlaps(acrd1, subj, inward1, acrd2, cutter, RTOL, L, ps_min, contact, dbg); //cutter is filtered
      
      if (contact != 0 && nb_faces2 == 0) // e2 was just in contact
      {
        return 0;
      }

      //gather points
      K_FLD::FloatArray crd(acrd1.array());
      E_Int nb_pts1 = crd.cols();
      crd.pushBack(acrd2.array());
      
      // go to a triangle view 
      DELAUNAY::Triangulator dt;
      err = subj.triangulate(dt, acrd1.array());
      if (err) 
      {
#ifdef DEBUG_CLIPPER
        std::cout << "ph_clipper::isolated_clip : error triangulating subject" << std::endl;
        return err;
#endif
      }

      err = cutter.triangulate(dt, acrd2.array());
      if (err) 
      {
#ifdef DEBUG_CLIPPER
        std::cout << "ph_clipper::isolated_clip : error triangulating cutter" << std::endl;
        /*
        E_Int faultyPG = err-1;
        std::vector<E_Int> toprocess;
        toprocess.push_back(faultyPG);
        medith::write("errorPG18", cutter.m_crd, cutter.m_pgs, &toprocess);*/
        return err;
#endif
      }

      std::vector<E_Int> ancPG1, ancPG2;
      subj.get_triangle_oids(ancPG1);
      cutter.get_triangle_oids(ancPG2);

      E_Int nb_tris1 = subj.nb_tris();
      E_Int nb_tris2 = cutter.nb_tris();

      K_FLD::IntArray cT3;
      K_FLD::FloatArray normalsT3(3, nb_tris1 + nb_tris2, NUGA::FLOAT_MAX);
      E_Int T[3];
      E_Float nT3[3];
      for (E_Int i=0; i < nb_tris1; ++i)
      {
        subj.triangle(i, T); //watchme : base ?
        cT3.pushBack(T, T+3);

        //normals
        K_MESH::Triangle::normal(acrd1.array(), &T[0], nT3);
        E_Float l2 = ::sqrt(nT3[0] * nT3[0] + nT3[1] * nT3[1] + nT3[2] * nT3[2]);
        if (::fabs(l2 - 1.) >= EPSILON) // DEGEN
          K_MESH::Polygon::normal<acrd_t, 3>(acrd1, subj.m_pgs.get_facets_ptr(ancPG1[i]), subj.m_pgs.stride(ancPG1[i]), 1, nT3);
        normalsT3(0, i) = nT3[0];
        normalsT3(1, i) = nT3[1];
        normalsT3(2, i) = nT3[2];
      }
      
      
      for (E_Int i=0; i < nb_tris2; ++i)
      {
        cutter.triangle(i, T); //watchme : base ?
 
        assert (ancPG2[i] < cutter.m_pgs.size());

        //normals
        //std::cout << "coord sz for cutter : " << acrd2.array().cols() << std::endl;
        //std::cout << "triangle " << i << " : " << T[0] << "/" << T[1] << "/" << T[2] << std::endl;
        assert (T[0] < acrd2.array().cols());
        assert (T[1] < acrd2.array().cols());
        assert (T[2] < acrd2.array().cols());
        K_MESH::Triangle::normal(acrd2.array(), T, nT3);
        E_Float l2 = ::sqrt(nT3[0] * nT3[0] + nT3[1] * nT3[1] + nT3[2] * nT3[2]);
        if (::fabs(l2 - 1.) >= EPSILON) // DEGEN
          K_MESH::Polygon::normal<acrd_t, 3>(acrd2, cutter.m_pgs.get_facets_ptr(ancPG2[i]), cutter.m_pgs.stride(ancPG2[i]), 1, nT3);

        assert ((i + nb_tris1) < normalsT3.cols());
        normalsT3(0, i + nb_tris1) = nT3[0];
        normalsT3(1, i + nb_tris1) = nT3[1];
        normalsT3(2, i + nb_tris1) = nT3[2];

        // shiftand append after normal computation
        T[0] += nb_pts1;
        T[1] += nb_pts1;
        T[2] += nb_pts1;
        cT3.pushBack(T, T + 3);
      }
      
      // type transmission
      std::vector<E_Int> type(subj.nb_tris() + cutter.nb_tris(), ANY);
      if (!subj._pgs->_type.empty())
        for (E_Int i=0; i<nb_tris1; ++i)
          type[i] = subj._pgs->_type[ancPG1[i]];
      if (!cutter._pgs->_type.empty())
        for (E_Int i=0; i<nb_tris2; ++i)
          type[nb_tris1 + i] = cutter._pgs->_type[ancPG2[i]];
      
      //medith::write("triangles.mesh", crd, cT3, "TRI");
#ifdef DEBUG_CLIPPER
      if (dbg){
        medith::write("triangles", crd, cT3, "TRI");
        //std::cout << crd << std::endl;
        //std::cout << cT3 << std::endl;
      }
#endif

      // compute an overall abstol
      E_Float min_d, max_d, abstol;
      NUGA::MeshTool::computeMinMaxEdgeSqrLength<3>(crd, cT3, min_d, max_d);
      abstol= ::sqrt(min_d) * RTOL;

      // conformize this cloud
      std::vector<E_Int> ancT3;
      TRI_Conformizer<3> conformizer(true/* keep track of nodes history*/);
      conformizer._split_swap_afterwards = false;
      
      //conformizer._brute_force=true;

#ifdef DEBUG_CLIPPER
      conformizer._silent_errors = false;
#else
      conformizer._silent_errors = true;
#endif
      
      E_Int nb_tris0 = cT3.cols();
      //medith::write("/home/slandier/projects/xcelln/CRM2/conferr.mesh", crd, cT3, "TRI");
//      std::cout << "abstol : " << abstol << std::endl;
//      std::cout << "nb_tris1 : " << nb_tris1 << std::endl;
//      std::cout << "iso : 9" << std::endl;
      err = conformizer.run(crd, cT3, ancT3, nullptr/*&priority*/, abstol, nb_tris1, 1 /*one iter only*/);
      if (err)
      {
#ifdef DEBUG_CLIPPER
        std::cout << "ph_clipper::isolated_clip : error conformizer" << std::endl;
#endif
        return err;
      }

      if (cT3.cols() == nb_tris0) return 0 ;// no intersections => fully visible or hidden

      true_clip = true;

#ifdef DEBUG_CLIPPER
      if (dbg)
        medith::write("conformized", crd, cT3, "TRI", 0, &ancT3);
#endif

      // recompute normals or take ancestor is not doable
      K_FLD::FloatArray new_normalsT3(3, cT3.cols(), NUGA::FLOAT_MAX);
      for (E_Int i = 0; i < cT3.cols(); ++i)
      {
        K_MESH::Triangle::normal(crd, cT3.col(i), nT3);
         E_Float l2 = ::sqrt(nT3[0] * nT3[0] + nT3[1] * nT3[1] + nT3[2] * nT3[2]);
        if (::fabs(l2 - 1.) >= EPSILON) // DEGEN
        {
          nT3[0] = normalsT3(0, ancT3[i]);
          nT3[1] = normalsT3(1, ancT3[i]);
          nT3[2] = normalsT3(2, ancT3[i]);
        }
        new_normalsT3(0, i) = nT3[0];
        new_normalsT3(1, i) = nT3[1];
        new_normalsT3(2, i) = nT3[2];
      }
      normalsT3 = std::move(new_normalsT3);

      // duplicates
      {
        E_Int nbt3o = cT3.cols();
        std::vector<E_Int> dupids;
        NUGA::MeshTool::removeDuplicated(cT3, dupids, false/*strict orient*/);

        if (cT3.cols() < nbt3o) // update ancT3
        {
          K_CONNECT::unchanged pred(dupids);
          K_CONNECT::IdTool::compress(ancT3, pred);
          K_CONNECT::IdTool::compress(normalsT3, pred);
        }
      }

      // update type
      {
        E_Int nbt3 = cT3.cols();
        std::vector<E_Int> new_type(nbt3, ANY);
        for (E_Int i=0; i <nbt3; ++i)
          new_type[i] = type[ancT3[i]];
        type = new_type;
      }

      // keep relevant pieces
      K_FLD::IntArray neighbors, neighbors_cpy;
      err = NUGA::EltAlgo<K_MESH::Triangle>::getNeighbours (cT3, neighbors, false/*means put somthing on manifolds to distinguish them from free edges*/);
      neighbors_cpy = neighbors; //save it as we are going to erase non manifold info to split but we need it afterwards
      // do the coloring to see non manifoldnesses and non-connexities
      //so add the cuts in the graph
      for (E_Int i=0; i < neighbors.cols(); ++i)
        for (E_Int j=0; j < 3; ++j)
          if (neighbors(j,i) == NON_MANIFOLD_COL) neighbors(j,i) = IDX_NONE;
      
      std::vector<E_Int> nonmfld_bits;
      NUGA::EltAlgo<K_MESH::Triangle>::coloring_pure (neighbors, nonmfld_bits);
      
      //
      E_Int nb_nonmfld_bits = *std::max_element(ALL(nonmfld_bits))+1;
      //
      STACK_ARRAY(bool, nb_nonmfld_bits, good_col);
      for (E_Int i=0; i < nb_nonmfld_bits; ++i)good_col[i]=true;
      // now use neighbors_cpy to check those with free edges to discard them (burning)
      for (E_Int i=0; i < neighbors_cpy.cols(); ++i)
        for (E_Int j=0; j < 3; ++j)
          if (neighbors_cpy(j,i) == IDX_NONE) good_col[nonmfld_bits[i]]=false;
        
      
      E_Int nb_t3 = cT3.cols();
      STACK_ARRAY(bool, nb_t3, keepT3);
      for (E_Int i=0; i < nb_t3; ++i)
        keepT3[i]= good_col[nonmfld_bits[i]] ? true : false;

      // for (size_t i=0; i < 111; ++i)std::cout << keepT3[i] << std::endl;
      
      K_CONNECT::keep2<bool> pred(keepT3.get(), nb_t3);
      std::vector<E_Int> nids;
      K_CONNECT::IdTool::compress(cT3, pred, nids);

      if (cT3.cols() == 0)
      {
        result.clear();
        return 0;
      }

      K_CONNECT::IdTool::compress(type, pred);
      K_CONNECT::IdTool::compress(normalsT3, pred);

      // update also neighbors_cpy
      K_CONNECT::IdTool::compress(neighbors, pred); //contains IDX_NONE at non-manfold edges
      K_FLD::IntArray::changeIndices(neighbors, nids);

      //std::cout << neighbors << std::endl;
#ifdef DEBUG_CLIPPER
      if (dbg){
        K_CONNECT::IdTool::compress(ancT3, pred);
        medith::write("reduced", crd, cT3, "TRI", 0, &ancT3);
      }
#endif
    
      // update neighbors at non-manifold edges
      
      using acnt_t = K_FLD::ArrayAccessor<K_FLD::IntArray>;
      using algoT3 = NUGA::EltAlgo<K_MESH::Triangle>;
      algoT3::BoundToEltType noE_to_oTs; // non oriented edge to oriented triangles (1 to n).
      //E_Int maxN = 
      acnt_t acT3(cT3);
      algoT3::getBoundToElements(acT3, noE_to_oTs);
      nb_t3 = cT3.cols();
      std::vector<bool> reversed;
      Vector_t<std::pair<E_Float, E_Int> > palmares;
      E_Float ERRORVAL(2.*NUGA::PI);
      for (auto it = noE_to_oTs.begin(); it != noE_to_oTs.end(); ++it)
      {
        const K_MESH::Triangle::boundary_type& E = it->first;
        Vector_t<E_Int>& T3s = it->second;
        
        assert (T3s.size() > 1); //should have been burned. hypothesis : do not consider open surface bits for output
        //if (T3s.size() == 2) continue;//fixme : neighbors must be updated at PREVIOUS and CURRENT non-manifold edges only
        
        reversed.clear();
        reversed.resize(T3s.size(), false);
        palmares.clear();
        
        E_Int E0 = E.node(0);
        E_Int E1 = E.node(1);
        
        E_Int K0 = T3s[0]; //K0 : the ref
        //std::cout << "K0 : " << cT3(0,K0) << "/" << cT3(1,K0) << "/" << cT3(2,K0) << std::endl;
        palmares.push_back(std::make_pair(-1., 0));//fixme : 0. is not enought o guarantee firstness ?
        reversed[0] = true; // for consistent logic on oi and oip1 below
        
        E_Int* pS = cT3.col(K0); //T0 : the ref
        E_Int i0 = K_MESH::Triangle::getLocalNodeId(pS, E0);
        if (*(pS+(i0+1)%3) != E1)
          std::swap(E0,E1); // Now E0E1 is oriented as in T0
        
        //const E_Float* norm0 = normals.col(K0);
        
        for (size_t j=1; j < T3s.size(); ++j)
        {
          E_Int& Kj = T3s[j];
          //std::cout << "Kj : " << cT3(0, Kj) << "/" << cT3(1, Kj) << "/" << cT3(2, Kj) << std::endl;
          //finding out Ap and use the right orientation for Ki.
          E_Int* pS = cT3.col(Kj);
          E_Int i0 = K_MESH::Triangle::getLocalNodeId(pS, E0);
          if (*(pS+(i0+1)%3) == E1) // the side to consider is the opposite
            reversed[j]=true;
          
          assert(Kj > -1 && Kj < normalsT3.cols());
          E_Float normj[] = {normalsT3(0, Kj), normalsT3(1, Kj), normalsT3(2,Kj)};
          
          if (reversed[j])
          {
            normj[0] = -normj[0];
            normj[1] = -normj[1];
            normj[2] = -normj[2]; 
          }
          
          assert(K0 > -1 && K0 < normalsT3.cols());
          E_Float q = NUGA::angle_measure(normalsT3.col(K0), normj, crd.col(E0), crd.col(E1));
          if (q == ERRORVAL)
          {
#if defined (DEBUG_CLIPPER)
            std::cout << "ph_clipper::isolated_clip : ERROR at edge E0E1 : " << E0 << "/" << E1 << std::endl;
            std::cout << "ph_clipper::isolated_clip : The conformizer missed some intersections there." << std::endl;
#endif
            err = 1;
            break;
          }
    
          palmares.push_back(std::make_pair(q, j));          
        }
        
        std::sort(ALL(palmares));
        E_Int sz = (E_Int)palmares.size();
        
#if defined (DEBUG_CLIPPER)
        if (err)  
        {  
          K_FLD::IntArray sorted_cnt;
          Vector_t<E_Int> colors;
          Vector_t<bool> keep(cT3.cols(), false);
          sorted_cnt.reserve(3, T3s.size());
          colors.resize(T3s.size(), 1);

          std::vector<E_Int> PGs;
          for (size_t i = 0; i < T3s.size(); ++i)
          {
            sorted_cnt.pushBack(cT3.col(T3s[i]), cT3.col(T3s[i])+3);
            keep[T3s[i]]=true;
            std::cout << "Triangle : " << T3s[i] << std::endl;
            if (i > 0)colors[i] = 0;
          }

          {
            std::ostringstream o;
            o << "sorted_on_edge_" << E0 << "_" << E1 << ".mesh";
            medith::write(o.str().c_str(), crd, sorted_cnt, "TRI", 0, &colors);
          }

          {
            std::ostringstream o;
            o << "Wsorted_on_edge_" << E0 << "_" << E1 << ".mesh";
            //TRI_debug::write_wired(o.str().c_str(), crd, cT3, normals, 0, &keep,true);
          }
        }
#endif
        if (err) return err;
        
        for (E_Int i=0; i< sz; ++i)
        {
          E_Int j = palmares[i].second;
          E_Int jp1 = palmares[(i+1)%sz].second;
          E_Int Ki = T3s[j];
          E_Int Kip1 = T3s[jp1];
          
          bool oi = reversed[j];
          bool oip1 = reversed[jp1];
          
          if (oi != true || oip1 != false) continue;
          
          E_Int Li = THIRD_NODE_POS(cT3.col(Ki), E0, E1);
          E_Int Lip1 = THIRD_NODE_POS(cT3.col(Kip1), E0, E1);
          
          neighbors(Li, Ki) = Kip1;
          neighbors(Lip1, Kip1) = Ki;     
        }
      }
      
      //std::cout << neighbors << std::endl;
      
      // cells detection
      std::vector<E_Int> T3_to_PHT3;
      NUGA::EltAlgo<K_MESH::Triangle>::coloring_pure(neighbors, T3_to_PHT3);
      E_Int nb_bits = *std::max_element(ALL(T3_to_PHT3))+1;
      
//#ifdef DEBUG_CLIPPER
//      std::vector<bool> kp(cT3.cols(), false);
//      kp[26]=true;
//      medith::write("k26.mesh", crd, cT3, "TRI", &kp/*, &colors*/);
//      TRI_debug::coloring_frames(crd, cT3, neighbors, 0);
//#endif
 
#ifdef DEBUG_CLIPPER
//      std::map<E_Int, Vector_t<E_Int> > PHT3s;
//      E_Int sz1(cT3.getSize());
//      for (E_Int i = 0; i < sz1; ++i)
//        PHT3s[T3_to_PHT3[i]].push_back(i);
//      
//      for (auto i = PHT3s.begin(); i != PHT3s.end(); ++i)
//        NGDBG::draw_PHT3(crd, cT3, PHT3s, i->first);
#endif
      
      // remove open bits
      //reuse good_col
      for (E_Int i=0; i < nb_bits; ++i)good_col[i]=true;
      // now use neighbors to check those with free edges to discard them 
      for (E_Int i=0; i < neighbors.cols(); ++i)
        for (E_Int j=0; j < 3; ++j)
          if (neighbors(j,i) == IDX_NONE) good_col[T3_to_PHT3[i]]=false;

      nb_t3 = cT3.cols();
      //reuse keepT3
      for (E_Int i=0; i < nb_t3; ++i)
        keepT3[i]= good_col[T3_to_PHT3[i]] ? true : false;

      {
        K_CONNECT::keep2<bool> pred(keepT3.get(), nb_t3);
        K_CONNECT::IdTool::compress(cT3, pred, nids);
        K_CONNECT::IdTool::compress(type, pred);
      }

      nb_t3 = cT3.cols();
      if (nb_t3 == 0)
      {
        result.clear();
        return 0;
      }

      //std::set<E_Int> all_types(ALL(type));
      //if (all_types.size() == 1) // error : only colors of one operand at this stage means wrong conformizing because other colors have been burnt
      {
        //medith::write("cutT3", crd, cT3, "TRI");
        //std::cout << "error clipping" << std::endl;
        //std::cout << "min_d / max_d : " << ::sqrt(min_d) << "/" << ::sqrt(max_d) << std::endl;
        //return 1;
        //
      }
      
      // build the output ph

      result.m_crd = crd;
      //if (dbg) std::cout << crd << std::endl;
      //std::cout << cT3 << std::endl;

      const std::vector<E_Int>& xpoids = conformizer.get_node_history();
      result.poids = xpoids;
      result.poids.resize(crd.cols(), IDX_NONE);//fixme : ?????

#ifdef DEBUG_CLIPPER
      medith::write("cutT3", crd, cT3, "TRI");
#endif

      ngon_unit::convert_fixed_stride_to_ngon_unit(cT3, 1, result.m_pgs);

      assert(type.size() == (size_t)result.m_pgs.size());
      result.m_pgs._type = type;

      E_Int nb_pgs = result.m_pgs.size();
      K_CONNECT::IdTool::init_inc(result.m_faces, nb_pgs, 1);
      result._nb_faces = nb_pgs;
      result.plug();

      return err;
    }
    
  } //INTERSECT
}   // NUGA

#endif /* NUGA_INTERSECT_HXX */

