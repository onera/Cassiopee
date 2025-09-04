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

#ifndef __NGON_T_HXX__
#define	__NGON_T_HXX__

#include "ngon_unit.h"
#include "openMP.h"
#include <iostream>
#include <map>
#include "Nuga/include/Edge.h"
#include "Nuga/include/Polygon.h"
#include "Nuga/include/Triangle.h"
#include "Nuga/include/Quadrangle.h"
#include "Nuga/include/Polyhedron.h"
#include "Nuga/include/polyhedron.hxx"
#include "Nuga/include/IdTool.h"
#include "Nuga/include/EltAlgo.h"
#include "Nuga/include/BARSplitter.h"
#include "Nuga/include/MeshTool.h"
#include "Nuga/include/FittingBox.h"
#include "Nuga/include/macros.h"
#include "Nuga/include/BbTree.h"
#include <limits.h>
#include <unordered_map>

#ifdef DEBUG_NGON_T
#include "IO/io.h"
#include "Nuga/include/NGON_debug.h"
#include "Nuga/include/TRI_debug.h"
#define NGDBG NGON_debug<K_FLD::FloatArray, Connectivity_t>
#include "Nuga/include/medit.hxx"
#endif
#ifdef FLAG_STEP
#include "Nuga/include/chrono.h"
#endif

#define STARTFACE 0
#define STARTELTS(cn) (2+cn[1])
#define CONVEXITY_TOL 1.e-8 //to be consistent with what is implemented

///
template <typename Connectivity_t>
struct ngon_t
{
  enum ePathoPG { PATHO_PG_NONE=0, SPIKE, HAT, DEGEN_PG, PG_DELAUNAY_FAILURE = dDELAUNAY_FAILURE};
  enum ePathoPH { PATHO_PH_NONE=dPATHO_PH_NONE, CENTROID_NOT_STAR =dCENTROID_NOT_STAR, 
                  ISO_BARY_NOT_STAR =dISO_BARY_NOT_STAR, OPEN_PHS = dOPEN_PHS, CONCAVITY_TO_SPLIT = dCONCAVITY_TO_SPLIT, PH_DELAUNAY_FAILURE = dDELAUNAY_FAILURE, PH_DEGEN = 999};
  enum eExtrudeStrategy { CST_ABS=0, CST_REL_MEAN, CST_REL_MIN, VAR_REL_MEAN, VAR_REL_MIN};

  enum eGEODIM { ERROR=-3, UNSET = -2, MIXED = -1, LINEIC=0, SURFACIC_CASSIOPEE = 1, SURFACIC = 2, VOLUMIC = 3};
  
  using unit_type = ngon_unit;
 
  ///
  ngon_t(const Connectivity_t& cNGON) :PGs(cNGON.begin() + STARTFACE), PHs(cNGON.begin() + STARTELTS(cNGON)){ updateFacets();}
  ///
  ngon_t(const ngon_t& ngt) :PGs(ngt.PGs), PHs(ngt.PHs){ updateFacets(); }

  /// WARNING : need a call to clean_connectivity
  template <typename ELT>
  static void convert(const K_FLD::IntArray& cnt,ngon_t& ng)
  {
    E_Int nbe = cnt.cols();
    
    using btype = typename ELT::boundary_type;
    
    E_Int molec[ELT::NB_BOUNDS];
    
    E_Int pg_count(1);
    
    for (E_Int i=0; i < nbe; ++i)
    {
      ELT e(cnt.col(i), 1); //convert to 1-based
      
      for (E_Int j=0; j< ELT::NB_BOUNDS; ++j)
      {
        btype b;
        e.getBoundary(j, b);
        
        ng.PGs.add(btype::NB_NODES, b.nodes());
        molec[j] = pg_count++;
      }
      
      ng.PHs.add(ELT::NB_BOUNDS, molec);
    }
    ng.PGs.updateFacets();
    ng.PHs.updateFacets();
  }
  /// 
  ngon_t(const E_Int* cNGon, E_Int sz1, E_Int nb_pgs, const E_Int* cNFace, E_Int sz2, E_Int nb_phs) :PGs(cNGon, sz1, nb_pgs), PHs(cNFace, sz2, nb_phs){ PGs.updateFacets(); PHs.updateFacets(); }
  ///
  ngon_t& operator=(const Connectivity_t& cNGON)
  {
    PGs = ngon_unit(cNGON.begin() + STARTFACE); PGs.updateFacets();
    PHs = ngon_unit(cNGON.begin() + STARTELTS(cNGON)); PHs.updateFacets();
    return *this;
  }
  ///
  ngon_t(const ngon_unit& pg, const ngon_unit& ph) :PGs(pg), PHs(ph){ PGs.updateFacets(); PHs.updateFacets(); }
  /// move version
  ///
  ngon_t(ngon_unit&& pg, ngon_unit&& ph) :PGs(std::move(pg)), PHs(std::move(ph)) { PGs.updateFacets(); PHs.updateFacets(); }
  ///
  ngon_t(const ngon_unit& pgs, bool one_ph_for_all=false):PGs(pgs)//false means one ph per pg.
  {
    PGs.updateFacets();
    E_Int sz= PGs.size();
        
    if (one_ph_for_all)
    {
      PHs._NGON.reserve(sz+3);
      PHs._NGON.resize(2,0);
      PHs._NGON.push_back(sz);
      for (E_Int i = 0; i < sz ; ++i) PHs._NGON.push_back(i+1);
      PHs._NGON[0] = 1;
      PHs._type.resize(1, INITIAL_SKIN);//only one PH => skin type
      PGs._type.resize(sz, INITIAL_SKIN);//only one PH => skin type
    }
    else
    {
      PHs._NGON.resize(2*sz+2, 1);
      for (E_Int i=0; i <sz ; ++i) PHs._NGON[2+2*i+1]=i+1;
      PHs._NGON[0] = sz;
    }
         
    PHs._NGON[1] = PHs._NGON.size()-2;
    PHs.updateFacets();
  }
  
  ///
  ngon_t(){}

  void updateFacets() const { PGs.updateFacets(); PHs.updateFacets(); }

  // for HEXA only
  static void create_ngon_from_connectivity(ngon_t& ngon, const K_FLD::IntArray& cnt)
  {
    auto& PHs = ngon.PHs._NGON;
    auto& PGs = ngon.PGs._NGON;

    // stride
    E_Int nfpc = 6;
    E_Int nvpc = 8;
    E_Int nvpf = 4;

    PGs.push_back(0);
    PGs.push_back(0);

    // number of cells
    E_Int ncells = cnt.cols();

    PHs.resize(2,0);
    PHs[0] = ncells;

    E_Int size = ncells * (nfpc + 1);
    PHs[1] = size;

    PHs.resize(size+2);

    // maps faces (polygons) to a unique key
    std::unordered_map<K_MESH::Quadrangle,uint32_t,Quadrangle_Hash> pg_map;

    // loop over the elements and create bottom, top, left, right, front, back
    E_Int pos(2);
    E_Int face_id(0);

    E_Int nodes[8];
      
    K_MESH::Quadrangle face;

    for (int i = 0; i < ncells; i++)
    {
      E_Int celli_pos = pos + i*(nfpc+1);

      PHs[celli_pos] = 6;
      E_Int local_face = 1;

      for (E_Int j = 0; j < nvpc; j++) nodes[j] = cnt(j,i)+1;

      // bottom
      face.setNodes(nodes[0], nodes[1], nodes[2], nodes[3]);
      ngon.insert_face(face, pg_map, ngon, celli_pos, local_face, face_id);

      // top
      face.setNodes(nodes[4], nodes[5], nodes[6], nodes[7]);
      ngon.insert_face(face, pg_map, ngon, celli_pos, local_face, face_id);

      // front
      face.setNodes(nodes[0], nodes[1], nodes[5], nodes[4]);
      ngon.insert_face(face, pg_map, ngon, celli_pos, local_face, face_id);

      // back
      face.setNodes(nodes[3], nodes[2], nodes[6], nodes[7]);
      ngon.insert_face(face, pg_map, ngon, celli_pos, local_face, face_id);

      // left
      face.setNodes(nodes[0], nodes[3], nodes[7], nodes[4]);
      ngon.insert_face(face, pg_map, ngon, celli_pos, local_face, face_id);

      // right
      face.setNodes(nodes[1], nodes[2], nodes[6], nodes[5]);
      ngon.insert_face(face, pg_map, ngon, celli_pos, local_face, face_id);
    }

    PGs[0] = face_id;
    PGs[1] = face_id * (1 + nvpf);

    ngon.PGs.updateFacets();
    ngon.PHs.updateFacets();
    
  }

  bool insert_face(K_MESH::Quadrangle& face, std::unordered_map<K_MESH::Quadrangle,uint32_t,Quadrangle_Hash>& pg_map, ngon_t& ngon, const E_Int celli, E_Int& local_face, E_Int& face_id)
  {
      bool inserted = false;

      auto& PHs = ngon.PHs._NGON;
      auto& PGs = ngon.PGs._NGON;

      // check if face is in map
      auto search = pg_map.find(face);
      if (search != pg_map.end()) // found
      {
          // retrieve its id
          E_Int id = search->second;
          // insert it
          PHs[celli+local_face] = id;
          inserted = false;
      }
      else
      {
          // first time making this face
          face_id++;
          PHs[celli+local_face] = face_id;
          // insert it in PGs with its nodes
          PGs.push_back(4);
          PGs.push_back(face.nodes()[0]);
          PGs.push_back(face.nodes()[1]);
          PGs.push_back(face.nodes()[2]);
          PGs.push_back(face.nodes()[3]);
          pg_map.emplace(face, face_id);
          inserted = true;
      }
      local_face++;
      return inserted;
  }

  static eGEODIM get_ngon_geodim(const K_FLD::IntArray& cnt)
  {
    ngon_t ng(cnt);

    E_Int nb_pgs(ng.PGs.size());
    E_Int nb_phs(ng.PHs.size());

    E_Int min_s(INT_MAX), max_s(0);

    for (E_Int i = 0; i < nb_pgs; ++i)
    {
      E_Int s = ng.PGs.stride(i);
      min_s = std::min(min_s, s);
      max_s = std::max(max_s, s);
    }

    if (min_s <= 0) // stride 0 => error
      return eGEODIM::ERROR;

    // => min_s & max_s >= 1

    if (max_s == 1)
      return eGEODIM::LINEIC;

    // => min_s >= 1 & max_s >= 2
    
    if (min_s == 1)
      return eGEODIM::MIXED;
    
    // => min_s >=2 & max_s >= 2
    
    if (max_s == 2)
      return eGEODIM::SURFACIC_CASSIOPEE;

    // => min_s >= 2 & max_s >= 3

    if (min_s == 2)
      return eGEODIM::MIXED;

    // => min_s >= 3 & max_s >= 3 == > SURFACIC or VOLUMIC ?

    if (ng.PHs.size() == 1 && ng.PHs.stride(0) == ng.PGs.size()) // SURFACIC with one PH for all
      return eGEODIM::SURFACIC;

    // compute max stride for PHs
    max_s = 0;
    for (E_Int i = 0; i < nb_phs; ++i)
    {
      E_Int s = ng.PHs.stride(i);
      max_s = std::max(max_s, s);
    }

    if (ng.PHs.size() == ng.PGs.size())
    {
      if (max_s == 1) // SURFACIC with one ph per face 
        return eGEODIM::SURFACIC;
      //else                      fixme : not sure how to handle this
      //  return eGEODIM::ERROR;
    }

    return eGEODIM::VOLUMIC;
  }

  bool is_consistent(E_Int max_node_id) const { 
    
    if (! attributes_are_consistent()) return false;
       
    E_Int maxPGid = PHs.get_facets_max_id();
    E_Int nb_pgs =  PGs.size();
    if (maxPGid > nb_pgs) return false;
    
    E_Int maxNodid = PGs.get_facets_max_id();
    if (maxNodid > max_node_id) return false;
    
    return true;
  }
  
  bool attributes_are_consistent() const
  {
    if (! PGs.attributes_are_consistent()) return false;
    if (! PHs.attributes_are_consistent()) return false;
    return true;
  }
  
  ///
  void addPH(const ngon_t& ngi, E_Int PHi, bool extractFaces = false)
  {
    if (PGs.size()) PGs.updateFacets();
    if (PHs.size()) PHs.updateFacets();
    
    if (extractFaces)
    {
      E_Int npgid = PGs.size()+1;
      E_Int nb_pgs = ngi.PHs.stride(PHi);
    
      Vector_t<E_Int> molecPH;
        
      const E_Int* pPGi = ngi.PHs.get_facets_ptr(PHi);
      for (E_Int i = 0; i <nb_pgs; ++i, ++npgid)
      {
        E_Int PGi = *(pPGi+i)-1;
        PGs.add(ngi.PGs.stride(PGi), ngi.PGs.get_facets_ptr(PGi));
        if (!ngi.PGs._type.empty())
          PGs._type.push_back(ngi.PGs._type[PGi]);
        if (ngi.PGs._ancEs.cols() != 0)
          PGs._ancEs.pushBack(ngi.PGs._ancEs.col(PGi), ngi.PGs._ancEs.col(PGi)+2);

        molecPH.push_back(npgid);
      }  
      PHs.add(molecPH.size(), &molecPH[0]);
    }
    else //warning : the caller must ensure to add these PGs outside  
      PHs.add(ngi.PHs.stride(PHi), ngi.PHs.get_facets_ptr(PHi));
    
    if (!ngi.PHs._type.empty())
      PHs._type.push_back(ngi.PHs._type[PHi]);
    if (ngi.PHs._ancEs.cols() != 0)
      PHs._ancEs.pushBack(ngi.PHs._ancEs.col(PHi), ngi.PHs._ancEs.col(PHi) + 2);
  }
  
  ///
  void export_to_array (Connectivity_t& c) const { 
    c.clear();
    E_Int ng_pgs = PGs._NGON.size();
    E_Int ng_phs = PHs._NGON.size();
    E_Int sum = ng_pgs + ng_phs;

    if (sum < 0)
    {
      std::cout << "export_to_array error : ngon too big to be exported (int32 issue)" << std::endl;
      return;
    }

    c.reserve(1, sum);

    if (PGs.size())
      c.pushBack(PGs._NGON);
    else
      c.resize((E_Int)1, (E_Int)2, (E_Int)0);

    if (PHs.size())
      c.pushBack(PHs._NGON);
    else
    {
      const E_Int& sizeFN = c[1];
      c.resize((E_Int)1, E_Int(sizeFN+4), (E_Int)0);
    }
  }
  
  // Converts a SURFACE formatted as Edges/PGs into a ngon_unit (PHs is cleared) representing the PGs
  // ORIENTATION IS PRESERVED IF POSSIBLE : consitent with the first PG
  E_Int export_surfacic_view(Connectivity_t& c)
  {
    c.clear();
    
    PGs.updateFacets();
    PHs.updateFacets();

    // Surfacic ngon ?
    for (E_Int i=0; (i<PGs.size()); ++i){
      if (PGs.stride(i) != 2)
        return 1;
    }
    
    c.resize((E_Int)1, (E_Int)2, (E_Int)0);
     
    E_Int nb_polys=0;
    std::vector<E_Int> molec,snodes;
    K_FLD::IntArray cB;
        
    ngon_unit& polys = PHs;
    ngon_unit& edges = PGs;
    
    for (E_Int i=0; (i<polys.size()); ++i){
            
      E_Int nb_edges = polys.stride(i);

      if (nb_edges < 3) continue; //degen

      const E_Int* pEs = polys.get_facets_ptr(i);

      bool edges_are_sorted = true;
      snodes.clear();

      for (E_Int j=0; (j<nb_edges-1); ++j) //assume first edge are sorted
      {
        E_Int Ej = pEs[j]-1;
        E_Int Ejp1 = pEs[j+1]-1;
        const E_Int* pNj = edges.get_facets_ptr(Ej);
        const E_Int* pNjp1 = edges.get_facets_ptr(Ejp1);

        edges_are_sorted = (pNj[0] == pNjp1[0] || pNj[0] == pNjp1[1] || pNj[1] == pNjp1[0] || pNj[1] == pNjp1[1]); //consecutive edge sharing a node
        if (!edges_are_sorted) break;

        E_Int ej[] = {pNj[0], pNj[1]};
        E_Int ejp1[] = {pNjp1[0], pNjp1[1]};

        if (ej[0] == ejp1[0])
        {
          //  ej[1]    ej[0]
          //    +--------+
          //             +----------+
          //           ejp1[0]    ejp1[1]
          std::swap(ej[0], ej[1]);
        }
        else if (ej[1] == ejp1[1])
        {
          //  ej[0]    ej[1]
          //    +--------+
          //             +----------+
          //           ejp1[1]    ejp1[0]
          std::swap(ejp1[0], ejp1[1]);
        }
        else if (ej[0] == ejp1[1])
        {
          //  ej[1]    ej[0]
          //    +--------+
          //             +----------+
          //           ejp1[1]    ejp1[0]
          std::swap(ej[0], ej[1]);
          std::swap(ejp1[0], ejp1[1]);
        }
        if (j==0)
        {
          snodes.push_back(ej[0]);
          // std::cout << "s0 : " << ej[0] << std::endl;
        }

        assert (ej[1] == ejp1[0]);
        snodes.push_back(ej[1]);

        // std::cout << "s1 : " << ej[1] << std::endl;
      }

      if (!edges_are_sorted)
      {
        cB.clear();
        for (E_Int j=0; j < nb_edges; ++j)
        {
          E_Int Ej = pEs[j]-1;
          const E_Int* pN = edges.get_facets_ptr(Ej);
          cB.pushBack(pN, pN+2);
        }
        // sort the nodes
        BARSplitter::getSortedNodes(cB, snodes);
      }
      
      molec.clear();
      molec.push_back(nb_edges);// stride : nb_edges = nb nodes
      
      if (snodes.size() != (size_t)nb_edges) //degen
      {
        //std::cout << "degen : " << snodes.size() << "/" << nb_edges << std::endl;
        continue;
      }
      
      ++nb_polys;
      
      molec.insert(molec.end(), snodes.begin(), snodes.end());
      
      //std::cout << molec.size() << std::endl;
      
      c.pushBack(molec);
    }
    
    c[0]=nb_polys;
    c[1]=c.getSize()-2;

    ngon_unit pgs(&c[0]);
    reorient_connex_PGs(pgs, false/*i.e. trust first PG orientation*/);
    
    c.clear();
    ngon_t<Connectivity_t> ngo(pgs, true);
    ngo.export_to_array(c);

    
    return 0;
    
  }
  
  // Converts a SURFACE formatted as Edges/PGs into a ngon_unit (PHs is cleared) representing the PGs
  // WARNING : inconsistent orientation of the surface upon exit
  // Morse version
  void export_surfacic_view(Connectivity_t& FN, Connectivity_t& xFN)
  {
    std::vector<E_Int> fn, xfn;

    const ngon_unit& polys = PHs;
    const ngon_unit& edges = PGs;
    
    edges.updateFacets();
    polys.updateFacets();
    
    E_Int nb_edges = edges.size();
    E_Int nb_pgs   = polys.size();

    // Surfacic ngon ?
    for (E_Int i=0; (i<nb_edges); ++i){
      if (edges.stride(i) != 2)
        return;
    }
    
    xfn.resize(nb_pgs+1);
    xfn[0]=0;
    
    //
    std::vector<E_Int> snodes;
    K_FLD::IntArray cB;

    //
    for (E_Int i=0; i<nb_pgs; ++i){
            
      E_Int nb_es = polys.stride(i);
      const E_Int* pEs = polys.get_facets_ptr(i);
      
      cB.clear();
      for (size_t j=0; j < nb_es; ++j)
      {
        E_Int Ej = pEs[j]-1;
        const E_Int* pN = edges.get_facets_ptr(Ej);
        cB.pushBack(pN, pN+2);
      }
                  
      // sort the nodes
      BARSplitter::getSortedNodes(cB, snodes);
      
      //if (snodes.size() != nb_es) //degen
        //continue;

      fn.insert(fn.end(), snodes.begin(), snodes.end());
      xfn[i+1]=fn.size();
    }
    
    FN = Connectivity_t(fn.size(), 1, &fn[0]);
    xFN = Connectivity_t(xfn.size(), 1, &xfn[0]);
  }
  
  /// reverse function : from surfacic FN to EFN
  static void export_surfacic_FN_to_EFN(const ngon_unit& PGS, Connectivity_t& c)
  {
    ngon_t ngo;

    ngo.PHs = PGS; //same structure
    
    E_Int E[2], count(1);
    E_Int nb_pgs = PGS.size();
    
    for (E_Int i=0; i < nb_pgs; ++i)
    {
      const E_Int* nodes = PGS.get_facets_ptr(i);
      E_Int nb_nodes = PGS.stride(i);
      
      E_Int* edges = ngo.PHs.get_facets_ptr(i);

      for (E_Int n=0; n < nb_nodes; ++n)
      {
        E[0] = *(nodes+n);
        E[1] = *(nodes+(n+1)%nb_nodes);
        
        ngo.PGs.add(2, E);
        *(edges+n)=count++;
      }
    }
    
    ngo.remove_duplicated_edges();
    
    Vector_t<E_Int> phnids, pgnids;
    ngo.remove_unreferenced_pgs(pgnids, phnids);
    
    ngo.export_to_array(c);
  }
  
  /// reverse function : from surfacic FN to EFN
  static void export_surfacic_FN_to_EFN(const ngon_unit& PGS, ngon_t& ngo)
  {
    ngo.clear();

    ngo.PHs = PGS; //same structure
    
    E_Int E[2], count(1);
    E_Int nb_pgs = PGS.size();
    
    for (E_Int i=0; i < nb_pgs; ++i)
    {
      const E_Int* nodes = PGS.get_facets_ptr(i);
      E_Int nb_nodes = PGS.stride(i);
      
      E_Int* edges = ngo.PHs.get_facets_ptr(i);

      for (E_Int n=0; n < nb_nodes; ++n)
      {
        E[0] = *(nodes+n);
        E[1] = *(nodes+(n+1)%nb_nodes);
        
        ngo.PGs.add(2, E);
        *(edges+n)=count++;
      }
    }
    
    ngo.remove_duplicated_edges();
    
    Vector_t<E_Int> phnids, pgnids;
    ngo.remove_unreferenced_pgs(pgnids, phnids);
  }

  ///
  inline bool empty() const {return (PGs._NGON.empty() || PHs._NGON.empty());}
  ///
  inline void clear() {PGs.clear(); PHs.clear();}
  ///
  inline const E_Int* begin() const {return PGs.begin();}
    
  ///
  static void select_phs(const ngon_t& NGin, const Vector_t<bool>& keep, Vector_t<E_Int>& new_pg_ids, ngon_t& NGout)
  {
    //WARNING : change the PG numerotation (even if keep is always true)
    
    new_pg_ids.clear();
    NGout.clear();
    
    Vector_t<E_Int> stacked_pg_ids;
    __get_selected_ids(NGin.PHs, keep, new_pg_ids, stacked_pg_ids);
    if (stacked_pg_ids.empty()) return;
    
    __create_selected_phs(NGin.PHs, keep, new_pg_ids, NGout.PHs);
    __create_selected_pgs(NGin.PGs, stacked_pg_ids, NGout.PGs);
  }
  
  ///
  static void split_phs(const ngon_t& NGin, const Vector_t<bool>& flag, ngon_t& cng_ok, ngon_t& cng_nok)
  {
    //WARNING : change the PG numerotation (even if flag is always true)
    
    ngon_t tmp(NGin); //use tmp because cng_ok or cng_nok could be NGin
    Vector_t<bool> nflag(flag); //copy at the beginnning because flag might change after compress (if it is NGin._external)
    K_CONNECT::IdTool::negative(nflag);
    Vector_t<E_Int> new_pg_ids;
    select_phs(tmp, flag, new_pg_ids, cng_ok);
    select_phs(tmp, nflag, new_pg_ids, cng_nok);
  }
  
  ///
  void append(const ngon_t& NGin)
  {     
    if (NGin.empty())
      return;
    
    E_Int shiftv = PGs.size(); //to apply to the PHin
    
    // PGs    
    PGs.append(NGin.PGs);
    // PHs
    ngon_unit tmp(NGin.PHs);
    tmp.shift(shiftv);
    PHs.append(tmp);
  }
  
  ///
  E_Int flag_externals(E_Int FLAG)
  {
    if (PHs.size() == PGs.size()) // if one pg per ph : it's a surface so all PGs/PHs are skin ones.
    {
      PHs._type.resize(PHs.size(), FLAG);
      PGs._type.resize(PGs.size(), FLAG);
    }
    else
    {
      flag_external_pgs(FLAG);
      flag_external_phs(FLAG);
    }
  
    return 0;
  }
  
  ///
  E_Int flag_external_pgs(E_Int FLAG)
  {
    PHs.updateFacets();
    
    E_Int nb_faces = PGs.size();
    E_Int nb_elts = PHs.size();
    Vector_t<E_Int> face_count(nb_faces, 0);
  
    // Loop through the elements and increment face_count
    for (E_Int i = 0; i < nb_elts; ++i)
    {
      E_Int nb_facets = PHs.stride(i);
      for (E_Int j = 0; j < nb_facets; ++j)
      {
        const E_Int& Fi = PHs.get_facet(i,j);
        //std::cout << "Fi : " << Fi << std::endl;
         ++face_count[Fi-1];
      }
    }
  
    // External faces are those with a count equal to 1.
    PGs._type.clear();
    PGs._type.resize(nb_faces, INNER);
    for (E_Int i = 0; i < nb_faces; ++i)
      PGs._type[i] = (face_count[i] == 1) ? FLAG : INNER;
  
    return 0;
  }
  
  ///
  E_Int flag_external_phs(E_Int FLAG)
  {
    PHs.updateFacets();
    
    E_Int nb_facets, nb_elts(PHs.size());
    
    // External elements are those with at least one external facets
    PHs._type.clear();
    PHs._type.resize(nb_elts, INNER);
    for (E_Int i = 0; i < nb_elts; ++i)
    {
      nb_facets = PHs.stride(i);
      for (E_Int j = 0; j < nb_facets; ++j)
      {
        const E_Int& pgi = PHs.get_facet(i,j);
        if (PGs._type[pgi - 1] == FLAG)
        {
          PHs._type[i] = FLAG;
          break;
        }
      }
    }
  
    return 0;
  }
  
  ///
  template <typename T>
  void flag_facets_of_elts(const std::vector<E_Int>& elts, E_Int index_start, std::vector<T> & flag) const
  {  
    PHs.updateFacets();
    PGs.updateFacets();
  
    //E_Int nb_elts = PHs.size();
  
    flag.clear();
    flag.resize(PGs.size(), 0);
  
    E_Int nb_in_set(elts.size());
    // 
    for (E_Int i = 0; i < nb_in_set; ++i)
    {
      E_Int ei =  elts[i];
      E_Int nb_facets = PHs.stride(ei);
      const E_Int* facets = PHs.get_facets_ptr(ei);
    
      for (E_Int n=0; n < nb_facets; ++n)
        flag[*(facets+n) - index_start] = 1;
    }
  }
  
  template <typename T>
  void flag_facets_of_elts(E_Int index_start, std::vector<T> & flag) const
  {  
    PHs.updateFacets();
    PGs.updateFacets();
  
    E_Int nb_elts = PHs.size();
  
    flag.clear();
    flag.resize(PGs.size(), 0);

    // 
    for (E_Int ei = 0; ei < nb_elts; ++ei)
    {
      E_Int nb_facets = PHs.stride(ei);
      const E_Int* facets = PHs.get_facets_ptr(ei);
    
      for (E_Int n=0; n < nb_facets; ++n)
        flag[*(facets+n) - index_start] = 1;
    }
  }
  
  ///
  void type_facets_as_elts_type(E_Int ignore_type, E_Int shft, E_Int index_start)
  {  
    PHs.updateFacets();
    PGs.updateFacets();
  
    E_Int nb_elts = PHs.size();
  
    PGs._type.resize(PGs.size(), ignore_type);
  
    // 
    for (E_Int i = 0; i < nb_elts; ++i)
    {
      const E_Int& typ = PHs._type[i];
      if (typ == ignore_type) continue;
      
      E_Int nb_facets = PHs.stride(i);
      const E_Int* facets = PHs.get_facets_ptr(i);
    
      for (E_Int n=0; n < nb_facets; ++n)
        PGs._type[*(facets+n) - index_start] = typ + shft;
    }
  }
  
  ///
  E_Int remove_pgs(const Vector_t<E_Int>& toremove)
  {
    Vector_t<E_Int> pgnids, phnids;
    return remove_pgs(toremove, pgnids, phnids);
  }
  ///
  E_Int remove_pgs(const Vector_t<E_Int>& toremove, Vector_t<E_Int>& pgnids, Vector_t<E_Int>& phnids)
  {
    // WARNING : index start dependant : 0
    if (toremove.empty())
      return 0;
    //
    /*std::cout << "stats before remove_entities" << std::endl;
    std::cout << "nb pgs : " << PGs.size() << std::endl;
    std::cout << "nb phs : " << PHs.size() << std::endl;
    std::cout <<"max face id " << PHs.get_facets_max_id() << std::endl;
    std::cout << "//////////////////////////////////" << std::endl;*/
        
    pgnids.clear();
    PGs.remove_entities(toremove, pgnids); //0-based nids
    PGs.updateFacets();
    
    
    //std::cout << "nids size " << nids.size() << std::endl;

    PHs.remove_facets(pgnids, phnids); //0-based nids
    
    /*std::cout << "stats after remove_entities" << std::endl;
    std::cout << "nb pgs : " << PGs.size() << std::endl;
    std::cout << "nb phs : " << PHs.size() << std::endl;
    std::cout <<"max face id  " << PHs.get_facets_max_id() << std::endl;
    std::cout << "//////////////////////////////////" << std::endl;*/
    
    //std::cout << "remove_pgs : 3 " << std::endl;
    return 0;
  }
  
  ///
  E_Int build_ph_neighborhood(ngon_unit& neighbor, std::vector<bool>& wall)
  {
    neighbor.clear();
    PHs.updateFacets();
    neighbor=PHs; //same molecules as the connectivity
    neighbor.reset_facets();
    E_Int maxPG = PHs.get_facets_max_id();
    
    Vector_t<E_Int> neigh(maxPG, IDX_NONE);
    E_Int nb_elts(PHs.size()), nb_facets, fid;
    E_Int count(0); // 2-pass required to fill both sides of the neighbor info
    while (count++ != 2)
    {
      for (E_Int eid = 0; eid < nb_elts; ++eid)
      {
        nb_facets = PHs.stride(eid);
        for (E_Int j = 0; j < nb_facets; ++j)
        {
          fid = PHs.get_facet(eid,j)-1;
          if (neigh[fid] != IDX_NONE) //2nd pass
            neighbor.get_facet(eid,j) = neigh[fid];
          if (!wall[fid]) neigh[fid] = eid;
        }
      }
    }
    return 0;
  }
  
  ///
  E_Int build_ph_neighborhood(ngon_unit& neighbor) const
  {
    neighbor.clear();
    
    if (PHs._NGON.empty())
      return 1;
    
    PHs.updateFacets();
    neighbor = PHs; //same molecules as the connectivity
    neighbor.reset_facets();
    E_Int maxPG = PHs.get_facets_max_id();
    if (maxPG < 0) return 1;

    Vector_t<E_Int> neigh(maxPG, IDX_NONE);
    E_Int nb_elts(PHs.size()), nb_facets, fid;
    E_Int count(0); // 2-pass required to fill both sides of the neighbor info
    while (count++ != 2)
    {
      for (E_Int eid = 0; eid < nb_elts; ++eid)
      {
        nb_facets = PHs.stride(eid);
        for (E_Int j = 0; j < nb_facets; ++j)
        {
          fid = PHs.get_facet(eid, j)-1;
          E_Int& eidn = neigh[fid];
          E_Int& Kn = neighbor.get_facet(eid, j);
          if (eidn != IDX_NONE && eidn != eid)
            Kn = eidn;
          neigh[fid] = eid;
        }
      }
    }
    return 0;
  }
 
  ///
  static void refine_pgs(const std::map<K_MESH::NO_Edge, Vector_t<E_Int> >& edge_to_refined_edge, ngon_unit& PGs)
  {
    ngon_unit refinedPGs;
    Vector_t<E_Int> pg_molec;
    
    for (E_Int i = 0; i < PGs.size(); ++i)
    {
      refine_pg(PGs.get_facets_ptr(i), PGs.stride(i), edge_to_refined_edge, pg_molec);
      refinedPGs.add(pg_molec);
    }
    refinedPGs._type = PGs._type;  //work around to not loose them
    refinedPGs._ancEs = PGs._ancEs;//work around to not loose them
    PGs = refinedPGs; //dirty
  }
  
  ///
  template<typename IntCont>
  static bool refine_pg
  (const E_Int* pg, E_Int sz, const std::map<K_MESH::NO_Edge, IntCont >& edge_to_refined_edge, Vector_t<E_Int>& pg_molec)
  {
    E_Int Ni, Nj;
    K_MESH::NO_Edge Ei;
    
    pg_molec.clear();
    pg_molec.resize(1, 0);
    for (E_Int i = 0; i < sz; ++i)
    {
      Ni = *(pg+i);
      Nj = *(pg + (i+1)%sz);
      Ei.setNodes(Ni, Nj);
      
      auto it = edge_to_refined_edge.find(Ei);
      if (it != edge_to_refined_edge.end())
      {
        const auto& nodes = it->second;
        if (Ni == Ei.node(0)) //same orientation
          for (size_t j=0; j< nodes.size()-1; ++j)
            pg_molec.push_back(nodes[j]);
        else
          for (E_Int j=nodes.size()-1; j>0; --j)
            pg_molec.push_back(nodes[j]);
      }
      else
      {
        pg_molec.push_back(Ni);
      }
    }
    
    pg_molec[0]=pg_molec.size()-1;
    return (sz < pg_molec[0]);
  }

  /// Remove hatty PGs (chinese hats when 3 nodes) from a mesh => open mesh to fix with close_phs
  static void remove_hatty_PGs(ngon_t& ngio, const K_FLD::FloatArray& crd, double ARTOL)
  {
    ngio.PGs.updateFacets();

    E_Int nb_pgs = ngio.PGs.size();
    Vector_t<E_Int> hatties;
    for (E_Int i = 0; i < nb_pgs; ++i)
    {
      const E_Int* nodes = ngio.PGs.get_facets_ptr(i);
      E_Int nnodes = ngio.PGs.stride(i);

      E_Int is(IDX_NONE), ie(IDX_NONE);
      if (!K_MESH::Polygon::is_hatty(crd, nodes, nnodes, 1, is, ie, ARTOL)) continue;
      hatties.push_back(i);
    }

#ifdef DEBUG_NGON_T
    NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PGs("hatties", crd, ngio.PGs, hatties, false);
#endif

    Vector_t<E_Int> pgnids, phnids;
    ngio.remove_pgs(hatties, pgnids, phnids); //sync the PHs accordingly
  }
  
  /// Close PHs having free edges (GEOM algo).
  static bool close_phs(ngon_t& ngio, const K_FLD::FloatArray& crd, std::vector<E_Int>* phs_to_process = nullptr)
  {
    ngio.PGs.updateFacets();
    ngio.PHs.updateFacets();

    bool had_open_cells{ false };

    using dint_t = std::deque<E_Int>;

    std::vector < std::pair<E_Float, E_Int>> sorter;
    std::map<K_MESH::NO_Edge, dint_t > edge_to_refine_nodesTMP, edge_to_refine_nodes;
    using edge_t = std::set<K_MESH::NO_Edge>;
    edge_t free_edges;
    using vedge_t = std::vector<K_MESH::NO_Edge>;
    vedge_t vfree_edges;
    std::vector<dint_t> chains;
    std::vector<E_Int> edgecol;

    E_Int nb_pgs = ngio.PGs.size();
    E_Int nb_phs(ngio.PHs.size());

    std::vector<bool> processPG, processPH;

    processPH.resize(nb_phs, (phs_to_process == nullptr)); // process all if nothing upon entry
    
    if (phs_to_process) // just set true for them
    {
      for (size_t k = 0; k < phs_to_process->size(); ++k)
        processPH[(*phs_to_process)[k]] = true;
    }

    bool carry_on(true);
    while (carry_on)
    {
      carry_on = false;
      edge_to_refine_nodes.clear();

      // Build the refined edges       
      for (E_Int i = 0; i < nb_phs; ++i)
      {
        if (!processPH[i]) continue;

        // closed ?
        if (K_MESH::Polyhedron<0>::is_closed(ngio.PGs, ngio.PHs.get_facets_ptr(i), ngio.PHs.stride(i), free_edges)) continue;

        carry_on = true;
        had_open_cells = true;

#ifdef DEBUG_NGON_T
        //if (i == 169804) medith::write("PH.mesh", crd, ngio, i);
#endif

        BARSplitter::split_eset_into_manifold_chains<edge_t>(free_edges, chains, edgecol);

        // now we have openLines and closed loops, so convert closed ones into 2 openLines based on worst angle nodes
        E_Int nchains = (E_Int)chains.size(); // initial because following sliiting will append it
        for (E_Int c = 0; c < nchains; ++c)
        {
          E_Int nnodes = chains[c].size();
          if (nnodes == 2) continue; // taken into account later when adding the edge made of edge chains ends to vfree_edges
                                     // get the edges assigned to that color
          vfree_edges.clear();
          E_Int j = 0;
          for (auto ite = free_edges.begin(); ite != free_edges.end(); ++ite, ++j) if (edgecol[j] == c)vfree_edges.push_back(*ite);

          bool is_closed_loop = (*chains[c].begin() == *chains[c].rbegin());

          if (is_closed_loop)
          {
            chains[c].pop_back(); // remove redundant end
            --nnodes;

            E_Int is(IDX_NONE), ie(IDX_NONE);
            if (!K_MESH::Polygon::is_hatty(crd, chains[c].begin(), nnodes, 1, is, ie)) continue;
            // put spiky nodes at the right place (start & end). following works beacause is < ie
            if (is != 0) std::swap(chains[c][0], chains[c][is]);
            if (ie != nnodes - 1) std::swap(chains[c][nnodes - 1], chains[c][ie]);
            // sort (GEOMETRIC)
            E_Int Ne(chains[c][nnodes - 1]);
            NUGA::MeshTool::reorder_nodes_on_edge<std::deque<E_Int>, 3>(crd, chains[c], 1, sorter);
            if (chains[c][nnodes - 1] != Ne) continue; // not handled : spiky nodes are not the extrema ones.
          }

          if (!is_closed_loop) vfree_edges.push_back(K_MESH::NO_Edge(chains[c][0], chains[c][nnodes - 1]));

          BARSplitter::compute_refinement<vedge_t, dint_t>(vfree_edges, chains[c], edge_to_refine_nodesTMP);

          // append it to global
          for (auto& itLOC : edge_to_refine_nodesTMP)
          {
            auto ii = edge_to_refine_nodes.find(itLOC.first);
            if (ii == edge_to_refine_nodes.end())// new
              edge_to_refine_nodes.insert(std::make_pair(itLOC.first, itLOC.second)); //set false to tell "no reorder required"
            else // append
              ii->second.insert(ii->second.end(), ALL(itLOC.second));
          }
        }
      }

      if (!carry_on) break;

      // complete refinement with ends and eventually sort in between if more than 1 point
      for (auto& e : edge_to_refine_nodes)
      {
        e.second.push_front(e.first.node(0));
        if (e.second.size() > 2)
          NUGA::MeshTool::reorder_nodes_on_edge<std::deque<E_Int>, 3>(crd, e.second, 1, sorter);
        e.second.push_back(e.first.node(1));
      }

      // Refine the PGs ( <=> close the PHs )
      Vector_t<E_Int> pg_molec;
      ngon_unit refinedPGs;

      processPG.clear();
      processPG.resize(nb_pgs, false);

      for (E_Int PGi = 0; PGi < nb_pgs; ++PGi)
      {
        processPG[PGi] = refine_pg(ngio.PGs.get_facets_ptr(PGi), ngio.PGs.stride(PGi), edge_to_refine_nodes, pg_molec);
        refinedPGs.add(pg_molec);
      }

      // update PGs ngon unit
      refinedPGs._type = ngio.PGs._type;  // hack to preserve flags (externality)
      refinedPGs._ancEs = ngio.PGs._ancEs;// hack
      ngio.PGs = refinedPGs;
      ngio.PGs.updateFacets();

      // build PH list for next iter
      ngio.flag_PHs_having_PGs(processPG, processPH);
    }

    return had_open_cells;
  }
  
  /// Change the node indice to reference the same when duplicated exist
  E_Int join_phs(const K_FLD::FloatArray& coord, E_Float tolerance = EPSILON, bool do_omp=false)
  {
    if (PGs.size() == 0) return 0;
    
    PGs.updateFacets();
    
    K_FLD::ArrayAccessor<K_FLD::FloatArray> ca(coord);
    Vector_t<E_Int> nids;
    E_Int nb_merges = ::merge(ca, tolerance, nids, do_omp);
    if (nb_merges)
      PGs.change_indices(nids);
    return nb_merges;
  }

  /// Change the node indice to reference the same when duplicated exist (Relative tolerance version)
  E_Int join_phs(const K_FLD::FloatArray& coord, const std::vector<E_Float>&nodal_metric2,  E_Float RTOL, bool do_omp=false)
  {
    if (PGs.size() == 0) return 0;

    PGs.updateFacets();

    K_FLD::ArrayAccessor<K_FLD::FloatArray> ca(coord);
    Vector_t<E_Int> nids;
    E_Int nb_merges = ::merge(ca, nodal_metric2, RTOL, nids, do_omp);
    if (nb_merges)
      PGs.change_indices(nids);
    return nb_merges;
  }
  
  /// Change the node indice to reference the same when duplicated exist (FldArrayF)
  void join_phs(const K_FLD::FldArrayF& coord, E_Int px, E_Int py, E_Int pz, E_Float tolerance = EPSILON, bool do_omp=false)
  {
    if (PGs.size() == 0)
      return;
    
    PGs.updateFacets();
    
    K_FLD::ArrayAccessor<K_FLD::FldArrayF> ca(coord, px, py, pz);
    Vector_t<E_Int> nids;
    E_Int nb_merges = ::merge(ca, tolerance, nids, do_omp);
    if (nb_merges)
      PGs.change_indices(nids, true/*zero based*/);
  }

  void join_phs(const K_FLD::FldArrayF& coord, E_Int px, E_Int py, E_Int pz, const std::vector<E_Float>&nodal_metric2, E_Float RTOL, bool do_omp=false)
  {
    if (PGs.size() == 0)
      return;

    PGs.updateFacets();

    K_FLD::ArrayAccessor<K_FLD::FldArrayF> ca(coord, px, py, pz);
    Vector_t<E_Int> nids;
    E_Int nb_merges = ::merge(ca, nodal_metric2, RTOL, nids, do_omp);
    if (nb_merges)
      PGs.change_indices(nids, true/*zero based*/);
  }

  /*E_Int join_phs2(const K_FLD::FldArrayF& coord, E_Int px, E_Int py, E_Int pz, const std::vector<E_Float>&nodal_metric2, E_Float RTOL)
  {
    if (PGs.size() == 0) return 0;

    PGs.updateFacets();

    flag_external_pgs(INITIAL_SKIN);

    std::vector<E_Int> anodes;

    {
      // get the exterior nodes
      E_Int nb_pgs = PGs.size();
      for (E_Int i = 0; i < nb_pgs; ++i)
      {
        if (PGs._type[i] != 1) continue;
        const E_Int* nodes = PGs.get_facets_ptr(i);
        E_Int nb_nodes = PGs.stride(i);
        anodes.insert(anodes.end(), nodes, nodes + nb_nodes);
      }
      
      std::set<E_Int> attaching_nodes;
      attaching_nodes.insert(ALL(anodes)); // make the list unique
      anodes.clear();
      anodes.insert(anodes.end(), ALL(attaching_nodes)); // put it back in a vector
      K_CONNECT::IdTool::shift(anodes, -1); // 0-based required for __merge
    }

    K_FLD::ArrayAccessor<K_FLD::FldArrayF> ca(coord, px, py, pz);

    Vector_t<E_Int> nids;
    E_Int nb_merges = ::__merge(ca, nodal_metric2, RTOL, anodes, anodes, nids);

    if (nb_merges)
      PGs.change_indices(nids);
    return nb_merges;
  }*/
  /*   
  ///
  E_Int unique_nodes_ph(E_Int PHi, Vector_t<E_Int>& indices, bool zerobased=false) const
  {
    PHs.updateFacets();
    PGs.updateFacets();
    
    indices.clear();
    
    std::set<E_Int> tmp;
    
    E_Int nb_pgs(PHs.stride(PHi)), nb_nodes;
   
    for (E_Int j = 0; j < nb_pgs; ++j)
    {
      const E_Int& pgi = PHs.get_facet(PHi,j);
      nb_nodes = PGs.stride(pgi-1);
      
      for (E_Int k = 0; k < nb_nodes; ++k)
        tmp.insert(PGs.get_facet(pgi-1,k));
    }
    
    indices.insert(indices.end(), tmp.begin(), tmp.end());
    if (zerobased)
      ngon_unit::shift(indices, -1);
    
    return 0;
  }*/


  E_Int nodes_ph(E_Int PHi, Vector_t<E_Int>& indices, bool zerobased = false) const
  {
    PHs.updateFacets();
    PGs.updateFacets();

    indices.clear();

    E_Int nb_pgs(PHs.stride(PHi)), nb_nodes;

    for (E_Int j = 0; j < nb_pgs; ++j)
    {
      const E_Int& pgi = PHs.get_facet(PHi, j);
      nb_nodes = PGs.stride(pgi - 1);

      for (E_Int k = 0; k < nb_nodes; ++k)
        indices.push_back(PGs.get_facet(pgi - 1, k));
    }

    if (zerobased)
      ngon_unit::shift(indices, (E_Int)-1);

    return 0;
  }

  ///
  //E_Int unique_nodes_phs(Vector_t<E_Int>& indices, bool zerobased=false) const
  //{
  //  PHs.updateFacets();
  //  indices.clear();
  //  
  //  Vector_t<E_Int> inds;
  //  std::set<E_Int> tmp;
  //  
  //  E_Int nb_facets, nb_elts(PHs.size());
  //  
  //  for (E_Int i = 0; i < nb_elts; ++i)
  //  {
  //    nodes_ph(i, inds, false); //false to do the shift only once
  //    tmp.insert(inds.begin(), inds.end());
  //  }
  //  
  //  indices.insert(indices.end(), tmp.begin(), tmp.end());
  //  if (zerobased)
  //    ngon_unit::shift(indices, -1);
  //  
  //  return 0;
  //}
  
  ///
  E_Int unique_nodes_ph(E_Int PHi, Vector_t<E_Int>& indices, bool zerobased=false) const
  {
    PHs.updateFacets();
    indices.clear();
    
    Vector_t<E_Int> inds;
    std::set<E_Int> tmp;
    
    nodes_ph(PHi, inds, false); //false to do the shift only once
    tmp.insert(inds.begin(), inds.end());
   
    indices.insert(indices.end(), tmp.begin(), tmp.end());
    if (zerobased)
      ngon_unit::shift(indices, E_Int(-1));
    
    return 0;
  }
  
  /// triangulate specified PGs
  template <typename TriangulatorType>
  static E_Int triangulate_pgs
    (const ngon_unit& PGs, const K_FLD::FloatArray& coord, K_FLD::IntArray& connectT3, Vector_t<E_Int>& colors, bool do_not_shuffle, bool improve_quality, const Vector_t<bool>* process=0) 
  {
    connectT3.clear();
    colors.clear();
    
    TriangulatorType dt;

    PGs.updateFacets();

    E_Int sz = PGs.size(), err(0);
    E_Int i=0;
    for (; (i<sz) && !err ; ++i)
    {
      if (process && (*process)[i] == false)
        continue;

#ifdef DEBUG_TRIANGULATOR
      if (i == 15138)
        dt.dbg_enabled=true;
      else
        dt.dbg_enabled=false;
#endif

      err = K_MESH::Polygon::triangulate(dt, coord, PGs.get_facets_ptr(i), PGs.stride(i), 1/*index start*/, connectT3, do_not_shuffle, improve_quality);
      if (err)
      {
        K_MESH::Polygon::cvx_triangulate(coord, PGs.get_facets_ptr(i), PGs.stride(i), 0/*ibest*/, 1/*index start*/, connectT3);
        err = 0;
      }
      
      colors.resize(connectT3.getSize(), i);
    }
    
    if (err) { //i-1 because have been incremented when exiting the loop
      std::cout << "triangulate_pgs error at PG : " << i - 1 << " (" << PGs.stride(i - 1) << " nodes)" << std::endl;
#ifdef DEBUG_TRIANGULATOR
      std::ostringstream o;
      o << "./pg_" << i - 1 << ".mesh";
      medith::write(o.str().c_str(), coord, PGs.get_facets_ptr(i - 1), PGs.stride(i - 1), 1);
#endif
    }

    return err;
  }

  /// star-triangulate specified PGs
  static E_Int triangulate_pgs
    (const ngon_unit& PGs, K_FLD::FloatArray& coord, K_FLD::IntArray& connectT3, Vector_t<E_Int>& colors, const Vector_t<bool>* process=0) 
  {
    connectT3.clear();
    colors.clear();

    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
    

    PGs.updateFacets();

    E_Float P[3];

    E_Int sz = PGs.size();

    for (E_Int j=0; j < sz; ++j)
    {
      if (process && (*process)[j] == false)
        continue;

      coord.pushBack(P,P+3);
      E_Int C = coord.cols()-1;

      //K_MESH::Polygon::centroid<3>(coord, PGs.get_facets_ptr(j), PGs.stride(j), 1, coord.col(C));
      K_MESH::Polygon::iso_barycenter(coord, PGs.get_facets_ptr(j), PGs.stride(j), 1, coord.col(C));

      acrd_t acrd(coord);
      K_MESH::Polygon::star_triangulate(acrd, PGs.get_facets_ptr(j), PGs.stride(j), 1, C, connectT3);

      colors.resize(connectT3.getSize(), j);
    }

    return 0;
  }

  ///
  template <typename Coordinate_t>
  static void compact_to_used_nodes (ngon_unit& PGs, Coordinate_t& coord)
  {
    Vector_t<E_Int> newIDs;
    compact_to_used_nodes(PGs, coord, newIDs);
  }

  ///
  template <typename Coordinate_t>
  static void compact_to_used_nodes(ngon_unit& PGs, Coordinate_t& coord, std::vector<E_Int>& nids)
  {
    nids.clear();

    Vector_t<E_Int> unic_nodes;
    PGs.unique_indices(unic_nodes);

    Vector_t<bool> flag(coord.getSize(), false);

    size_t sz(unic_nodes.size());
    for (size_t i = 0; i < sz; ++i) flag[unic_nodes[i] - 1] = true; // indices start at 1.

    E_Int nb_removed = Coordinate_t::compact(coord, flag, nids);
    if (nb_removed == 0) return;

    PGs.change_indices(nids);
  }
  
  ///
  static K_FLD::FloatArray compress_to_used_nodes (ngon_unit& PGs, const K_FLD::FloatArray& coord, std::vector<E_Int>& nids)
  {
    Vector_t<E_Int> unic_nodes;
    PGs.unique_indices(unic_nodes);

    K_CONNECT::IdTool::shift(unic_nodes, -1);

    K_FLD::FloatArray lcrd = K_CONNECT::IdTool::compress_(coord, unic_nodes, 0/*idx start*/);
     
    if (lcrd.cols() == coord.cols())
    {
      K_CONNECT::IdTool::init_inc(nids, coord.cols());
      return coord;
    }

    K_CONNECT::IdTool::reverse_indirection(unic_nodes, nids);
    nids.resize(coord.cols(), IDX_NONE);
  
    PGs.change_indices(nids); 
    return lcrd;
  }
   
  ///
  template <typename CoordAccType>
  bool remove_duplicated_pgs (const CoordAccType& coord, Vector_t<E_Int>& pgnids, bool do_omp)
  {
    //detect the matching and replace the ids with the first one found
    bool found = replace_duplicated_pgs(coord,pgnids, do_omp); 
    if (found)
      PHs.remove_duplicated(); //clean
    return found;
  }
  
  bool remove_duplicated_edges()
  {
    bool found = replace_duplicated_edges();
    if (found)
      PHs.remove_duplicated();//clean
    return found;
  }
  
  bool remove_duplicated_nodes()
  {
    bool found = replace_duplicated_nodes();
    if (found)
      PHs.remove_duplicated();//clean
    return found;
  }

  template <typename CoordAccType>
  bool detect_duplicated_pgs(const CoordAccType& coord, std::vector<E_Int>& nids, bool check_node_ids)
  {
    nids.clear();

    if (PGs.size() == 0)
      return false;

    PGs.updateFacets();
    E_Int nb_pgs(PGs.size());

    E_Float G[3];
    K_FLD::FloatArray barys;
    barys.reserve(3, PGs.size());

    Vector_t<E_Int> ISO(PGs.size(), IDX_NONE), PGI;

    // Loop over PH rather than directly over PGs to consider only used PGs
    E_Int count{ 0 };
    for (E_Int i = 0; i < nb_pgs; ++i)
    {
      E_Int s = PGs.stride(i);
      K_MESH::Polygon PGi(PGs.get_facets_ptr(i), s, -1/*to make it zero based*/);
      PGi.iso_barycenter<CoordAccType, 3>(coord, G);
      barys.pushBack(G, G + 3);

      ISO[i] = count++;
    }

    K_CONNECT::IdTool::reverse_indirection(ISO, PGI);

    // Merge
    K_FLD::ArrayAccessor<K_FLD::FloatArray> cab(barys);
    E_Int nmerges = ::merge(cab, EPSILON, nids);

    if (!nmerges) //no duplicated faces.
      return false;

    if (!check_node_ids)
      return (nmerges > 0);

    // check if matching isoG means really identical PGs
    bool found = false;
    for (size_t i = 0; i < nids.size(); ++i)
    {
      if (nids[i] == (E_Int)i)
        continue;

      const E_Int &pgi1 = PGI[i];
      const E_Int &pgi2 = PGI[nids[i]];

      if (PGs.stride(pgi1) != PGs.stride(pgi2)) // fast treatment : not the same number of nodes
      {
        nids[i] = i;
        continue;
      }
      if (!K_CONNECT::IdTool::equal(PGs.get_facets_ptr(pgi1), PGs.get_facets_ptr(pgi2), PGs.stride(pgi1), true/*permut*/, false/*strict orient*/))
        nids[i] = i;
      else
        found = true;
    }

    return found;
  }

  ///
  bool collapse_micro_edge(const K_FLD::FloatArray& crd, double edge_ratio, double Lmax, std::vector<E_Int>& nids)
  {
    bool has_changed=false;

    if (edge_ratio > 1.) edge_ratio = 1./edge_ratio;
    if (edge_ratio <= 0.) return false;

    double edge_ratio2 = edge_ratio*edge_ratio;

    nids.clear();
    K_CONNECT::IdTool::init_inc(nids, crd.cols()); 

    K_FLD::FloatArray L;
    NUGA::MeshTool::computeIncidentEdgesSqrLengths(crd, PGs, L);

    E_Int npgs = PGs.size();
    E_Int n_bad_nodes;

    // rule : collapse an edge iff this edge is small compared to some surrounding edges attached at both ends
    for (E_Int i=0; i < npgs; ++i)
    {
      n_bad_nodes = 0;
      const E_Int* pnodes = PGs.get_facets_ptr(i);
      int nnodes = PGs.stride(i);

      for (int n=0; n < nnodes; ++n)
      {
        E_Int Ni   = pnodes[n]-1;
        E_Int Nip1 = pnodes[(n+1) % nnodes]-1;

        double NiNj[3];
        NUGA::diff<3>(crd.col(Ni), crd.col(Nip1), NiNj);
        double Lref = ::sqrt(NUGA::sqrNorm<3>(NiNj));

        const double& emin2 = L(0, Ni);
        const double& emax2 = L(1, Ni);
        if (emin2 < edge_ratio2*emax2) ++n_bad_nodes;

        const double& emin21 = L(0, Nip1);
        const double& emax21 = L(1, Nip1);
        if (emin21 < edge_ratio2*emax21) ++n_bad_nodes;

        double emin = ::sqrt(std::min(emin2, emin21));
        if (Lmax > 0. && emin > Lmax) continue; // consider only edges under Lmax (if valid value)

        bool small_edge = ::fabs(emin - Lref) < 1.e-6 * emin ; // NiNj is (or very near) the smallest incident edge

        if (n_bad_nodes == 2 && small_edge)
        {
          nids[std::min(Ni, Nip1)] = std::max(Ni, Nip1); //X-interface preserving policy : assign max, most chances to be an X point 
          has_changed = true;
        }
      }
    }

    // update the pointers to point to the leaves
    for (size_t i =0; i < nids.size(); ++i)
    {
      E_Int Fi=nids[i];
      while (Fi != nids[Fi])Fi=nids[Fi];
      nids[i]=Fi;
    }

    return has_changed;

  }

  ///
  template <typename CoordAccType>
  bool replace_duplicated_pgs (const CoordAccType& coord, Vector_t<E_Int>& pgnids, bool do_omp)
  { 
    if (PHs.size()==0 || PGs.size() == 0)
      return false;

    pgnids.clear(); 
    pgnids.resize(PGs.size(), IDX_NONE);
    
    PHs.updateFacets();
    PGs.updateFacets();
    
    E_Int nb_pgs, nb_phs(PHs.size()), pgi, count(0);
    E_Float G[3];
    K_FLD::FloatArray barys;
    barys.reserve(3, PGs.size());
        
    Vector_t<E_Int> ISO(PHs.get_facets_max_id(), IDX_NONE), PGI;
    
    // Loop over PH rather than directly over PGs to consider only used PGs
    for (E_Int i = 0; i < nb_phs; ++i)
    {
      nb_pgs = PHs.stride(i);
      for (E_Int j = 0; j < nb_pgs; ++j)
      {
        pgi = PHs.get_facet(i, j)-1; 
        if (ISO[pgi] != IDX_NONE) //already stored
          continue;
        E_Int s = PGs.stride(pgi);
        K_MESH::Polygon PGi(PGs.get_facets_ptr(pgi), s, -1/*to make it zero based*/);
        PGi.iso_barycenter<CoordAccType, 3>(coord, G);
        barys.pushBack(G, G+3);
        
        ISO[pgi]=count++;
      }
    }
    
    K_CONNECT::IdTool::reverse_indirection(ISO, PGI);
    
    // Merge
    Vector_t<E_Int> nids;
    K_FLD::ArrayAccessor<K_FLD::FloatArray> cab(barys);
    E_Int nmerges = ::merge(cab, EPSILON, nids, do_omp);

    if (!nmerges) //no duplicated faces.
    {
      pgnids.clear(); // to get empty pgnids
      return false;
    }
    
    // check if matching isoG means really identical PGs
    bool found = false;
    for (size_t i=0; i < nids.size(); ++i)
    {
      if (nids[i] == (E_Int)i)
        continue;
      
      const E_Int &pgi1 = PGI[i];
      const E_Int &pgi2 = PGI[nids[i]];
      
      if (PGs.stride(pgi1) != PGs.stride(pgi2)) // fast treatment : not the same number of nodes
      {
        nids[i]=i;
        continue;
      }
      if (!K_CONNECT::IdTool::equal(PGs.get_facets_ptr(pgi1), PGs.get_facets_ptr(pgi2), PGs.stride(pgi1), true/*permut*/, false/*strict orient*/))
        nids[i]=i;
      else
        found=true;
    }
    
    if (!found)
    {
      pgnids.clear(); // to get empty pgnids
      return false;
    }
    
    K_CONNECT::IdTool::propagate(nids, ISO);
    
    for (E_Int i = 0; i < nb_phs; ++i)
    {
      nb_pgs = PHs.stride(i);
      for (E_Int j = 0; j < nb_pgs; ++j)
      {
        E_Int & pgi = PHs.get_facet(i, j); 
        E_Int npgid = PGI[ISO[pgi-1]]+1; 
	pgnids[pgi-1]=npgid-1;
        pgi = npgid;
      }
    }
    return true;
  }
  
  /// surfacic version of replace_duplicated_pgs : here PHs are polygones and PGs are edges..
  bool replace_duplicated_edges()
  {
    if (PHs.size()*PGs.size() == 0)
      return false;
    
    PHs.updateFacets();
    PGs.updateFacets();
    
    E_Int nb_edges(PGs.size()); // nb_polys(PHs.size());
    std::map<K_MESH::NO_Edge, E_Int> edge_to_id;
    Vector_t<E_Int> nids(nb_edges, IDX_NONE);
    K_MESH::NO_Edge E;
    bool found = false;
    for (E_Int i = 0; i < nb_edges; ++i)
    {
      const E_Int* ptr = PGs.get_facets_ptr(i);
      E.setNodes(ptr);
      if (edge_to_id.find(E) != edge_to_id.end())
      {
        nids[i] = edge_to_id[E];
        found = true;
      }
      else
      {
        edge_to_id[E]=i;
        nids[i]=i;
      }
    }
    
    if (!found)
      return false;
    
    PHs.change_indices(nids);
    return true;
  }
  
  /// lineic version of replace_duplicated_pgs : here PHs are edges and PGs are nodes..
  bool replace_duplicated_nodes()
  {
    if (PHs.size()*PGs.size() == 0)
      return false;
    
    PHs.updateFacets();
    PGs.updateFacets();
    
    E_Int Ni, nb_nodes(PGs.size());
    Vector_t<E_Int> nids(nb_nodes, IDX_NONE);
    Vector_t<E_Int> pgid(PGs.get_facets_max_id(), IDX_NONE);
    bool found = false;
    for (E_Int i = 0; i < nb_nodes; ++i)
    {
      const E_Int* ptr = PGs.get_facets_ptr(i);
      Ni = *ptr-1;
      if (pgid[Ni] != IDX_NONE)
      {
        nids[i] = pgid[Ni];
        found = true;
      }
      else
        nids[i]=pgid[Ni]=i;
    }
    
    if (!found)
      return false;
    
    PHs.change_indices(nids);
    return true;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// DETECTION /////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  ///
  template <typename CoordAccType>
  E_Int detect_phs_with_same_centroid (const CoordAccType& coord, Vector_t<E_Int>& nids)
  {
    if (PHs.size()*PGs.size() == 0)
      return false;
    
    PHs.updateFacets();
    PGs.updateFacets();
    
    E_Int nb_pgs, nb_phs(PHs.size()), pgi;
    E_Float g[3], G[3], K;
    K_FLD::FloatArray barys;
    barys.reserve(3, nb_phs);
    
    //
    for (E_Int i = 0; i < nb_phs; ++i)
    {
      nb_pgs = PHs.stride(i);
      G[0]=G[1]=G[2]=0.;
      for (E_Int j = 0; j < nb_pgs; ++j)
      {
        pgi = PHs.get_facet(i, j)-1; 
        E_Int s = PGs.stride(pgi);
        K_MESH::Polygon PGi(PGs.get_facets_ptr(pgi), s, -1/*to make it zero based*/);
        PGi.iso_barycenter<CoordAccType, 3>(coord, g);
        
        for (size_t k=0; k < 3; ++k)
          G[k] += g[k];
      }
      K = 1./(E_Float)nb_pgs;
      for (size_t k=0; k < 3; ++k)
        G[k] *= K;
        
      barys.pushBack(G, G+3);
    }
    
    // Merge
    nids.clear();
    K_FLD::ArrayAccessor<K_FLD::FloatArray> cab(barys);
 
    return ::merge(cab, EPSILON, nids, true /*do_omp*/);
  }

  ///
  static E_Int detect_concave_PHs(const K_FLD::FloatArray& crd, const ngon_t& ng, const ngon_unit& orient, E_Float concavity_ratio, E_Float convexity_ratio, std::vector<E_Int>& PHlist)
  {
    bool f;
    E_Int nb_phs = ng.PHs.size();

    PHlist.clear();
    
    for (E_Int PHi = 0; PHi < nb_phs; ++PHi)
    {
      //std::cout << PHi << " is .." << std::endl;
      E_Int err = K_MESH::Polyhedron<UNKNOWN>::is_concave(crd, ng.PGs, ng.PHs.get_facets_ptr(PHi), ng.PHs.stride(PHi), false/*i.e. closed PH*/, orient.get_facets_ptr(PHi), f, concavity_ratio);

      if (err /*|| f*/)
      {
        std::cout << " faulty PH : " << PHi << " . cannot find out its concavity status." << std::endl;
#ifdef DEBUG_NGON_T
        std::ostringstream o;
        if (err)
        o << "faulty_PH_" << PHi << ".plt";
        else
          o << "concave_PH_" << PHi << ".plt";
        NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PH(o.str().c_str(), crd, ng, PHi);
#endif
        continue;
      }
      
      if (f) PHlist.push_back(PHi);
    }

    return 0;
  }

  ///
  static E_Int detect_concave_PGs
  (const K_FLD::FloatArray& crd, const ngon_unit& pgs, E_Float convexity_tol, std::vector<bool>& flag)
  {
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
    acrd_t acrd(crd);

    E_Int nb_pgs = pgs.size();

    flag.clear();
    flag.resize(nb_pgs, false);
    E_Float normal[3];
    for (E_Int PGi = 0; PGi < nb_pgs; ++PGi)
    {
      const E_Int* nodes = pgs.get_facets_ptr(PGi);
      E_Int nb_nodes = pgs.stride(PGi);
      K_MESH::Polygon::normal<acrd_t, 3>(acrd, nodes, nb_nodes, 1, normal);
      flag[PGi] = !K_MESH::Polygon::is_convex(crd, nodes, nb_nodes, 1, normal, convexity_tol);
    }

    return 0;
  }
  
  ///
  template<typename TriangulatorType>
  static bool detect_pathological_PHs(const K_FLD::FloatArray& crd, ngon_t& ng, const ngon_unit& orient, std::vector<ePathoPH>& PHtype)
  {
    TriangulatorType dt;
    
    E_Int nb_phs = ng.PHs.size();
    
    PHtype.clear();
    PHtype.resize(nb_phs, PATHO_PH_NONE);
    
    bool has_problems(false);
    E_Int open(0), showstopper(0), to_split(0);

#ifndef DEBUG_NGON_T
#pragma omp parallel for shared(has_problems) private(dt)
#endif
    for (E_Int PHi=0; PHi < nb_phs; ++PHi)
    {
      const E_Int* first_pg = ng.PHs.get_facets_ptr(PHi);
      E_Int nb_pgs = ng.PHs.stride(PHi);
      
      PHtype[PHi] = (ePathoPH)K_MESH::Polyhedron<UNKNOWN>::is_pathological(dt, crd, ng.PGs, first_pg, nb_pgs, orient.get_facets_ptr(PHi));
            
      if (PHtype[PHi] == dOPEN_PHS) //this should be a bug
      {
        has_problems = true;
        std::cout << "OPEN PHS : " << PHi << std::endl;
        ++open;
        continue;
        
#ifdef DEBUG_NGON_T
        std::cout << PHi << " : is open !." << std::endl;
#endif
      }
      
      if (PHtype[PHi] == PATHO_PH_NONE ) // OK
        continue;
      else if (PHtype[PHi] == CONCAVITY_TO_SPLIT) // FLAG IT FOR SPLIT
      {
        has_problems = true;
        ++to_split;
      }
      else if (PHtype[PHi] == PH_DELAUNAY_FAILURE) /*if (nb_pgs > 6)*/ /* fixme : temporary hack to avoid detecting highly anisotropic prisms*/
      {
        has_problems = true;
        std::cout << "DELAUNAY_FAILURE : " << PHi << std::endl;
        ++showstopper;
      }
      else if (PHtype[PHi] == PH_DEGEN)
      {
        has_problems = true;
        std::cout << "SPIKE OR HAT OR DEGEN PG : " << PHi << std::endl;
        ++showstopper;
      }
      else
      {
        //todo
      }
      
      //if (concave) std::cout << "concave " << std::endl;
 
#ifdef DEBUG_NGON_T
      if (PHtype[PHi] == PATHO_PH_NONE) continue;

//      ngon_unit phs;
//      phs.add(ng.PHs.stride(PHi), ng.PHs.get_facets_ptr(PHi));
//      ngon_t ngo;
//      ngo.PGs=ng.PGs;
//      ngo.PHs= phs;
//      std::vector<E_Int> pgnids, phnids;
//      ngo.remove_unreferenced_pgs(pgnids, phnids);
//      K_FLD::IntArray cnto;
//      ngo.export_to_array(cnto);
//      medith::write("patho.plt", crd, cnto, "NGON");
//    
//      DELAUNAY::Triangulator t;
//      E_Float v, centroid[3];
//      K_MESH::Polyhedron<UNKNOWN>::metrics2<DELAUNAY::Triangulator>(t, crd, ng.PGs, ng.PHs.get_facets_ptr(PHi), ng.PHs.stride(PHi), v, centroid);
//          
//      E_Float isoG[3];
//      typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
//      acrd_t acrd(crd);
//      std::vector<E_Int> phnodes;
//      K_MESH::Polyhedron<UNKNOWN>::unique_nodes(ng.PGs, ng.PHs.get_facets_ptr(PHi), ng.PHs.stride(PHi), phnodes);
//      K_MESH::Polyhedron<UNKNOWN>::iso_barycenter(acrd, &phnodes[0], phnodes.size(), 1, isoG);
//        
//      K_FLD::IntArray cn(2,1,0);
//      cn(1,0)=1;
//      K_FLD::FloatArray crd;
//      crd.pushBack(isoG, isoG+3);
//      crd.pushBack(centroid, centroid+3);
//
//      medith::write("baryTocentroid.mesh", crd, cn, "BAR");    
#endif
    }
    
    if (has_problems)
    {
      std::cout << "Open cells : " << open << std::endl;
      std::cout << "Showstoppers : " << showstopper << std::endl;
      std::cout << "NON-star to split : " << to_split << std::endl;
    }
    
    return has_problems;
  }

  ///
  template <typename TriangulatorType>
  static E_Int detect_uncomputable_pgs(const K_FLD::FloatArray& crd, const ngon_unit& pgs, std::vector<ePathoPG> &flagPG)
  {
    E_Int errcount(0), normal_err_hat_count(0), normal_err_spike_count(0), normal_err_unknown_count(0), delaunay_err_count(0);
    flagPG.clear();
    flagPG.resize(pgs.size(), PATHO_PG_NONE);
  
    E_Int err = 0;
    TriangulatorType dt;
    K_FLD::IntArray cT3;
    
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray > acrd_t;
    acrd_t acrd(crd);
    E_Float W[3];
    
    
#ifndef DEBUG_NGON_T
#pragma omp parallel for private(err, dt, cT3, W) reduction(+:normal_err_hat_count, normal_err_spike_count, normal_err_unknown_count, delaunay_err_count)
#endif
    for (size_t i = 0; i < flagPG.size(); ++i)
    {
      const E_Int * nodes = pgs.get_facets_ptr(i);
      E_Int nb_nodes = pgs.stride(i);
      
      // test normal computation
      K_MESH::Polygon::normal<acrd_t, 3>(acrd, nodes, nb_nodes, 1, W);
      E_Float l2 = ::sqrt(W[0]*W[0]+W[1]*W[1]+W[2]*W[2]);
      
      if (::fabs(l2 - 1.) >= EPSILON) // NORMAL CALCULATION FAILED
      {
        /*E_Float Lmin=NUGA::FLOAT_MAX;
        E_Float Lmax = -1;
        for (E_Int n=0; n < nb_nodes; ++n)
        {
          E_Int ni = *(nodes+n)-1;
          E_Int nj = *(nodes+(n+1)%nb_nodes)-1;
          E_Float d2 = NUGA::sqrDistance(crd.col(ni), crd.col(nj), 3);
          Lmin = std::min(Lmin, d2);
          Lmax = std::max(Lmax, d2);
        }
        
        Lmin = ::sqrt(Lmin);
        Lmax = ::sqrt(Lmax);*/

        if (nb_nodes == 3)
        {
          E_Int nspecial;
          K_MESH::Triangle::eDegenType type = K_MESH::Triangle::degen_type(crd, nodes[0]-1, nodes[1]-1, nodes[2]-1, 1.e-12, 0., nspecial);
          
          if (type == K_MESH::Triangle::HAT)
          {
            flagPG[i]=HAT;
            ++normal_err_hat_count;

#ifdef DEBUG_NGON_T
            std::cout << "NORMAL FAILURE : HAT: " << i << std::endl;
        /*ngon_unit pg;
        pg.add(pgs.stride(i), pgs.get_facets_ptr(i));
        ngon_t ngo(pg);
        K_FLD::IntArray cnto;
        ngo.export_to_array(cnto);
        medith::write("hat.plt", crd, cnto, "NGON");*/
#endif
          }
          else if (type == K_MESH::Triangle::SPIKE)
          {
            flagPG[i]=SPIKE;
            ++normal_err_spike_count;

#ifdef DEBUG_NGON_T
            std::cout << "NORMAL FAILURE : SPIKE: " << i << std::endl;
        /*ngon_unit pg;
        pg.add(pgs.stride(i), pgs.get_facets_ptr(i));
        ngon_t ngo(pg);
        K_FLD::IntArray cnto;
        ngo.export_to_array(cnto);
        medith::write("spike.plt", crd, cnto, "NGON");*/
#endif
          }
          else if (type == K_MESH::Triangle::SMALL)
          {
            flagPG[i]=SPIKE;//fixme
             ++normal_err_spike_count;
          }
        }
        else
        {
          //todo
          ++normal_err_unknown_count;
          std::cout << "NORMAL FAILURE : UNKNOWN: " << i << std::endl;
          flagPG[i] = DEGEN_PG;
        }
      }
      
      //test triangulation
      if (flagPG[i] == PATHO_PG_NONE) 
      {
        err = K_MESH::Polygon::triangulate(dt, crd, nodes, nb_nodes, 1, cT3);
        flagPG[i] = (err != 0) ? PG_DELAUNAY_FAILURE : PATHO_PG_NONE;
        if (err) ++delaunay_err_count;
      }
      
      if (flagPG[i] != PATHO_PG_NONE)
      {
#ifdef DEBUG_NGON_T
        ngon_unit pg;
        pg.add(pgs.stride(i), pgs.get_facets_ptr(i));
        ngon_t ngo(pg);
        K_FLD::IntArray cnto;
        ngo.export_to_array(cnto);
        medith::write("fPG.plt", crd, cnto, "NGON");
#endif
      }
    }
    
    errcount = normal_err_hat_count + normal_err_spike_count + normal_err_unknown_count + delaunay_err_count;
    if (errcount){
      std::cout << errcount << " uncomputable surfaces were detected : " << std::endl;
      std::cout << "  -- Normal failure    : " << normal_err_hat_count + normal_err_spike_count + normal_err_unknown_count << std::endl;
      std::cout << "  ----- hat       : " << normal_err_hat_count << std::endl;
      std::cout << "  ----- spike     : " << normal_err_spike_count << std::endl;
      std::cout << "  ----- degen PG  : " << normal_err_unknown_count << std::endl;
      std::cout << "  -- Delaunay failure  : " << delaunay_err_count << std::endl;
    }
    else
      std::cout << "OK : all polygons are computables." << std::endl;

  return errcount;
}
  
///
template <typename TriangulatorType>
static void detect_uncomputable_phs(const K_FLD::FloatArray& crd, const ngon_t& ngi, Vector_t<E_Int>& PHlist)
{
  ngi.PGs.updateFacets();
  ngi.PHs.updateFacets();
  
  PHlist.clear();
  Vector_t<ePathoPG> flagPGs;
  detect_uncomputable_pgs<TriangulatorType>(crd, ngi.PGs, flagPGs);
  
  E_Int nb_phs = ngi.PHs.size();
  for (E_Int i=0; i < nb_phs; ++i)
  {
    const E_Int* pgs = ngi.PHs.get_facets_ptr(i);
    E_Int nb_pgs = ngi.PHs.stride(i);
    
    bool caught=false;
    for (E_Int n=0; (n<nb_pgs) && !caught; ++n)
      caught = (flagPGs[*(pgs+n)-1] == true);
    
    if (caught) PHlist.push_back(i);
  }
}

///
template <typename TriangulatorType>
static E_Int detect_bad_volumes(const K_FLD::FloatArray& crd, const ngon_t& ngi, const ngon_unit& neighborsi, E_Float vmin, E_Float vratio, Vector_t<E_Int>& PHlist)
{
  ngi.PGs.updateFacets();
  ngi.PHs.updateFacets();
  
  PHlist.clear();

  //logic with varatio i [0,1]
  if (vratio > 1.) vratio = 1./vratio;
  
  E_Int nb_phs = ngi.PHs.size();
 
  std::vector<E_Float> vols;
  ngon_t::volumes<TriangulatorType>(crd, ngi, vols, false/*i.e. not all cvx*/, true/* ! new algo*/);
  
  using pair_t = std::pair<E_Float, E_Int>;
  std::vector<pair_t> v_to_id;

  E_Float vi, vj;
  
  for (E_Int i=0; i < nb_phs; ++i)
  {
    vi = vols[i];
    if (vi < vmin)
    {
      //std::cout << " cell " << i << " has volume " << vi << " which is less than " << vmin << std::endl;
      v_to_id.push_back(std::make_pair(vi, i));
      continue;
    }
      
    if (vratio > 0.)
    {

      E_Int nb_neighs = neighborsi.stride(i);
      bool found = false;
      for (E_Int n = 0; n < nb_neighs && !found; ++n)
      {
        E_Int j = neighborsi.get_facet(i, n);
        if (j == IDX_NONE)
          continue;
        vj = vols[j];

        if (vi< vratio*vj) //OK because vration in [0.,1.]
        {
          found = true;
          //std::cout << " cell " << i << " has volume " << vi << " which is " << vratio << " times smaller than cell " << j << " with vol " << vj << std::endl;
        }
      }
      if (found)
      {
        v_to_id.push_back(std::make_pair(vi, i));
        continue;
      }
    }
  }

  std::stable_sort(v_to_id.begin(), v_to_id.end(), 
  [](const pair_t& a, const pair_t& b)
  {
    double tol = (std::max(a.first, b.first) < 1.e-15) ? 1.e-30 : 1.e-15;
    return (a.first < (b.first + tol));
  }
  );

  for (size_t i = 0; i < v_to_id.size(); ++i)
  {
    //std::cout << "small cell detected : " << v_to_id[i].second << " with volume : " << v_to_id[i].first << std::endl;
    PHlist.push_back(v_to_id[i].second);
  }

  
  return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// EXTRACTION ///////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
///
template <typename TriangulatorType>
static E_Int extract_uncomputables
(const K_FLD::FloatArray& crd, const ngon_t& ngi, E_Int neigh_level,
 ngon_unit& upgs, ngon_t& uphs, ngon_t& uphs_wv1, ngon_t& remaining)
{
  //detect cells with uncomputable polygons
  Vector_t<E_Int> PHlist;
  ngi.PGs.updateFacets();
  ngi.PHs.updateFacets();
  
  upgs.clear();
  uphs.clear();
  uphs_wv1.clear();
  remaining.clear();

  Vector_t<ePathoPG> flagPGs;
  E_Int nb_bad_pgs = detect_uncomputable_pgs<TriangulatorType>(crd, ngi.PGs, flagPGs);

  if (nb_bad_pgs == 0)
    return 0;
  else
  {
    Vector_t<E_Int> indices;
    for (size_t i=0; i < flagPGs.size(); ++i) if (flagPGs[i] != PATHO_PG_NONE) indices.push_back(i);
    Vector_t<E_Int> oids;
    ngi.PGs.extract(indices, upgs, oids);
  }
  
  // get uncomputable phs
  E_Int nb_phs = ngi.PHs.size();
  for (E_Int i=0; i < nb_phs; ++i)
  {
    const E_Int* pgs = ngi.PHs.get_facets_ptr(i);
    E_Int nb_pgs = ngi.PHs.stride(i);
    
    bool caught=false;
    for (E_Int n=0; (n<nb_pgs) && !caught; ++n)
      caught = (flagPGs[*(pgs+n)-1] != PATHO_PG_NONE);
    
    if (caught) PHlist.push_back(i);
  }
     
  std::vector<bool> keep(ngi.PHs.size(), false);
  for (size_t i=0; i < PHlist.size(); ++i) keep[PHlist[i]] = true;
  
  {
    Vector_t<E_Int> ngpids;
    ngi.select_phs(ngi, keep, ngpids, uphs);
  }

  // extend with second neighborhood and separate from non-involved polyhedra
  for (E_Int j=0; j< neigh_level; ++j)
    flag_neighbors(ngi, keep);
    
  split_phs(ngi, keep, uphs_wv1, remaining);
  
  return 0;
}

///
template <typename TriangulatorType>
static E_Int extract_pathological_PHs
(const K_FLD::FloatArray& crd, ngon_t& ngi, E_Int neigh_level, std::vector<ngon_t>& phsv, ngon_t& neigh_phs, ngon_t& remaining_phs)
{
  //detect non-centroid-star-shaped cells
  
  ngi.PGs.updateFacets();
  ngi.PHs.updateFacets();
  
  phsv.clear();
  neigh_phs.clear();
  remaining_phs.clear();

  //std::cout << "extract patho : build orient..." << std::endl;
  ngon_unit orienti;
  E_Int err = build_orientation_ngu<TriangulatorType>(crd, ngi, orienti);
  if (err)
    return err;

  Vector_t<ePathoPH> PHtypes;
  //std::cout << "extract patho : detect_pathological_PHs..." << std::endl;
  /*bool has_non_star = */detect_pathological_PHs<TriangulatorType>(crd, ngi, orienti, PHtypes);
  
  // OVERRIDE WITH UNCOMPUTABLES PGS (and therefore PHs) missed in detect_pathological_PHs (because some test have suceed before getting to the delaunay)
  //std::cout << "extract patho : detect_uncomputable_pgs..." << std::endl;
  E_Int nb_bad_pgs = 0;
  {
    Vector_t<ePathoPG> flagPGs;
    nb_bad_pgs = detect_uncomputable_pgs<TriangulatorType>(crd, ngi.PGs, flagPGs);
    // get uncomputable phs
    //std::cout << "nb uncomp pgs : " << nb_bad_pgs << std::endl;
    if (nb_bad_pgs)
    {
      E_Int nb_phs = ngi.PHs.size();
      for (E_Int i=0; i < nb_phs; ++i)
      {
       const E_Int* pgs = ngi.PHs.get_facets_ptr(i);
       E_Int nb_pgs = ngi.PHs.stride(i);
    
       for (E_Int n=0; (n<nb_pgs); ++n){
         E_Int flag = flagPGs[*(pgs+n)-1];
         
         if (flag == PATHO_PG_NONE) continue;
         else if (flag == PG_DELAUNAY_FAILURE) { 
           if (PHtypes[i] != PATHO_PH_NONE) std::cout << "PH : " << i << " pathology switched from " << PHtypes[i] << " to " << PH_DELAUNAY_FAILURE << std::endl;
           PHtypes[i] = PH_DELAUNAY_FAILURE; break;
         }
         else /*SPIKE HAT DEGEN_PG*/ { 
           if (PHtypes[i] != PATHO_PH_NONE) std::cout << "PH : " << i << " pathology switched from " << PHtypes[i] << " to " << PH_DEGEN << " (PG patho : " << flag << ")" << std::endl;
           PHtypes[i] = PH_DEGEN; break;
         }
       }
      }
    }
  }
  
  phsv.resize(4); //currently four type
  
  E_Int nb_phs = ngi.PHs.size();
  
  
  //std::cout << "detect_pathological_PHs : " << has_non_star << std::endl;
 
  //OPEN CELL
  {
    std::vector<bool> keep(nb_phs, false);
    for (E_Int i=0; i < nb_phs; ++i) keep[i] = (PHtypes[i] == OPEN_PHS);
    Vector_t<E_Int> ngpids;
    select_phs(ngi, keep, ngpids, phsv[0]);
  }
    
  //DEGEN CELLS ( containing degen PGs : might be fixable by PG collapse)
  {
    std::vector<bool> keep(nb_phs, false);
    for (E_Int i=0; i < nb_phs; ++i) keep[i] = (PHtypes[i] == PH_DEGEN);
    Vector_t<E_Int> ngpids;
    select_phs(ngi, keep, ngpids, phsv[1]);
  }
  
  //DELAUNAY_FAILURE ( containing degen PGs : might be fixable by PG collapse)
  {
    std::vector<bool> keep(nb_phs, false);
    for (E_Int i=0; i < nb_phs; ++i) keep[i] = (PHtypes[i] == PH_DELAUNAY_FAILURE);
    Vector_t<E_Int> ngpids;
    select_phs(ngi, keep, ngpids, phsv[2]);
  }
  
  //CONCAVITY_TO_SPLIT (might be fixable by split)
  {
    std::vector<bool> keep(nb_phs, false);
    for (E_Int i=0; i < nb_phs; ++i) keep[i] = (PHtypes[i] == CONCAVITY_TO_SPLIT);
    Vector_t<E_Int> ngpids;
    select_phs(ngi, keep, ngpids, phsv[3]);
  }
  
  // GOOD CELLS
//  Vector_t<bool> keep(nb_phs, false);
//  for (E_Int i=0; i < nb_phs; ++i) keep[i] = (PHtypes[i] == PATHO_PH_NONE);
//  Vector_t<E_Int> ngpids;
//  select_phs(ngi, keep, ngpids, remaining_phs);
//  
//  if (neigh_level)
//  {
//    Vector_t<bool> keep(nb_phs, false);
//    for (E_Int i=0; i < nb_phs; ++i) keep[i] = (PHtypes[i] != PATHO_PH_NONE);
//    for (E_Int j=1; j< neigh_level; ++j)
//      flag_neighbors(ngi, keep, false/*append*/);
//    Vector_t<E_Int> ngpids;
//    select_phs(ngi, keep, ngpids, neigh_phs);
//  }

  return 0;
}

///
static E_Int extract_n_outer_layers (const K_FLD::FloatArray& crd, ngon_t& ngi, E_Int N, ngon_t& outers, ngon_t& remaining, bool discard_external)
{
  //
  E_Int nb_phs = ngi.PHs.size();
  std::vector<bool> keep(nb_phs, false);
  // detect external PGs.
  ngi.flag_external_pgs(INITIAL_SKIN);
  if (discard_external) discard_by_box(crd, ngi, false/*i.e we want to discard external*/);
  //
  std::set<E_Int> attaching_nodes;
  {
    std::vector<E_Int> anodes;
    // initialize attaching nodes with external PGs nodes
    E_Int nb_pgs = ngi.PGs.size();
    for (E_Int i = 0; i < nb_pgs; ++i)
    {
      if (ngi.PGs._type[i] != 1) continue;
      const E_Int* nodes = ngi.PGs.get_facets_ptr(i);
      E_Int nb_nodes = ngi.PGs.stride(i);
      anodes.insert(anodes.end(), nodes, nodes + nb_nodes);
    }
    //K_CONNECT::IdTool::shift(anodes, -1);
    attaching_nodes.insert(ALL(anodes));
  }

  E_Int n = 1;
  std::vector<bool> kp;
  std::vector<E_Int> vtmp1, vtmp1bis, vtmp2, vtmp3;
  E_Int noutercell(0);

  //
  do
  {
  	//std::cout << "nb of attaching_nodes : " << attaching_nodes.size() << std::endl;
    // flag attached elements to current nodes set
    E_Int nb_elts = get_elements_having_nodes(ngi, attaching_nodes, keep/*to skip already in*/, kp);
    if (nb_elts == 0) break;
    //std::cout << "nb of attached element : "  << nb_elts << std::endl;

    // append flags
    for (E_Int i = 0; i < nb_phs; ++i) keep[i] = keep[i] || kp[i];

    // get new nodes set
    vtmp1.clear();
    
    for (E_Int i = 0; i < nb_phs; ++i)
    {
      if (!kp[i]) continue;
      ++noutercell;

      ngi.nodes_ph(i, vtmp1bis, false);
      //E_Int sz0 = vtmp1.size();
      vtmp1.insert(vtmp1.end(), ALL(vtmp1bis));
      //std::cout << "added nodes : " << toto.size()-sz0 << std::endl;
    }
    vtmp2.clear();
    vtmp3.clear();

    vtmp2.insert(vtmp2.end(), ALL(attaching_nodes));
    std::sort(ALL(vtmp1));

    std::set_difference(ALL(vtmp1), ALL(vtmp2), std::back_inserter(vtmp3)); // new set = all PH attached nodes minus current set
    //std::cout << "number of nodes in the new front : " << vtmp3.size() << std::endl;

    if (vtmp3.empty()) break;
    
    attaching_nodes.clear();
    attaching_nodes.insert(ALL(vtmp3));

  } while (n++ < N);

  if (noutercell == nb_phs)
  	std::cout << "N is too big : everything is in the outer zone." << std::endl;

  //
  split_phs(ngi, keep, outers, remaining);
  return 0;
}

///
static E_Int get_elements_having_nodes(const ngon_t& ngi, const std::set<E_Int>& tnodes/*starting at 1*/, const std::vector<bool>& skip, std::vector<bool>& flag)
{
  E_Int nb_phs = ngi.PHs.size();
  //
  flag.clear();
  flag.resize(nb_phs, false);

  E_Int count(0);
  
  for (E_Int i = 0; i < nb_phs; ++i)
  {
    if (skip[i]) continue;

    const E_Int* pgs = ngi.PHs.get_facets_ptr(i);
    E_Int nb_pgs = ngi.PHs.stride(i);

    bool is_attached = false;

    for (E_Int p = 0; p < nb_pgs; ++p)
    {
      E_Int PGi = *(pgs + p) - 1;

      const E_Int* nodes = ngi.PGs.get_facets_ptr(PGi);
      E_Int stride = ngi.PGs.stride(PGi);
      //
      for (E_Int j = 0; j < stride; ++j)
      {
        E_Int Ni = *(nodes + j);
        if (tnodes.find(Ni) != tnodes.end())
        {
          is_attached = true;
          break;
        }
      }

      if (is_attached) break;
    }

    flag[i] = is_attached;
    if (is_attached)++count;
  }
  return count;
}

/// 
static E_Int discard_by_box(const K_FLD::FloatArray& coord, ngon_t& ngi, bool holes)
{
  assert (ngi.PHs.size() < ngi.PGs.size()); //i.e. not one ph per pg beacuse we need at least a closed PH

#ifdef FLAG_STEP
  NUGA::chrono c;
  c.start();
#endif
  
  // 1. Get the skin PGs
  Vector_t<E_Int> oids;
  ngon_unit pg_ext;
  ngi.PGs.extract_of_type (INITIAL_SKIN, pg_ext, oids);
  if (pg_ext.size() == 0)
    return 0; //should be only the case where One of the mesh contains completely the other

#ifdef FLAG_STEP
  std::cout << "__discard_holes : extract_of_type : " << c.elapsed() << std::endl;
  c.start();
#endif
    
  // 2. Build the neighbourhood for skin PGs
  ngon_unit neighbors;
  K_MESH::Polygon::build_pg_neighborhood(pg_ext, neighbors);

#ifdef FLAG_STEP
  std::cout << "__discard_holes : build_pg_neighborhood : " << c.elapsed() << std::endl;
  c.start();
#endif
  
  // Color to get connex parts
  E_Int nb_connex = 1;
  Vector_t<E_Int> colors;
  if (ngi.PHs.size() > 1)
  {
    NUGA::EltAlgo<K_MESH::Polygon>::coloring (neighbors, colors);
    nb_connex = 1+*std::max_element(colors.begin(), colors.end());
  }

#ifdef FLAG_STEP
  std::cout << "__discard_holes : coloring : " << c.elapsed() << std::endl;
  c.start();
#endif
  
  if (nb_connex == 1) //no holes
    return 0;
  
  // Find out which is the external (the one with the biggest bbox)
  Vector_t<K_SEARCH::BBox3D> bbox_per_color(nb_connex);
  //init boxes
  for (E_Int i=0; i < nb_connex; ++i)
  {
    bbox_per_color[i].maxB[0]=bbox_per_color[i].maxB[1]=bbox_per_color[i].maxB[2]=-NUGA::FLOAT_MAX;
    bbox_per_color[i].minB[0]=bbox_per_color[i].minB[1]=bbox_per_color[i].minB[2]=NUGA::FLOAT_MAX;
  }
  
  Vector_t<E_Int> nodes;
  
  E_Int nb_pgex=colors.size();
  K_SEARCH::BBox3D box;
  for (E_Int i=0; i < nb_pgex; ++i)
  {
    const E_Int& s = pg_ext.stride(i);
    const E_Int* pN = pg_ext.get_facets_ptr(i);
    
    nodes.clear();
    for (E_Int j = 0; j < s; ++j, ++pN)
      nodes.push_back((*pN)-1);//indices convention : start at 1 

    box.compute(coord, nodes);
    
    K_SEARCH::BBox3D& b = bbox_per_color[colors[i]];
    
    for (size_t j = 0; j < 3; ++j) 
    {
      b.minB[j] = (b.minB[j] > box.minB[j]) ? box.minB[j] : b.minB[j];
      b.maxB[j] = (b.maxB[j] < box.maxB[j]) ? box.maxB[j] : b.maxB[j];
    }
  }
  // And the hole color are..
  std::set<E_Int> hole_colors;
  for (E_Int i = 0; i < nb_connex; ++i)
  {
    for (E_Int j = i + 1; j < nb_connex; ++j)
    {
      E_Int I=i;
      E_Int J=j;
      if (!holes) std::swap(I,J); //i.e. discard external
      
      if (K_SEARCH::BbTree3D::box1IsIncludedinbox2(&bbox_per_color[i], &bbox_per_color[j], EPSILON))
        hole_colors.insert(I);
      else if (K_SEARCH::BbTree3D::box1IsIncludedinbox2(&bbox_per_color[j], &bbox_per_color[i], EPSILON))
        hole_colors.insert(J);
    }
  }

  // Reset flag for holes
  for (E_Int i=0; i < nb_pgex; ++i)
  {
    if (hole_colors.find(colors[i]) != hole_colors.end())
      ngi.PGs._type[oids[i]]=INNER;
  }
  ngi.flag_external_phs(INITIAL_SKIN);//update consistently PHs flags
  
  return 0;
} 

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////:: FOR VISU /////////////////////////////////////////////////////////////////////

///
template <typename TriangulatorType>
static E_Int stats_bad_volumes
(const K_FLD::FloatArray& crd, const ngon_t& ngi, const ngon_unit& neighborsi, 
 const std::vector<E_Float>&vols, E_Float vmin, Vector_t<E_Float>& aspect_ratio)
{
  ngi.PGs.updateFacets();
  ngi.PHs.updateFacets();
  
  E_Int nb_phs = ngi.PHs.size();
  
  aspect_ratio.clear();
  aspect_ratio.resize(nb_phs, 0.);
 
  
  E_Float vi, vj;
  
  for (E_Int i=0; i < nb_phs; ++i)
  {
    vi = vols[i];
    //if (vi < EPSILON) continue;
    if (vi < vmin) continue;  
    //std::cout << "vi : " << vi << std::endl;

    E_Int nb_neighs = neighborsi.stride(i);
    
    E_Float ar = NUGA::FLOAT_MAX;
    
    for (E_Int n=0; n < nb_neighs; ++n)
    {
      E_Int j = neighborsi.get_facet(i,n);
      if (j == IDX_NONE)
        continue;
      
      vj = vols[j];
      //std::cout << "vj : " << vj << std::endl;
      if (vj < vmin)
      {
        ar = 0.;
        break;
      }
      
      if (vj < vi)
        ar = std::min(vj / vi, ar);
      else
        ar = std::min(vi / vj, ar);   
    }
    
    aspect_ratio[i] = ar;
  }
  
  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
///
E_Int remove_degenerated_pgs(E_Int ngon_dim, Vector_t<E_Int>& pgnids, Vector_t<E_Int>& phnids)
{
  Vector_t<E_Int> pgindices;
  
  PGs.get_degenerated(ngon_dim, pgindices); // Degenrated edge(<2) or Polygons (<3) => topological test

  remove_pgs(pgindices, pgnids, phnids); //sync the PHs accordingly
  PGs.updateFacets();//fixme : required ?
  return pgindices.size();
}
  
///
E_Int remove_unreferenced_pgs(Vector_t<E_Int>& pgnids, Vector_t<E_Int>& phnids)
{
  PHs.updateFacets();
  PGs.updateFacets();

#ifdef DEBUG_NGON_T
  assert(attributes_are_consistent());
#endif
  
  Vector_t<E_Int> visited(PGs.size(),0);
  for (E_Int i=0; i < PHs.size(); ++i)
    for (E_Int j=0; j < PHs.stride(i); ++j)
      ++visited[PHs.get_facet(i,j)-1];
 
  Vector_t<E_Int> unused;
  for (size_t i=0; i < visited.size(); ++i)if (!visited[i])unused.push_back(i);
  remove_pgs(unused, pgnids, phnids);
  return (E_Int)unused.size();
}
  
  // "PRIVATE" METHODS
  ///
  static void __get_selected_ids(const ngon_unit& ngon_in, const Vector_t<bool>& keep,
                                 Vector_t<E_Int>& new_ids, Vector_t<E_Int>& stack_ids)
 {
    //WARNING : change the numerotation
    
    size_t sz0 = ngon_in.size();
    assert(sz0 == keep.size());
  
    ngon_in.updateFacets();
  
    size_t old_facets_sz = ngon_in.get_facets_max_id();
    
    new_ids.resize(old_facets_sz+1, IDX_NONE);
    stack_ids.clear();
  
    size_t count(0);
    E_Int pos, sz, Fij;
  
    for (size_t i = 0; i < sz0; ++i)
    {
      if (!keep[i])
        continue;
    
      pos = ngon_in._facet[i];
      sz = ngon_in._NGON[pos];
    
      for (E_Int j = 0; j < sz; ++j)
      {
        Fij = ngon_in.get_facet(i,j);
        if (new_ids[Fij] == IDX_NONE)
        {
          new_ids[Fij]=++count;
          stack_ids.push_back(Fij);
        }
      }
    }
  }

  ///
  static void __create_selected_phs
  (const ngon_unit& ngon_in, const Vector_t<bool>& keep,
   const Vector_t<E_Int>& nids, ngon_unit& ngon_out)
  {
    size_t ksz(keep.size());
    E_Int pos, sz, Fij;
  
    ngon_out.clear();
    ngon_out._NGON.resize(2, 0);
  
    ngon_in.updateFacets();//refresh in case
  
    for (size_t i = 0; i < ksz; ++i)
    {
      if (!keep[i]) continue;
    
      ++ngon_out._NGON[0];
    
      pos = ngon_in._facet[i];
      sz = ngon_in._NGON[pos];
      ngon_out._NGON.push_back(sz);

      for (E_Int j = 0; j < sz; ++j)
      {
        Fij = ngon_in.get_facet(i,j);
        ngon_out._NGON.push_back(nids[Fij]);
      } 
    }
  
    ngon_out._NGON[1]=ngon_out._NGON.size()-2;
    
    //transfer externality
    ngon_out.compact_attributes(ngon_in, keep);
    
  }

  ///
  static void __create_selected_pgs
  (const ngon_unit& ngon_in, const Vector_t<E_Int>& stack_ids, ngon_unit& ngon_out)
  {
    size_t ksz(stack_ids.size());
    E_Int pos,sz, Fij;
  
    ksz = stack_ids.size();
  
    ngon_out.clear(); 
    ngon_out._NGON.resize(2, 0);
    ngon_out._NGON[0]=ksz;
  
    ngon_in.updateFacets();//refresh in case
  
    bool tranfer_type = !ngon_in._type.empty();
    bool transfer_anc = ngon_in._ancEs.cols();
  
    for (size_t i = 0; i < ksz; ++i)
    {
      const E_Int & Fi = stack_ids[i];
    
      pos = ngon_in._facet[Fi-1];
      sz = ngon_in._NGON[pos];
      ngon_out._NGON.push_back(sz);
      //transfer externality
      if (tranfer_type)
        ngon_out._type.push_back(ngon_in._type[Fi-1]);
      if (transfer_anc)
        ngon_out._ancEs.pushBack(ngon_in._ancEs.col(Fi - 1), ngon_in._ancEs.col(Fi - 1) + 2);
    
      for (E_Int j = 0; j < sz; ++j)
      {
        Fij = ngon_in.get_facet(Fi-1,j);
        ngon_out._NGON.push_back(Fij);
      }
    }
  
    ngon_out._NGON[1]=ngon_out._NGON.size()-2;
  }
  
  ///
  E_Int keep_PHs_having_PGs(const Vector_t<E_Int>& PGs_tokeep)
  {
    // index start dependant
    Vector_t<bool> keep(PGs.size(), false);
    for (size_t i = 0; i < PGs_tokeep.size(); ++i) keep[PGs_tokeep[i]]=true;

    keep_PHs_having_PGs(keep);

    return 0;
  }
  
  ///
  E_Int keep_PHs_having_PGs(const Vector_t<bool>& keep)
  {
    // index start dependant
    ngon_unit keptPHs;
    PHs.updateFacets();
    bool kp;
    E_Int sz, *pN;
    for (E_Int i = 0; i < PHs.size(); ++i)
    {
      sz = PHs.stride(i);
      pN = PHs.get_facets_ptr(i);
      kp = false;
      for (E_Int j = 0; (j < sz) && !kp; ++j)
        kp = keep[*(pN + j) - 1];

      if (kp)
        keptPHs.__add(PHs, i);
    }
  
    PHs=keptPHs;
    PHs.updateFacets();

    return 0;
  }
  
  ///
  E_Int flag_PHs_having_PGs(const Vector_t<bool>& keepPG, Vector_t<bool>& keepPH) const 
  {
    // index start dependant
    
    PHs.updateFacets();
    
    E_Int sz(PHs.size()), nb_pgs;
    
    keepPH.clear();
    keepPH.resize(sz, false);
    
    for (E_Int i = 0; i < sz; ++i)
    {
      nb_pgs = PHs.stride(i);
      const E_Int* pN = PHs.get_facets_ptr(i);
 
      for (E_Int j = 0; (j < nb_pgs); ++j)
      {
        if (keepPG[*(pN + j) - 1])
        {
          keepPH[i]=true;
          break;
        }
      }
    }
  
    return 0;
  }
  
  ///
  void get_PHs_having_PGs(const std::set<E_Int>& pgids, E_Int idx_start, std::set<E_Int>& PHlist) const 
  {
    E_Int sz(PHs.size());
    
    PHlist.clear();
        
    for (E_Int i = 0; i < sz; ++i)
    {
      E_Int nb_pgs = PHs.stride(i);
      const E_Int* faces = PHs.get_facets_ptr(i);
 
      for (E_Int j = 0; (j < nb_pgs); ++j)
      {
        E_Int PGi = *(faces + j);
        PGi = (idx_start == 0) ? PGi - 1 : PGi;
        if (pgids.find(PGi) != pgids.end()) PHlist.insert(i);
      }
    }
  }
  
  ///
  E_Int remove_baffles()
  {
    PHs.updateFacets();
    
    //std::cout << "baffle 1" << std::endl;
    PHs.updateFacets();
    //std::cout << "baffle 2" << std::endl;
    std::map<K_MESH::NO_Edge, E_Int> w_map;
    E_Int max_stride = PHs.get_facets_max_id();
    //std::cout << "baffle 3" << std::endl;
    //std::cout << "max_stride : " << max_stride << std::endl;
    std::vector<E_Int> buffer_flag_pgs(max_stride), molecule(max_stride);
    //std::cout << "asking for size" << std::endl;
    //std::cout << "size : " << PHs.size() << std::endl;
    E_Int nb_phs(PHs.size());
    //std::cout << "baffle 4" << std::endl;
    ngon_unit new_phs;
    E_Int nb_cells_w_baffles{0};

    for (E_Int i=0; i < nb_phs; ++i)
    {
      const E_Int* first_pg = PHs.get_facets_ptr(i);
      E_Int nb_pgs = PHs.stride(i);
      
      if (max_stride < nb_pgs) std::cout << "max stride error for PH : " << i << " . nb_pgs/max_stride : " << nb_pgs << "/" << max_stride << std::endl;
      //std::cout << i << std::endl;
      bool has_baff = K_MESH::Polyhedron<UNKNOWN>::has_baffles(PGs, first_pg, nb_pgs, w_map, &buffer_flag_pgs[0]);

      //std::cout << "has baf : " << has_baff << std::endl;
      
      if (!has_baff)
        new_phs.add(nb_pgs, first_pg);
      else
      {
        ++nb_cells_w_baffles;
        E_Int new_nb_pgs=0;
        for (E_Int j = 0; j < nb_pgs; ++j)
        {
          if (buffer_flag_pgs[j] == IDX_NONE) continue;
          molecule[new_nb_pgs++] = *(first_pg + j);
        }
        if (new_nb_pgs) 
        {
          new_phs.add(new_nb_pgs, &molecule[0]);
          //std::cout << "old/new nb pgs : " << nb_pgs << "/" << new_nb_pgs << std::endl;
        }
      }
    }
    
    //std::cout << "baffle 5" << std::endl;
    
    PHs=new_phs;
    PHs.updateFacets();

    Vector_t<E_Int> pgnids, phnids;
    remove_unreferenced_pgs(pgnids, phnids);
    
    //std::cout << "baffle 6" << std::endl;
    
    return nb_cells_w_baffles;
  }

  ///
  void get_facets_connexion_nb(Vector_t<E_Int>& nfconnects)
  {
    nfconnects.clear();
    E_Int nb_facets = PGs.size();
    E_Int nb_elts = PHs.size();

    nfconnects.resize(nb_facets, 0);

    for (E_Int i = 0; i < nb_elts; ++i)
    {
      const E_Int * faces = PHs.get_facets_ptr(i);
      E_Int strid = PHs.stride(i);

      for (E_Int n = 0; n < strid; ++n)
      {
        E_Int PGi = *(faces + n) - 1;
        ++nfconnects[PGi];
      }
    }
  }
  
  /// Warning : The coordinates are not cleaned, only the connectivity.
  static E_Int clean_connectivity
  (ngon_t& NG, const K_FLD::FloatArray& f, E_Int ngon_dim, E_Float tolerance, bool remove_dup_phs, bool do_omp,
   std::vector<E_Int>* pgnids=nullptr, std::vector<E_Int>* phnids=nullptr,
   std::vector<E_Float>* Lmin2=nullptr)
  {   
    bool histo = false;
    if (pgnids!=nullptr) histo = true;
    
    // Init pgnids and phnids 
    if (histo)
    {
      K_CONNECT::IdTool::init_inc(*pgnids, NG.PGs.size()); 
      K_CONNECT::IdTool::init_inc(*phnids, NG.PHs.size());
    }
    
    E_Int nb_phs0(NG.PHs.size()), nb_pgs0(NG.PGs.size());
    
    if (nb_pgs0 == 0) // fast return
      return 0;
    
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
    acrd_t fcA(f);
  
    NG.PGs.updateFacets();
    NG.PHs.updateFacets();
  
    // type of NGON: volumic , surfacic or lineic?
    if (ngon_dim != 1 && ngon_dim != 2 && ngon_dim != 3) // not a kwnown value
    {
      for (E_Int i=0; i<nb_pgs0; ++i)
        ngon_dim = std::max(ngon_dim, NG.PGs.stride(i));
      
      if (ngon_dim != 1 && ngon_dim !=2)
        ngon_dim=3;//volumic
    }

    // 1- Referencement unique pour les noeuds confondus par kdtree (les doublons seront supprimes a la fin)
    // si la tolerance est :
    // * positive : absolue
    // * nulle => pas de merge
    // * negative => relative
    E_Int nb_merges = 0;
    if (tolerance >= 0.)
    {
      nb_merges = NG.join_phs(f, tolerance, do_omp);
    }
    else //if (tolerance < 0.)
    {
      //std::cout << "RELATIVE TOL" << std::endl;
      std::vector<double> nodal_metric2;
      if (Lmin2 == nullptr)
      {
        //std::cout << "computing Lmin2" << std::endl;
        NUGA::MeshTool::computeNodalDistance2<K_FLD::FloatArray, ngon_unit>(f, NG.PGs, nodal_metric2);
        Lmin2 = &nodal_metric2;
      }

      //std::cout << "Limn2 size vs crd : " << Lmin2->size() << "/" << f.cols() << std::endl;
      
      double RTOL = -tolerance;
      nb_merges = NG.join_phs(f, *Lmin2, RTOL, do_omp);
    }

    // 2- Elimination des faces degenerees
    Vector_t<E_Int> pgnidstmp, phnidstmp; // required to update the history (PG/PH)
    if (ngon_dim != 1)
    {
      //E_Int nb_degen_faces = 
      NG.remove_degenerated_pgs(ngon_dim, pgnidstmp, phnidstmp);

      //Propagation de pgnidstmp/phnidstmp
      if ((histo) and (not pgnidstmp.empty())) K_CONNECT::IdTool::propagate(pgnidstmp, *pgnids);
      if ((histo) and (not phnidstmp.empty())) K_CONNECT::IdTool::propagate(phnidstmp, *phnids);

      //E_Int nb_consec_changes = 
      NG.PGs.remove_consecutive_duplicated(); //removing duplicated facets : compact representation
    }
    
    // 3- Faces confondues : identification et suppression des references.
    /*bool has_dups = false;*/
    if (ngon_dim == 3) //volumic
    {
      // /*has_dups = */NG.remove_duplicated_pgs(fcA, pgnidstmp);
      /*bool has_dups =*/NG.remove_duplicated_pgs(fcA, pgnidstmp, do_omp);
      
      //Propagation de pgnidstmp/phnidstmp
      if ((histo) and (not pgnidstmp.empty())) K_CONNECT::IdTool::propagate(pgnidstmp, *pgnids);
      if ((histo) and (not phnidstmp.empty())) K_CONNECT::IdTool::propagate(phnidstmp, *phnids);
    }
        
    else if (ngon_dim == 2) //surfacic
      /*has_dups = */NG.remove_duplicated_edges();
    else // lineic
      /*has_dups = */NG.remove_duplicated_nodes();

    // remove duplicated references to PGs within each elements
    /*E_Int nb_phs_dups = */NG.PHs.remove_duplicated(); //fixme : redundant ?
    
    // 4a- Elimination des elts degeneres
    Vector_t<E_Int> toremove;
    E_Int min_nb_facets = ngon_dim + 1;
    NG.PHs.get_degenerated(min_nb_facets, toremove);

    // 5- Elimination des elts doubles : WARNING : do not care of multiple occ in toremove as remove_entities handles it.
    if (remove_dup_phs)
    {
      std::vector<E_Int> duphnids;
      NG.detect_phs_with_same_centroid(f, duphnids);
      for (size_t k = 0; k < duphnids.size(); ++k)
      {
        if (duphnids[k] != (E_Int)k)
          toremove.push_back(k);
      }
    }
      
    /*E_Int nb_removed_phs = */NG.PHs.remove_entities(toremove, phnidstmp);
    if ((histo) and (not phnidstmp.empty())) K_CONNECT::IdTool::propagate(phnidstmp, *phnids);
 
    // 6- Suppression des faces non referencees
    /*E_Int nb_unrefs = */NG.remove_unreferenced_pgs(pgnidstmp, phnidstmp); //Important, sinon tecplot en ASCII (.tp) n'aime pas.
    
    //Propagation de pgnidstmp/phnidstmp
    if ((histo) and (not pgnidstmp.empty())) K_CONNECT::IdTool::propagate(pgnidstmp, *pgnids);
    if ((histo) and (not phnidstmp.empty())) K_CONNECT::IdTool::propagate(phnidstmp, *phnids);
      
    // 7- Compression du ngon aux seuls noeuds utilises
    //ngon_t::compact_to_used_nodes(NG.PGs, f);
    
    E_Int nb_phs1 = NG.PHs.size();
    E_Int nb_pgs1 = NG.PGs.size();

#ifdef DEBUG_NGON_T
  assert (NG.is_consistent(f.cols()));
#endif
    
    return nb_merges + (nb_phs0-nb_phs1) + (nb_pgs0-nb_pgs1);
  }
  
  static void surface_extrema(const ngon_unit& PGs, const K_FLD::FloatArray& crd, E_Float& smin, E_Int& imin, E_Float& smax, E_Int& imax)
  {
    //
    E_Int nb_pgs = PGs.size();
    E_Int i, id;
    E_Float s;
    
    smin = NUGA::FLOAT_MAX;
    imin = IDX_NONE;
    smax = -1.;
    imax = IDX_NONE;
    
    E_Int nb_max_threads = __NUMTHREADS__;
    
    std::vector<E_Int> im(nb_max_threads, IDX_NONE);
    std::vector<E_Float> sm(nb_max_threads, NUGA::FLOAT_MAX);
    std::vector<E_Int> iM(nb_max_threads, IDX_NONE);
    std::vector<E_Float> sM(nb_max_threads, -1.);

#pragma omp parallel shared(sm, im, PGs, crd, nb_pgs) private (i, s, id)
    {
      id = __CURRENT_THREAD__;
    
#pragma omp for 
      for (i=0; i < nb_pgs; ++i)
      {
        s = K_MESH::Polygon::surface<K_FLD::FloatArray, 3>(crd, PGs.get_facets_ptr(i), PGs.stride(i), 1);
        //std::cout << "s : " << s << std::endl;
        if (s < sm[id])
        {
          sm[id]=s;
          im[id]=i;
        }
        if (s > sM[id])
        {
          sM[id]=s;
          iM[id]=i;
        }
      }
    }
    
    for (E_Int i=0; i < nb_max_threads; ++i)
    {
      if (sm[i] < smin)
      {
        imin = im[i];
        smin = sm[i];
      }
      if (sM[i] > smax)
      {
        imax = iM[i];
        smax = sM[i];
      }
    }
  }
  
  ///
  template <typename Triangulator_t>
  static void volume_extrema(ngon_t& ngi, const K_FLD::FloatArray& crd, E_Float& vmin, E_Int& imin, E_Float& vmax, E_Int& imax)
  {
    vmin = NUGA::FLOAT_MAX;
    imin = imax = IDX_NONE;
    vmax = -1.;

    //
    E_Int nb_phs = ngi.PHs.size();
    E_Int nb_pgs = ngi.PGs.size(); 
    
    //std::cout << "build orient..." << std::endl;
    
    // give orient for each PG within each PH
    ngon_unit orienti;
    E_Int err = build_orientation_ngu<Triangulator_t>(crd, ngi, orienti);
    if (err)
      return;
    
    //std::cout << "Triangulate..." << std::endl;
    
    // Triangulate
    
    std::vector<E_Int> xT3(nb_pgs+1);
    // find out size and xranges
    E_Int szT3(0);
     xT3[0] = 0;
    for (E_Int i=0; i < nb_pgs; ++i)
    {
      szT3 += (ngi.PGs.stride(i)-2) ;
      xT3[i+1] = szT3;
    }
    
    K_FLD::IntArray connectT3(3, szT3);
    Triangulator_t t;
    
#ifndef DEBUG_NGON_T
#pragma omp parallel for private (t) shared(connectT3, xT3)
#endif
    for (E_Int i=0; i < nb_pgs; ++i)
    {
      const E_Int* nodes = ngi.PGs.get_facets_ptr(i);
      E_Int nb_nodes = ngi.PGs.stride(i);
      K_MESH::Polygon::triangulate_inplace(t, crd, nodes, nb_nodes, 1, connectT3.col(xT3[i]));//WARNING : connectT3 is filled IN PLACE (not cleared upon entry)
    }
  
    //
    E_Float G[3], v;
    K_FLD::IntArray cT3;
    E_Int T[3], id, i, f;
    
    E_Int nb_max_threads = __NUMTHREADS__;
    
    std::vector<E_Int> im(nb_max_threads, IDX_NONE);
    std::vector<E_Float> vm(nb_max_threads, NUGA::FLOAT_MAX);
    std::vector<E_Int> iM(nb_max_threads, IDX_NONE);
    std::vector<E_Float> vM(nb_max_threads, -1.);
    
    //std::cout << "Volumes..." << std::endl;
    
#pragma omp parallel shared(vm, im, crd, ngi, orienti, connectT3) private(G, v, cT3, T, id, i, f)
   {
      id = __CURRENT_THREAD__;
    
#pragma omp for 
      for (i=0; i < nb_phs; ++i)
      {
        // Gather tha appropriate T3s with the correct orientation
        cT3.clear();
        E_Int nb_f = ngi.PHs.stride(i);
        const E_Int* fs = ngi.PHs.get_facets_ptr(i);
        const E_Int* ori = orienti.get_facets_ptr(i);
      
        for (f=0; f < nb_f; ++f)
        {
          E_Int PGi = *(fs + f) -1;
          E_Int o = *(ori+f);
          for (E_Int k=xT3[PGi]; k < xT3[PGi+1]; ++k)
          {
            T[0] = connectT3(0,k);
            T[1] = connectT3(1,k);
            T[2] = connectT3(2,k);
            if (o == -1) std::swap(T[1], T[2]);
            cT3.pushBack(T, T+3);
          }
        }
      
        K_MESH::Polyhedron<UNKNOWN>::metrics(crd, cT3, v, G);
        v = ::fabs(v);
        if (v < vm[id])
        {
          vm[id]=v;
          im[id]=i;
        }
        if (v > vM[id])
        {
          vM[id]=v;
          iM[id]=i;
        }
      }
   }

   for (E_Int i=0; i < nb_max_threads; ++i)
   {
     if (vm[i] < vmin)
     {
       imin = im[i];
       vmin = vm[i];
     }
     if (vm[i] > vmax)
     {
       imax = im[i];
       vmax = vm[i];
     }
   }
  }

  static void edge_length_extrema(ngon_unit& PGs, const K_FLD::FloatArray& crd, E_Float& Lmin, E_Int& imin, E_Float& Lmax, E_Int& imax)
  {
    Lmin = NUGA::FLOAT_MAX;
    Lmax = -1.;
    imin=imax=IDX_NONE;

    E_Int nb_pgs = PGs.size(), nb_nodes;
    E_Float L2;

    for (E_Int i = 0; i < nb_pgs; ++i)
    {
      nb_nodes = PGs.stride(i);
      const E_Int* pNi = PGs.get_facets_ptr(i);

      for (E_Int j = 0; j < nb_nodes; ++j)
      {
        L2 = NUGA::sqrDistance(crd.col(*(pNi + j) - 1), crd.col(*(pNi + (j + 1) % nb_nodes) - 1), 3);
        imin = (L2 < Lmin) ? i : imin;
        Lmin = (L2 < Lmin) ? L2 : Lmin;
        
        imax = (Lmax < L2) ? i : imax;
        Lmax = (Lmax < L2) ? L2 : Lmax;
      }
    }

    Lmin = ::sqrt(Lmin);
    Lmax = ::sqrt(Lmax);
  }

  struct link_t{
      link_t():count(-1){}
      E_Int nodes[2], count;
      bool add(E_Int n){ 
        if (count == -1)
        {
          nodes[++count] = n;
          return true;
        }
        else if (count == 0)
        {
          if (n != nodes[0])
            nodes[++count] = n;
          return true;
        }
        else // count = 1
        {
          if (n == nodes[0] || n == nodes[1])
            return true;
        }
        return false;//non manifold
      }
    };
  
  static E_Int flag_manifold_nodes(const ngon_unit& PGs, const K_FLD::FloatArray& crd, Vector_t<bool>& mnfld_nodes)
  {
    PGs.updateFacets();

    Vector_t<link_t> links(crd.cols());
    
    mnfld_nodes.resize(crd.cols(), true);
    E_Int nb_nodes;
    for (E_Int i = 0; i < PGs.size(); ++i)
    {
      nb_nodes = PGs.stride(i);
      const E_Int* nodes = PGs.get_facets_ptr(i);

      for (E_Int n = 0; n < nb_nodes; ++n)
      {
        E_Int Ni = *(nodes + n) - 1;
        E_Int Nip1 = *(nodes + (n + 1) % nb_nodes) - 1;

        mnfld_nodes[Ni] = mnfld_nodes[Ni] && links[Ni].add(Nip1);
        mnfld_nodes[Nip1] = mnfld_nodes[Nip1] && links[Nip1].add(Ni);
      }
    }

    return 0;

  }

  /// symplify PGs
  static void simplify_pgs(ngon_unit& PGs, const K_FLD::FloatArray& crd)
  {
    Vector_t<bool> man_nodes;

    flag_manifold_nodes(PGs, crd, man_nodes);

#ifdef DEBUG_NGON_T
    E_Int count{ 0 };
#endif
    Vector_t<E_Int> pgnids, node_ids(crd.cols());//0-based ids
    for (E_Int i = 0; i < crd.cols(); ++i)
    {
      node_ids[i] = (man_nodes[i]) ? IDX_NONE : i;
      //std::cout << "manifold ? : " << man_nodes[i] << " . so new val is : " << node_ids[i] << std::endl;
#ifdef DEBUG_NGON_T
      if (node_ids[i] == IDX_NONE)++count;
#endif
    }

#ifdef DEBUG_NGON_T
    std::cout << " nb of removed nodes : " << count << std::endl;

    count = 0;
    for (size_t u = 0; u < ng.PGs.size(); ++u)
      count += PGs.stride(u);

    std::cout << "cumulated stride before cleaning : " << count << std::endl;
#endif

    //E_Int nb_remove = PGs.remove_facets(node_ids, pgnids, 2);//a surface must have at least 3 nodes

#ifdef DEBUG_NGON_T 
    count = 0;
    for (size_t u = 0; u < PGs.size(); ++u)
      count += PGs.stride(u);

    std::cout << "cumulated stride after cleaning : " << count << std::endl;
#endif
  }
  
 
  /// Clean superfluous nodes on PGs and sync PHs
  static void simplify_pgs (ngon_t& ng, const K_FLD::FloatArray& crd, bool process_externals = true)
  {
    Vector_t<bool> man_nodes;

    flag_manifold_nodes(ng.PGs, crd, man_nodes);

    // unflag external nodes
    if (!process_externals)
    {
      //std::cout << "unflag external nodes" << std::endl;
      ng.PHs._type.clear();
      ng.PGs._type.clear();
      ng.flag_externals(1);

      for (E_Int i = 0; i < ng.PGs.size(); ++i)
      {
        if (ng.PGs._type[i] != 1) continue;
        
        E_Int nnodes = ng.PGs.stride(i);
        const E_Int* nodes = ng.PGs.get_facets_ptr(i);

        for (E_Int j=0; j < nnodes; ++j)
        {
          E_Int Ni = nodes[j] - 1;
          man_nodes[Ni] = false; //reset
        }
      }
    }

#ifdef DEBUG_NGON_T
    E_Int count{0};
#endif
    Vector_t<E_Int> pgnids, node_ids(crd.cols());//0-based ids
    for (E_Int i = 0; i < crd.cols(); ++i)
    {
      node_ids[i] = (man_nodes[i]) ? IDX_NONE : i;
      //std::cout << "manifold ? : " << man_nodes[i] << " . so new val is : " << node_ids[i] << std::endl;
#ifdef DEBUG_NGON_T
      if (node_ids[i] == IDX_NONE)++count;
#endif
    }

#ifdef DEBUG_NGON_T
    std::cout << " nb of removed nodes : " << count << std::endl;

    count=0;
    for (size_t u=0; u < ng.PGs.size(); ++u)
      count += ng.PGs.stride(u);

    std::cout << "cumulated stride before cleaning : " << count << std::endl;
#endif

    E_Int nb_remove = ng.PGs.remove_facets(node_ids, pgnids, 2);//a surface must have at least 3 nodes

    if (nb_remove) //update PHs accordingly
    {
      Vector_t<E_Int> phnids;
      //E_Int nb_remove = 
      ng.PHs.remove_facets(pgnids, phnids);
      //std::cout << "nb PG removed : " << nb_remove << std::endl;
    }

#ifdef DEBUG_NGON_T 
    count=0;
    for (size_t u=0; u < ng.PGs.size(); ++u)
      count += ng.PGs.stride(u);

    std::cout << "cumulated stride after cleaning : " << count << std::endl;
#endif
  }
  
  ///
  template <typename TriangulatorType>
  static E_Int __reorient_skin
(const K_FLD::FloatArray& coord, const ngon_t& NGZ, const ngon_unit& pg_ext, const Vector_t<E_Int>& pg_ext_to_wPG, const ngon_unit& neighbors, Vector_t<E_Int>& orient)
{
  E_Int refPG(IDX_NONE), ori(IDX_NONE);
  E_Int err = __set_ref_PGs_for_orient<TriangulatorType>(coord, NGZ, refPG,ori);
      
  if (err)
    return err;
  
  if (refPG == IDX_NONE) // the mesh might be corrupted, anyway it must not happen
    return 1;
  //assert (refPG != IDX_NONE);

  assert (ori != IDX_NONE);
  
  //find the refPG
  E_Int iref = IDX_NONE;
  for (size_t i=0; i < pg_ext_to_wPG.size(); ++i)
  {
    if (pg_ext_to_wPG[i] == refPG)
    {
      iref = i;
      break;
    }
  }
  assert (iref != IDX_NONE);
  refPG=iref;
  
#ifdef DEBUG_NGON_T 
  /*K_FLD::IntArray tmpE;
  for (std::set<K_MESH::NO_Edge>::const_iterator it = walls.begin(); it != walls.end(); ++it)
    tmpE.pushBack(it->begin(), it->begin()+2);
  tmpE.shift(-1);
  std::cout << tmpE << std::endl;
  medith::write("walls.mesh", coord, tmpE, "BAR");*/
#endif
    
  orient.resize(pg_ext.size(), 1);
  if (ori == -1)
    orient[iref]=-1;
  NUGA::EltAlgo<K_MESH::Polygon>::reversi_connex(pg_ext, neighbors, refPG, orient);
   
  return 0;
}
  
///
template <typename TriangulatorType>
static E_Int __set_ref_PGs_for_orient(const K_FLD::FloatArray& coord, const ngon_t& ng, E_Int& PGref, E_Int& orient)
{
  E_Int err(0);
  K_FLD::ArrayAccessor<K_FLD::FloatArray> ac(coord);
  
  orient=IDX_NONE;
  PGref=IDX_NONE;
  
#ifdef DEBUG_NGON_T
  assert(ng.is_consistent(coord.cols()));
#endif

  if ((E_Int)ng.PHs._type.size() < ng.PHs.size())
  {
    std::cout << "__set_ref_PGs_for_orient : input PHs does not have type set (looking for INITIAL_SKIN or CONNEXION_SKIN)" << std::endl;
    return 1;
  }
  
  for (E_Int PHi = 0; (PHi < ng.PHs.size()) && (PGref == IDX_NONE) && !err; ++PHi)
    if ((ng.PHs._type[PHi] == INITIAL_SKIN) || (ng.PHs._type[PHi] == CONNEXION_SKIN))
      err = __find_reference_orientation_triangle<TriangulatorType>(ng, coord, PHi, PGref, orient);
  
  return err;
}

///
template <typename TriangulatorType>
static E_Int __find_reference_orientation_triangle
(const ngon_t& ng, const K_FLD::FloatArray& coord, E_Int PHi, E_Int& PGgoal, E_Int& orient)
{
  PGgoal = IDX_NONE;
  orient = 0;
  
#ifdef DEBUG_NGON_T
  K_FLD::ArrayAccessor<K_FLD::FloatArray> ac(coord);
  //NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::write("ng_for_ref_orient.plt", ac, ng);
  //NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PGT3s(coord, ng.PGs);
  //NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::highlight_PH<DELAUNAY::Triangulator>(ng, coord, PHi);
#endif
  
  // WARNING : assume here PHi is external
  E_Int nb_pg = ng.PHs.stride(PHi);
  const E_Int* ptr = ng.PHs.get_facets_ptr(PHi);
  
  Vector_t<E_Int> indices;
  indices.insert(indices.end(), ptr, ptr+nb_pg);
  for (size_t i = 0; i < indices.size(); ++i)indices[i] -= 1;//indcies start at one
  
  ngon_unit pgs;
  Vector_t<E_Int> oIds;
  ng.PGs.extract(indices, pgs, oIds);
  
  K_FLD::IntArray connectT3;
  Vector_t<E_Int> T3_to_nPG;
  E_Int err = triangulate_pgs<TriangulatorType>(pgs, coord, connectT3, T3_to_nPG, true, false);
  if (err)
    return err;
  
#ifdef DEBUG_NGON_T
  //medith::write("PHref.mesh", coord, connectT3, "TRI", 0, &T3_to_nPG);
#endif
  
  {
  K_FLD::IntArray::const_iterator pS;
  E_Float Ni[3], P0[3], P1[3];
  E_Int dirp, dirm, PGi, PGj, PGprocessed;
  bool intersect;
  size_t sz = T3_to_nPG.size();
  for (size_t Ti = 0; Ti < sz; ++Ti)
  {
    // test only T3s belonging to skin PGs.
    PGi = oIds[T3_to_nPG[Ti]];
    if (!IS_EXTERNAL(ng.PGs._type[PGi]))
      continue;

 #ifdef DEBUG_NGON_T   
    /*std::vector<E_Int> colors;
    colors.resize(connectT3.cols(), 0);
    colors[Ti]=1;
    TRI_debug::write_wired("refT3wire.mesh", coord, connectT3, true, &colors);*/
#endif   
    
    pS = connectT3.col(Ti);
    K_MESH::Triangle::isoG(coord, pS, P0); // P0 is the center of Ti
    K_MESH::Triangle::normal(coord, pS, Ni); // Ni is the normal to Ti
    NUGA::sum<3>(P0, Ni, P1); //P0P1 == Ni
    
    dirp=dirm=0;
    PGprocessed = IDX_NONE;
    
    E_Float lambda, UV[2], min_d;
    E_Bool  coincident, parallel;
    
    for (size_t Tj = 0; Tj < sz; ++Tj)
    {
      PGj = oIds[T3_to_nPG[Tj]];
      if ((PGj == PGi) || (PGj ==  PGprocessed))
        continue;
      
      PGprocessed = IDX_NONE;
      const E_Int& N1 = connectT3(0,Tj);
      const E_Int& N2 = connectT3(1,Tj);
      const E_Int& N3 = connectT3(2,Tj);
      
      const E_Float* q0 = coord.col(N1);
      const E_Float* q1 = coord.col(N2);
      const E_Float* q2 = coord.col(N3);
      
      K_MESH::Triangle::planeLineMinDistance<3>(q0, q1, q2, P0, P1, EPSILON, true/*tol_is_absolute*/, lambda, UV, parallel, coincident, min_d);
      //intersect = K_MESH::Triangle::intersect<3>(q0, q1, q2, P0, P1, EPSILON, true, u0, u1, tx, overlap);
      //intersect |= (u0 != NUGA::FLOAT_MAX) && (u1 == NUGA::FLOAT_MAX);
      intersect= (((1. - UV[0] -  UV[1]) >= -EPSILON)  && (UV[0] >= -EPSILON) && (UV[1] >= -EPSILON));
      
      if (coincident)
        continue;
      else if (intersect)
      {
        if (lambda > EPSILON)++dirp;
        else if (lambda < -EPSILON)++dirm;
        PGprocessed = PGj;
      }
    }
    
    if (((dirp+dirm)%2 == 1))
    {
      if (dirp%2 == 0)orient = 1;
      else orient = -1;
      PGgoal = PGi;
      break;
    }
  }
  }
  return err;
}

///
template <typename TriangulatorType>
static E_Int reorient_skins(const TriangulatorType& t, const K_FLD::FloatArray& coord, ngon_t& wNG, bool& has_been_reversed)
{
  if (wNG.PHs.size()) assert(wNG.PHs.size() < wNG.PGs.size()); //i.e. not one ph per pg beacuse we need at least a closed PH

  has_been_reversed = false;

#ifdef FLAG_STEP
  NUGA::chrono c;
  c.start();
#endif

  // 1. Get the skin (and connexion skin) PGs
  Vector_t<E_Int> oids;
  ngon_unit pg_ext;
  wNG.PGs.extract_of_type(INITIAL_SKIN, pg_ext, oids);

  if (pg_ext.size() == 0)
    return 0;

#ifdef FLAG_STEP
  std::cout << "__reorient_skins : extract_of_type : " << c.elapsed() << std::endl;
  c.start();
#endif

#ifdef DEBUG_NGON_T
  {
    ngon_t ng(pg_ext, true);
    K_FLD::IntArray cnto;
    ng.export_to_array(cnto);
    medith::write("pgski1.plt", coord, cnto, "NGON");
    
    K_FLD::IntArray connectT3;
    Vector_t<E_Int> T3_to_nPG;
    E_Int err = ngon_t::triangulate_pgs<DELAUNAY::Triangulator>(pg_ext, coord, connectT3, T3_to_nPG, true, false);
    TRI_debug::write_wired("before_orient.mesh", coord, connectT3, true);
  }
#endif

  // 2. Build the neighbourhood for skin PGs
  ngon_unit neighbors;
  K_MESH::Polygon::build_pg_neighborhood(pg_ext, neighbors);
  
#ifdef FLAG_STEP
  std::cout << "__reorient_skins : build_pg_neighborhood : " << c.elapsed() << std::endl;
  std::cout << "reorient_skins : cut non-manifolds in graph..." << std::endl;
  c.start();
#endif

  // update it to treat non-manifoldness as void neighbor.
  {
    // count edges sharing number
    std::map<K_MESH::NO_Edge, E_Int> edge_to_count;
    std::map<K_MESH::NO_Edge, E_Int>::iterator it;
    K_MESH::NO_Edge E;
    for (E_Int i = 0; i < pg_ext.size(); ++i)
    {
      E_Int* pN = pg_ext.get_facets_ptr(i);
      E_Int nb_nodes = pg_ext.stride(i);

      for (E_Int n = 0; n < nb_nodes; ++n)
      {
        E_Int ni = *(pN + n);
        E_Int nj = *(pN + (n + 1) % nb_nodes);
        E.setNodes(ni, nj);
        it = edge_to_count.find(E);
        if (it == edge_to_count.end())
          edge_to_count.insert(std::make_pair(E, 1));
        else
          ++it->second;
      }
    }
    // traverse the pgs and unset the non manifold neighbors
    for (E_Int i = 0; i < pg_ext.size(); ++i)
    {
      E_Int* pN = pg_ext.get_facets_ptr(i);
      E_Int nb_nodes = pg_ext.stride(i);

      E_Int* pKn = neighbors.get_facets_ptr(i);

      for (E_Int n = 0; n < nb_nodes; ++n)
      {
        E_Int ni = *(pN + n);
        E_Int nj = *(pN + (n + 1) % nb_nodes);
        E.setNodes(ni, nj);

        if (edge_to_count[E] != 2)
          *(pKn + n) = IDX_NONE;
      }
    }
  }

#ifdef FLAG_STEP
  std::cout << "__reorient_skins : cut non-manifolds in graph : " << c.elapsed() << std::endl;
  std::cout << "reorient_skins : coloring..." << std::endl;
  c.start();
#endif

  // Color to get connex parts
  E_Int nb_connex = 1;
  Vector_t<E_Int> colors;
  if (wNG.PHs.size() > 1)
  {
    NUGA::EltAlgo<K_MESH::Polygon>::coloring(neighbors, colors);
    nb_connex = 1 + *std::max_element(colors.begin(), colors.end());
  }

#ifdef FLAG_STEP
  std::cout << "__reorient_skins : coloring : " << c.elapsed() << std::endl;
  std::cout << "NB CONNEX FOUND : " << nb_connex << std::endl;
  c.start();
#endif

  E_Int err = 0;
  Vector_t<E_Int> orient;
  if (nb_connex > 1)
  {
    Vector_t<bool> keep;
    ngon_t NGZ;
    for (E_Int z = 0; z < nb_connex; ++z)
    {
#ifdef FLAG_STEP
      c.start();
#endif

      NGZ = wNG;
      NGZ.PHs._type = wNG.PHs._type;
      NGZ.PGs._type = wNG.PGs._type;

      keep.resize(NGZ.PGs.size(), false);
      for (size_t i = 0; i < colors.size(); ++i) keep[oids[i]] = (colors[i] == z);
      NGZ.keep_PHs_having_PGs(keep);

#ifdef FLAG_STEP
      std::cout << "__reorient_skins : prepare connex part number : " << z << " " << c.elapsed() << std::endl;
      c.start();
#endif

      err = __reorient_skin<TriangulatorType>(coord, NGZ, pg_ext, oids, neighbors, orient); //reorient a connex part NGZ of wNG. pg_ext are related to wNG
      if (err)
      {
        std::cout << "could not orient properly zone " << z << std::endl;
        
#ifdef DEBUG_NGON_T
        std::ostringstream o;
        o << "reorient_error_zone_" << z << ".tp";
        K_FLD::IntArray cnto;
        NGZ.export_to_array(cnto);
        medith::write(o.str().c_str(), coord, cnto, "NGON");
#endif
        return err;
      }

#ifdef FLAG_STEP
      std::cout << "__reorient_skins : __reorient_skin number : " << z << " " << c.elapsed() << std::endl;
      c.start();
#endif
    }
  }
  else //one connex zone
  {
    err = __reorient_skin<TriangulatorType>(coord, wNG, pg_ext, oids, neighbors, orient); //reorient a connex part NGZ of wNG. pg_ext are related to wNG

    if (err)
    {
      std::cout << "could not orient properly zone 0" << std::endl;

#ifdef DEBUG_NGON_T
      std::ostringstream o;
      o << "reorient_error_zone_0.tp";
      K_FLD::IntArray cnto;
      wNG.export_to_array(cnto);
      medith::write(o.str().c_str(), coord, cnto, "NGON");
#endif
      return err;
    }

#ifdef FLAG_STEP
    std::cout << "__reorient_skins : __reorient_skin : " << c.elapsed() << std::endl;
    c.start();
#endif
  }

  //Apply new orientation
  for (E_Int i = 0; i < pg_ext.size(); ++i)
  {
    if (orient[i] == -1)
    {
      E_Int PGi = oids[i];
      E_Int s = wNG.PGs.stride(PGi);
      E_Int* p = wNG.PGs.get_facets_ptr(PGi);
      std::reverse(p, p + s);
      has_been_reversed = true;
    }
  }

#ifdef FLAG_STEP
  std::cout << "__reorient_skins : apply orientation  : " << c.elapsed() << std::endl;
#endif

#ifdef DEBUG_NGON_T
//  {
//    K_FLD::IntArray connectT3;
//    Vector_t<E_Int> T3_to_nPG;
//        
//    Vector_t<E_Int> oids, oids_cnx;
//    ngon_unit pg_ext, pg_cnx;
//    wNG.PGs.extract_of_type(INITIAL_SKIN, pg_ext, oids);
//    wNG.PGs.extract_of_type(CONNEXION_SKIN, pg_cnx, oids_cnx);
//
//    oids.insert(oids.end(), oids_cnx.begin(), oids_cnx.end());
//    pg_ext.append(pg_cnx);
//    
//    E_Int err = ngon_t::triangulate_pgs<DELAUNAY::Triangulator>(pg_ext, coord, connectT3, T3_to_nPG, true, false);
//    TRI_debug::write_wired("oriented.mesh", coord, connectT3, true);
//  }
#endif

  return 0;
}

///
static E_Int reorient_connex_PGs(ngon_unit& PGs, bool reverse_first)
{
  // WARNING : use the first Polygon as the reference. Caller can reverse it if required.
  
  PGs.updateFacets();
  
  Vector_t<E_Int> orient(PGs.size(), 1);
  if (reverse_first)
    orient[0]=-1;
  
  ngon_unit neighbors;
  K_MESH::Polygon::build_pg_neighborhood(PGs, neighbors);
  
  NUGA::EltAlgo<K_MESH::Polygon>::reversi_connex(PGs, neighbors, 0/*reference PG*/, orient);
  
  //Apply new orientation
  for (E_Int PGi = 0; PGi < PGs.size(); ++PGi)
  {
    if (orient[PGi] == -1)
    {
      E_Int s = PGs.stride(PGi);
      E_Int* p = PGs.get_facets_ptr(PGi);
      std::reverse(p, p + s);
    }
  }
  
  return 0;
}

static E_Int reorient_connex_PGs(ngon_unit& PGs, bool reverse_first, Vector_t<E_Int>& orient)
{
  // WARNING : use the first Polygon as the reference. Caller can reverse it if required.
   
  std::cout << "passed list : " << PGs.size() << std::endl;
  PGs.updateFacets();
  
  orient.clear(); 
  orient.resize(PGs.size(), 1);
  
  if (reverse_first)
    orient[0]=-1;
  
  ngon_unit neighbors;
  K_MESH::Polygon::build_pg_neighborhood(PGs, neighbors);
  
  NUGA::EltAlgo<K_MESH::Polygon>::reversi_connex(PGs, neighbors, 0/*reference PG*/, orient);
  
  return 0;
}
  
/// coloring-like algo
void
build_F2E(const ngon_unit& neighbors, K_FLD::IntArray& F2E) const
{

  // WARNING : assuming that external PGs have been propoerly oriented toward exterior
  E_Int NB_PHS(PHs.size());
  E_Int NB_PGS(PGs.size());

  F2E.clear();
  F2E.resize(2, NB_PGS, IDX_NONE);

  Vector_t<E_Int> exPH(NB_PHS, IDX_NONE);

  //init for external elements : setting the left
  E_Int nb_phs = neighbors.size();
  for (E_Int PHi = 0; PHi < nb_phs; ++PHi)
  {
    const E_Int* pGi = PHs.get_facets_ptr(PHi);
    const E_Int& nb_neigh = PHs.stride(PHi);
    const E_Int* pKn = neighbors.get_facets_ptr(PHi);

    for (E_Int i = 0; (i<nb_neigh); ++i)
      if (*(pKn + i) == IDX_NONE){
        //std::cout << "PHi/PGi : " << PHi << "/" << *(pGi + i) - 1 << std::endl;
        F2E(0, *(pGi + i) - 1) = PHi;  // left PH
        exPH[PHi] = *(pGi + i); //store i-th external PG for this external PH
        break;
      }
  }

  E_Int PHseed(0), PH, iref, PGref;
  Vector_t<bool> processed(NB_PHS, false);
  ngon_unit PGneighbors, pgs;
  Vector_t<E_Int> oids, orient, cpool;
  bool revers;

  while (1) // as many iter as connex bits
  {
    while ((PHseed < NB_PHS) && ((processed[PHseed] == true) || (exPH[PHseed] == IDX_NONE))) ++PHseed; // take first external element
    if (NB_PHS - 1 < PHseed)
      return;

    cpool.push_back(PHseed);

    while (!cpool.empty())
    {
      PH = cpool.back();
      cpool.pop_back();

      if (processed[PH])
        continue;

      processed[PH] = true;

      const E_Int* pGi = PHs.get_facets_ptr(PH);
      const E_Int& nb_neigh = PHs.stride(PH);
      const E_Int* pKn = neighbors.get_facets_ptr(PH);

      // build PG neighborhood for that given PH
      PGs.extract(pGi, nb_neigh, pgs, oids);
      K_MESH::Polygon::build_pg_neighborhood(pgs, PGneighbors);

      orient.clear();
      orient.resize(nb_neigh, 1);
      revers = false;
      PGref = exPH[PH]; //satrt at 1, can be negative
      if (PGref < 0)
      {
        revers = true;
        PGref = -PGref;
      }

      //found it in oids
      iref = IDX_NONE;
      for (size_t i = 0; i < oids.size(); ++i)
        if (PGref == (oids[i] + 1)){ iref = i; break; }
      assert(iref != IDX_NONE);
      if (revers)
        orient[iref] = -1;

      NUGA::EltAlgo<K_MESH::Polygon>::reversi_connex(pgs, PGneighbors, iref, orient);

      for (E_Int i = 0; (i<nb_neigh); ++i)
      {
        E_Int PGi = *(pGi + i);
        const E_Int& PHn = *(pKn + i);

        F2E(0, PGi - 1) = PH;
        F2E(1, PGi - 1) = PHn;

        if (PHn == IDX_NONE)
          continue;

        exPH[PHn] = -PGi; //starts at 1 so logic with negative value is OK

        if (orient[i] == -1){ //swap
          std::swap(F2E(0, PGi - 1), F2E(1, PGi - 1));
          exPH[PHn] = PGi;
        }

        if (!processed[PHn])
          cpool.push_back(PHn);
      }
    }
  }
}

// non-oriented F2E is formatted as cassiopee : non-interleaved, 0 for none, start index at 1
void
build_noF2E(K_FLD::FldArrayI& F2E) const
{
  // rule : first fill the left, then the right => upon exit, at least F2E is oriented for exterior PGs
  PGs.updateFacets();
  PHs.updateFacets();

  E_Int NB_PHS(PHs.size());
  E_Int NB_PGS(PGs.size());
  
  //std::cout << "build_noF2E : phs/pgs : " << NB_PHS << "/" << NB_PGS << std::endl;

  F2E.clear();
  F2E.resize(2, NB_PGS, 0);
    
  for (E_Int i = 0; i < NB_PHS; ++i)
  {
    E_Int stride = PHs.stride(i);
    const E_Int* pPGi = PHs.get_facets_ptr(i);
    
    for (E_Int n = 0; n < stride; ++n)
    {
      E_Int PGi = *(pPGi + n) - 1;
      E_Int k= (F2E(PGi, 1) == 0) ? 1 : 2;
      
      F2E(PGi, k) = i+1;
    }
  }
}

void
build_noF2E(K_FLD::IntArray& F2E) const
{
  // rule : first fill the left, then the right => upon exit, at least F2E is oriented for exterior PGs
  PGs.updateFacets();
  PHs.updateFacets();

  E_Int NB_PHS(PHs.size());
  E_Int NB_PGS(PGs.size());
  
  //std::cout << "build_noF2E : phs/pgs : " << NB_PHS << "/" << NB_PGS << std::endl;

  F2E.clear();
  F2E.resize(2, NB_PGS, IDX_NONE);
    
  for (E_Int i = 0; i < NB_PHS; ++i)
  {
    E_Int stride = PHs.stride(i);
    const E_Int* pPGi = PHs.get_facets_ptr(i);
    
    for (E_Int n = 0; n < stride; ++n)
    {
      E_Int PGi = *(pPGi + n) - 1;
      E_Int k= (F2E(0, PGi) == IDX_NONE) ? 0 : 1;
      
      F2E(k, PGi) = i;
    }
  }
}

template <typename TriangulatorType>
static E_Int build_orientation_ngu(const K_FLD::FloatArray& crd, ngon_t& ng, ngon_unit& orient)
{
  // WARNING : type are lost : fixme
  
  ng.flag_externals(1);
  TriangulatorType dt;
  bool has_been_reversed;
  E_Int err = ngon_t::reorient_skins(dt, crd, ng, has_been_reversed);
  if (err)
    return err;

  K_FLD::IntArray F2E2;
  ngon_unit neighbors;
  ng.build_ph_neighborhood(neighbors);
  ng.build_F2E(neighbors, F2E2);

  orient = ng.PHs;

  E_Int nb_phs = ng.PHs.size();

  for (E_Int PHi = 0; PHi < nb_phs; ++PHi)
  {
    E_Int nb_pgs = ng.PHs.stride(PHi);
    E_Int* pPGi = ng.PHs.get_facets_ptr(PHi);

    for (E_Int j = 0; j < nb_pgs; ++j)
    {
      E_Int PGi = *(pPGi + j) - 1;
      if (F2E2(0, PGi) == PHi)
        orient.get_facet(PHi, j) = 1;
      else
        orient.get_facet(PHi, j) = -1;
      assert(F2E2(0, PGi) == PHi || F2E2(1, PGi) == PHi);
    }
  }
  return 0;
}

static E_Int flag_neighbors(const ngon_t& ng, Vector_t<bool>& flag, bool overwrite=false)
{
  ng.PHs.updateFacets();
  ng.PGs.updateFacets();
  
  E_Int nb_phs = ng.PHs.size();
  if (nb_phs == 0)
    return 0;
  
#ifdef DEBUG_NGON_T
  std::cout << "number of phs : " << nb_phs << std::endl;
  std::cout << "flag size : " << flag.size() << std::endl;
#endif 
  if (flag.size() != (size_t)nb_phs)
    return 1;
  
  Vector_t<bool> new_flags(nb_phs, false);
  
  ngon_unit neighbors;
  ng.build_ph_neighborhood(neighbors);

  for (E_Int i = 0; i < nb_phs; ++i)
  {
    if (!flag[i])
      continue;
    new_flags[i]=overwrite ? false : true;
    
    E_Int nb_neighs = neighbors.stride(i);
    
    for (E_Int n=0; n < nb_neighs; ++n)
    {
      E_Int j = neighbors.get_facet(i,n);
      if (j == IDX_NONE)
        continue;
      new_flags[j]=true;
    }
  }
  
  flag = new_flags;

  return 0;
}

///
template<typename TriangulatorType>
static E_Int volumes (const K_FLD::FloatArray& crd, const ngon_t& ng, std::vector<E_Float>& vols, bool all_pgs_cvx, bool new_algo)
{
  ng.PGs.updateFacets();
  ng.PHs.updateFacets();
  E_Int nb_phs = ng.PHs.size();
  //std::cout << "nb of pHs : " << nb_phs << std::endl;

  vols.clear();
  vols.resize(nb_phs, -1.);

  E_Float Gdum[3];
  TriangulatorType dt;

  E_Float v;
  E_Int err, errcount(0);

  
  if (!new_algo)
  {
#ifndef DEBUG_NGON_T
#pragma omp parallel for private(err, dt, v, Gdum) reduction(+:errcount)
#endif
    for (E_Int i = 0; i < nb_phs; ++i){
      err = K_MESH::Polyhedron<UNKNOWN>::metrics<TriangulatorType>(dt, crd, ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i), v, Gdum);
      if (!err) vols[i] = v;
      else ++errcount;
    }
  }
  else
  {
#ifndef DEBUG_NGON_T
#pragma omp parallel for private(err, dt, v, Gdum) reduction(+:errcount)
#endif
    for (E_Int i = 0; i < nb_phs; ++i){
      err = K_MESH::Polyhedron<UNKNOWN>::metrics2<TriangulatorType>(dt, crd, ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i), v, Gdum, all_pgs_cvx);
      v = ::fabs(v);
      if (!err) vols[i] = v;
      else ++errcount;
    }
  }
 
  if (errcount) std::cout << "nb of volume error : " << errcount << std::endl;

  return 0;
}

///
template<typename TriangulatorType>
static E_Int volume (const K_FLD::FloatArray& crd, const ngon_t& ng, E_Float& total_vol, std::vector<E_Float>* coefs = 0, bool new_algo=true)
{
  total_vol = 0.;
  std::vector<E_Float> vols;

  E_Int err = volumes<TriangulatorType>(crd, ng, vols, new_algo);
  
  if (err){
    std::cout << "total volume compuattion failed" << std::endl;
    return err;
  }
  
  if (coefs && coefs->size() != vols.size())
  {
    std::cout << "total volume computation failed : mismatch between nb of entities and coefs size" << std::endl;
    return 1;
  }
  
  if (!coefs)
    for (E_Int i = 0; i < vols.size(); ++i)total_vol += vols[i];
  else
    for (E_Int i = 0; i < vols.size(); ++i)total_vol += (*coefs)[i]*vols[i];
    
  
  return 0;
}

///
static E_Int centroids(const ngon_unit& PGS, const K_FLD::FloatArray& crd, K_FLD::FloatArray& centroids)
{
  PGS.updateFacets();

  E_Int nb_pgs = PGS.size();
  //std::cout << "nb of pgs : " << nb_pgs << std::endl;

  centroids.clear();
  centroids.resize(3, nb_pgs, 0.);

  for (E_Int i = 0; i < nb_pgs; ++i)
    //K_MESH::Polygon::centroid<3>(crd, PGS.get_facets_ptr(i), PGS.stride(i), 1, centroids.col(i));
    K_MESH::Polygon::iso_barycenter(crd, PGS.get_facets_ptr(i), PGS.stride(i), 1, centroids.col(i));

  return 0;
}

///
template<typename TriangulatorType>
static E_Int centroids(const ngon_t& ng, const K_FLD::FloatArray& crd, K_FLD::FloatArray& centroids, bool new_algo=true)
{
  ng.PGs.updateFacets();
  ng.PHs.updateFacets();

  E_Int nb_phs = ng.PHs.size();
  //std::cout << "nb of pHs : " << nb_phs << std::endl;

  centroids.clear();
  centroids.resize(3, nb_phs, 0.);

  E_Float v;
  E_Int err, errcount=0;
  std::vector<E_Float> vols(nb_phs);
  TriangulatorType dt;

  if (!new_algo)
  {
#ifndef DEBUG_NGON_T
#pragma omp parallel for private(err, dt, v) reduction(+:errcount)
#endif
    for (E_Int i = 0; i < nb_phs; ++i){
      err = K_MESH::Polyhedron<UNKNOWN>::metrics<TriangulatorType>(dt, crd, ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i), v, centroids.col(i));
      if (!err) vols[i] = v;
      else ++errcount;
    }
  }
  else
  {
#ifndef DEBUG_NGON_T
#pragma omp parallel for private(err, dt, v) reduction(+:errcount)
#endif
    for (E_Int i = 0; i < nb_phs; ++i){
      err = K_MESH::Polyhedron<UNKNOWN>::metrics2<TriangulatorType>(dt, crd, ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i), v, centroids.col(i), false);
      v = ::fabs(v);
      if (!err) vols[i] = v;
      else ++errcount;
    }
  }

  return 0;
}

///
static E_Int check_phs_closure(const ngon_t& ng)
{
  E_Int nb_phs = ng.PHs.size();
  NUGA::non_oriented_edge_set_type edges;
  //
  for (E_Int i = 0; i < nb_phs; ++i)
  {
    bool is_closed = K_MESH::Polyhedron<UNKNOWN>::is_closed(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i), edges);
    if (!is_closed)
    {
      std::cout << "Cell " << i << " is not closed !" << std::endl;
      return 1;
    }
  }
  std::cout << "check_phs_closure : OK, all cells are closed" << std::endl;
  return 0;
}


#ifdef DEBUG_NGON_T
// Writes a single zone into a file (element type must be specified to distinguish QUAD/TETRA, it is not ambiguous otherwise.)
static E_Int write
(const char* filename, const K_FLD::FloatArray& crd, const ngon_t& ng, const std::vector<bool>* mask=0, const std::vector<E_Int>* colors=0)
{
  K_FLD::IntArray cnt;
  ng.export_to_array(cnt);
  return medith::write(filename, crd, cnt, "NGON", mask, colors);
}
#endif

///
/*E_Int keep_PHs(const Vector_t<E_Int>& PHs_tokeep)
{
  ngon_unit keptPHs;
  for (size_t i = 0; i < PHs_tokeep.size(); ++i)
    keptPHs.__add(PHs, PHs_tokeep[i]);

  PHs=keptPHs;
  PHs.updateFacets();

  return 0;
}*/
  
///
E_Int remove_phs(const std::set<E_Int>& PHslist)
{
  if (PHslist.empty()) return 0;
  
  ngon_unit keptPHs;
 
  size_t sz = PHs.size();
  for (size_t i = 0; i < sz; ++i)
  {
    if (PHslist.find(i) == PHslist.end())
      keptPHs.__add(PHs, i);
  }
      
  PHs=keptPHs;
  PHs.updateFacets();
  
  Vector_t<E_Int> pgnids, phnids;
  remove_unreferenced_pgs(pgnids, phnids);
    
  return 0;
}
  
///
/*E_Int reorient_PGs(const Vector_t<E_Int>& PGs_toreverse)
{
  E_Int sz1(PGs_toreverse.size()), sz2;
  for (size_t i = 0; i < sz1; ++i)
  {
    const E_Int& PGi = PGs_toreverse[i];
    sz2 = PGs.get_facets_nb(PGi) >> 1; //half size
    //do the reverse by swapping
    E_Int* ptr = PGs.get_facets_ptr(PGi);
    for (size_t j = 0; j < sz2; ++j)
      std::swap(ptr[j], ptr[sz2-1-j]);
  }
  return 0;
}*/

  static E_Float toto_surface(const std::vector<E_Int>& indices, const K_FLD::FloatArray& crd)
  {
    K_SEARCH::BBox3D box;

    box.compute(crd, indices);

    E_Float dx = box.maxB[0] - box.minB[0];
    E_Float dy = box.maxB[1] - box.minB[1];
    E_Float dz = box.maxB[2] - box.minB[2];
    
    return std::max(dx*dy, std::max(dx*dz, dy*dz));

  }

  static E_Int update_PGs_flags_for_aggregating(ngon_t& ng, E_Float vmin, const Vector_t<E_Float>& vs, const K_FLD::FloatArray& crd, ngon_unit& neighbors)
  {
    E_Int nb_mates = 0;

    E_Int nb_phs = ng.PHs.size();
    Vector_t<E_Float> surfaces_per_K(nb_phs, 0.), surfaces(ng.PGs.size(), 0.);
    Vector_t<E_Int> indices;
    Vector_t<E_Int> attractor(nb_phs, IDX_NONE);

    std::vector<E_Int>& flags = ng.PGs._type;
    
    //flags.resize(ng.PGs.size(), 0);

    // aggregate K to the best mate (meaning one way : K "moves" into best Kn)
    for (E_Int K = 0; K < nb_phs; ++K)
    {
      if (attractor[K] == K) continue; // cannot move since it's an attractor
      if (vs[K] >= vmin) continue;

      const E_Int* pKn = neighbors.get_facets_ptr(K);
      const E_Int* pFn = ng.PHs.get_facets_ptr(K);
      E_Int nb_pgs = ng.PHs.stride(K);

      for (E_Int j = 0; j < nb_pgs; ++j)
      {
        const E_Int& Kn = *(pKn + j);
        E_Int Fn = *(pFn + j) - 1;

        if (Kn == IDX_NONE) continue;

        // append to the surface between K and Kn
        E_Float & sn = surfaces[Fn];
        if (sn == 0.) //compute it. fixme better crterium as it can be really 0.
        {
          indices.clear();
          indices.insert(indices.end(), ng.PGs.get_facets_ptr(Fn), ng.PGs.get_facets_ptr(Fn) + ng.PGs.stride(Fn));
          K_CONNECT::IdTool::shift(indices, -1);
          sn = toto_surface(indices, crd);
        }
        E_Int attracK = (attractor[Kn] == IDX_NONE)? Kn : attractor[Kn];
        surfaces_per_K[attracK] += sn;
      }
      
      // now get the best
      E_Int jmax(IDX_NONE);
      E_Float smax = -NUGA::FLOAT_MAX;
      for (E_Int j = 0; j < nb_pgs; ++j)
      {
        const E_Int& Kn = *(pKn + j);
        //E_Int Fn = *(pFn + j) - 1;

        if (Kn == IDX_NONE) continue;
        E_Int attracK = (attractor[Kn] == IDX_NONE)? Kn : attractor[Kn];
        if (surfaces_per_K[attracK] > smax)
        {
          jmax = j;
          smax = surfaces_per_K[attracK];
        }
        surfaces_per_K[attracK] = 0.;//reinit for next i.
      }

      if (jmax == IDX_NONE) // shoud this happen ?
        continue;

      // disable all the common faces with this attractor AND those shared by neighbor having the same attractor
      E_Int Kmax = neighbors.get_facet(K, jmax);
      
      attractor[Kmax] = Kmax; // set as attractor
      attractor[K] = Kmax;
      ++nb_mates;
    }

    // disable the face inside aggregates
    for (E_Int K = 0; K < nb_phs; ++K)
    {
      if (attractor[K] == IDX_NONE) continue; // not involved

      const E_Int* pKn = neighbors.get_facets_ptr(K);
      const E_Int* pFn = ng.PHs.get_facets_ptr(K);
      E_Int nb_pgs = ng.PHs.stride(K);

      for (E_Int j = 0; j < nb_pgs; ++j)
      {
        const E_Int& Kn = *(pKn + j);
        if (Kn == IDX_NONE) continue;
        if (attractor[K] == attractor[Kn])
        {
          const E_Int& Fn = *(pFn + j) - 1;
          flags[Fn] = 1;
        }
      }
    }

    return nb_mates;
  }

  /// 
  static E_Int check_neigh_consistency(const ngon_t& ng, const ngon_unit& neighbors)
  {
    E_Int nb_neighs(0);
    for (size_t K = 0; K < ng.PHs.size(); ++K)
    {
      E_Int nb_pgs = ng.PHs.stride(K);
      const E_Int* pKn = neighbors.get_facets_ptr(K);
      //const E_Int* pFn = ng.PHs.get_facets_ptr(K);

      for (size_t j = 0; j < nb_pgs; ++j)
      {
        const E_Int& Kn = *(pKn + j);
        //E_Int Fn = *(pFn + j) - 1;
        assert(K != Kn);
        if (Kn == IDX_NONE)
        {
          ++nb_neighs;
          //assert(ng.PGs._type[Fn] == 0);
        }
      }
    }
    return nb_neighs;
  }
  /// Warning : The coordinates are not cleaned, only the connectivity.
  static E_Int aggregate_phs(ngon_t& ng, ngon_unit& neighbors, Vector_t<E_Int>& nids)
  {
    E_Int nb_phs(ng.PHs.size());
    ngon_unit new_neigh, new_phs;
    Vector_t<E_Int> processed(nb_phs, false), molec, neigh_molec;
    std::vector<E_Int>& flags = ng.PGs._type;
    nids.resize(nb_phs);
    E_Int nnid = 0;
    std::set<E_Int> pool;
        
    for (E_Int i = 0; i < nb_phs; ++i)
    {
      if (processed[i])
        continue;

      E_Int nb_pgs = ng.PHs.stride(i);
      bool toprocess = false;
      //const E_Int* pKn = neighbors.get_facets_ptr(i);
      const E_Int* pFn = ng.PHs.get_facets_ptr(i);
      for (size_t j = 0; (j < nb_pgs) && !toprocess; ++j)
      {
        E_Int Fn = *(pFn + j) - 1;
        toprocess = (flags[Fn] == 1);
      }

      if (!toprocess)
      {
        new_phs.add(nb_pgs, ng.PHs.get_facets_ptr(i));
        new_neigh.add(nb_pgs, neighbors.get_facets_ptr(i));
        nids[i] = nnid++;
        if (new_neigh.size() > 1)
          new_neigh.get_facet(1, 0);
        continue; // cannot add to a separate ngon_t as doesnt mean that this will not be an aggregate attractor
      }
      pool.clear();
      pool.insert(i);
      molec.resize(1);
      neigh_molec.resize(1);

      while (!pool.empty())
      {
        E_Int K = *pool.begin();
        pool.erase(K);
        processed[K] = true;
        nids[K] = nnid;
        const E_Int* pKn = neighbors.get_facets_ptr(K);
        const E_Int* pFn = ng.PHs.get_facets_ptr(K);
        
        nb_pgs = ng.PHs.stride(K);
        for (size_t j = 0; j < nb_pgs; ++j)
        {
          const E_Int& Kn = *(pKn + j);
          E_Int Fn = *(pFn + j) -1;
          assert(K != Kn);
          assert(!((Kn == IDX_NONE) && (flags[Fn] == 1)));
          
          if (Kn == IDX_NONE || flags[Fn] == 0) // including case for externals
          {
            const E_Int& Fn = *(pFn + j);
            molec.push_back(Fn);
            neigh_molec.push_back(Kn);
          }
          else if (flags[Fn] == 1 && !processed[Kn])//to be aggregated
            pool.insert(Kn);
        }
      }
      molec[0] = molec.size() - 1;
      neigh_molec[0] = neigh_molec.size() - 1;
      assert(molec[0] * neigh_molec[0]);
      new_phs.add(molec);
      new_neigh.add(neigh_molec);
      ++nnid;
    }

    //update the data
    ng.PHs = new_phs;
    neighbors = new_neigh;
    
    // propagate new ids into neighbors
    for (size_t i = 0; i < neighbors.size(); ++i)
    {
      for (size_t j = 0; j < neighbors.stride(i); ++j)
      {
        E_Int& Kn = neighbors.get_facet(i, j);
        //assert(i != Kn);
        if (Kn == IDX_NONE) continue;
        Kn = nids[Kn];
      }
        
    }
    
    ng.PHs.updateFacets();
    neighbors.updateFacets();

    return 0;

  }

  static void shared_faces_list(const ngon_t& ng, E_Int PHi, E_Int PHj, std::vector<E_Int>& shared_pgs)
  {
    shared_pgs.clear();
    std::set<E_Int> tmp;
    const E_Int* begin = ng.PHs.get_facets_ptr(PHi); 
    tmp.insert(begin, begin+ng.PHs.stride(PHi));
  
    const E_Int* ptr = ng.PHs.get_facets_ptr(PHj); 
    for (E_Int i=0; i < ng.PHs.stride(PHj); ++i)
    {
      if (!tmp.insert(*(ptr+i)).second)//already in
        shared_pgs.push_back(*(ptr+i));
    }
  }

  static void edge_shell(const ngon_t& ng, E_Int N0, E_Int N1, E_Int PHseed, const ngon_unit& neighbors, Vector_t<E_Int>& PHlist)
  {
    PHlist.clear();

    Vector_t<E_Int> pool;
    std::map<E_Int, bool> processed;

    pool.push_back(PHseed);

    while (!pool.empty())
    {
      E_Int PH = pool.back();
      pool.pop_back();

      processed[PH] = true;
      PHlist.push_back(PH);

      E_Int nb_pgs = ng.PHs.stride(PH);
      const E_Int* pgs = ng.PHs.get_facets_ptr(PH);
      const E_Int* neighs = neighbors.get_facets_ptr(PH);

      for (E_Int i = 0; i < nb_pgs; ++i)
      {
        E_Int PGi = *(pgs + i) - 1;
        E_Int nb_nodes = ng.PGs.stride(PGi);
        const E_Int* nodes = ng.PGs.get_facets_ptr(PGi);
        

        // does PGi contain (N0,N1) ?
        bool is_in = false;
        for (E_Int n = 0; n < nb_nodes; ++n)
        {
          E_Int Ni = *(nodes + n);
          E_Int Nj = *(nodes + (n + 1) % nb_nodes);

          if ( (Ni == N0 && Nj == N1) || (Nj == N0 && Ni == N1) )
          {
            is_in = true;
            break;
          }
        }
        if (!is_in)
          continue;

        // is the neighbor unprocesed ?
        E_Int PHn = *(neighs + i);// neighbors.get_facet(PH, i);
        if (PHn == IDX_NONE)
          continue;
        if (processed.find(PHn) != processed.end())
          continue;

        pool.push_back(PHn);
      }
    }
  }

  static void node_shell(const ngon_t& ng, E_Int N0, E_Int PHseed, const ngon_unit& neighbors, Vector_t<E_Int>& PHlist)
  {
    PHlist.clear();

    Vector_t<E_Int> pool;
    std::map<E_Int, bool> processed;

    pool.push_back(PHseed);

    while (!pool.empty())
    {
      E_Int PH = pool.back();
      pool.pop_back();

      processed[PH] = true;
      PHlist.push_back(PH);

      E_Int nb_pgs = ng.PHs.stride(PH);
      const E_Int* pgs = ng.PHs.get_facets_ptr(PH);
      const E_Int* neighs = neighbors.get_facets_ptr(PH);

      for (E_Int i = 0; i < nb_pgs; ++i)
      {
        E_Int PGi = *(pgs + i) - 1;
        E_Int nb_nodes = ng.PGs.stride(PGi);
        const E_Int* nodes = ng.PGs.get_facets_ptr(PGi);


        // does PGi contain N0 ?
        bool is_in = false;
        for (E_Int n = 0; n < nb_nodes; ++n)
        {
          if (*(nodes + n) == N0)
          {
            is_in = true;
            break;
          }
        }
        if (!is_in)
          continue;

        // is the neighbor unprocesed ?
        E_Int PHn = *(neighs + i);// neighbors.get_facet(PH, i);
        if (PHn == IDX_NONE)
          continue;
        if (processed.find(PHn) != processed.end())
          continue;

        pool.push_back(PHn);
      }
    }
  }

///
static void ph_shell(const ngon_t& ng, E_Int PHi, const ngon_unit& neighbors, Vector_t<E_Int>& shellPHs, Vector_t<E_Int>&boundPGs, Vector_t<bool>& wprocess)
{
  wprocess.clear();
  shellPHs.clear();
  boundPGs.clear();

  wprocess.resize(ng.PHs.size(), false);
  std::set<E_Int> unodes;
  std::vector<E_Int> uniqnodes;
  
  K_MESH::Polyhedron<0>::unique_nodes(ng.PGs, ng.PHs.get_facets_ptr(PHi), ng.PHs.stride(PHi), uniqnodes);
  unodes.insert(ALL(uniqnodes));

  std::vector<E_Int> pool;
  std::set<E_Int> shell;
  
  // init : put PHi and its first neighborhood in the shell

  pool.push_back(PHi);
  shell.insert(PHi);
  wprocess[PHi] = true;

  //
  const E_Int * neighs = neighbors.get_facets_ptr(PHi);
  for (E_Int n=0; n < neighbors.stride(PHi); ++n)
  {
    if (neighs[n] != IDX_NONE) 
    {
      pool.push_back(neighs[n]);
      wprocess[neighs[n]] = true;
      shell.insert(neighs[n]);
    }
  }

  while(!pool.empty())
  {
    E_Int Ki = pool.back();
    pool.pop_back();

    E_Int nneighs = neighbors.stride(Ki);
    const E_Int * neighs = neighbors.get_facets_ptr(Ki);
    const E_Int* faces = ng.PHs.get_facets_ptr(Ki);

    for (E_Int n=0; n < nneighs; ++n)
    {
      E_Int Kn = neighs[n];

      if (Kn == IDX_NONE) continue;
      if (wprocess[Kn] == true) continue;

      // check that shared PG has a node of PHi
      const E_Int* PGnodes = ng.PGs.get_facets_ptr(faces[n]-1);
      E_Int nnodes = ng.PGs.stride(faces[n]-1);

      bool connected=false;
      for (E_Int u=0; (u < nnodes) && !connected; ++u)
        connected = (unodes.find(PGnodes[u]) != unodes.end());
      if (!connected) continue; 

      wprocess[Kn]=true;
      pool.push_back(Kn);
      shell.insert(Kn);
    }

  }

  shellPHs.insert(shellPHs.end(), ALL(shell));

  // shell PGs are unicly referenced PGs in shell PHs
  std::set<E_Int> sboundPGs;
  for (size_t i=0; i < shellPHs.size(); ++i)
  {
    E_Int PHi = shellPHs[i];
    const E_Int * faces = ng.PHs.get_facets_ptr(PHi);
    E_Int nfaces = ng.PHs.stride(PHi);

    for (E_Int f=0; f < nfaces; ++f)
      if (!sboundPGs.insert(faces[f]).second)//already in
        sboundPGs.erase(faces[f]);
  }

  boundPGs.insert(boundPGs.end(), ALL(sboundPGs));
  //K_CONNECT::IdTool::shift(boundPGs, -1); //0-based

}
  
///
static E_Int extrude_faces
(K_FLD::FloatArray& coord, ngon_t& wNG, const Vector_t<E_Int>& PGlist, E_Float height_factor, bool build_cells, 
 eExtrudeStrategy strategy = CST_ABS, E_Int smooth_iters = 0, Vector_t<E_Int>* topPGlist = 0)
{
  
  if (PGlist.empty()) return 0;
  
#ifdef DEBUG_NGON_T
  E_Int minid = *std::min_element(PGlist.begin(), PGlist.end());
  E_Int maxid = *std::max_element(PGlist.begin(), PGlist.end());
  E_Int nb_pgs = wNG.PGs.size();
  assert (minid >= 0);
  assert (maxid < nb_pgs);

  if (strategy == CST_REL_MEAN) std::cout      << "CST_REL_MEAN" << std::endl;
  else if (strategy == CST_REL_MIN) std::cout  << "CST_REL_MIN" << std::endl; 
  else if (strategy == VAR_REL_MEAN) std::cout << "VAR_REL_MEAN" << std::endl; 
  else if (strategy == VAR_REL_MIN)  std::cout << "VAR_REL_MIN" << std::endl;

  std::cout << "input factor : " << height_factor << std::endl;


#endif
  
  if (height_factor <= 0.) return 1;
  
  // 1. compute nodal normals and sizes

  ngon_unit ghost_pgs;
  Vector_t<E_Int> obids; //initial bottom ids  : for indirection when building ghost cells at the end
  wNG.PGs.extract(PGlist, ghost_pgs, obids);

  K_FLD::FloatArray normals;
  E_Int nb_new_points = NUGA::MeshTool::computeNodeNormals(coord, ghost_pgs, normals, smooth_iters);
  
  // Fix uncomputed normals by averaging

  // Compute layer's height
  std::vector<E_Float> heightv(coord.cols(), height_factor);
  if (strategy != CST_ABS)
  {
    E_Float Lcomp(0.);
    height_factor =  std::min(1., height_factor); // 100% max
    K_FLD::FloatArray L;
    NUGA::MeshTool::computeIncidentEdgesSqrLengths(coord, ghost_pgs, L);

    E_Int count(0);
    for (E_Int i=0; i < L.cols(); ++i)
    {
      if (L(0,i) == NUGA::FLOAT_MAX) continue;
      if (L(1,i) == -NUGA::FLOAT_MAX) continue;
      ++count;
      if (strategy == CST_REL_MEAN) Lcomp += 0.5 *(::sqrt(L(0,i)) + ::sqrt(L(1,i)));           // GLOBAL MEAN
      else if (strategy == CST_REL_MIN) Lcomp += ::sqrt(L(0,i));                               // GLOBAL MEAN OF MINS
      else if (strategy == VAR_REL_MEAN) heightv[i] *= 0.5 *(::sqrt(L(0,i)) + ::sqrt(L(1,i))); // LOCAL MEAN
      else if (strategy == VAR_REL_MIN) heightv[i] *= ::sqrt(L(0,i));                          // LOCAL MIN
    }
    if (count) Lcomp /= count;

    if (strategy == CST_REL_MEAN || strategy == CST_REL_MIN)
      for (size_t i=0; i < heightv.size(); ++i)heightv[i] *= Lcomp;

  }
  
  E_Int nb_points = coord.cols();
  E_Int nb_normals(normals.cols());
  assert (nb_normals <= nb_points);
  
  // 2. extrude or create and append smoothly the ghost cell without requiring a clean_conenctivity
  
  if (!build_cells)
  {
    for (E_Int i = 0; i < nb_normals; ++i)
    {
      const E_Float* norm = normals.col(i);
      if (norm[0] == NUGA::FLOAT_MAX) continue;
      E_Float* p = coord.col(i);
      NUGA::sum<3>(1., p, heightv[i], norm, p); // move point
    }
    return 0;
  }
  

  // a. conversion NGON3D -> NGON2D
  ngon_t ng2D;
  ngon_t::export_surfacic_FN_to_EFN(ghost_pgs, ng2D);
  // sync with wNG for future bulks appending
  E_Int shft = wNG.PGs.size();
  ng2D.PHs.shift(shft);
  // b. image points
  coord.resize(3, nb_points + nb_new_points);
  Vector_t<E_Int> img(coord.cols(), IDX_NONE);
  E_Int nid(nb_points);
  
  for (E_Int i = 0; i < nb_normals; ++i)
  {
    const E_Float* norm = normals.col(i);
    if (norm[0] == NUGA::FLOAT_MAX) continue;
    const E_Float* p = coord.col(i);
    //E_Float Lmax = ::sqrt(L(1, i));
    //E_Float Lmin = ::sqrt(L(0, i));
    //E_Float Lref = FACTOR*Lmean;//std::min(FACTOR*Lmin, Lmax);
    NUGA::sum<3>(1., p, heightv[i], norm, coord.col(nid)); // creating image point
    img[i] = nid++;
  }
  
  //ng2D.flag_external_pgs(INITIAL_SKIN);
  
#ifdef DEBUG_NGON_T
  ngon_unit bulks;
#endif

  // c. Q4 bulkheads
  E_Int nb_bulks = ng2D.PGs.size(); // for each edge, a bulk
  E_Int bulk[4];
  for (E_Int i = 0; i < nb_bulks; ++i)
  {
    const E_Int* nodes = ng2D.PGs.get_facets_ptr(i);

    bulk[0]=nodes[0];
    bulk[1]=nodes[1];
    
#ifdef DEBUG_NGON_T
    assert ( (img[nodes[0]-1] != IDX_NONE) && (img[nodes[1]-1] != IDX_NONE) );
#endif
    
    bulk[2] = img[nodes[1]-1] + 1;
    bulk[3] = img[nodes[0]-1] + 1;
    
    wNG.PGs.add(4, &bulk[0]); // this PG id is sync with ng2D (due to previous shift)
    
    //wNG.PGs._type.push_back(ng2D.PGs._type[i]);
    
#ifdef DEBUG_NGON_T
    bulks.add(4, &bulk[0]);   // this PG id is sync with ng2D (due to previous shift)
#endif
  }
  
  wNG.PGs._type.resize(wNG.PGs.size(), 0);
  wNG.PGs._ancEs.resize(2, wNG.PGs.size(), IDX_NONE);
  
  // d. TOPs
  ngon_unit tops = ghost_pgs;
  tops.change_indices(img);
  shft = wNG.PGs.size(); // for top indirection when building ghost cells at the end
  
  tops._type.resize(ghost_pgs.size(), INITIAL_SKIN);
  tops._ancEs.resize(2, ghost_pgs.size(), IDX_NONE);
  
  if (topPGlist)
  {
    E_Int PGb = wNG.PGs.size();    
    for (E_Int k=0; k < tops.size(); ++k)
      topPGlist->push_back(PGb+k);
  }
  
  wNG.PGs.append(tops);
#ifdef DEBUG_NGON_T
  assert(wNG.PGs.attributes_are_consistent());
#endif
  // e. ghost PHs
  E_Int nb_cells = ghost_pgs.size();
  Vector_t<E_Int> buff;
  //
  for (E_Int i = 0 ; i < nb_cells; ++i)
  {
    const E_Int* pgs = ng2D.PHs.get_facets_ptr(i);
    E_Int nb_pgs = ng2D.PHs.stride(i);
    
    buff.clear();
    buff.insert(buff.end(), pgs, pgs + nb_pgs); // adding bulkheads
    buff.push_back(obids[i] + 1);               // adding bottom
    buff.push_back(i + shft + 1);               // adding top
    
    wNG.PHs.add(buff.size()/*nb_pgs + 2*/, &buff[0]);
  }
  //std::cout << "size before adding ghost : " << wNG.PHs._type.size() << std::endl;
  wNG.PHs._type.resize(wNG.PHs.size(), INITIAL_SKIN);
  wNG.PHs._ancEs.resize(2, wNG.PHs.size(), IDX_NONE);
  //std::cout << "nb of added ghost : " << nb_cells << std::endl;
  
#ifdef DEBUG_NGON_T
  {
    ngon_t just_tops(tops, false);
    K_FLD::IntArray cnto;
    just_tops.export_to_array(cnto);
    medith::write("tops.plt", coord, cnto, "NGON");
  }
  
  {
    ngon_t just_bulks(bulks, false);
    just_bulks.PGs.updateFacets();
    just_bulks.PHs.updateFacets();
    K_FLD::IntArray cnto;
    just_bulks.export_to_array(cnto);
    medith::write("bulks.plt", coord, cnto, "NGON");
  }
  
  {
    K_FLD::IntArray cnto;
    wNG.export_to_array(cnto);
    medith::write("ghost.plt", coord, cnto, "NGON");
  }
  
#endif

  wNG.PGs.updateFacets();
  wNG.PHs.updateFacets();

  return 0;
}

static E_Int extrude_revol_faces
(K_FLD::FloatArray& coord, ngon_t& wNG, const Vector_t<E_Int>& PGlist, const E_Float* ax, const E_Float* ax_pt,
  E_Float angle, Vector_t<E_Int>* topPGlist = 0)
{
  if (PGlist.empty()) return 0;

  //E_Float axis[] = {0., 0., 1.};
  E_Float* axis_pt = nullptr;
  
#ifdef DEBUG_NGON_T
  E_Int minid = *std::min_element(PGlist.begin(), PGlist.end());
  E_Int maxid = *std::max_element(PGlist.begin(), PGlist.end());
  E_Int nb_pgs = wNG.PGs.size();
  assert (minid >= 0);
  assert (maxid < nb_pgs);

  std::cout << "input angle : " << height_factor << std::endl;
#endif

  coord.pushBack(ax_pt, ax_pt+3); // to embark it in the transfo

  K_FLD::FloatArray P(3,3), iP(3,3);
  NUGA::computeAFrame(ax, P);
  iP = P;
  K_FLD::FloatArray::inverse3(iP);
  NUGA::transform(coord, iP);// Now we are in the reference cylindrical coordinate system.

  axis_pt = coord.col(coord.cols()-1);
  
  // 1.

  ngon_unit ghost_pgs;
  Vector_t<E_Int> obids; //initial bottom ids  : for indirection when building ghost cells at the end
  wNG.PGs.extract(PGlist, ghost_pgs, obids);

  E_Float x0 = axis_pt[0];
  E_Float y0 = axis_pt[1];
  std::vector<E_Float> angles, radius;
  E_Int nb_new_points = NUGA::MeshTool::computeNodeRadiusAndAngles(coord, ghost_pgs, x0, y0, radius, angles);

  // a. conversion NGON3D -> NGON2D
  ngon_t ng2D;
  ngon_t::export_surfacic_FN_to_EFN(ghost_pgs, ng2D);

  // sync with wNG for future bulks appending
  E_Int shft = wNG.PGs.size();
  ng2D.PHs.shift(shft);

  // b. image points
  E_Int nb_points = coord.cols();
  coord.resize(3, nb_points + nb_new_points);
  Vector_t<E_Int> img(coord.cols(), IDX_NONE);
  E_Int nid(nb_points);
  
  for (E_Int i = 0; i < nb_points; ++i)
  {
    E_Float& a = angles[i];
    if (a == NUGA::FLOAT_MAX) continue; // not a point to be rotated
    a += angle;

    const E_Float* p = coord.col(i);

    E_Float* newp = coord.col(nid);

    newp[0] = radius[i] * ::cos(a);
    newp[1] = radius[i] * ::sin(a);
    newp[2] = p[2];

    img[i] = nid++;
  }

  NUGA::transform(coord, P); // back to original ref frame  
  
#ifdef DEBUG_NGON_T
  ngon_unit bulks;
#endif

  // c. Q4 bulkheads
  E_Int nb_bulks = ng2D.PGs.size(); // for each edge, a bulk
  std::cout << "nb bulks : " << nb_bulks << std::endl;
  E_Int bulk[4];
  for (E_Int i = 0; i < nb_bulks; ++i)
  {
    const E_Int* nodes = ng2D.PGs.get_facets_ptr(i);
    E_Int nnodes = ng2D.PGs.stride(i);

    if (nnodes != 2)
    {

    //std::cout << "bulk : " << i << "nnode : " << nnodes << std::endl;
    //std::cout << "nodes[0] : " << nodes[0] << std::endl;
    //std::cout << "nodes[1] : " << nodes[1] << std::endl;
    return 1;
  }

    bulk[0]=nodes[0];
    bulk[1]=nodes[1];
    
#ifdef DEBUG_NGON_T
    assert ( (img[nodes[0]-1] != IDX_NONE) && (img[nodes[1]-1] != IDX_NONE) );
#endif
    
    bulk[2] = img[nodes[1]-1] + 1;
    bulk[3] = img[nodes[0]-1] + 1;
    
    wNG.PGs.add(4, &bulk[0]); // this PG id is sync with ng2D (due to previous shift)
    
    //wNG.PGs._type.push_back(ng2D.PGs._type[i]);
    
#ifdef DEBUG_NGON_T
    bulks.add(4, &bulk[0]);   // this PG id is sync with ng2D (due to previous shift)
#endif
  }
  
  wNG.PGs._type.resize(wNG.PGs.size(), 0);
  wNG.PGs._ancEs.resize(2, wNG.PGs.size(), IDX_NONE);
  
  // d. TOPs
  ngon_unit tops = ghost_pgs;
  tops.change_indices(img);
  shft = wNG.PGs.size(); // for top indirection when building ghost cells at the end
  
  tops._type.resize(ghost_pgs.size(), INITIAL_SKIN);
  tops._ancEs.resize(2, ghost_pgs.size(), IDX_NONE);
  
  if (topPGlist)
  {
    E_Int PGb = wNG.PGs.size();    
    for (E_Int k=0; k < tops.size(); ++k)
      topPGlist->push_back(PGb+k);
  }
  
  wNG.PGs.append(tops);
#ifdef DEBUG_NGON_T
  assert(wNG.PGs.attributes_are_consistent());
#endif

  // e. PHs
  E_Int nb_cells = ghost_pgs.size();
  Vector_t<E_Int> buff;

  //
  for (E_Int i = 0 ; i < nb_cells; ++i)
  {
    const E_Int* pgs = ng2D.PHs.get_facets_ptr(i);
    E_Int nb_pgs = ng2D.PHs.stride(i);
    
    buff.clear();
    buff.insert(buff.end(), pgs, pgs + nb_pgs); // adding bulkheads
    buff.push_back(obids[i] + 1);               // adding bottom
    buff.push_back(i + shft + 1);               // adding top
    
    wNG.PHs.add(buff.size()/*nb_pgs + 2*/, &buff[0]);
  }

  //std::cout << "size before adding ghost : " << wNG.PHs._type.size() << std::endl;
  wNG.PHs._type.resize(wNG.PHs.size(), INITIAL_SKIN);
  wNG.PHs._ancEs.resize(2, wNG.PHs.size(), IDX_NONE);
  //std::cout << "nb of added ghost : " << nb_cells << std::endl;
  
#ifdef DEBUG_NGON_T
  {
    ngon_t just_tops(tops, false);
    K_FLD::IntArray cnto;
    just_tops.export_to_array(cnto);
    medith::write("tops.plt", coord, cnto, "NGON");
  }
  
  {
    ngon_t just_bulks(bulks, false);
    just_bulks.PGs.updateFacets();
    just_bulks.PHs.updateFacets();
    K_FLD::IntArray cnto;
    just_bulks.export_to_array(cnto);
    medith::write("bulks.plt", coord, cnto, "NGON");
  }
  
  {
    K_FLD::IntArray cnto;
    wNG.export_to_array(cnto);
    medith::write("ghost.plt", coord, cnto, "NGON");
  }
  
#endif

  wNG.PGs.updateFacets();
  wNG.PHs.updateFacets();

  return 0;
}

///
static E_Int add_flat_ghosts(ngon_t& wNG, const Vector_t<E_Int>& PGlist, bool create_tops)
{
  
  if (PGlist.empty()) return 0;
  
#ifdef DEBUG_NGON_T
  E_Int minid = *std::min_element(PGlist.begin(), PGlist.end());
  E_Int maxid = *std::max_element(PGlist.begin(), PGlist.end());
  E_Int nb_pgs = wNG.PGs.size();
  assert (minid >= 0);
  assert (maxid < nb_pgs);
#endif
  
  // 1. extract ghost bottom PGs

  ngon_unit ghost_pgs;
  Vector_t<E_Int> obids; //initial bottom ids  : for indirection when building ghost cells at the end
  wNG.PGs.extract(PGlist, ghost_pgs, obids);

  // 2. add flat ghost cells without requiring a clean_conenctivity
  
  // a. conversion NGON3D -> NGON2D
  ngon_t ng2D;
  ngon_t::export_surfacic_FN_to_EFN(ghost_pgs, ng2D);
  // sync with wNG for future bulks appending
  E_Int shft = wNG.PGs.size();
  ng2D.PHs.shift(shft);
  
#ifdef DEBUG_NGON_T
  ngon_unit bulks;
#endif

  // c. Q4 bulkheads
  E_Int nb_bulks = ng2D.PGs.size(); // for each edge, a bulk
  E_Int bulk[4];
  for (E_Int i = 0; i < nb_bulks; ++i)
  {
    const E_Int* nodes = ng2D.PGs.get_facets_ptr(i);

    bulk[0] = nodes[0];
    bulk[1] = nodes[1];
    bulk[2] = nodes[1];
    bulk[3] = nodes[0];
    
    wNG.PGs.add(4, &bulk[0]); // this PG id is sync with ng2D (due to previous shift)
    
    //wNG.PGs._type.push_back(ng2D.PGs._type[i]);
    
#ifdef DEBUG_NGON_T
    bulks.add(4, &bulk[0]);   // this PG id is sync with ng2D (due to previous shift)
#endif
  }
  
  E_Int top_shift=wNG.PGs.size();
  if (create_tops)
  {
    for (size_t i=0; i < PGlist.size(); ++i)
    {
      const E_Int* faces = wNG.PGs.get_facets_ptr(PGlist[i]);
      E_Int str = wNG.PGs.stride(PGlist[i]);
      wNG.PGs.add(str, faces);
    }
  }
  
  wNG.PGs._type.resize(wNG.PGs.size(), 99/*PG_GHOST*/);
  wNG.PGs._ancEs.resize(2, wNG.PGs.size(), IDX_NONE);
  
#ifdef DEBUG_NGON_T
  assert(wNG.PGs.attributes_are_consistent());
#endif
  // e. ghost PHs
  E_Int nb_cells = ghost_pgs.size();
  Vector_t<E_Int> buff;
  //
  for (E_Int i = 0 ; i < nb_cells; ++i)
  {
    const E_Int* pgs = ng2D.PHs.get_facets_ptr(i);
    E_Int nb_pgs = ng2D.PHs.stride(i);
    
    buff.clear();
    buff.insert(buff.end(), pgs, pgs + nb_pgs); // adding bulkheads
    buff.push_back(obids[i] + 1);               // adding bottom
    if (!create_tops)
      buff.push_back(obids[i] + 1);               // adding "top"
    else
      buff.push_back(i+top_shift+1);
    
    wNG.PHs.add(buff.size()/*nb_pgs + 2*/, &buff[0]);
  }
  //std::cout << "size before adding ghost : " << wNG.PHs._type.size() << std::endl;
  wNG.PHs._type.resize(wNG.PHs.size(), 99/*PH GHOST*/);
  wNG.PHs._ancEs.resize(2, wNG.PHs.size(), IDX_NONE);
  //std::cout << "nb of added ghost : " << nb_cells << std::endl;
  
  wNG.PGs.updateFacets();
  wNG.PHs.updateFacets();
  
#ifdef DEBUG_NGON_T  
  {
    ngon_t just_bulks(bulks, false);
    just_bulks.PGs.updateFacets();
    just_bulks.PHs.updateFacets();
    K_FLD::IntArray cnto;
    just_bulks.export_to_array(cnto);
    medith::write("bulks.plt", coord, cnto, "NGON");
  }
  
  {
    K_FLD::IntArray cnto;
    wNG.export_to_array(cnto);
    medith::write("ghost.plt", coord, cnto, "NGON");
  }
#endif

  return 0;
}

///
static void get_pgs_with_non_manifold_edges(const ngon_unit& PGs, std::set<E_Int>& pgsid)
{
  pgsid.clear();
  
  E_Int nb_pgs = PGs.size();
  
  //get non-manifold egdes
  std::set<K_MESH::NO_Edge> edges;
  {
    std::map<K_MESH::NO_Edge, E_Int> edge_to_fcount;
    K_MESH::NO_Edge noE;
    //
    for (E_Int i=0; i < nb_pgs; ++i)
    {
      const E_Int* nodes = PGs.get_facets_ptr(i);
      E_Int nb_nodes = PGs.stride(i);

      for (E_Int n=0; n < nb_nodes; ++n)
      {
        noE.setNodes(*(nodes+n), *(nodes+(n+1)%nb_nodes));

        auto it = edge_to_fcount.find(noE);
        if (it == edge_to_fcount.end())
          edge_to_fcount[noE] = 1;
        else ++it->second;
      }
    }

    for (auto &i : edge_to_fcount)
      if (i.second > 2)edges.insert(i.first);   
  }
  
  // get pg ids having these edges
  K_MESH::NO_Edge noE;
  for (E_Int i=0; i < nb_pgs; ++i)
  {
    const E_Int* nodes = PGs.get_facets_ptr(i);
    E_Int nb_nodes = PGs.stride(i);

    for (E_Int n=0; n < nb_nodes; ++n)
    {
      noE.setNodes(*(nodes+n), *(nodes+(n+1)%nb_nodes));
      if (edges.find(noE) != edges.end())
      {
        pgsid.insert(i+1); break;
      }
    }
  }
  
}

///
static E_Int discard_holes_by_box(const K_FLD::FloatArray& coord, ngon_t& wNG)
{
  assert(wNG.PHs.size() < wNG.PGs.size()); //i.e. not one ph per pg beacuse we need at least a closed PH

  // 1. Get the skin PGs
  Vector_t<E_Int> oids;
  ngon_unit pg_ext;
  wNG.PGs.extract_of_type(INITIAL_SKIN, pg_ext, oids);
  if (pg_ext.size() == 0)
    return 0; //should be only the case where One of the mesh contains completely the other

  // 2. Build the neighbourhood for skin PGs
  ngon_unit neighbors;
  K_MESH::Polygon::build_pg_neighborhood(pg_ext, neighbors);

  // Color to get connex parts
  E_Int nb_connex = 1;
  Vector_t<E_Int> colors;
  if (wNG.PHs.size() > 1)
  {
    NUGA::EltAlgo<K_MESH::Polygon>::coloring(neighbors, colors);
    nb_connex = 1 + *std::max_element(colors.begin(), colors.end());
  }

  if (nb_connex == 1) //no holes
    return 0;

  // Find out which is the external (the one with the biggest bbox)
  Vector_t<K_SEARCH::BBox3D> bbox_per_color(nb_connex);
  //init boxes
  for (E_Int i = 0; i < nb_connex; ++i)
  {
    bbox_per_color[i].maxB[0] = bbox_per_color[i].maxB[1] = bbox_per_color[i].maxB[2] = -NUGA::FLOAT_MAX;
    bbox_per_color[i].minB[0] = bbox_per_color[i].minB[1] = bbox_per_color[i].minB[2] = NUGA::FLOAT_MAX;
  }

  Vector_t<E_Int> nodes;

  E_Int nb_pgex = colors.size();
  K_SEARCH::BBox3D box;
  for (E_Int i = 0; i < nb_pgex; ++i)
  {
    const E_Int& s = pg_ext.stride(i);
    const E_Int* pN = pg_ext.get_facets_ptr(i);

    nodes.clear();
    for (E_Int j = 0; j < s; ++j, ++pN)
      nodes.push_back((*pN) - 1);//indices convention : start at 1 

    box.compute(coord, nodes);

    K_SEARCH::BBox3D& b = bbox_per_color[colors[i]];

    for (size_t j = 0; j < 3; ++j)
    {
      b.minB[j] = (b.minB[j] > box.minB[j]) ? box.minB[j] : b.minB[j];
      b.maxB[j] = (b.maxB[j] < box.maxB[j]) ? box.maxB[j] : b.maxB[j];
    }
  }
  // And the hole color are..
  std::set<E_Int> hole_colors;
  for (E_Int i = 0; i < nb_connex; ++i)
  {
    for (E_Int j = i + 1; j < nb_connex; ++j)
    {
      if (K_SEARCH::BbTree3D::box1IsIncludedinbox2(&bbox_per_color[i], &bbox_per_color[j], EPSILON))
        hole_colors.insert(i);
      else if (K_SEARCH::BbTree3D::box1IsIncludedinbox2(&bbox_per_color[j], &bbox_per_color[i], EPSILON))
        hole_colors.insert(j);
    }
  }

  // Reset flag for holes
  // Hole-cell bug fix : rather than creating an extra color for reseted PGs (INNER doesn't work in case of a single layer solid), we use CONNEXION_SKIN :
  // it prevents to take them into account in the workingPGs AND allow to discard the big parasite corresponfding to holes
  for (E_Int i = 0; i < nb_pgex; ++i)
  {
    if (hole_colors.find(colors[i]) != hole_colors.end())
      wNG.PGs._type[oids[i]] = CONNEXION_SKIN;
  }
  wNG.flag_external_phs(INITIAL_SKIN);//update consistently PHs flags

  return 0;
}

template <typename TriangulatorType>
static int validate_moves_by_fluxes
(
  std::vector<E_Int>& nids,
  const K_FLD::FloatArray& crd,
  ngon_t& ngio,
  const ngon_unit& neighborsi,
  const std::vector<E_Int>& PHlist
)
{
  bool has_moves = false;
  // identify moves as neg vals, then validate by setting them positive
  for (size_t n = 0; n < nids.size(); ++n)
  {
    if (nids[n] != E_Int(n))
    {
      nids[n] = -nids[n] - 1;
      has_moves=true;
    }
  }

  if (!has_moves) return 0;

  ngon_unit orient;
  //E_Int err = 
  build_orientation_ngu<TriangulatorType>(crd, ngio, orient);//fixme hpc : should be deduced from the input PH orientation

  std::vector<double> fluxes0(ngio.PHs.size(), NUGA::FLOAT_MAX);
  std::vector<double> vols0  (ngio.PHs.size(), NUGA::FLOAT_MAX);

  //computes initial flux at cells
  for (size_t i = 0; i < PHlist.size(); ++i)
  {
    E_Int PHi = PHlist[i];
    K_MESH::Polyhedron<0> PH0(ngio, PHi);
    double fluxvec[3];
    PH0.flux(crd, orient.get_facets_ptr(PHi), fluxvec);

    E_Float f = ::sqrt(NUGA::sqrNorm<3>(fluxvec));
    E_Float s = PH0.surface(crd);
    f /= s;

    fluxes0[PHi] = f;
  }


  Vector_t<bool> keep, wprocessed/*externalized to not reallocate it each time*/;
  Vector_t<E_Int> to_remove;
  Vector_t<E_Int> shellPHs, boundPGs, pgnids, nidsshell;
  ngon_t ngshell;
  K_CONNECT::IdTool::init_inc(nidsshell, crd.cols());
  TriangulatorType dt;

  for (size_t i = 0; i < PHlist.size(); ++i) //fixme : shell might contain more than one bad ph => process by group of bads ?
  {
    E_Int PHi = PHlist[i];
    //std::cout << "PHi : " << PHi << std::endl;
    

    // apply all node moves for PHi (fixme : shell might contain more than one bad ph => process by group of bads ?)
    const E_Int* faces = ngio.PHs.get_facets_ptr(PHi);
    int nfaces = ngio.PHs.stride(PHi);

    has_moves = false;

    for (int f = 0; f < nfaces; ++f)
    {
      E_Int Fi = faces[f] - 1;
      const E_Int* nodes = ngio.PGs.get_facets_ptr(Fi);
      int nnodes = ngio.PGs.stride(Fi);

      for (E_Int n = 0; n < nnodes; ++n)
      {
        E_Int Ni = nodes[n] - 1;
        if (nids[Ni] >= 0) nidsshell[Ni] = nids[Ni];
        else {
          nidsshell[Ni] = -nids[Ni] - 1;
          has_moves = true;
        }
      }
    }

    if (!has_moves) continue; //static shell

    //std::cout << "compute shell for " << PHlist[i] << std::endl;
    shellPHs.clear();
    boundPGs.clear();
    keep.clear();
    keep.resize(ngio.PHs.size(), false);

    ph_shell(ngio, PHi, neighborsi, shellPHs, boundPGs, wprocessed);

    double maxflux = fluxes0[PHi];
    double minvol  = NUGA::FLOAT_MAX;

    for (size_t u = 0; u < shellPHs.size(); ++u)
    {
      E_Int PHn = shellPHs[u];
      maxflux = std::max(maxflux, fluxes0[PHn]);

      auto & v = vols0[PHn];
      if (v == NUGA::FLOAT_MAX)  //compite it
      {
        K_MESH::Polyhedron<0> PH0(ngio, PHn);
        PH0.volume<TriangulatorType>(crd, orient.get_facets_ptr(PHn), v, dt);
      }

      minvol = std::min(minvol, v);
    }

    // extract shell
    keep[PHi] = true;
    for (size_t u = 0; u < shellPHs.size(); ++u)keep[shellPHs[u]] = true;

    select_phs(ngio, keep, pgnids, ngshell);

    ngshell.PGs.change_indices(nidsshell);
    clean_connectivity(ngshell, crd, 3/*ngon_dim*/, 0./*tolerance*/, false/*remove_dup_phs*/, false/*do_omp*/);

    ngon_unit orientshell;
    build_orientation_ngu<TriangulatorType>(crd, ngshell, orientshell);

    // compute new max flux : must decrease to validate
    double newmaxflux = -1.;
    double newminvol = NUGA::FLOAT_MAX;
    TriangulatorType dt;
    for (E_Int u = 0; u < ngshell.PHs.size(); ++u)
    {
      K_MESH::Polyhedron<0> PH0(ngshell, u);
      double fluxvec[3];
      PH0.flux(crd, orientshell.get_facets_ptr(u), fluxvec);

      E_Float f = ::sqrt(NUGA::sqrNorm<3>(fluxvec));
      E_Float s = PH0.surface(crd);
      f /= s;

      newmaxflux = std::max(newmaxflux, f);
      double v;
      PH0.volume<TriangulatorType>(crd, orientshell.get_facets_ptr(u), v, dt);
      if (v < newminvol) newminvol = v;
    }

    if ( (newmaxflux <= maxflux) || (newmaxflux < ZERO_M && newminvol > minvol) )
    {
      // improvement => validate moves
      for (E_Int f = 0; f < nfaces; ++f)
      {
        E_Int Fi = faces[f] - 1;
        const E_Int* nodes = ngio.PGs.get_facets_ptr(Fi);
        int nnodes = ngio.PGs.stride(Fi);

        for (E_Int n = 0; n < nnodes; ++n)
        {
          E_Int Ni = nodes[n] - 1;
          if (nids[Ni] >= 0) continue;
          nids[Ni] = -nids[Ni] - 1;
        }
      }
    }
    else
    {
      //std::cout << "current collapse strategy not good for PH : " << PHi << " => skipped" << std::endl;
    }
  }

  // discard non-valids
  for (size_t n = 0; n < nids.size(); ++n) if (nids[n] < 0)nids[n] = n;

  E_Int nb_valid_moves = 0;
  for (size_t n = 0; n < nids.size(); ++n)
  {
    if (nids[n] != E_Int(n))
    {
      ++nb_valid_moves;
    }
  }
  return nb_valid_moves;
}
    
  ngon_unit PGs;
  ngon_unit PHs;
  
};

#endif	/* NGON_T_HXX */

