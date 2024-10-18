/*    
    Copyright 2013-2024 Onera.

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

#ifndef __TRI_CONFORMIZER_CXX__
#define __TRI_CONFORMIZER_CXX__

//#define DEBUG_EXTRACT
//#define DEBUG_LIGHT

#include "Nuga/include/TRI_Conformizer.h"
#include "Nuga/include/merge.h"
#include "Nuga/include/ContourSplitter.h"
#include "Nuga/include/GeomAlgo.h"
#include "Nuga/include/SwapperT3.h"
#include "Nuga/include/Triangulator.h"
#include "Nuga/include/defs.h"

#if (defined DEBUG_TRI_CONFORMIZER) /*|| (defined DEBUG_EXTRACT)*/
#include "IO/io.h"
#endif
#ifdef DEBUG_TRI_CONFORMIZER
#include <iostream>
#include <iomanip>
#include "TRI_debug.h"
#include "chrono.h"
static bool xtest=false;
#endif
#ifdef DEBUG_EXTRACT
#include <fstream>
#include <sstream>
#endif
#if defined(DEBUG_MESHER) || defined(DEBUG_LIGHT)
#include "Nuga/include/medit.hxx"
#endif

#ifdef FLAG_STEP
#include "Nuga/include/chrono.h"
#endif

#define OLD_STYLE

namespace NUGA
{

#define IS_STRICTLY_IN_SEG(val, tol) ((val > tol) && (val < (1. - tol)))
#define IS_IN_SEG(val, tol) ((val > -tol) && (val < (1. + tol)))

///
template<short DIM>
TRI_Conformizer<DIM>::TRI_Conformizer(bool wnh) : Conformizer<DIM, K_MESH::Triangle>(wnh)
{
  _P.resize(3,3);
  _iP.resize(3,3);
  parent_type::_tolerance = 0.;
  parent_type::_iter = 0;
}

///
template<short DIM>
void TRI_Conformizer<DIM>::__set_tolerances(E_Float Lmin, E_Float Lmax, E_Float  user_tolerance)
{
  // min edge length (MEL)
  // We the previous code, passing inputol=EPSILON was in fact like passing 0.5e-3*mel;
  
  //  the following are chosen :
  // - max allowed is 5% of MEL
  E_Float tol_max = std::max(MAX_PERCENT*Lmin, ZERO_MACHINE);
  // - min allowed is ZERO MACHINE
  E_Float tol_min = std::max(ZERO_MACHINE*Lmax, ZERO_MACHINE);

#ifdef DEBUG_TRI_CONFORMIZER
  if (tol_min > tol_max)
  {
    std::cout << "WARNING : ABSOLUTE TOLERANCE CAN NOT BE SET PROPERLY GLOBALLY FOR THIS CASE :" << std::endl;
    std::cout << "Lmin/Lmax : " << Lmin << "/" << Lmax << std::endl;
    std::cout << "RATIO MAX : " << MAX_PERCENT << std::endl;
  }
#endif
  
  if (user_tolerance < 0.) //input relative value
    user_tolerance = - user_tolerance * Lmin;
  
  if (user_tolerance > 0.)
  {
    parent_type::_tolerance = std::min(user_tolerance, tol_max);
    parent_type::_tolerance = std::max(parent_type::_tolerance, tol_min);
  }
  else // == 0 : defaulted to maximum
    parent_type::_tolerance = tol_max;
  
  parent_type::_tol_x = tol_min;
  parent_type::_tol_clean = 1.e-9;
}

///
template <short DIM>
E_Int
TRI_Conformizer<DIM>::__split_Elements
(const K_FLD::FloatArray& pos, K_FLD::IntArray & connect,
 NUGA::bool_vector_type& xc,
 NUGA::int_vector_type& ancestors)
{
  K_FLD::IntArray::const_iterator pS;
  K_FLD::IntArray& connectIn(connect), connectOut;
  NUGA::int_vector_type ancOut;
  NUGA::bool_vector_type xcOut;

  NUGA::int_vector_type& colors = parent_type::_colors; //size as new connect
  NUGA::int_vector_type& xr = parent_type::_xr; //sized as old connect+1(one pass end) : give the beginning of the split for each ancestor
  E_Int Ci = 0;
  colors.clear();
  xr.resize(connect.cols()+1, 0);

  DELAUNAY::MeshData data;
  DELAUNAY::MesherMode mode;
  mode.mesh_mode = mode.TRIANGULATION_MODE;
  mode.remove_holes = false;  // option to speed up.
  mode.ignore_coincident_nodes = true;
  
  DELAUNAY::T3Mesher<E_Float> mesher(mode);
  
  NUGA::int_vector_type gnids, lnids;
  if (mode.ignore_coincident_nodes) K_CONNECT::IdTool::init_inc(gnids, pos.cols());
  
#ifdef DEBUG_TRI_CONFORMIZER
  int err_count(0);
#endif
#ifdef DEBUG_EXTRACT
  std::ostringstream xfname;
  xfname << /*MIO::wdir*/"./" << "_extract_boxes.txt" ;
  std::ofstream extract_file (xfname.str().c_str());
#endif

#ifdef FLAG_STEP
  NUGA::chrono c;
#endif
  
#ifdef DEBUG_TRI_CONFORMIZER
  drawT3(pos, connect, zTi, true);
#endif

  K_FLD::IntArray ci, ci2;
  K_FLD::FloatArray pi;
  std::set<K_MESH::Edge> sci;
  std::vector<E_Int> revIDs;

  NUGA::int_vector_type hnodes;
  E_Int x[2], err;
  E_Int ret;
  // for swapping
  K_MESH::NO_Edge E;
  std::set<K_MESH::NO_Edge> origE;
  std::vector<std::pair<E_Int, E_Int> > swapE;

  for (E_Int i = 0; i < connectIn.cols(); ++i)
  {
    pS = connectIn.col(i);
    xr[i] = connectOut.cols();
    T3& tri = parent_type::_elements[i];

#ifdef DEBUG_TRI_CONFORMIZER
//    if (ancestors[i] == zTi)
//      drawT3(pos, connect, i);
#endif
    
#ifdef FLAG_STEP
    c.start();  
#endif

#ifdef OLD_STYLE
    ret = __get_connectB2(pos, connectIn, tri, _edges, sci, hnodes);
#else
    ret = __get_mesh_data(pos, connectIn, tri, _edges, pi, ci2, revIDs);
#endif

#ifdef FLAG_STEP
    tcon += c.elapsed();
    c.start();  
#endif

#ifdef OLD_STYLE
    if ((ret == 2) || ((ret == 0) && (sci.size() == 3) && hnodes.size() == 3)) // new nodes/edges are matching exactly the initial triangle.
#else
    if (ret == 2)
#endif
    {
      connectOut.pushBack(pS, pS+3);
      ancOut.push_back(ancestors[i]);
      xcOut.push_back(false);
      colors.push_back(Ci++);
      
#ifdef DEBUG_CONFORMIZER
      NUGA::ConformizerRoot::split_fastdiscard_counter++;
#endif

#ifdef FLAG_STEP
      tnot += c.elapsed();
#endif
    }
    else if (ret == 1) // degenerated
    {
#ifdef DEBUG_TRI_CONFORMIZER
      //drawT3(pos, connectIn, i);
      NUGA::ConformizerRoot::degen_counter++;    
//      {
//          K_FLD::FloatArray coord1(pos);
//          K_FLD::IntArray tmp;
//          tmp.pushBack(connectIn.col(i), connectIn.col(i)+3);
//          std::vector<E_Int> newIDs;
//          DELAUNAY::MeshTool::compact_to_mesh(coord1, tmp, newIDs);   
//          std::ostringstream o;
//          o << "TRI_degen_compact_" <<parent_type::_iter << "_" << i << ".mesh";
//          medith::write(o.str().c_str(), coord1, tmp, "TRI");
//        }
#endif
      //return ret;
    }
    else // need a split.
    {
#ifdef DEBUG_TRI_CONFORMIZER
      
      NUGA::ConformizerRoot::split_counter++;
      
      //if (parent_type::_iter == 10)
        //drawT3(pos, connect, tri.Si);
#endif

#ifdef OLD_STYLE
#ifdef FLAG_STEP
      c.start();  
#endif
      ci.clear();

      // connectivity.
      for (std::set<K_MESH::Edge>::iterator it = sci.begin(); it != sci.end(); ++it)
      {
        x[0] = (*it).node(0);
        x[1] = (*it).node(1);
        if (x[0] == x[1]) continue;
        ci.pushBack(x, x+2);
      }
      
#ifdef DEBUG_TRI_CONFORMIZER
    //if (ancestors[i]==401)//parent_type::_iter == 10 || xtest
      //medith::write("conto.mesh", pos, ci, "BAR");
#endif

      // coordinates
      // compact
      if (hnodes.empty())
        __compact_to_mesh(pos, ci, pi, ci2, revIDs);
      else
        __compact_to_mesh(pos, ci, pi, ci2, revIDs, &hnodes);
      
#ifdef DEBUG_TRI_CONFORMIZER
    //if (ancestors[i]==401) //parent_type::_iter == 10 || xtest
    //{
    //  medith::write("contoci.mesh", pi, ci2, "BAR");
    //  std::cout << "hard nodes ? " << hnodes.size() << std::endl;
    //}
#endif

#ifdef FLAG_STEP
      tcomp += c.elapsed();
      c.start();
#endif

      // transform (in the coord. sys of the triangle)
      __transform(pos.col(*pS), pos.col(*(pS+1)), pos.col(*(pS+2)), pi);
      pi.resize(2, pi.cols());

#ifdef FLAG_STEP
      ttra += c.elapsed();
      c.start();
#endif

#endif
#ifdef DEBUG_MESHER
      //if (ancestors[i]==1167)//parent_type::_iter == 10 || xtest
        //medith::write("contot.mesh", pi, ci2, "BAR");
#endif
      // Now mesh
      data.clear(); //reset containers
      data.pos = &pi;
      data.connectB = &ci2;
#ifdef OLD_STYLE
      data.hardNodes = hnodes;
#endif
      
#ifdef DEBUG_TRI_CONFORMIZER
#ifdef DEBUG_MESHER
      mesher.dbg_flag=false;
      //if (/*parent_type::_iter == 10 || xtest*/ancestors[i]==Ti) DBG_MESHER_FLAG = 1;
      //if (i==59259) mesher.dbg_flag=true;
      /*if (i==59259 || i == 59281)
        {
        std::ostringstream o;
        o << "contot_" << i << ".mesh";
        medith::write(o.str().c_str(), pi, ci2, "toto");
        //continue;
      }*/
#endif
#endif
      // iterative for robustness : try shuffling the input data, then try without forcing edges (hasardous)
      mesher.seed_random(3*i);
      err = __iterative_run (mesher, pi, ci2, hnodes, data, lnids, false/*i.e. try to force all edge*/, true/*i.e silent also last it*/);
      if (err != 0)
      {
        //try again with a normalized contour
        K_SEARCH::BBox2D box;
        box.compute(pi);
        double dX = box.maxB[0] - box.minB[0];
        double dY = box.maxB[1] - box.minB[1];

        /*{
          K_FLD::FloatArray tmp(pi);
          tmp.resize(3, tmp.cols(), 0.);
          medith::write("contour_init.mesh", tmp, ci2, "BAR");
        }*/

        for (int u = 0; u < pi.cols(); ++u)
        {
          pi(0, u) = (pi(0, u) - box.minB[0]) / dX;
          pi(1, u) = (pi(1, u) - box.minB[1]) / dY;
        }

        /*{
          K_FLD::FloatArray tmp(pi);
          tmp.resize(3, tmp.cols(), 0.);
          medith::write("contour_norma.mesh", tmp, ci2, "BAR");
        }*/

        err = __iterative_run(mesher, pi, ci2, hnodes, data, lnids, false/*i.e. try to force all edge*/, true/*i.e silent also last it*/);
        //if (err == 0)
        //  std::cout << " found " << std::endl;
      }
      if (err != 0)
        err = __iterative_run(mesher, pi, ci2, hnodes, data, lnids, true/*i.e. ignore unforceable edges*/, mode.silent_errors/*i.e output last it error eventually*/);
      
      if (err != 0)
      {
#ifdef DEBUG_LIGHT
        std::ostringstream o;
        o << "TRIConformizer_err_" << i << ".mesh";
        medith::write(o.str().c_str(), pi, ci2, "BAR");
        o.str("");
        K_FLD::FloatArray coord1(pos);
        std::vector<E_Int> newIDs;
        NUGA::MeshTool::compact_to_mesh(coord1, ci2, newIDs);
        o << "TRIConformizer_err3D_" << i << ".mesh";
        medith::write(o.str().c_str(), coord1, ci2, "BAR");
#endif
        
#ifdef DEBUG_TRI_CONFORMIZER
        drawT3(pos, connectIn, i, true);
#endif
        {
#if (defined DEBUG_TRI_CONFORMIZER) || (defined DEBUG_EXTRACT)
          K_FLD::FloatArray coord1(pos);
          std::vector<E_Int> newIDs;
          NUGA::MeshTool::compact_to_mesh(coord1, ci, newIDs);
#endif
#ifdef DEBUG_TRI_CONFORMIZER
          {
            std::ostringstream o;
            o << "TRI_mesher_err_compact_" <<parent_type::_iter << "_" <<  i << ".mesh";
            medith::write(o.str().c_str(), coord1, ci, "BAR");
          }
#endif
#ifdef DEBUG_EXTRACT
          K_SEARCH::BBox3D box;
          box.compute(K_FLD::ArrayAccessor<K_FLD::FloatArray >(coord1));

          extract_file << "boxes.push_back(K_SEARCH::BBox3D());" << std::endl;
          extract_file << "boxes[boxes.size()-1].minB[0]=" << box.minB[0] << ";"<< std::endl;
          extract_file << "boxes[boxes.size()-1].minB[1]=" << box.minB[1] << ";"<< std::endl;
          extract_file << "boxes[boxes.size()-1].minB[2]=" << box.minB[2] << ";"<< std::endl;
          extract_file << "boxes[boxes.size()-1].maxB[0]=" << box.maxB[0] << ";"<< std::endl;
          extract_file << "boxes[boxes.size()-1].maxB[1]=" << box.maxB[1] << ";"<< std::endl;
          extract_file << "boxes[boxes.size()-1].maxB[2]=" << box.maxB[2] << ";"<< std::endl << std::endl;
#endif
        }
#ifdef DEBUG_TRI_CONFORMIZER
        medith::write("contot.mesh", pi, ci2, "BAR");
        ++err_count;
#endif

#ifdef DEBUG_EXTRACT
        std::cout << "extract box for " << i << std::endl;
        continue;
#else
        return err;
#endif  
      }

#ifdef FLAG_STEP
      trun += c.elapsed();
#endif

#ifdef DEBUG_TRI_CONFORMIZER
      if (i == zTi)
      {
        K_FLD::FloatArray c(*data.pos);
        c.resize(3, data.pos->cols(), 0.);
        medith::write("mesh0.mesh", c, data.connectM, "TRI");
      }
#endif   
   
      // Do some swapping to avoid degen T3s
      __improve_triangulation_quality(tri, revIDs, swapE, origE, data);
      
#ifdef DEBUG_TRI_CONFORMIZER
      if (i == zTi)
      {
        K_FLD::FloatArray c(*data.pos);
        c.resize(3, data.pos->cols(), 0.);
        medith::write("meshSwap.mesh", c, data.connectM, "TRI");
      }
#endif  
     
      // Get back to initial ids
      K_FLD::IntArray::changeIndices(data.connectM, revIDs);
      
      if (mode.ignore_coincident_nodes && data.unsync_nodes)
      {
        for (size_t n = 0; n < lnids.size(); ++n)
        {
          if ((size_t)lnids[n] != n) gnids[revIDs[n]] = gnids[revIDs[lnids[n]]]; //we used gnids in the rhs to propagate the moves
        }
      }

      // Append upon exit.
      connectOut.pushBack(data.connectM);
      ancOut.resize(ancOut.size() + data.connectM.cols(), ancestors[i]);
      xcOut.resize(xcOut.size() + data.connectM.cols(), true);
      for (E_Int j = 0; j < data.connectM.cols(); ++j)
        colors.push_back(data.colors[j] + Ci);
      Ci += 1 + *std::max_element(data.colors.begin(), data.colors.end());
      
//#ifdef DEBUG_CONFORMIZER
//      if (detect_duplis_and_baffles(connectOut))
//      {
//        std::cout << "oddities : after Mesher : FOUND !!!! for element " << tri.Si << std::endl;
//        return 1;
//      }
//#endif
      
    }
  }
  
  xr[connect.cols()] = connectOut.cols();
  
  if (mode.ignore_coincident_nodes)
    parent_type::__clean(gnids, connectOut, ancOut, &xcOut);

  connect = connectOut;
  ancestors = ancOut;
  xc = xcOut;
  
#ifdef DEBUG_EXTRACT
  extract_file.close();
#endif
  
#ifdef DEBUG_TRI_CONFORMIZER
#ifdef FLAG_STEP
  //if (NUGA::chrono::verbose > 1)
  {
#ifdef DEBUG_CONFORMIZER
    std::cout << "split : nb of quicly discarded : " << NUGA::ConformizerRoot::split_fastdiscard_counter << std::endl;
    std::cout << "split : nb of done : " << NUGA::ConformizerRoot::split_counter << std::endl;
    std::cout << "split : nb of degen discarded : " << NUGA::ConformizerRoot::degen_counter << std::endl;

    std::cout << "split : nb of resulting bits in more : " << connect.cols() - connectIn.cols() << std::endl;
    std::cout << "split : NB ERRORS : " << err_count << std::endl;
#endif
    std::cout << "TIMES : "<< std::endl;
    std::cout << "false X : " << tnot << std::endl;
    std::cout << "prepare contour : " << tcon << std::endl;
//    std::cout << "pepare contour (form 0) : " << tcon0 << std::endl;
//    std::cout << "pepare contour (form 1) : " << tcon1<< std::endl;
//    std::cout << "pepare contour (form 2) : " << tcon2 << std::endl;
//    std::cout << "pepare contour (compact) : " << tcomp << std::endl;
//    std::cout << "pepare contour (transform) : " << ttra << std::endl;
#ifdef E_TIME
    typedef DELAUNAY::Mesher<E_Float, DELAUNAY::VarMetric<E_Float> > mesher_t;
    typedef DELAUNAY::Kernel<E_Float> kernel_t;
    std::cout << "mesher : init                           : " << mesher_t::tinit << std::endl;
    std::cout << "mesher : tria : cav time                : " << kernel_t::cavity_time << std::endl;
    std::cout << "mesher : tria :     init cav time       : " << kernel_t::init_cavity_time << std::endl;
    std::cout << "mesher : tria :             base time   : " << kernel_t::base_time << std::endl;
    std::cout << "mesher : tria :             append time : " << kernel_t::append_time << std::endl;
    std::cout << "mesher : tria :     fix cav time        : " << kernel_t::fix_cavity_time << std::endl;
    std::cout << "mesher : tria :     sort bound time     : " << kernel_t::sorting_bound_time << std::endl;
    std::cout << "mesher : tria : remesh cavity           : " << kernel_t::remesh_time << std::endl;
    std::cout << "mesher : tria : invalidate prev elts    : " << kernel_t::inval_time << std::endl;
    std::cout << "mesher : tria : TOTAL                   : " << mesher_t::ttria << std::endl;
    std::cout << "mesher : restore boundaries             : " << mesher_t::tbound << std::endl;
    std::cout << "mesher : set colors                     : " << mesher_t::tcolors << std::endl;
    std::cout << "mesher : finalize                       : " << mesher_t::tfinalize << std::endl;
    std::cout << "mesher : TOTAL                          : " << trun << std::endl;
#endif
  }
#endif
#endif
  
  return 0;
}

///
template <short DIM>
E_Int
TRI_Conformizer<DIM>::__iterative_run
(DELAUNAY::T3Mesher<E_Float>& mesher, K_FLD::FloatArray& crd, K_FLD::IntArray& cB, 
 NUGA::int_vector_type& hnodes, DELAUNAY::MeshData& data, std::vector<E_Int>& lnids,
 bool ignore_unforceable_edge, bool silent_last_iter)
{
  E_Int err(0), nb_pts(crd.cols());
  
  mesher.mode.silent_errors = silent_last_iter;
  mesher.mode.ignore_unforceable_edges = ignore_unforceable_edge;
    
  lnids.clear();
  if (mesher.mode.ignore_coincident_nodes) //and eventually unforceable edge
    K_CONNECT::IdTool::init_inc(lnids, nb_pts);
 
  // Multiple attempts : do not shuffle first, and then do shuffle
  E_Int railing = -1;
  E_Int nb_attemps = std::min(cB.cols(), (E_Int)6);//to test all RTOL2 values in range [1.e-8, 1.e-3]
  E_Int k = 9;

  while (railing++ < nb_attemps)
  {
    //if (err) std::cout << "atempt " << railing << " return error : " << err << std::endl;
    
    k = std::max((E_Int)3, k-1);
    E_Float RTOL2 = ::pow(10., -k);
    
    data.clear(); //reset containers
    mesher.clear();
    
    mesher.mode.do_not_shuffle=(railing == 0); // i.e. shuffle every time but the first

#ifdef DEBUG_TRI_CONFORMIZER
    if (!silent_last_iter && (railing == nb_attemps))mesher.mode.silent_errors = false; // only listen the last iter
#endif
    
    data.pos = &crd;
    data.connectB = &cB;
#ifdef OLD_STYLE
    data.hardNodes = hnodes;
#endif
    
    err = mesher.run(data);

    crd.resize(2, nb_pts);//discard box nodes
    
    // process errors and update cB, hNodes and lnids
    if (mesher.mode.ignore_unforceable_edges)
    {
      if (!mesher._edge_errors.empty())//has errors
      {
        // get (lambda,d2) for each xedge node regarding current faulty hard edge
        for (size_t i=0; i < mesher._edge_errors.size(); ++i)
        {
          DELAUNAY::edge_error_t& edge_err = mesher._edge_errors[i];
          E_Int& Ni = edge_err.Ni;
          E_Int& Nj = edge_err.Nj;
          std::set<E_Int>& nodes = edge_err.nodes;
          
          E_Float L2 = NUGA::sqrDistance(crd.col(Ni), crd.col(Nj), 2);
//        
          std::vector<std::pair<E_Float, E_Int>> sNodes;
          
          for (auto& N : nodes)
          {
            if (N >= crd.cols()) continue; // box nodes
            
            E_Float dNNi2 = NUGA::sqrDistance(crd.col(N), crd.col(Ni), 2);
            E_Float dNNj2 = NUGA::sqrDistance(crd.col(N), crd.col(Nj), 2);
            //std::cout <<data.hnids.size() << std::endl;
            if (dNNi2 < RTOL2*L2)
            {
              data.hnids[MAX(N, Ni)] = MIN(N,Ni);  //lower id priorization
            }
            else if (dNNj2 < RTOL2*L2)
            {
              data.hnids[MAX(N, Nj)] = MIN(N,Nj);  //lower id priorization
            }
            else // do we split NiNj ?
            {
              E_Float lambda;
              E_Float h2 = K_MESH::Edge::linePointMinDistance2<DIM>(crd.col(Ni),crd.col(Nj), crd.col(N), lambda);
              if (h2 >= RTOL2*L2 || lambda < 0. || lambda > 1.) continue;
              sNodes.push_back(std::make_pair(lambda, N));
            }
          }
          
          size_t sz = sNodes.size();
          if (sz > 1) std::sort(ALL(sNodes));
          
          if (sz)
          {
            data.hardEdges.erase(K_MESH::NO_Edge(Ni,Nj)); //remove edge to refine
            data.hardEdges.insert(K_MESH::NO_Edge(Ni,sNodes[0].second)); //first bit
            data.hardEdges.insert(K_MESH::NO_Edge(sNodes[sz-1].second, Nj)); //last bit
            //in-between bits
            for (size_t i=0; i < sz-1; ++i)
              data.hardEdges.insert(K_MESH::NO_Edge(sNodes[i].second,sNodes[i+1].second));
          }
        }
        
        data.sync_hards();//reflect hnids update to hard edges and hard node
        
        // now rebuild cB from current data
        cB.clear();
        for (auto& he : data.hardEdges)
          cB.pushBack(he.begin(), he.end());

#ifdef OLD_STYLE
        hnodes.clear();
        for (size_t i=0; i < data.hardNodes.size(); ++i)
        {
          hnodes.push_back(data.hardNodes[i]);
        }
#endif
        
        //
        err = 1;
      }
    }
    if (mesher.mode.ignore_coincident_nodes)
    {
      for (E_Int i = 0; i < nb_pts; ++i)
      {
        if ((data.hnids[i] != i) && (lnids[i] == i))
        {
          lnids[i] = data.hnids[i];
        }
      }
    }
    
    if (err == 0) return 0;
  }
  
  
  return err;
}

///
template <short DIM>
void
TRI_Conformizer<DIM>::__improve_triangulation_quality(const T3& tri, const std::vector<E_Int>& revIDs,
                                                      std::vector<std::pair<E_Int, E_Int> >& wpair_set,
                                                      std::set<K_MESH::NO_Edge>& wedge_set,
                                                      DELAUNAY::MeshData& data)
{
  // swap ? if original edge are present in the triangulation whereas they got splitting points
  wpair_set.clear();
  wedge_set.clear();
  
  for (size_t t=0; t < 3; ++t)
  {
    size_t sz11 = _edges[tri.edges[t]].size() ;
    if (sz11> 2) // has splitting points
      wedge_set.insert(K_MESH::NO_Edge(_edges[tri.edges[t]][0], _edges[tri.edges[t]][sz11-1])); //insert the original edge
  }
  if (!wedge_set.empty())
  {
    K_FLD::IntArray::const_iterator pK;
    K_MESH::NO_Edge E;
    std::set<K_MESH::NO_Edge>::const_iterator itE;
    //
    for (E_Int k=0; k < data.connectM.cols(); ++k)
    {
      pK = data.connectM.col(k);
      for (size_t n=0; n < 3; ++n)
      {
        E.setNodes(revIDs[*(pK+n)], revIDs[*(pK+(n+1)%3)]);
        itE = wedge_set.find(E);
        if (itE != wedge_set.end())
        {
          wpair_set.push_back(std::make_pair(k,(n+2)%3));
          wedge_set.erase(itE);
        }
      } 
    }
  }
  
  if (!wpair_set.empty())
    NUGA::EltAlgo<K_MESH::Triangle>::fast_swap_edges(wpair_set, data.connectM, data.neighbors);
  
  // Also swap poor quality triangles (ITERATIVE)
  E_Int railing = data.connectM.cols();
  while (--railing)
  {
    wpair_set.clear();
      
    // also swap poor quality triangles
    NUGA::GeomAlgo<K_MESH::Triangle>::get_swapE(*data.pos, data.connectM, data.neighbors, data.hardEdges, parent_type::_tolerance, wpair_set); //warning : quality<2> doesnt give the same result as quality<3>. need to convert to tolerance criterium.

    if (!wpair_set.empty())
    {
      NUGA::EltAlgo<K_MESH::Triangle>::fast_swap_edges(wpair_set, data.connectM, data.neighbors);
#ifdef DEBUG_TRI_CONFORMIZER
      //std::ostringstream o;
      //K_FLD::FloatArray c(*data.pos);
      //c.resize(3, data.pos->cols(), 0.);
      //o << "swap_" << data.connectM.cols() - railing << ".mesh";
      //medith::write(o.str().c_str(), c, data.connectM, "TRI");
#endif
    }
    else
      break;
  }
}

///
template <short DIM>
E_Int
TRI_Conformizer<DIM>::__intersect
(K_FLD::FloatArray& pos, const K_FLD::IntArray & connect, T3& t1, T3& t2, E_Float tol)
{
  E_Int  i, r(0);
  E_Int ret(0);
  
//  if (t1.Si == 19494 || t2.Si == 19494)
//  {
//      std::cout << "caught" << std::endl;
//  }
  
//  //////////////////////////////////////////////////////////////////////////////////
//  const E_Float* P1 = pos.col(connect(0,t1.Si));
//  const E_Float* Q1 = pos.col(connect(1,t1.Si)); 
//  const E_Float* R1 = pos.col(connect(2,t1.Si));
//  const E_Float* P2 = pos.col(connect(0,t2.Si));
//  const E_Float* Q2 = pos.col(connect(1,t2.Si));
//  const E_Float* R2 = pos.col(connect(2,t2.Si));
//  
//  if (!K_MESH::Triangle::fast_intersectT3<3>(P1, Q1, R1, P2, Q2, R2))
//  {
//#ifdef DEBUG_TRI_CONFORMIZER
//    ++NUGA::ConformizerRoot::fastdiscard_counter;
//#endif
//    return false;
//  }
//      
//  //////////////////////////////////////////////////////////////////////////////////
  
  // Fast return : Triangles far from each others  
  if (__fast_discard(pos, connect.col(t2.id), connect.col(t1.id), tol))
  {
#ifdef DEBUG_TRI_CONFORMIZER
    ++NUGA::ConformizerRoot::fastdiscard_counter;
#endif
    return 0;
  }
    
  if (__fast_discard(pos, connect.col(t1.id), connect.col(t2.id), tol))
  {
#ifdef DEBUG_TRI_CONFORMIZER
    ++NUGA::ConformizerRoot::fastdiscard_counter;
#endif
    return 0;
  }
  
  _wnodes_vec.clear();

  E_Int pidx[2];
  //bool insideT1(false), insideT2(false);
  E_Bool coplanarE;
  E_Int nb_coplanarE(0);
  bool are_coplanar=false;
  
// T2->T1
  for (i = 0; i < 3; ++i)
  {
    r = __intersect(pos, connect, t1, _edges, t2.edges[i], tol, _wnodes_vec, pidx, coplanarE);
    nb_coplanarE = (r==2) ? nb_coplanarE+1 : nb_coplanarE;
    //insideT1 |= (r==3);
    ret |=r;
  }
  
  if (nb_coplanarE > 1)
  {
    are_coplanar=true; // overlap might require a second pass as each triangle split is done separtely so triangulation can be different
  }
  nb_coplanarE=0;
  for (i = 0; i < 3; ++i)
  {
    r = __intersect(pos, connect, t2, _edges, t1.edges[i], tol, _wnodes_vec, pidx, coplanarE);
    nb_coplanarE = (r==2) ? nb_coplanarE+1 : nb_coplanarE;
    ret |=r;
    //insideT2 |= (r==3);
  }
  if (nb_coplanarE > 1)
  {
    are_coplanar=true;; // overlap might require a second pass as each triangle split is done separtely so triangulation can be different
  }
  
  if (are_coplanar)
    ret=2;

  size_t sz = _wnodes_vec.size();
  
  if (sz==0 )
    return ret;
  
  if (sz >= 2)
  {
    _wnodes_set.clear();
    _wnodes_vecclean.clear();
    for (size_t i=0; i< sz; ++i)
        if (_wnodes_set.insert(_wnodes_vec[i]).second)
          _wnodes_vecclean.push_back(_wnodes_vec[i]);
    _wnodes_vec=_wnodes_vecclean;
    sz = _wnodes_vec.size();
  }

  if (sz == 2)
  {
    E_Int E[2];
    E[0] = _wnodes_vec[0];
    E[1] = _wnodes_vec[1];
#ifdef DEBUG_TRI_CONFORMIZER
    if (t1.id == zTi || t2.id == zTi)
      std::cout << "caught" << std::endl;
#endif    
    t1.Xedges.pushBack(E, E+2);
    t2.Xedges.pushBack(E, E+2);
  }
  else if (sz > 2)
  {
#ifdef DEBUG_TRI_CONFORMIZER
    /*if (sz > 2)
    {
      std::vector<T3> elts;
      elts.push_back(t1);
      elts.push_back(t2);
      std::ostringstream o;
      o << "i_" <<t1.Si << "_" << t2.Si << ".mesh";
      this->drawElements(o.str().c_str(), "fmt_mesh", pos, connect, elts);
      drawT3(pos, connect, t1.Si, true);
      drawT3(pos, connect, t2.Si, true);
    }*/
#endif     
    parent_type::_sorterFI.clear();
    parent_type::_sorterFI.reserve(sz);
    E_Float P0P1[DIM], P0Pi[DIM], d;
    const E_Float* P0 = pos.col(_wnodes_vec[0]);
    const E_Float* P1 = pos.col(_wnodes_vec[1]);
    NUGA::diff<DIM>(P1, P0, P0P1);

    parent_type::_sorterFI.push_back(std::make_pair(0., _wnodes_vec[0]));
    parent_type::_sorterFI.push_back(std::make_pair(1., _wnodes_vec[1]));

    for (size_t i = 2; i < sz; ++i)
    {
      const E_Float* Pi = pos.col(_wnodes_vec[i]);
      NUGA::diff<DIM>(Pi, P0, P0Pi);
      d = NUGA::dot<DIM>(P0Pi, P0P1);
      parent_type::_sorterFI.push_back(std::make_pair(d, _wnodes_vec[i]));
    }

    sz = parent_type::_sorterFI.size();
    if (sz > 2) std::sort(parent_type::_sorterFI.begin(), parent_type::_sorterFI.end());

    E_Int E[2];

#ifdef DEBUG_TRI_CONFORMIZER
    /*if (sz > 2)
    {
      std::vector<T3> elts;
      elts.push_back(t1);
      elts.push_back(t2);
      std::ostringstream o;
      o << "i_" <<t1.Si << "_" << t2.Si << ".mesh";
      this->drawElements(o.str().c_str(), "fmt_mesh", pos, connect, elts);
      drawT3(pos, connect, t1.Si, false);
      drawT3(pos, connect, t2.Si, false);
    }*/
#endif     
    for (size_t i = 0; i < sz-1; ++i)
    {
      E[0] = parent_type::_sorterFI[i].second;
      E[1] = parent_type::_sorterFI[i+1].second;
#ifdef DEBUG_TRI_CONFORMIZER
      if (t1.id == zTi || t2.id == zTi)
        std::cout << "caught" << std::endl;
#endif
      
      t1.Xedges.pushBack(E, E+2);
      t2.Xedges.pushBack(E, E+2);
    }
  }
  else if (sz == 1) // ???
  {
    E_Int E[2];
    E[0] = _wnodes_vec[0];
    E[1] = _wnodes_vec[0];
#ifdef DEBUG_TRI_CONFORMIZER
    if (t1.id == zTi || t2.id == zTi)
      std::cout << "caught" << std::endl;
#endif    
    t1.Xedges.pushBack(E, E+2);
    t2.Xedges.pushBack(E, E+2);
#ifdef DEBUG_TRI_CONFORMIZER
    /*std::vector<T3> elts;
    elts.push_back(t1);
    elts.push_back(t2);
    std::ostringstream o;
    o << "X_" << t1.Si << "_" << t2.Si << ".mesh";
    this->drawElements(o.str().c_str(), "fmt_mesh", pos, connect, elts);*/
#endif
    //return false;
  }

  return ret;
}

inline bool getBoundary(const E_Int* t0, const E_Int* t1, E_Int& i, E_Int& j)
{
  for (i=0; i < 3; ++i)
  {
    const E_Int& t00=*(t0+i);
    const E_Int& t01=*(t0+(i+1)%3);
      
    for (j=0; j < 3; ++j)
    {
      const E_Int& t10=*(t1+j);
      const E_Int& t11=*(t1+(j+1)%3);
      
      if ((t00==t10) && (t01==t11))
        return true;
      if ((t00==t11) && (t01==t10))
        return true;
    }
  }
  
  i=j=IDX_NONE;
  return false;
}
    

///
template <short DIM>
bool
TRI_Conformizer<DIM>::__fast_discard
(const K_FLD::FloatArray& pos, const E_Int* T0, const E_Int* T1, E_Float tol)
{
  
  const E_Float* P0 = pos.col(T1[0]);
  const E_Float* P1 = pos.col(T1[1]);
  const E_Float* P2 = pos.col(T1[2]);
  const E_Float* Q0 = pos.col(T0[0]);
  const E_Float* Q1 = pos.col(T0[1]);
  const E_Float* Q2 = pos.col(T0[2]);
  
  // 1. Check if points are on the same side of the t's plane
  
  NUGA::diff<DIM>(P1, P0, _U1);
  NUGA::diff<DIM>(P2, P0, _U2);
  NUGA::crossProduct<DIM>(_U1,_U2,_U3);
  NUGA::normalize<DIM>(_U3);
        
  bool is_far[] = {false,false, false};
  E_Float h0(0.), h1(0.), h2(0.);
  NUGA::diff<DIM>(Q0, P0, _U1);
  h0 = NUGA::dot<DIM>(_U3, _U1);
  is_far[0] = (h0 >= tol) || (h0 <= -tol);
    
  NUGA::diff<DIM>(Q1, P0, _U1);
  h1 = NUGA::dot<DIM>(_U3, _U1);
  is_far[1] = (h1 >= tol) || (h1 <= -tol);
    
  NUGA::diff<DIM>(Q2, P0, _U1);
  h2 = NUGA::dot<DIM>(_U3, _U1);
  is_far[2] = (h2 >= tol) || (h2 <= -tol);
    
  E_Int s[3];
  s[0]=SIGN(h0);
  s[1]=SIGN(h1);
  s[2]=SIGN(h2);
    
  if (is_far[0] && is_far[1] && is_far[2] && (s[0] == s[1]) && (s[0]==s[2]))
    return true;
  
  bool overlapping = (!is_far[0] && !is_far[1] && !is_far[2]);
  E_Int shn=IDX_NONE;//shared node's rank
  //E_Int shn1=IDX_NONE;
  
       if (T0[0] == T1[0]) {shn=0; /*shn1=0;*/}
  else if (T0[0] == T1[1]) {shn=0; /*shn1=1;*/}
  else if (T0[0] == T1[2]) {shn=0; /*shn1=2;*/}
  else if (T0[1] == T1[0]) {shn=1; /*shn1=0;*/}
  else if (T0[1] == T1[1]) {shn=1; /*shn1=1;*/}
  else if (T0[1] == T1[2]) {shn=1; /*shn1=2;*/}
  else if (T0[2] == T1[0]) {shn=2; /*shn1=0;*/}
  else if (T0[2] == T1[1]) {shn=2; /*shn1=1;*/}
  else if (T0[2] == T1[2]) {shn=2; /*shn1=2;*/}
    
  // 2. Check if they are not overlapping but they share a point and the 2 others are far and in the same side => discard
  
  if (!overlapping) //nor coplanar
  {
    if (shn != IDX_NONE)
    {
      if ( (s[(shn+1)%3] == s[(shn+2)%3]) && (is_far[(shn+1)%3] && is_far[(shn+2)%3]) )
        return true;  
    
      // 3. Check if they are not overlapping but they share an edge => discard
      E_Int i,j;    
      if (getBoundary(T0,T1, i,j)) //true if an edge is shared and therefore i,j are valued
        return true;
      
      // 4. Check if they share a node S, another one is on the plane P => check if SP intersects T1
// WARNING : COMMENTED BECAUSE NOT CONSISTENT WITH REAL ALGO : some X are missed here whereas they exist and the real algo do the job correctly
//      const E_Float* Pplane=0;
//      E_Float tmp[3];
//      if (!is_far[(shn+1)%3]) // one node shared and one node on the plane
//        Pplane=pos.col(T0[(shn+1)%3]);
//      else if (!is_far[(shn+2)%3]) // one node shared and one node on the plane
//        Pplane=pos.col(T0[(shn+2)%3]);
//      else if (s[(shn+1)%3] == -s[(shn+2)%3])
//      {
//        // Compute X between Line E0E1 and Plane (Shared, n)== (T0[shn], _U3)
//        const E_Float* E0=pos.col(T0[(shn+1)%3]);
//        const E_Float* E1=pos.col(T0[(shn+2)%3]);
//        NUGA::diff<DIM>(E1, E0, _U1); // Line
//        NUGA::diff<DIM>(pos.col(T0[shn]), E0, _U2); // vector : shared node to one edge vertex
//        
//        E_Float s=NUGA::dot<DIM>(_U3, _U1);
//        if (SIGN(s)==0)
//          return false;
//        s=1./s;
//        s*=NUGA::dot<DIM>(_U3, _U2);
//        
//        NUGA::sum<3>(s, _U1, E0, tmp);
//        Pplane=&tmp[0];
//      }
//      
//      if (Pplane == 0)
//        return false;
//      
//      NUGA::diff<DIM>(pos.col(T1[(shn1+1)%3]), pos.col(T1[shn1]), _U1);
//      NUGA::diff<DIM>(Pplane, pos.col(T1[shn1]), _U2);
//      NUGA::crossProduct<DIM>(_U1,_U2,_U3);
//      
//      if (SIGN(_U3[0]) == -1)
//        return true;
//      
//      NUGA::diff<DIM>(pos.col(T1[(shn1+2)%3]), pos.col(T1[shn1]), _U1);
//      NUGA::crossProduct<DIM>(_U2,_U1,_U3);
//      
//      if (SIGN(_U3[0]) == -1)
//        return true;  
    }
  }
  else //overlapping or just coplanar : COMMENTED FOR NGON BOOLEAN DOINT ONE PASS ONLY (very few case like that occur at first iter)
  {
    //if sharing an edge and opposite nodes are on each side => just coplanar
    if (shn != IDX_NONE) //they must share at least a node
    {
      E_Int i,j;    
      if (getBoundary(T0,T1, i,j)) //true if an edge is shared and therefore i,j are valued
      {
        E_Float shE[3], U[3], n0[3], n1[3];
        E_Int n0op = T0[(i+2)%3];
        E_Int n1op = T1[(j+2)%3];
        
        NUGA::diff<3>(pos.col(T0[(i+1)%3]), pos.col(T0[i]), shE);
        
        NUGA::diff<3>(pos.col(n0op), pos.col(T0[i]), U);
        NUGA::crossProduct<3>(shE,U,n0);
        
        NUGA::diff<3>(pos.col(n1op), pos.col(T0[i]), U);
        NUGA::crossProduct<3>(shE,U,n1);
        
        if (NUGA::dot<3>(n0, n1) < 0.)
          return true;
      }
    }   
  }
  
  return false;
}

///
template <short DIM>
E_Bool
TRI_Conformizer<DIM>::__intersect
(K_FLD::FloatArray& pos, const K_FLD::IntArray & connect, T3& t, edge_container_type& edges, E_Int idxE,
 E_Float tol, std::vector<E_Int>& nodes, E_Int* pidx, E_Bool& coplanar)
{
  std::vector<E_Int>& e = edges[idxE];
  E_Float u0[2], E[DIM], eps(/*100.*EPSILON*/tol), IP[DIM];
  E_Bool  overlap, intersect;
  E_Int e0(e[0]), e1(e[1]), Ni, k, tx[2];  
  K_FLD::IntArray::const_iterator pS = connect.col(t.id);

  if (*pS == e0 && *(pS+1) == e1)
    return false;
  if (*pS == e1 && *(pS+1) == e0)
    return false;

  if (*(pS+1) == e0 && *(pS+2) == e1)
    return false;
  if (*(pS+1) == e1 && *(pS+2) == e0)
    return false;

  if (*pS == e0 && *(pS+2) == e1)
    return false;
  if (*pS == e1 && *(pS+2) == e0)
     return false;
  
  intersect = K_MESH::Triangle::intersect<3>
    (pos, *pS, *(pS+1), *(pS+2), e[0], e[1], tol, true, u0[0], u0[1], tx, overlap, coplanar);
  
  if (!intersect)
  {
    // these 3 if are added to not miss some trace with an existing summit as an end
    if (tx[0] == 3 && ((*(pS+1)==e0)|| (*(pS+1)==e1))) //false X : they just share anode
      nodes.push_back(*(pS+1));
    if (tx[0] == 5 && ((*pS==e0)|| (*pS==e1))) //false X : they just share anode
      nodes.push_back(*pS);
    if (tx[0] == 6 && ((*(pS+2)==e0)|| (*(pS+2)==e1))) //false X : they just share anode
      nodes.push_back(*(pS+2));
    return false;
  }
    

  bool one_single_x_point = (u0[1] == NUGA::FLOAT_MAX) || (::fabs(u0[1] - u0[0]) < eps);
  bool share_a_node = (*pS == e0 || *(pS+1) == e0 || *(pS+2) == e0) || (*pS == e1 || *(pS+1) == e1 || *(pS+2) == e1);
  bool one_inside = ((u0[0] > eps) && (u0[0] < 1.-eps));
  
  if (share_a_node && one_single_x_point && !one_inside) //not a real intersection
    return false;
  
  bool add2E, add2T;
  intersect=false;

  NUGA::diff<DIM>(pos.col(e1), pos.col(e0), E);

  for (E_Int n = 0; n < 2; ++n)
  {
    if (u0[n] == NUGA::FLOAT_MAX)
      continue;
    
    add2E = (u0[n] >= eps) && (u0[n] <= (1. - eps));
    add2T = (tx[n] == 0);
    
    Ni = IDX_NONE;
    if (add2E)
    {
      Ni = pos.cols();//id of the next point to create
      e.push_back(Ni); // push back the node in the edge's node list
    }
    else if (u0[n] <eps)
      Ni = e0;
    else if (u0[n] > 1.-eps)
      Ni = e1; 
    
    assert (Ni != IDX_NONE);
    
    // add its mate to the appropropiate triangle edge.
    if (tx[n] == 1 || tx[n] == 2 || tx[n] == 4)
    {
      E_Int M0, M1;
      if (tx[n]==1)
      {
        M0 = *pS;
        M1 = *(pS+1);
      }
      else if (tx[n]==2)
      {
        M0 = *(pS+1);
        M1 = *(pS+2);
      }
      else // tx==4
      {
        M0 = *(pS+2);
        M1 = *pS;
      }
      
      for (E_Int j = 0; j < 3; ++j)
      {
        E_Int& idx = t.edges[j];
        if ((edges[idx][0] == M0 && edges[idx][1] == M1) || 
          (edges[idx][1] == M0 && edges[idx][0] == M1))
        {
          pidx[n]=j;
          edges[idx].push_back(Ni);/*add2T=true*/;break;
        }
      }
    }
    
    if (add2E) // Create new point Ni
    {
      for (k = 0; k < DIM; ++k)
        IP[k] = pos(k, e0) + u0[n]*E[k];
      pos.pushBack(IP, IP+3);
      intersect = 1;
    }
    else if (add2T)
      intersect = 3;
    
    if (!overlap)
      nodes.push_back(Ni);
  }
  
  if (overlap)
  {
    t.edges.push_back(idxE);
    //drawTandE(pos, pS, e0, e1);
    intersect=2;
  }

  return intersect;
}

///
template <short DIM>
void
TRI_Conformizer<DIM>::__update_data
(const K_FLD::FloatArray& coord, const K_FLD::IntArray& dum_connect, const std::vector<E_Int>& newIDs)
{
  size_t i, j, NB_EDGES(_edges.size()), nb_tris, nb_nodes;
  
  if (!newIDs.empty())
  {

    for (i = 0; i < NB_EDGES; ++i)
    {
      std::vector<E_Int>& E = _edges[i];
      nb_nodes = E.size();

      for (j = 0; j < nb_nodes; ++j)
        E[j] = newIDs[E[j]];
    }

    nb_tris = parent_type::_elements.size();
    for (i = 0; i < nb_tris; ++i)
    {
      T3& t = parent_type::_elements[i];
      K_FLD::IntArray::changeIndices(t.Xedges, newIDs);
    }
  }

  // Tidy up the edge nodes: sort them and remove duplicated ones.
  __tidy_edges(coord, _edges);
}

///
template <short DIM>
E_Int
TRI_Conformizer<DIM>::__get_mesh_data
(const K_FLD::FloatArray & pos,
 const K_FLD::IntArray & connect,
 const T3& t, edge_container_type& Edges,
 K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, std::vector<E_Int>& oids)
{
  size_t nb_edges = t.edges.size();
  
  cnt.clear();
  crd.clear();
  
  if ((nb_edges == 3) && (t.Xedges.cols() == 0)) // Fast return : clean triangle.
  {
    if ((Edges[t.edges[0]].size() == 2) && (Edges[t.edges[1]].size() == 2) && (Edges[t.edges[2]].size() == 2))
      return 2;
  }

#ifdef FLAG_STEP
  NUGA::chrono c;
  c.start();
#endif
  K_FLD::IntArray cnt0, cnt1;
  for (size_t e = 0; e < 3; ++e)
  {
    const std::vector<E_Int>& Ei = Edges[t.edges[e]];
    
    for (size_t n=0; n < Ei.size()-1; ++n)
    {
      cnt0.pushBack(&Ei[n], &Ei[n]+2);
    }
  }
  E_Int n0 = cnt0.cols();
  for (size_t e = 3; e < nb_edges; ++e)
  {
    const std::vector<E_Int>& Ei = Edges[t.edges[e]];
    
    for (size_t n=0; n < Ei.size()-1; ++n)
    {
      cnt0.pushBack(&Ei[n], &Ei[n]+2);
    }
  }
  
  //medith::write("cnt0.mesh", pos, cnt0, "BAR");
  
  // X edges
  K_FLD::IntArray::const_iterator pS;
  E_Int nb_x_edges = t.Xedges.cols();
  for (E_Int e = 0; e < nb_x_edges; ++e)
  {
    pS = t.Xedges.col(e);
    cnt0.pushBack(pS, pS+2);
  }
  
  std::vector<E_Int> unodes;
  cnt0.uniqueVals(unodes);
  if (unodes.size() == 3)
    return 2; //clean triangle
  
  //medith::write("cnt0x.mesh", pos, cnt0, "BAR");
  
  NUGA::MeshTool::compact_to_mesh(pos, cnt0, crd, cnt1, oids);
    
  //medith::write("cnt1.mesh", crd, cnt1, "BAR");
  
  // transform (in the coord. sys of the triangle)
  pS = connect.col(t.id);
  __transform(pos.col(*pS), pos.col(*(pS+1)), pos.col(*(pS+2)), crd);
  
  bool do_the_reject_test=false;
  std::vector<E_Int> keep;
  if (n0 < cnt0.cols())
  {
    keep.resize(crd.cols(), -1);
    for (E_Int i=0; i < n0; ++i)
    {
      keep[cnt1(0,i)]=keep[cnt1(1,i)]=1; //always keep contour points
    }
  
    for (E_Int i=n0; i < cnt1.cols(); ++i)
    {
      if (keep[cnt1(0,i)] != 1 || keep[cnt1(1,i)] != 1)
      {
        do_the_reject_test=true; break;
      }
    }
  }
  
  if (do_the_reject_test)
  {
    K_FLD::FloatArray crdT3; //hack to get the triangle summit in the local frame
    crdT3.reserve(3,3);
    crdT3.pushBack(pos.col(*pS), pos.col(*pS)+3);
    crdT3.pushBack(pos.col(*(pS+1)), pos.col(*(pS+1))+3);
    crdT3.pushBack(pos.col(*(pS+2)), pos.col(*(pS+2))+3);
  
    __transform(pos.col(*pS), pos.col(*(pS+1)), pos.col(*(pS+2)), crdT3);
  
    //std::cout << crdT3 << std::endl;
  
    E_Float h0(crdT3(2,0)), h;
    for (size_t i=0; i < crd.cols(); ++i) //reject any node not falling on the triangle's plane
    {
      if (keep[i] != -1)
        continue;
      h=::fabs(crd(2,i)-h0);
      keep[i] = (h < EPSILON) ? -1 : 0; //reject (set to 0) only if not on the plane.)
    }
  
    crd.resize(2, crd.cols()); //now in 2D
    crdT3.resize(2, crdT3.cols()); //now in 2D

    //
    for (size_t i=0; i < crd.cols(); ++i)
    {
      if (keep[i] != -1)
        continue;
      keep[i]=1;
      for(size_t n=0; (n<3) && keep[i]; ++n)
      {
        E_Float s = K_MESH::Triangle::surface<2>(crd.col(i), crdT3.col(n), crdT3.col((n+1)%3));
        keep[i] &= (s >= -EPSILON);
      }
    }
    
    // now remove any bit that has a rejected node or degenerated
    for (size_t i=0; i< cnt1.cols(); ++i)
    {
      if (keep[cnt1(0,i)]==1 && keep[cnt1(1,i)]==1)
      {
        if (cnt1(0,i) != cnt1(1,i))
          cnt.pushBack(cnt1.col(i), cnt1.col(i)+2);
      }
    }
  }
  else
  {
    crd.resize(2, crd.cols()); //now in 2D
    
    // now remove any bit that is degenerated
    for (size_t i=0; i< cnt1.cols(); ++i)
    {
      if (cnt1(0,i) != cnt1(1,i))
        cnt.pushBack(cnt1.col(i), cnt1.col(i)+2);
    }
  }
  
  if (cnt.cols() < 3) return 1;
  if (cnt.cols() == 3) return 2;
  
  cnt.uniqueVals(unodes);
  if (unodes.size() == 3)
    return 2; //clean triangle
  
#ifdef FLAG_STEP
  
#endif
  
  return 0;
}


///
template <short DIM>
E_Int
TRI_Conformizer<DIM>::__get_connectB2
(const K_FLD::FloatArray & pos,
 const K_FLD::IntArray & connect,
 const T3& t, edge_container_type& Edges,
 std::set<K_MESH::Edge>& hBO, std::vector<E_Int>& hNodes)
{
  size_t nb_edges = t.edges.size();
  
  if ((nb_edges == 3) && (t.Xedges.cols() == 0)) // Fast return : clean triangle.
  {
    if ((Edges[t.edges[0]].size() == 2) && (Edges[t.edges[1]].size() == 2) && (Edges[t.edges[2]].size() == 2))
      return 2;
  }

  hNodes.clear();
  hBO.clear();
  _wnodes_set.clear();

#ifdef FLAG_STEP
  NUGA::chrono c;
  c.start();
#endif

  // Contour edges.
  bool degenerated = __get_B0_edges(connect, t, Edges, hBO, _wnodes_set);

#ifdef FLAG_STEP
  tcon0 += c.elapsed();
  c.start();
#endif

  if (degenerated)
    return 1;

  // Overlapping edges (add only edges bits that are inside the triangle).
  __get_Inner_edges(pos, connect, t, Edges, hBO, _wnodes_set);

#ifdef FLAG_STEP
  tcon1 += c.elapsed();
  c.start();
#endif

  // Imprinted edges.
  __get_Imprint_edges(t, hBO, _wnodes_set);

#ifdef FLAG_STEP
  tcon2 += c.elapsed();
  c.start();
#endif

  for (std::set<E_Int>::const_iterator it = _wnodes_set.begin(); it != _wnodes_set.end(); ++it)
    hNodes.push_back(*it);

  return 0;
}

///
template <short DIM>
inline bool
TRI_Conformizer<DIM>::__get_B0_edges
(const K_FLD::IntArray & connect,
 const T3& t, edge_container_type& Edges,
 std::set<K_MESH::Edge>& hBO, std::set<E_Int>& Nodes0)
{
  E_Int             n, start, sz0, end, s;
  E_Bool            same_orient(false);
  size_t            e, NB_BOUND(3);
  K_MESH::Triangle  tri(connect.col(t.id));

  std::set<K_MESH::NO_Edge> no_tmp;
  
  for (e = 0; e < NB_BOUND; ++e)
  {
    std::vector<E_Int>& Ei = Edges[t.edges[e]];

    start = 0;
    sz0 = Ei.size();
    end = sz0-1;

    while ((start < end) && (Ei[start] != tri.node(0) && Ei[start] != tri.node(1) && Ei[start] != tri.node(2)))++start;
    while ((end > start) && (Ei[end] != tri.node(0) && Ei[end] != tri.node(1) && Ei[end] != tri.node(2)))--end;

    E_Int & N0 = Ei[start];
    E_Int & N1 = Ei[end];

    K_MESH::Triangle::getOrientation(tri, N0, N1, same_orient);
    //assert (err == 0);//fixme : handle properly the errors.

    if (!same_orient)
    {
      std::reverse(Ei.begin(), Ei.end());
      s = start;
      start = sz0 - end - 1;
      end = sz0 - s - 1;
    }

    for (n = start; n < end; ++n)
    {
      if (Ei[n] != Ei[n+1])
      {
        if (no_tmp.insert(K_MESH::NO_Edge(Ei[n], Ei[n+1])).second)
          hBO.insert(K_MESH::Edge(Ei[n], Ei[n+1]));
        else
        {
          hBO.erase(K_MESH::Edge(Ei[n], Ei[n+1]));
          hBO.erase(K_MESH::Edge(Ei[n+1], Ei[n]));
        }
      }
    }

    Nodes0.insert(Ei.begin() + start, Ei.begin() + end +1);
  }

  return (hBO.empty()); // degneration = empty

}

///
template <short DIM>
inline void
TRI_Conformizer<DIM>::__get_Inner_edges
(const K_FLD::FloatArray & pos,
 const K_FLD::IntArray & connect,
 const T3& t, edge_container_type& Edges,
 std::set<K_MESH::Edge>& hBO, std::set<E_Int>& Nodes0)
{
  E_Int             n, start, sz0, end;
  size_t            e, NB_BOUND(3);
  K_MESH::Triangle  tri(connect.col(t.id));
  size_t            nb_edges = t.edges.size();
  
  for (e = NB_BOUND; e < nb_edges; ++e)
  {
    const std::vector<E_Int>& Ei = Edges[t.edges[e]];

    start = 0;
    sz0 = Ei.size();
    end = sz0-1;

    if (!is_inside(connect, t, pos, Ei[start], EPSILON))
      while ((start < end) && (_wnodes_set.find(Ei[start]) == _wnodes_set.end()))++start;
    if (!is_inside(connect, t, pos, Ei[end], EPSILON))
      while ((end > start) && (_wnodes_set.find(Ei[end]) == _wnodes_set.end()))--end;

    for (n = start; n < end; ++n)
    {
      if (Ei[n] == Ei[n+1])
        continue;
      if (hBO.find(K_MESH::Edge(Ei[n+1], Ei[n])) == hBO.end())//fixme
        hBO.insert(K_MESH::Edge(Ei[n], Ei[n+1]));
    }

    Nodes0.insert(Ei.begin() + start, Ei.begin() + end +1);
  }
}

///
template <short DIM>
inline void
TRI_Conformizer<DIM>::__get_Imprint_edges
(const T3& t, std::set<K_MESH::Edge>& hBO, std::set<E_Int>& Nodes0)
{
  // X edges
  K_FLD::IntArray::const_iterator pS;
  E_Int nb_x_edges = t.Xedges.cols();
  for (E_Int e = 0; e < nb_x_edges; ++e)
  {
    pS = t.Xedges.col(e);
    Nodes0.insert(*pS);
    Nodes0.insert(*(pS+1));

    if (*pS == *(pS+1))
      continue;

    hBO.insert(K_MESH::Edge(pS));
  }
}

///
template <short DIM>
E_Bool
TRI_Conformizer<DIM>::is_inside
(const K_FLD::IntArray& connect, const T3& t, const K_FLD::FloatArray& pos, E_Int Ni, E_Float tol)
{
  E_Bool inside = true;
  K_FLD::IntArray::const_iterator pS = connect.col(t.id);
  E_Float s0 = K_MESH::Triangle::surface<DIM>(pos.col(*pS), pos.col(*(pS+1)), pos.col(*(pS+2)));
  E_Float s = 0.;
  for (E_Int n = 0; (n < K_MESH::Triangle::NB_NODES) && inside; ++n)
  {
    s += K_MESH::Triangle::surface<DIM>(pos.col(Ni), pos.col(*(pS+n)), pos.col(*(pS+(n+1)%K_MESH::Triangle::NB_NODES)));
  }
  return (::fabs(s - s0) < tol);
}

template <short DIM>
void
TRI_Conformizer<DIM>::__tidy_edges
(const K_FLD::FloatArray& coord, edge_container_type& edges)
{
  size_t  NB_EDGES(edges.size());
 
  for (size_t i = 0; i < NB_EDGES; ++i)
    __tidy_edge(coord, edges[i]);
}

template <short DIM>
void
TRI_Conformizer<DIM>::__tidy_edge (const K_FLD::FloatArray& coord, std::vector<E_Int>& nodes)
{
  std::vector<E_Int> tmp;
  tmp.reserve(nodes.size());

  for (size_t i = 0; i < nodes.size(); ++i)
    if (nodes[i] != IDX_NONE)tmp.push_back(nodes[i]);

  nodes = tmp;

  size_t nb_nodes = nodes.size();
  
  if (nb_nodes == 3)
  {
      std::swap(nodes[1], nodes[2]);
      return;
  }
  else if (nb_nodes > 3)
    NUGA::MeshTool::reorder_nodes_on_edge<std::vector<E_Int>, DIM>(coord, nodes, 0, parent_type::_sorterFI);
}

///
template <short DIM>
E_Int
TRI_Conformizer<DIM>::__simplify_and_clean
(const K_FLD::FloatArray& pos, E_Float tolerance,
 K_FLD::IntArray& connect, NUGA::int_vector_type& ancestors, NUGA::bool_vector_type& xc)
{
  // Merge with priorization : Existing common nodes (ECN) > X common nodes (XCN) > Other intersected elements nodes (OXEN)
  
  NUGA::int_vector_type nodesX, nodesXC;
  
#ifdef FLAG_STEP
  NUGA::chrono c;
  c.start();
#endif

  // Get the intersecting elements' nodes.
  {
    NUGA::int_set_type for_unicity;
    K_FLD::IntArray::const_iterator pS;
    E_Int nb_tris(connect.cols()), stride(connect.rows());
    for (E_Int Si = 0; Si < nb_tris; ++Si)
    {
      if (xc[Si])
      {
        pS = connect.col(Si);
        for (E_Int n=0; n < stride; ++n)
          if (for_unicity.insert(*(pS+n)).second) //not already in
            nodesX.push_back(*(pS+n));
      }
    }
  }

#ifdef DEBUG_TRI_CONFORMIZER
  //TRI_debug::draw_connected_to_node_T3s(pos, connectX, 1109);
  //TRI_debug::draw_connected_to_node_T3s(pos, connectX, 1110);
  /*std::vector<bool> mask(connectX.cols(), false);
  for (size_t i=0; i < mask.size(); ++i)
    if ((ancX[i]==407) || (ancX[i]==424)  || (ancX[i]==1361) || (ancX[i]==1378))
      mask[i]=true;
  medith::write("before.mesh", pos, connectX, "TRI", &mask, &ancX);*/
#endif

  
  // Get the common nodes.
  {
    K_FLD::IntArray connectXC;
    //fixme : instead of entire connect, should be only elements providing from intersection from the beginning (only elements around the contour)
    // connectX only have newly intersction, so we need to keep a xc history..
    NUGA::MeshTool::getNonManifoldEdges(connect, connectXC); // Common edges.
    connectXC.uniqueVals(nodesXC);

#ifdef DEBUG_TRI_CONFORMIZER
    medith::write("commonE0.mesh", pos, connectXC, "BAR");
#endif
  }
  
  // Priorisation : Existing common nodes > X common nodes
  std::sort(nodesXC.begin(), nodesXC.end()); // to have : 

  // Priorisation complete : adding all the nodes at the end to have the OXEN with lowest priority.
  nodesXC.insert(nodesXC.end(), nodesX.begin(), nodesX.end());
    
  NUGA::int_vector_type nids;
  K_FLD::ArrayAccessor<K_FLD::FloatArray> posAcc(pos); 
  K_FLD::ArrayAccessor<K_FLD::IntArray> connectAcc(connect);   
  E_Int nb_merges = merge(posAcc, connectAcc, tolerance, nodesXC, nodesXC, nids);
  
#ifdef DEBUG_TRI_CONFORMIZER
  medith::write("beforem.mesh", pos, connect, "TRI");
#endif
  
  parent_type::__clean(nids, connect, ancestors, &xc);

#ifdef DEBUG_TRI_CONFORMIZER
  medith::write("afterm.mesh", pos, connect, "TRI");
#endif
  
  __process_duplicates(connect, ancestors, xc); // remove baffles due to duplication by merging.

  return nb_merges;
}

#ifdef DEBUG_TRI_CONFORMIZER
template <short DIM>
void
TRI_Conformizer<DIM>::draw(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, double tol, E_Int what, K_FLD::IntArray& out)
{
  // OK=0, HAT = 1, SPIKE = 2, SMALL = 4, BADQUAL = 5
  
  out.clear();

  double MINQUAL{ 1.e-3 };
  double tol2 = tol * tol;

  for (E_Int i = 0; i < connect.cols(); ++i)
  {
    K_FLD::IntArray::const_iterator pS = connect.col(i);

    E_Int N0 = *pS;
    E_Int N1 = *(pS + 1);
    E_Int N2 = *(pS + 2);

    E_Int ni;
    SwapperT3::eDegenType type = SwapperT3::degen_type2(pos, N0, N1, N2, tol2, 0., ni);

    bool skip = true;

    if (what == 5)
    {
      E_Float q = K_MESH::Triangle::qualityG<3>(pos.col(connect(0, i)), pos.col(connect(1, i)), pos.col(connect(2, i)));
      if (q < MINQUAL)
        skip = false;
    }
    else if (what == 1 && type == SwapperT3::eDegenType::HAT)
      skip = false;
    else if (what == 2 && type == SwapperT3::eDegenType::SPIKE)
      skip = false;
    else if (what == 4 && type == SwapperT3::eDegenType::SMALL)
      skip = false;

    if (skip) continue;
 
    out.pushBack(connect.col(i), connect.col(i) + 3);
  }
}
#endif
///
template <short DIM>
E_Int
TRI_Conformizer<DIM>::__simplify_and_clean2
(const K_FLD::FloatArray& pos, E_Float tolerance,
K_FLD::IntArray& connect, NUGA::int_vector_type& ancestors, NUGA::bool_vector_type& xc)
{
  //std::cout << "FIXING SPIKES" << std::endl;
  E_Int nb_elts0 = connect.cols(), nb_elts;
  //E_Int it(-1);
  bool has_degn = false;
  do
  {
    //std::cout << "iter : " << ++it << std::endl;
    nb_elts = connect.cols();
    std::vector<E_Int> cnids;
    SwapperT3::clean(pos, tolerance, connect, ancestors, cnids);//compress inside
    
    K_CONNECT::valid pred(cnids);
    K_CONNECT::IdTool::compress(parent_type::_colors, pred); //must be compress for below xr update and validity

    // xr must keep the same size to refer to ancestors : might have empty slots for totally vanished elts
    std::vector<E_Int> xr = parent_type::_xr;
    for (size_t i = 0; i < parent_type::_xr.size()-1; ++i)
    {
      E_Int v = 0;
      for (E_Int u = parent_type::_xr[i]; u < parent_type::_xr[i + 1]; ++u)
        if (cnids[i] != IDX_NONE)++v;
      xr[i + 1] = xr[i] + v;
    }
    parent_type::_xr = xr;

    //std::cout << "xr size : " << parent_type::_xr.size() << std::endl;
    //std::cout << "xr last : " << parent_type::_xr[parent_type::_xr.size() - 1] << std::endl;

    // SwapperT3::run changes connect sorting => need to recompute xr, colors, and ancestors
    std::map<E_Int, E_Int> col_to_anc;
    for (E_Int i = 0; i < connect.cols(); ++i)
      col_to_anc[parent_type::_colors[i]] = ancestors[i];
    
    SwapperT3::run(pos, tolerance, connect, parent_type::_colors);

    xc.resize(connect.cols(), true);

    std::map<E_Int, K_FLD::IntArray> col_to_elts;
    for (E_Int i = 0; i < connect.cols(); ++i)
      col_to_elts[parent_type::_colors[i]].pushBack(connect.col(i), connect.col(i) + 3);

    ancestors.clear();
    parent_type::_colors.clear();
    connect.clear();
    E_Int nbanc = parent_type::_xr.size();
    parent_type::_xr.clear();
    parent_type::_xr.resize(nbanc, 0);
 
    for (auto it : col_to_elts)
    {
      E_Int subcol = it.first;
      K_FLD::IntArray& cnt = it.second;
      E_Int anc = col_to_anc[subcol];

      //std::cout << "sub/anc/nbe : " << subcol << "/" << anc << "/" << cnt.cols() << std::endl;

      parent_type::_xr[anc+1] += cnt.cols(); // accum all sub colors size for a given anc

      for (E_Int k = 0; k < cnt.cols(); ++k)
      {
        connect.pushBack(cnt.col(k), cnt.col(k) + 3);
        parent_type::_colors.push_back(subcol);
        ancestors.push_back(anc);
      }
    }

    // xr : now go to accumulated storage as it should be
    for (size_t i = 0; i < parent_type::_xr.size() - 1; ++i)
      parent_type::_xr[i + 1] += parent_type::_xr[i];

    has_degn = SwapperT3::has_degen(pos, connect, tolerance*tolerance);
    
    //std::cout << "has degen ?" << has_degn << std::endl;
    //std::cout << "has cleaned ?" << (nb_elts != connect.cols()) << std::endl;
  } while (has_degn  && nb_elts != connect.cols());

  return (connect.cols() - nb_elts0);
}

///
template <short DIM>
void
TRI_Conformizer<DIM>::__prepare_data
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect)
{
  // Builds the BbTree
  parent_type::__prepare_data(pos, connect);

  // Edges
  algo_type::BoundToEltType                 E_to_T;
  algo_type::BoundToEltType::const_iterator itE;

  K_FLD::ArrayAccessor<K_FLD::IntArray> actv(connect);
  algo_type::getBoundToElements(actv, E_to_T);

  _edges.clear();
  _edges.resize(E_to_T.size());

  // Create edge objects and link them to triangles.
  E_Int ne = 0, Si, nb_tris;
  for (itE = E_to_T.begin(); itE != E_to_T.end(); ++itE, ++ne)
  {
    _edges[ne].push_back(itE->first.node(0));
    _edges[ne].push_back(itE->first.node(1));

    nb_tris = itE->second.size();
    for (E_Int i = 0; i < nb_tris; ++i)
    {
      Si = itE->second[i];
      parent_type::_elements[Si].edges.push_back(ne);
    }
  }

  _edges.resize(ne);//fixme : ???
}

///
template <short DIM>
void
TRI_Conformizer<DIM>::__compact_to_mesh
(const K_FLD::FloatArray& posIn, const K_FLD::IntArray& connectIn, 
 K_FLD::FloatArray& posOut, K_FLD::IntArray& connectOut, 
 std::vector<E_Int>& revIDs, std::vector<E_Int> * hard_nodes)
{
  K_FLD::FloatArray::const_iterator pS;
  size_t                            nb_nodes;
  std::map<E_Int, E_Int>            newIDs;

  revIDs.clear();
  posOut.clear();
  connectOut.clear();
  
  std::set<E_Int> uids;
  connectIn.uniqueVals(uids);
  if (hard_nodes)
    uids.insert(hard_nodes->begin(), hard_nodes->end());

  revIDs.insert(revIDs.end(), ALL(uids));

  nb_nodes = revIDs.size();
  posOut.reserve(3, nb_nodes);  

  for (size_t i = 0; i < nb_nodes; ++i)
  {
    pS = posIn.col(revIDs[i]);
    posOut.pushBack(pS, pS+3);
    newIDs[revIDs[i]] = i;
  }

  connectOut.resize(2, connectIn.cols());
  for (E_Int i = 0; i < connectIn.cols(); ++i)
  {
    connectOut(0, i) = newIDs[connectIn(0, i)];
    connectOut(1, i) = newIDs[connectIn(1, i)];
  }
  if (hard_nodes)
  {
    for (size_t i = 0; i < hard_nodes->size(); ++i)
      (*hard_nodes)[i] = newIDs[(*hard_nodes)[i]];
  }
}

/// WARNING ; transform but do not resize
template <short DIM>
void
TRI_Conformizer<DIM>::__transform
(const E_Float* P0, const E_Float* P1, const E_Float* P2, K_FLD::FloatArray& pos)
{
  // Compute a orthonormal frame having U along P0P1 and W along the normal.
  E_Float U[3], V[3], W[3];
  K_MESH::Triangle::normal(P0, P1, P2, W);
  NUGA::diff<3>(P1, P0, U);
  NUGA::normalize<3>(U);
  NUGA::crossProduct<3>(W, U, V);

  // Build the transformation matrix.
  for (E_Int i = 0; i < 3; ++i)
  {
    _P(i, 0) = U[i];
    _P(i, 1) = V[i];
    _P(i, 2) = W[i];
  }

  // Transform the working coordinates.
  K_FLD::FloatArray::inverse3(_iP = _P);
  __transform (pos, _iP);
}

///
template <short DIM>
void
TRI_Conformizer<DIM>::__transform
(K_FLD::FloatArray& pos, const K_FLD::FloatArray& t)
{
  K_FLD::FloatArray::iterator pN;
  E_Float Q[3];
  E_Int i, j, n;
  for (i = 0; i < pos.cols(); ++i)
  {
    pN = pos.col(i);

    for (j = 0; j < 3; ++j)
    {
      Q[j] = 0.;
      for (n = 0; n < 3; ++n)
        Q[j] += t(j, n) * (*(pN+n)) ;
    }

    for (j = 0; j < 3; ++j)
      pos(j, i) = Q[j];
  }
}

///
template <short DIM>
void
TRI_Conformizer<DIM>::__process_duplicates
(K_FLD::IntArray& connect, std::vector<E_Int>& ancestors, NUGA::bool_vector_type& xc)
{
  std::vector<E_Int> dupIds;
  bool has_duplis = NUGA::MeshTool::detectDuplicated(connect, dupIds, false /*strict orient*/);
  
  if (!has_duplis)
    return;
  
  std::set<K_MESH::NO_Edge> dupE;
  K_FLD::IntArray::const_iterator pS;
  for (size_t i=0; i < dupIds.size(); ++i)
  {
    if (dupIds[i] != (E_Int)i)
    {
      //std::cout << i << "(" << ancestors[i] << ") --> " << dupIds[i] << "(" << ancestors[dupIds[i]] << ")" << std::endl;
      pS = connect.col(i);
      dupE.insert(K_MESH::NO_Edge(*pS,*(pS+1)));
      dupE.insert(K_MESH::NO_Edge(*(pS+1),*(pS+2)));
      dupE.insert(K_MESH::NO_Edge(*(pS+2),*pS));
    }
  }
  
  algo_type::BoundToEltType E_to_T;
  algo_type::BoundToEltType::const_iterator itE;
  K_FLD::ArrayAccessor<K_FLD::IntArray> actv(connect);
  NUGA::EltAlgo<K_MESH::Triangle>::getBoundToElements(actv, E_to_T);
  
  // discard duplicates with a free edge by setting to -1 the new id.
  bool free_edge;
  for (std::set<K_MESH::NO_Edge>::const_iterator itE = dupE.begin(); itE != dupE.end(); ++itE)
  {
    free_edge = true;
    const std::vector<E_Int>& T3s = E_to_T[*itE];
    for (size_t i=1; (i<T3s.size()) && free_edge; ++i)
    {
      free_edge &= (dupIds[T3s[i]] == dupIds[T3s[0]]);
    }
    if (free_edge)
    {
      for (size_t i=0; i<T3s.size(); ++i)
      {
        //std::cout << "remove ; " << T3s[i] << std::endl;
        dupIds[T3s[i]] = -1;
      }
    }
  }
  
  for (size_t i=0; i < dupIds.size(); ++i)
    if (dupIds[i] != -1)dupIds[i]=i;
     
  // now remove any faulty duplicates and preserve one single copy for valid duplicates
  K_CONNECT::unchanged pred(dupIds);
  K_CONNECT::IdTool::compress(connect, pred);
  K_CONNECT::IdTool::compress(ancestors, pred);
  K_CONNECT::IdTool::compress(xc, pred);
}

#ifdef DEBUG_TRI_CONFORMIZER
///

inline E_Float dist2(const E_Float* pt, const E_Float* pt2)
{
  E_Float d2 = 0.;
  for (E_Int k = 0; k < 3; ++k)
    d2 += (pt2[k] - pt[k])*(pt2[k] - pt[k]);
  return d2;
}

template <short DIM>
void
TRI_Conformizer<DIM>::drawElements
(const char* fname, const char* filefmt, const K_FLD::FloatArray& coord, const K_FLD::IntArray& connect,
 const std::vector<T3> & elts, bool localid, std::vector<E_Int>* colors)
{
  K_FLD::IntArray cOut;
  E_Int Ei[2];
  
  K_FLD::IntArray::const_iterator pS;
  std::set<E_Int> nids;
  for (size_t i = 0; i < elts.size(); ++i)
  {
    const T3& t = elts[i];
    if (t.id >= connect.cols())
      continue;
    pS = connect.col(t.id);
    cOut.pushBack(pS, pS+3);
  }
  
  std::vector<E_Int> newIDs;
  K_FLD::FloatArray coord1(coord);
  if (localid)
    NUGA::MeshTool::compact_to_mesh(coord1, cOut, newIDs);
  
  //E_Float d01 = ::sqrt(dist2(coord1.col(0), coord1.col(1)));
  //E_Float d03 = ::sqrt(dist2(coord1.col(0), coord1.col(3)));
  //E_Float d13 = ::sqrt(dist2(coord1.col(1), coord1.col(3)));
  //std::cout << "distances : " << d01 << " " << d03 << " " << d13 << std::endl;
      
  medith::write(fname, coord1, cOut, "TRI", 0, colors);

}
#endif


#ifdef DEBUG_TRI_CONFORMIZER

///
template <short DIM>
void
TRI_Conformizer<DIM>::drawTandE
(const K_FLD::FloatArray& pos, K_FLD::IntArray::const_iterator pS, E_Int Ni, E_Int Nj)
{
  E_Int           E[2];
  K_FLD::IntArray toto;
  K_FLD::FloatArray p;

  E[0] = *pS; E[1] = *(pS+1);
  toto.pushBack(E, E+2);
  E[0] = *(pS+1); E[1] = *(pS+2);
  toto.pushBack(E, E+2);
  E[0] = *(pS+2); E[1] = *pS;
  toto.pushBack(E, E+2);
  E[0] = Ni; E[1] = Nj;
  toto.pushBack(E, E+2);
  p = pos;

  std::vector<E_Int> nids;
  NUGA::MeshTool::compact_to_mesh(p, toto, nids);
  medith::write("TandE.mesh", p, toto);

}

///
template <short DIM>
void
TRI_Conformizer<DIM>::drawT3(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, E_Int ith_elt, bool compact)
{
  
  if (ith_elt >= parent_type::_elements.size())
    return;
  
  const T3& t = parent_type::_elements[ith_elt];

  std::vector<E_Int> colors;
  E_Int Ei[2];
    
  K_FLD::IntArray cOut;

  //initial tri
  /*K_FLD::IntArray::const_iterator pS = connect.col(t.id);
  for (E_Int n = 0; n < 3; ++n)
  {
    Ei[0] = *(pS+n);
    Ei[1] = *(pS+(n+1)%3);
    cOut.pushBack(Ei, Ei+2);
  }
  colors.resize(cOut.cols(), 0);*/

  for (size_t e = 0; e < t.edges.size(); ++e)
  {
    std::vector<E_Int>& E = _edges[t.edges[e]];

    size_t sz = E.size();
    for (size_t k = 0; k < sz-1; ++k)
    {
      Ei[0] = E[k];
      Ei[1] = E[k+1];
      cOut.pushBack(Ei, Ei+2);
    }
  }
  colors.resize(cOut.cols(), 1);
  cOut.pushBack(t.Xedges);
  colors.resize(cOut.cols(), 2);

  std::ostringstream o;
  o << "T3_" << ith_elt << ".mesh";
  //std::cout << cOut << std::endl;
  
  if (compact)
  {
    std::vector<E_Int> newIds;
    K_FLD::FloatArray p(pos);
    NUGA::MeshTool::compact_to_mesh(p, cOut, newIds);
      
    medith::write(o.str().c_str(), p, cOut, "BAR", 0, &colors);
  }
  else
    medith::write(o.str().c_str(), pos, cOut, "BAR", 0, &colors);
  
}

///
template <short DIM>
void
TRI_Conformizer<DIM>::T3nodes
(const K_FLD::IntArray& connect, const T3& t, std::set<E_Int>& nodes)
{

  nodes.clear();
  K_FLD::IntArray::const_iterator pS = connect.col(t.id);

  for (E_Int n = 0; n < 3; ++n)
  {
    nodes.insert(*(pS+n));
    nodes.insert(*(pS+(n+1)%3));
  }

  for (size_t e = 0; e < t.edges.size(); ++e)
  {
    std::vector<E_Int>& E = _edges[t.edges[e]];
    size_t sz = E.size();
    for (size_t k = 0; k < sz-1; ++k)
    {
      nodes.insert(E[k]);
      nodes.insert(E[k+1]);
    }
  }

  std::vector<E_Int> tmp;
  t.Xedges.uniqueVals(tmp);

  nodes.insert(tmp.begin(), tmp.end());
}

///
template <short DIM>
bool
TRI_Conformizer<DIM>::detect_duplis_and_baffles(const K_FLD::IntArray& connect)
{
  return false;
  std::vector<E_Int> dupIds;
  bool has_duplis = NUGA::MeshTool::detectDuplicated(connect, dupIds, false /*strict orient*/);
  return has_duplis;
  
  algo_type::BoundToEltType E_to_T;
  algo_type::BoundToEltType::const_iterator itE;
  K_FLD::ArrayAccessor<K_FLD::IntArray> actv(connect);
  NUGA::EltAlgo<K_MESH::Triangle>::getBoundToElements(actv, E_to_T);
  bool has_free_edges = false;
  for (itE = E_to_T.begin(); (itE != E_to_T.end()) && !has_free_edges; ++itE)
  {
    has_free_edges = (itE->second.size() == 1);
    
    has_free_edges |= (itE->second.size() == 2) && (dupIds[itE->second[0]] == dupIds[itE->second[1]] );
      
  }
  
  return has_free_edges || has_duplis;
}

#endif


}

#endif
