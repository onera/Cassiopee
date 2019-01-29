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
//Author : Sâm Landier (sam.landier@onera.fr)

#include "Triangulator.h"
#include "T3Mesher.h"
#include "Connect/MeshTool.h"
#include "Connect/IdTool.h"
#include "Nuga/GapFixer/FittingBox.h"
#include "Connect/GeomAlgo.h"

#ifdef DEBUG_TRIANGULATOR
#include "IO/io.h"
#endif

#ifdef FLAG_STEP
#include "chrono.h"
#endif

#ifdef DEBUG_TRIANGULATOR
    bool DELAUNAY::Triangulator::dbg_enabled = false;
#endif

#define Vector_t std::vector
namespace DELAUNAY
{
  
#ifdef NETBEANSZ
inline 
#endif
E_Int Triangulator::run
(const K_FLD::FloatArray& coord, const E_Int* pNodes, E_Int nb_nodes, E_Int index_start, 
       K_FLD::IntArray& connectM, K_FLD::IntArray& neighbors, bool do_not_shuffle, bool improve_quality) const 
{
   //OVERWRITE mesh upon exit.
  #ifdef FLAG_STEP
  chrono c;
  c.start();
#endif
  
  K_FLD::IntArray connectE2, connectE2b;
  __set_connectE2(pNodes, nb_nodes, connectE2, index_start);
  
#ifdef FLAG_STEP
  tconnect +=c.elapsed();
  c.start();
#endif
  
  K_FLD::FloatArray Wpos;
  Vector_t<E_Int> oldIds;
  K_CONNECT::MeshTool::compact_to_mesh(coord, connectE2, Wpos, connectE2b, oldIds);
  
#ifdef FLAG_STEP
  tcompact +=c.elapsed();
  c.start();
#endif
  
  
  // Computes the fitting box coordinate system optimizing the view over the contour
  K_FLD::FloatArray P(3,3), iP(3,3);
  E_Float W[3];
  FittingBox::computeNormalToContour(Wpos, connectE2b, W);
  FittingBox::computeAFrame(W, P);
  iP = P;
  K_FLD::FloatArray::inverse3(iP);
  FittingBox::transform(Wpos, iP);// Now we are in the fitting coordinate system.
  
  Wpos.resize(2, Wpos.cols()); // Wpos is 2D now.
  
#ifdef FLAG_STEP
  ttra +=c.elapsed();
  c.start();
#endif
  
  // Triangulate the projected contour.
  DELAUNAY::MesherMode mode;
  mode.mesh_mode = mode.TRIANGULATION_MODE;
  mode.do_not_shuffle = do_not_shuffle; 
  mode.silent_errors = true;
  DELAUNAY::T3Mesher<E_Float> mesher(mode);
  DELAUNAY::MeshData data(Wpos, connectE2b);
  
#ifdef FLAG_STEP
  tcreate +=c.elapsed();
  c.start();
#endif
  
#if defined(DEBUG_TRIANGULATOR) && defined(DEBUG_MESHER)
  if (dbg_enabled)
    mesher.dbg_flag=true;
#endif
    
  E_Int err = mesher.run(data);
  
  if (!err && (data.connectM.cols() == 0) && (data.connectB->cols() != 0))
    err = 1;

  if (err)
  {
#ifdef DEBUG_TRIANGULATOR
      MIO::write("data0.mesh", coord, connectE2, "BAR");
      MIO::write("dataT.mesh", Wpos, connectE2b, "BAR");
#endif
      return err;
  }
  
  // In case we need to improve the output quality by swapping.
  if (improve_quality) // fixme : means specified (what to do for relative tol is given as negative ?)
  {
#ifdef DEBUG_TRIANGULATOR
    if (dbg_enabled)
    {
      K_FLD::FloatArray crd(*data.pos);
      crd.resize(3,crd.cols(), 0.);
      MIO::write("before_swap.mesh", crd, data.connectM, "TRI");
    }
#endif
    E_Float quality_tol = 1.e-6;// E_EPSILON;
    E_Int railing = data.connectM.cols();
    while (--railing)
    {
      _swapE.clear();
      
      // also swap poor quality triangles
      K_CONNECT::GeomAlgo<K_MESH::Triangle>::get_swapE(*data.pos, data.connectM, data.neighbors, data.hardEdges, quality_tol, _swapE); //warning : quality<2> doesnt give the same result as quality<3>. need to convert to tolerance criterium.

      if (!_swapE.empty())
        K_CONNECT::EltAlgo<K_MESH::Triangle>::fast_swap_edges(_swapE, data.connectM, data.neighbors);
      else
        break;
    }

#ifdef DEBUG_TRIANGULATOR
    if (dbg_enabled)
    {
      {
      K_FLD::FloatArray crd(*data.pos);
      crd.resize(3,crd.cols(), 0.);
      MIO::write("after_swap.mesh", crd, data.connectM, "TRI");
      }
    }
#endif
  }
  
#ifdef FLAG_STEP
  trun +=c.elapsed();
  c.start();
#endif
  
  K_FLD::IntArray::changeIndices(data.connectM, oldIds);//Back to original ids
  
  connectM  = data.connectM;
  neighbors = data.neighbors;

  return 0;
}
  


#ifdef NETBEANSZ
inline 
#endif
E_Int Triangulator::run
(const K_FLD::FloatArray& coord, const E_Int* pNodes, E_Int nb_nodes, E_Int index_start, K_FLD::IntArray& connectM, bool do_not_shuffle, bool improve_quality) const 
{
  //APPENDING mesh upon exit.
  K_FLD::IntArray cM, neighbors;
  E_Int err = run(coord, pNodes, nb_nodes, index_start, cM, neighbors, do_not_shuffle, improve_quality);
  
  if (!err)
    connectM.pushBack(cM);
  
  return err;
}

///
#ifdef NETBEANSZ
inline 
#endif
E_Int Triangulator::__set_connectE2
(const E_Int* pNodes, E_Int nb_nodes, K_FLD::IntArray& connectE2, E_Int index_start)
{
  connectE2.clear();
  connectE2.reserve(2, nb_nodes);
  
  E_Int E[2];
  const E_Int* p = pNodes;
  for (E_Int i = 0; i < nb_nodes-1; ++i, ++p)
  {
    E[0]=(*p)-index_start;
    E[1]=*(p+1)-index_start;
    connectE2.pushBack(&E[0], &E[0]+2);
  }
  
  E[0]=(*p)-index_start;
  E[1]=(*pNodes)-index_start;
  connectE2.pushBack(&E[0], &E[0]+2);
  
  return 0;
}

}
