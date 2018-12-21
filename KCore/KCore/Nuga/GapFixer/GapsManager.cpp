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

#include "GapsManager.h"
#include "PostNodeAssociator.h"
#include "Zipper.h"
#include "PatchMaker.h"
#include "Plaster.h"
#include "GapFixer.h"
#include "Connect/MeshTool.h"
#include "Connect/merge.h"
#include "Connect/ContourSplitter.h"
#include <iostream>

#ifdef DEBUG_GAPSMANAGER
#include <sstream>
#include "IO/io.h"
#include "Fld/DynArray.h"
#endif

#ifdef E_TIME
#include "Delaunay/chrono.h"
#endif

//std::vector<E_Int> GapsManager::_hard_nodes;

///
void
GapsManager::run
(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray*>& components,
 std::vector<K_FLD::FloatArray>& posFs, std::vector<K_FLD::IntArray>& connectFs,
 GapsManager::eCase pcase, E_Bool refine, eMode mode)
{
  posFs.clear();
  connectFs.clear();

  //
  if ((pcase == POST_NODAL) || (pcase == POST_CENTER))
    refine = false;

  //
  std::vector<E_Int> nmates;
  K_FLD::IntArray connectFixed, OneSurface1, OneSurface2, zipped;

#ifdef E_TIME
  DELAUNAY::chrono c, glob;
  c.start();
  glob.start();
#endif

#ifdef E_TIME
  std::cout << "associating nodes : " ;
  c.start();
#endif

  // Make node association.
  NodeAssociator* assoc = __buildAssociator(pcase);
  assoc->make_pairs(pos, components, nmates, OneSurface1);
  delete assoc;
  
#ifdef DEBUG_GAPSMANAGER
  K_CONVERTER::DynArrayIO::write("OneSurface1.mesh", pos, OneSurface1, "TRI");
#endif

#ifdef E_TIME
  std::cout << c.elapsed() << std::endl;
  c.start();
#endif

  OneSurface2 = OneSurface1;
  // Zip what can be zipped.
  if ((refine == false) /*&& (pcase == POST_CENTER)*/)
  {
    K_FLD::IntArray connectB;
    K_CONNECT::MeshTool::getBoundary(OneSurface1, connectB);
#ifdef DEBUG_GAPSMANAGER
    K_CONVERTER::DynArrayIO::write("OneSurface1B.mesh", pos, connectB, "BAR");
#endif
    
    Zipper::zip(connectB, nmates, zipped);
    OneSurface2.pushBack(zipped);
  }

#ifdef E_TIME
  std::cout << "zip " << c.elapsed() << std::endl;
  c.start();
#endif
  
#ifdef DEBUG_GAPSMANAGER
  K_CONVERTER::DynArrayIO::write("OneSurface2.mesh", pos, OneSurface2, "TRI");
#endif

  std::vector<K_FLD::IntArray> connectBout;
  PatchMaker::run(pos, OneSurface2, nmates,  K_CONST::E_PI_4, connectBout);

#ifdef E_TIME
  std::cout << "patch maker " << c.elapsed() << std::endl;
  c.start();
#endif

  // GAP FIXER
  Plaster p;
  std::vector<E_Int> colError ;
  K_FLD::FloatArray plaster, posG, posFixed;
  K_FLD::IntArray connectG, connectError;
  E_Int col(0), /* freecol(0),*/ ni, err, nb_cols;
  E_Int nb_contours = connectBout.size();
#ifdef DEBUG_GAPSMANAGER
  E_Int errcol(0);
#endif
  
  E_Int i0=-1;
  if (mode == PLANAR)
  {
    // flag biggest contour (external) to discard it.
    __get_external_contour(pos, connectBout, i0);
  }

  for (E_Int i = 0; i < nb_contours; ++i)
  {    
    if (i == i0) continue;
    
#ifdef DEBUG_GAPSMANAGER
    std::cout << "PATCH " << i << std::endl;
#endif
    // check closure roughly
    {
      std::vector<E_Int> unodes;
      connectBout[i].uniqueVals(unodes);
      if (unodes.size() != size_t(connectBout[i].cols()))
      {
#ifdef DEBUG_GAPSMANAGER
        std::cout << "GapsManager Warning : discarding this unclosed contour." << std::endl;
#endif
        continue;
      }
    }

    
#ifdef DEBUG_GAPSMANAGER
    // ease reading of connectivity to check unsane ones
    /*{
      K_FLD::FloatArray p(pos);
      K_FLD::IntArray c(connectBout[i]);
      std::vector<E_Int> new_IDs;
      K_CONNECT::MeshTool::compact_to_mesh(p, c, new_IDs);
      K_CONVERTER::DynArrayIO::write("toto.mesh", p, c, "BAR");
    }*/
#endif
    
    // Skip if it's a free contour.
    /*if (__isFreeContour(connectBout[i], nmates))
    {
    ++freecol;
    continue;
    }*/

    // Flip the contour orientation (to have a counter clockwise for the outer contour)
    nb_cols = connectBout[i].cols();
    for (E_Int c = 0; c < nb_cols; ++c)
      std::swap(connectBout[i](0,c), connectBout[i](1,c));
   
    err = p.make(pos, connectBout[i], /*_hard_nodes,*/ plaster, ni);
    
    if (err)
    {
      std::cout << "GapsManager Error : could not create a plaster for " << i << "-th contour." << std::endl;
#ifdef DEBUG_GAPSMANAGER
      std::ostringstream o;
      o << "plasterErr_" << i << ".mesh";
      K_CONVERTER::DynArrayIO::write(o.str().c_str(), pos, connectBout[i], "BAR");
#endif
    }
    
    if (!err)
    {
      err = GapFixer::run(plaster, ni, pos, connectBout[i], /*_hard_nodes,*/ posG, connectG, refine);
      
      if (err)
      {
      std::cout << "GapsManager Error : could not run GapFixer for " << i << "-th contour." << std::endl;
#ifdef DEBUG_GAPSMANAGER
      std::ostringstream o;
      o << "gapfixerErr_" << i << ".mesh";
      K_CONVERTER::DynArrayIO::write(o.str().c_str(), pos, connectBout[i], "BAR");
      connectError.pushBack(connectBout[i]);
      colError.resize(connectError.cols(), ++errcol);
#endif
      }
      else
      {
        ++col;
        connectG.shift(posFixed.cols());
        posFixed.pushBack(posG);
        connectFixed.pushBack(connectG);
      }
    }
  }
#ifdef E_TIME
  std::cout << "FIXING TIME : " << c.elapsed() << std::endl;
  c.start();
#endif
  
#ifdef DEBUG_GAPSMANAGER
  K_CONVERTER::DynArrayIO::write("connectError.mesh", pos, connectError, "BAR", 0, &colError);
#endif

  // Get back to original ids and assign colors for each connex new connectivity.
  K_FLD::FloatArray posTmp = pos;
  {   
    connectFixed.shift(posTmp.cols());
    assert (posFixed.rows() <= posTmp.rows());
    if (posFixed.rows() != posTmp.rows())
      posFixed.resize(posTmp.rows(), posFixed.cols());
    posTmp.pushBack(posFixed);

    connectFixed.pushBack(zipped); // Contains now all the new mesh.

    K_FLD::ArrayAccessor<K_FLD::FloatArray> pAcc(posTmp);
    std::vector<E_Int> new_Ids;
    ::merge(pAcc, E_EPSILON, new_Ids);
    K_FLD::IntArray::changeIndices(connectFixed, new_Ids);
  }

  // Fill the meshes for exit.
  {
    std::vector<K_FLD::IntArray> cOutS;
    K_CONT_DEF::non_oriented_edge_set_type dummyS;
    // Split the new mesh by connex bits, assign colors and eventually swap.
    ContourSplitter<K_MESH::Triangle, K_MESH::NO_Edge>::splitConnectivity(connectFixed, dummyS, cOutS);
    size_t nb_colors(cOutS.size());

    connectFs.resize(1+nb_colors);

    // Add the cleaned (without overlap) original mesh.
    connectFs[0] = OneSurface1;

    for (size_t i = 0; i < nb_colors; ++i)
    {
/*
      if (refine == false)// Do some edge swapping to improve the resulting mesh (post mode only)
      {
        //fixme ? doable ?
      }
*/
      connectFs[i+1] = cOutS[i];
    }
  }

  // Split the coordinates for exit.
  ContourSplitter<K_MESH::Triangle, K_MESH::NO_Edge>::splitCoordinates(connectFs, posTmp, posFs);

#ifdef E_TIME
  std::cout << "aggregate the connectivities " << c.elapsed() << std::endl;
  std::cout << "TOTAL TIME : " << glob.elapsed() << std::endl;
  std::cout << "SUCCESS : " << col << "/" << connectBout.size() << std::endl;
  std::cout << "FAILED  : " << errcol << "/" << connectBout.size() << std::endl;
  //std::cout << "FREE    : " << freecol << "/" << connectBout.size() << std::endl;
#endif
  //std::cout << "SUCCESS : " << col << "/" << connectBout.size() << std::endl;
  //std::cout << "FAILED  : " << errcol << "/" << connectBout.size() << std::endl;
}

NodeAssociator*
GapsManager::__buildAssociator(eCase pcase)
{
  switch (pcase)
  {
  case POST_NODAL:return new PostNodeAssociator(true);
  case POST_CENTER:return new PostNodeAssociator(false);
  default:return new NodeAssociator();
  }
}

bool
GapsManager::__isFreeContour
(const K_FLD::IntArray& connectB, const std::vector<E_Int> & nmates)
{
  E_Int countFree = 0;
  std::vector<E_Int> cnodes;

  connectB.uniqueVals(cnodes);
  for (size_t i = 0; i < cnodes.size(); ++i)
    if (nmates[cnodes[i]] == Zipper::FREE)
      ++countFree;

  return  ( (E_Float(countFree)/E_Float(cnodes.size())) > 0.8); // Free contour;
}

#define IS_IN_BOX(b1, b2) (b1.minB[0] <= b2.minB[0]) && (b1.minB[1] <= b2.minB[1]) && (b1.minB[2] <= b2.minB[2]) && \
                      (b2.maxB[0] <= b1.maxB[0]) && (b2.maxB[1] <= b1.maxB[1]) && (b2.maxB[2] <= b1.maxB[2])

///
void
GapsManager::__get_external_contour
(const K_FLD::FloatArray& coord, const std::vector<K_FLD::IntArray>& connectBs, E_Int& i0)
{
  K_SEARCH::BBox3D bbox, BBOX;
  
  BBOX.minB[0]=BBOX.minB[1]=BBOX.minB[2]=K_CONST::E_MAX_FLOAT;
  BBOX.maxB[0]=BBOX.maxB[1]=BBOX.maxB[2]=-K_CONST::E_MAX_FLOAT;
  
  std::vector<E_Int> nodes;
  
  for (size_t c = 0; c < connectBs.size(); ++c)
  {
    nodes.clear();
    connectBs[c].uniqueVals(nodes);
    bbox.compute(coord, nodes);
    
    if (IS_IN_BOX(bbox, BBOX))
    {
      BBOX.minB[0]=bbox.minB[0];BBOX.minB[1]=bbox.minB[1];BBOX.minB[2]=bbox.minB[2];
      BBOX.maxB[0]=bbox.maxB[0];BBOX.maxB[1]=bbox.maxB[1];BBOX.maxB[2]=bbox.maxB[2];
      
      i0=c;
    }
  }
}
