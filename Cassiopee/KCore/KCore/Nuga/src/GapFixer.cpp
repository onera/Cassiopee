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

#ifndef __DELAUNAY_GAP_FIXER_H__
#define __DELAUNAY_GAP_FIXER_H__

#include "Nuga/include/GapFixer.h"
#include "Nuga/include/SurfaceMesher.h"
#include "Nuga/include/UBSSurface.h"
#include "Nuga/include/Imprinter.h"
#include <memory>
#ifdef DEBUG_GAPFIXER
#include "IO/DynArrayIO.h"
#include <sstream>
#include "Nuga/include/EltAlgo.h"
#endif

///
GapFixer::GapFixer(){}
///
GapFixer::~GapFixer(){}


///
E_Int
GapFixer::run
(const K_FLD::FloatArray& posC, E_Int nj, const K_FLD::FloatArray& posB0,
 const K_FLD::IntArray& connectB0, K_FLD::FloatArray& posG,
 K_FLD::IntArray& connectG, E_Bool refine, const K_FLD::FloatArray* coordHP)
{
  connectG.clear();
  posG.clear();
  
  // Fast Return
  if (posC.cols() == 0)      return 0;
  if (posB0.cols() == 0)     return 0;
  if (connectB0.cols() == 0) return 0;
  if (nj == 0)               return 1;

  // Compact the data to work only on contour nodes.
  K_FLD::IntArray              cB0(connectB0);
  K_FLD::FloatArray            pB0(posB0);
  std::vector<E_Int>           new_IDs, hN;

  NUGA::MeshTool::compact_to_mesh(pB0, cB0, new_IDs);
  
  //add specified hard points.
  if (coordHP)
  {
    E_Int n0 = pB0.cols();
    pB0.pushBack(*coordHP);
    hN.resize(coordHP->cols());
    for (E_Int i=n0; i < pB0.cols(); ++i) hN[i-n0] = i;
  }
  
  // Flip the contour orientation (to have a counter clockwise for the outer contour)
  //for (E_Int c = 0; c < cB0.cols(); ++c)
    //std::swap(cB0(0,c), cB0(1,c));

  // Fast Return (just a triangle to fix)
  int_vector_type nodes;
  cB0.uniqueVals(nodes);
  if (nodes.size() == 3)
  {
    posG = pB0;
    connectG.resize(3,1);
    connectG(0,0) = 0; connectG(1,0) = 1; connectG(2,0) = 2;
    return 0;
  }
 
  // Build the spline surface.
  std::unique_ptr<UBSSurface> ubs(UBSSurface::buildUBSSurface(posC, nj));
  if (!ubs.get())
    return 1;

  // Build the triangulated spline surface.
  K_FLD::FloatArray               posST3;
  K_FLD::IntArray                 connectST3;
  
  ubs->triangulate(posST3, connectST3);

#ifdef DEBUG_GAPFIXER
  K_CONVERTER::DynArrayIO::write("triaSpline.mesh", posST3, connectST3, "TRI");
  K_FLD::FloatArray pp = pB0;
  pp.pushBack(posST3);
  K_CONVERTER::DynArrayIO::write("contourOverSpline.mesh", pp, cB0, "BAR");
#endif
  
  // Imprint the contour on the spline surface (get the (u,v) parameters).
  int_vector_type         cell_indices;
  DELAUNAY::Imprinter     imprinter(posST3, connectST3);
  K_FLD::FloatArray       posUV;

  imprinter.run(pB0, cB0, posUV, cell_indices, &hN);

  // Convert to global Q-interpolation.
  __convertToGlobalQInterp(pB0, nj, cell_indices, posUV);
  
#ifdef DEBUG_GAPFIXER
  K_CONVERTER::DynArrayIO::write("imprinted.mesh", posUV, cB0, "BAR");
#endif
 
  // Mesh in the parameters space
  DELAUNAY::SurfaceMeshData<UBSSurface> data(posUV, pB0, cB0, *ubs);
  DELAUNAY::SurfaceMesherMode           mode;
  if (refine == false)
    mode.mesh_mode = mode.TRIANGULATION_MODE;
  else
    mode.symmetrize = true;

  data.hardNodes = hN;
  mode.growth_ratio = 1.2;

  DELAUNAY::SurfaceMesher<UBSSurface> mesher(mode);
  mesher.seed_random(1);
  E_Int err = mesher.run (data);
  if (err || (data.connectM.cols() == 0)) return 1;

  posG = data.pos3D;
  connectG = data.connectM;
  
#ifdef DEBUG_GAPFIXER
  K_CONVERTER::DynArrayIO::write("param.mesh", *data.pos, data.connectM, "TRI");
#endif

  return err;
}

void
GapFixer::__convertToGlobalQInterp
(const K_FLD::FloatArray& posB0, E_Int nj,
 const NUGA::int_vector_type& cell_indices,
 K_FLD::FloatArray& posUV)
{
  E_Int     Ni, Ci, Qi, I, J;
  E_Int     n(nj+1);
  E_Int     COLS(posB0.cols());

  for (Ni = 0; Ni < COLS; ++Ni)
  {
    Ci = cell_indices[Ni];
    if (Ci % 2 == 1)
    {
      posUV(0,Ni) = 1. - posUV(0,Ni);
      posUV(1,Ni) = 1. - posUV(1,Ni);
    }

    Qi = Ci / 2;
    J  = Qi % n;
    I  = Qi / n;
    
    // Global parametrization.
    posUV(0,Ni) += J;
    posUV(1,Ni) += I;
  }
}

#endif
