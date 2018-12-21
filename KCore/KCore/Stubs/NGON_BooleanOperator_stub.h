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

#ifndef __NGON_BOOLEANOPERATOR_H__
#define	__NGON_BOOLEANOPERATOR_H__

//#define DEBUG_EXTRACT
//#define DEBUG_W_PYTHON_LAYER

#include "Fld/ngon_t.hxx"
#include "Search/BbTree.h"
#include "Search/ZoneExtractor.h"
#include "TRI_Conformizer.h"
#include "Connect/BARSplitter.h"
#include "Connect/merge.h"
#include "MeshElement/Polyhedron.h"
#include "Connect/MeshTool.h"
#ifdef FLAG_STEP
#include "chrono.h"
#endif
#include "Nuga/Boolean/Agglomerator.h"
// For the T3Mesher
#include "Nuga/Delaunay/Triangulator.h"
//////////////////

#include "Nuga/include/macros.h"
#include "Splitter.h"

#ifdef DEBUG_EXTRACT
#include "debug.h"
#endif
#ifdef DEBUG_W_PYTHON_LAYER
#include <fstream>
#include <sstream>
#endif
//#include "medit.hxx"

#if defined(DEBUG_BOOLEAN) && !defined(MANUAL)
#include "TRI_debug.h"
#include "NGON_debug.h"
#include <fstream>
//#include "TRI_BooleanOperator.h"
#include <iostream>
#endif

#define Vector_t std::vector
#define TEMPLATE_COORD_CONNECT template <typename Coordinate_t, typename Connectivity_t>
#define NGON_BOOLEAN_CLASS NGON_BooleanOperator<Coordinate_t,Connectivity_t>
#define NGON_BOOLEAN_FLD_SPECIALIZED NGON_BooleanOperator<K_FLD::FldArrayF,K_FLD::FldArrayI>
#define NGON_BOOLEAN_DYN_SPECIALIZED NGON_BooleanOperator<K_FLD::FloatArray,K_FLD::IntArray>
#define NGON_DBG NGON_debug<Coordinate_t, Connectivity_t>
#define TRI_DBG TRI_debug

#define ROUND(x) (((x<E_EPSILON) && (x>-E_EPSILON)) ? 0.: x)
#define robust_det4(a,b,c,d) K_FUNC::zzdet3(ROUND(d[0]-a[0]), ROUND(d[1]-a[1]), ROUND(d[2]-a[2]), ROUND(b[0]-d[0]), ROUND(b[1]-d[1]), ROUND(b[2]-d[2]), ROUND(c[0]-d[0]), ROUND(c[1]-d[1]), ROUND(c[2]-d[2]))

#define zABS(a) ((a < 0) ? -a : a)

#define BLANKED 0.
#define VISIBLE 1.

namespace NUGA
{

/// Generalized Polyhedra Meshes Boolean Operator
TEMPLATE_COORD_CONNECT
class NGON_BooleanOperator {
public : 
    enum eInterPolicy { SOLID_NONE, SOLID_RIGHT, SURFACE_RIGHT, BOTH_SURFACE};
    enum eMergePolicy {PRESERVE_LEFT, PRESERVE_RIGHT/*, PRESERVE_X*/};
    enum eAggregation {NONE=0, CONVEX=1, FULL=2};
    enum eOperation   {INTER, DIFF, UNION, MODIFIED_SOLID, DIFFSURF};
    
    enum eRetCode {EMPTY_X=-7, OK=0,ERROR=1, IMMERSED_1};
    
    typedef ngon_t<Connectivity_t>               ngon_type;    
    typedef K_FLD::ArrayAccessor<Coordinate_t>   ACoordinate_t;
    typedef K_FLD::ArrayAccessor<Connectivity_t> AConnectivity_t;
    
private:
      enum eZone {Z_1=0, Z_2=1, Z_IN=2, /*Z_ALL=3,*/ Z_NONE=E_IDX_NONE};

public:
  /// Constructor with the 2 input GPM M1 & M2.
  //For FldArrays
  NGON_BooleanOperator(const K_FLD::FldArrayF& pos1, E_Int px, E_Int py, E_Int pz, const K_FLD::FldArrayI& cNGON1,
                       const K_FLD::FldArrayF& pos2, E_Int px2, E_Int py2, E_Int pz2, const K_FLD::FldArrayI& cNGON2,
                       E_Float tolerance, eAggregation aggtype){}
  //For DynArrays
  NGON_BooleanOperator(const K_FLD::FloatArray& pos1, const K_FLD::IntArray& cNGON1,
                       const K_FLD::FloatArray& pos2, const K_FLD::IntArray& cNGON2,
                       E_Float tolerance, eAggregation aggtype){}
  /// Destructor.
  ~NGON_BooleanOperator(void){}
  
public:
  /// 
  E_Int Intersection(Coordinate_t& coord, Connectivity_t& connect, eInterPolicy XPol = SOLID_RIGHT, eMergePolicy MergePol = PRESERVE_RIGHT) {return 0;}
  ///
  E_Int Diff(Coordinate_t& coord, Connectivity_t& connect, eInterPolicy XPol = SOLID_RIGHT, eMergePolicy MergePol = PRESERVE_RIGHT) {return 0;}
  /// 
  E_Int Union(Coordinate_t& coord, Connectivity_t& connect, eInterPolicy XPol = SOLID_RIGHT, eMergePolicy MergePol = PRESERVE_RIGHT) {return 0;}
  ///
  E_Int Modified_Solid(Coordinate_t& coord, Connectivity_t& connect, eMergePolicy MergePol = PRESERVE_RIGHT) {return 0;}
  ///
  E_Int Diffsurf(Coordinate_t& coord, Connectivity_t& connect) {return 0;}
  
public:
  ///
  void passPGs(E_Int rk, const Vector_t<E_Int>& PGlist){}
  

  void setTriangulatorParams(bool do_no_shuff, bool improve_qual){}
  
  void setConformizerParams(bool split_swap_afterwards){}

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////// Functions based on mapping between input and output meshes //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /// Conservative Interpolation for Chimera
  ///rec_ids[i] gives the i-th receiver, its donnors range from don_ids[xdon[i]] to don_ids[xdon[1+1]-1].
  // Same for coefs : from don_coefs[xdon[i]] to don_coefs[xdon[1+1]-1].
  E_Int volume_coefficients
       (std::vector<E_Int>& rec_ids, std::vector<E_Int>& xdon, std::vector<E_Int>& don_ids, std::vector<E_Float>& don_coefs){return 0;}
  ///
  E_Int volume_and_centroid_coefficients
       (E_Int rec_op, std::vector<E_Int>& rec_ids, std::vector<E_Int>& xdon, std::vector<E_Int>& don_ids,
        std::vector<E_Float>& don_coefs, K_FLD::FloatArray& piece_centroids, K_FLD::FloatArray& crd, ngon_type& ngoper, std::vector<E_Int>& piece_ids) {return 0;}
  ///
  E_Int conservative_transfer
       (E_Int rec_op /*0 or 1*/, const K_FLD::FloatArray& fields_don, K_FLD::FloatArray& fields_rec, E_Int field = -1 /*-1 means all*/, const K_FLD::FloatArray* grads_don = NULL) {return 0;}
  
  E_Int chim_conservative_transfer(E_Int rec_op /*0 or 1*/, const K_FLD::FloatArray& fields_don, K_FLD::FloatArray& fields_rec, E_Int field = -1/*-1 means all*/){return 0;}
  ///
  E_Int total_mass(const K_FLD::FloatArray& crd, const K_FLD::IntArray& ng_cnt, const K_FLD::FloatArray& fields, const std::vector<E_Float>& vols, std::vector<E_Float>& masses){return 0;}
  
  ///
  E_Int XcellN (std::vector<E_Float>& cellN){return 0;}

  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

};


}

#endif	/* NGON_BOOLEANOPERATOR_H */

