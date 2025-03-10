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

#ifndef NUGA_CONSERVATIVE_CHIMERA_H
#define NUGA_CONSERVATIVE_CHIMERA_H

#include "Nuga/include/conservative_chimera.hxx"
#include "Nuga/include/Triangulator.h"
#include "Nuga/include/DynArray.h"
#include "Nuga/include/mesh_t.hxx"
#include "Nuga/include/supermesh.hxx"

namespace NUGA
{
  namespace P1_CONSERVATIVE
  {
    /**
     * Computes conservative chimera donnor coefficients and indices for receiver cells to be interpolated (with cellN == 2)
     * fldR : coordinates accessor of the receiver mesh
     * posRx, posRy, posRz : cooridnates starts in fldR
     * cnR : receiver NGON connectivity
     * fldD : accessor for coordinates of the receiver mesh
     * posDx, posDy, posDz : cooridnates starts in fldD
     * cntD : donnor NGON connectivity
     * cellND : celn nature field for the receiver mesh
     * dindices : donnor indices per receiver cell (ranged with xr)
     * dcoeffs : donnor coefficients per receiver cell (ranged with xr)
     * xr : range delimiter in dindices and dcoefs : the i-th receiver cell has info between xr[i] and xr[i+1]
     * roids : receiver global ids : roids[i] is the i-th receiver's global id
     */

    inline E_Int compute_chimera_coeffs
    (const K_FLD::FldArrayF& fldR, E_Int posRx, E_Int posRy, E_Int posRz,
      const K_FLD::FldArrayI& cnR,
      const K_FLD::FldArrayF& fldD, E_Int posDx, E_Int posDy, E_Int posDz,
      const K_FLD::FldArrayI& cnD,
      const Vector_t<E_Int>& cellNR,
      Vector_t<E_Int>& dindices, Vector_t<E_Float>& dcoeffs,
      Vector_t<E_Int>& xr, Vector_t<E_Int>& roids)
    {
      // HACK CONVERSION OF THE FLD INTO DYNS
      // TEMPORARY HACK COPY to pass to FloaArrays ///////
      K_FLD::FloatArray crdR(fldR, posRx, posRy, posRz);
      K_FLD::FloatArray crdD(fldD, posDx, posDy, posDz);
      K_FLD::IntArray cntR(cnR), cntD(cnD);
      /////////////////////////////////////////////////////

      typedef NUGA::P1_Conservative_Chimera<K_FLD::FloatArray, K_FLD::IntArray> chimera_t;

      return chimera_t::compute_coeffs<DELAUNAY::Triangulator>(crdR, cntR, crdD, cntD, cellNR, dindices, dcoeffs, xr, roids);
    }
  }

  namespace INTERPOL
  {
    enum eMode { CLOSEST_CENTROID, CONSERVATIVE_O1 };

    template <typename acrd_t, typename cnt_t>
    int compute_surf_coeffs_basic(const acrd_t& crdR,
      const cnt_t& cntR,
      const acrd_t& crdD,
      const cnt_t& cntD,
      std::vector<int>& dindices,
      std::vector<double>& dcoeffs,
      std::vector<int>& xr, bool do_omp = false)
    {
      dindices.clear();
      dcoeffs.clear();
      xr.clear();

      E_Int ret(0);
      ngon_unit pgsR(cntR.begin()), pgsD(cntD.begin());

      E_Int npgsD = pgsD.size();

      // 1.1 Compute donnor centroids
      K_FLD::FloatArray isoGD(3, npgsD);
#pragma omp parallel for if(do_omp)
      for (E_Int i = 0; i < npgsD; ++i)
        K_MESH::Polygon::iso_barycenter<acrd_t, 3>(crdD, pgsD.get_facets_ptr(i), pgsD.stride(i), 1, isoGD.col(i));
      // 1.2 Build donnor KdTree
      K_FLD::ArrayAccessor<K_FLD::FloatArray> aiso(isoGD);
      K_SEARCH::KdTree<K_FLD::FloatArray> treeD(aiso);

      // 2. Loop through receivers ang get closest
      E_Int npgsR = pgsR.size();

      dindices.resize(npgsR, IDX_NONE);

#pragma omp parallel for if(do_omp)
      for (E_Int i = 0; i < npgsR; ++i)
      {
        E_Float isoGR[3];
        K_MESH::Polygon::iso_barycenter<acrd_t, 3>(crdR, pgsR.get_facets_ptr(i), pgsR.stride(i), 1, isoGR);

        E_Float d2;
        dindices[i] = treeD.getClosest(isoGR, d2);
      }

      //
      dcoeffs.resize(npgsR, 1.);
      K_CONNECT::IdTool::init_inc(xr, npgsR + 1);

      return ret;
    }

    template <typename acrd_t, typename cnt_t>
    int compute_surf_coeffs_conservative_order1(
      const acrd_t& crdR,
      const cnt_t& cntR,
      const acrd_t& crdD,
      const cnt_t& cntD,
      std::vector<int>& dindices,
      std::vector<double>& dcoeffs,
      std::vector<int>& xr, bool do_omp = false)
    {
      dindices.clear();
      dcoeffs.clear();
      xr.clear();

      E_Int ret(0);
      ngon_unit pgsR(cntR.begin()), pgsD(cntD.begin());

      NUGA::pg_smesh_t m0, m1;

      if (typeid(acrd_t) == typeid(K_FLD::ArrayAccessor<K_FLD::FldArrayF>))
      {
        // HACK CONVERSION OF THE FLD INTO DYNS
        // TEMPORARY HACK COPY to pass to FloaArrays ///////
        K_FLD::FloatArray crdR1(crdR.array(), crdR.posX(0), crdR.posX(1), crdR.posX(2));
        K_FLD::FloatArray crdD1(crdD.array(), crdD.posX(0), crdD.posX(1), crdD.posX(2));
        m0.crd = crdR1;
        m1.crd = crdD1;
      }
      else if (typeid(acrd_t) == typeid(K_FLD::ArrayAccessor<K_FLD::FloatArray>))
      {
        m0.crd = crdR.array();
        m1.crd = crdD.array();
      }

      m0.cnt = pgsR;
      m1.cnt = pgsD;

      NUGA::interpol_coeffs_for_first(m0, m1, 1.e-6, dindices, dcoeffs, xr, do_omp);

      return ret;
    }

    /**
    * Computes donnor coefficients and indices for receiver cells
    * crdR : coordinates accessor of the receiver mesh
    * cntR : receiver surface NGON connectivity (Face to Node representation)
    * crdD : accessor for coordinates of the receiver mesh
    * cntD : donnor NGON connectivity (Face to Node representation)
    * mode : type of interpolation (only CLOSEST_CENTROID is available right now)
    * dindices : 0-based donnor indices per receiver cell (ranged with xr)
    * dcoeffs : donnor coefficients per receiver cell (ranged with xr)
    * xr : range delimiter in dindices and dcoefs : the i-th receiver cell has info between xr[i] and xr[i+1]
    * do_omp : activate Open MP or not
    */
    template <typename acrd_t, typename cnt_t>
    int compute_surf_coeffs(const acrd_t& crdR,
      const cnt_t& cntR,
      const acrd_t& crdD,
      const cnt_t& cntD,
      eMode mode,
      std::vector<int>& dindices,
      std::vector<double>& dcoeffs,
      std::vector<int>& xr,
      bool do_omp = false)
    {
      if (mode == CLOSEST_CENTROID)
        return compute_surf_coeffs_basic(crdR, cntR, crdD, cntD, dindices, dcoeffs, xr, do_omp);
      else if (mode == CONSERVATIVE_O1)
        return compute_surf_coeffs_conservative_order1(crdR, cntR, crdD, cntD, dindices, dcoeffs, xr, do_omp);
      else return 1;
    }

  }

}

#endif
