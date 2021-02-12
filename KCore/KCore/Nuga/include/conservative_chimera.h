/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_CONSERVATIVE_CHIMERA_H
#define NUGA_CONSERVATIVE_CHIMERA_H

#include "Nuga/include/conservative_chimera.hxx"
#include "Nuga/include/Triangulator.h"
#include "Nuga/include/DynArray.h"


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
    enum eMode { CLOSEST_CENTROID, CONSERVATIVE };
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
      //if (mode == CLOSEST_CENTROID)
      return compute_surf_coeffs_basic(crdR, cntR, crdD, cntD, dindices, dcoeffs, xr, do_omp);
      //else
    }

    template <typename acrd_t, typename cnt_t>
    int compute_surf_coeffs_basic(const acrd_t& crdR,
      const cnt_t& cntR,
      const acrd_t& crdD,
      const cnt_t& cntD,
      std::vector<int>& dindices,
      std::vector<double>& dcoeffs,
      std::vector<int>& xr, bool do_omp=false)
    {
      dindices.clear();
      dcoeffs.clear();
      xr.clear();

      int ret(0);
      ngon_unit pgsR(cntR.begin()), pgsD(cntD.begin());

      E_Int npgsD = pgsD.size();

      // 1.1 Compute donnor centroids
      K_FLD::FloatArray isoGD(3, npgsD);
#pragma omp parallel for if(do_omp)
      for (size_t i = 0; i < npgsD; ++i)
        K_MESH::Polygon::iso_barycenter<acrd_t, 3>(crdD, pgsD.get_facets_ptr(i), pgsD.stride(i), 1, isoGD.col(i));
      // 1.2 Build donnor KdTree
      K_FLD::ArrayAccessor<K_FLD::FloatArray> aiso(isoGD);
      K_SEARCH::KdTree<K_FLD::FloatArray> treeD(aiso);

      // 2. Loop through receivers ang get closest
      int npgsR = pgsR.size();

      dindices.resize(npgsR, IDX_NONE);

#pragma omp parallel for if(do_omp)
      for (int i = 0; i < npgsR; ++i)
      {
        E_Float isoGR[3];
        K_MESH::Polygon::iso_barycenter<acrd_t, 3>(crdR, pgsR.get_facets_ptr(i), pgsR.stride(i), 1, isoGR);

        double d2;
        dindices[i] = treeD.getClosest(isoGR, d2);
      }

      //
      dcoeffs.resize(npgsR, 1.);
      K_CONNECT::IdTool::init_inc(xr, npgsR + 1);

      return ret;
    }
  }
}

#endif
