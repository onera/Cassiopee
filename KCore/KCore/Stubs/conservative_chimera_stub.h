/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

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
    return 0;
  }
  }

  namespace INTERPOL
  {
    enum eMode { CLOSEST_CENTROID, CONSERVATIVE_O1 };
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
      return 0;
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
      return 0;
    }
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
    return 0;
  }
}

#endif
