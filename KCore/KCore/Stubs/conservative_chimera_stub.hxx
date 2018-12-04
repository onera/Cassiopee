/*
 
 
 
              NUGA 
 
 
 
 */

#ifndef NUGA_CONSERVATIVE_CHIMERA_HXX
#define NUGA_CONSERVATIVE_CHIMERA_HXX

#include "Nuga/Boolean/NGON_BooleanOperator.h"
#include "MeshElement/Polygon.h"

namespace NUGA
{
  static const E_Int INTERPOL = 2;
  static const E_Float RESET_VAL = K_CONST::E_MAX_FLOAT;
  static const E_Float COEF_THRESHOLD2 = E_EPSILON*E_EPSILON;
  
  template <typename crd_t, typename cnt_t>
  class P1_Conservative_Chimera
  {
    public:
      typedef ngon_t<cnt_t> ngon_type;
      typedef K_FLD::ArrayAccessor<crd_t> acrd_t;
    
    public:
        
  /**
   * Computes conservative chimera donnor coefficients and indices for receiver cells to be interpolated (with cellN == 2)
   * acrdR : accessor for coordinates of the donnor mesh
   * cntR : receiver connectivity
   * acrdD : accessor for coordinates of the receiver mesh
   * cntD : donnor connectivity
   * cellNR : celn nature field for the receiver mesh
   * dindices : donnor indices per receiver cell (ranged with xr)
   * dcoeffs : donnor coefficients per receiver cell (ranged with xr)
   * xr : range delimiter in dindices and dcoefs : the i-th receiver cell has info between xr[i] and xr[i+1]
   */
  
  template <typename TriangulatorType> 
  static E_Int compute_coeffs
  (const crd_t& crdR, const cnt_t& cntR,
   const crd_t& crdD, const cnt_t& cntD,
   const Vector_t<E_Int>& cellNR,
   Vector_t<E_Int>& dindices, Vector_t<E_Float>& dcoeffs,
   Vector_t<E_Int>& xr, Vector_t<E_Int>& roids)
  {
  	return 0;
  }
 
};

}

#endif
