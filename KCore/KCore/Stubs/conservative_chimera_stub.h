/*
 
 
 
              NUGA 
 
 
 
 */

#ifndef NUGA_CONSERVATIVE_CHIMERA_H
#define NUGA_CONSERVATIVE_CHIMERA_H

#include "Nuga/include/conservative_chimera.hxx"
#include "Nuga/Delaunay/Triangulator.h"


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
  
  E_Int compute_chimera_coeffs
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
}

#endif