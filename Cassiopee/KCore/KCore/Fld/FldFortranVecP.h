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
#ifndef _KCORE_FLD_FLDFORTRANVECP_H_
#define _KCORE_FLD_FLDFORTRANVECP_H_

extern "C"
{
  void k6fldgather_(const E_Int& neq, const E_Int& N, const E_Int& szRhs,
                    const E_Int* indic, const E_Float* rhs,
                    E_Float* lhs);
  void k6fldgather2_(const E_Int& neq, const E_Int& N,
                     const E_Int& szindic2, const E_Int& szRhs,
                     const E_Int* indic1, const E_Int* indic2, 
                     const E_Float* rhs,
                     E_Float* lhs);
  // It may be better to exchange lhs and rhs
  void k6setallvaluesf_(const E_Int& N, E_Float* lhs, const E_Float* rhs);
  void k6setallvaluesi_(const E_Int& N, E_Int* lhs, const E_Int* rhs);
  void k6setallvaluesi4_(const E_Int& N, int* lhs, const int* rhs);
  //  void setallvaluesfrompartf_  (const E_Int& N1, const E_Int& N2,
  //                                E_Float* x, E_Float* y);
  void k6setallvaluesatf_(const E_Int& N, E_Float* x, const E_Float& val);
  void k6setallvaluesati_(const E_Int& N, E_Int* x, const E_Int& val);
  void k6operatoregalf_(E_Float* x, const E_Int& first,
                      const E_Int& last, const E_Float& val);
  void k6operatorplusf_(const E_Int& N, E_Float* x, const E_Float* y);
  void k6operatorsubf_(const E_Int& N, E_Float* x, const E_Float* y);
  void k6fldopdivf_(const E_Int& sz, const E_Int& nfld, const E_Int rhsNfld,
                    const  E_Float* lhs, E_Float* rhs);

  void k6multfractionf_(const E_Int& length,
                        E_Float* array,
                        const E_Float* dt,
                        const E_Float* vol);
  void k6multden_(const E_Int& length, const E_Int& nFld,
                const E_Float* vol, E_Float* array);
  void k6multnum_(const E_Int& length, const E_Int& nFld,
                  const E_Float* dt,  E_Float* array);
  void k6fldlib2a_(E_Float* tab, const E_Float* tab0 ,const E_Float* tab1,
                   const E_Float& c0, const E_Float& c1 ,
                   const E_Int& first , const E_Int& last);
  void
  k6conv2double_(const E_Int& size, const E_Float*, double*);
  void
  k6conv2simple_(const E_Int& size, const E_Float*, float*);
  void
  k6multmatvecf_(const E_Int& lenghti, const E_Int& lenghtj, E_Float* x,
                 const E_Float* array);
}
#endif

// ===== Fld/FldFortranVecP.h === Last line ===
