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
# ifndef _INTERP_INTERPDATA_H
#define _INTERP_INTERPDATA_H

# include "Fld/FldArray.h"
# include <list>

using namespace K_FLD;

namespace K_INTERP 
{
// ============================================================================
// Remarque : en non structure tetra, les noeuds de l adt correspondent
// aux elements et non aux sommets des elts 
//=============================================================================
class InterpData
{
public: 
  
  enum InterpolationType
  {
    O2CF, O2ABC, O3ABC, O5ABC, O4ABC
  };
  
  ///+ 1- Constructor / Destructor
  
  /* 1- Destructor. */
  virtual ~InterpData();

   InterpData();  
   InterpData(E_Int topology, E_Float xmin, E_Float ymin, E_Float zmin);

  public:
    virtual
    short searchInterpolationCellUnstruct(E_Float* xt, E_Float* yt, E_Float* zt,
                                          FldArrayI& connectEV,
                                          E_Float x, E_Float y, E_Float z,
                                          E_Int& noelt, FldArrayF& cf)=0;
    virtual
    short searchInterpolationCellStruct(E_Int ni, E_Int nj, E_Int nk, 
                                        E_Float* xl, E_Float* yl, E_Float* zl,
                                        E_Float x, E_Float y, E_Float z,
                                        E_Int& ic, E_Int& jc, E_Int& kc,
                                        FldArrayF& cf)=0;

    virtual
    short searchExtrapolationCellUnstruct(E_Float* xt, E_Float* yt, E_Float* zt,
                                          E_Float* cellNp,
                                          FldArrayI& connectEV,
                                          E_Float x, E_Float y, E_Float z,
                                          E_Int& noet, FldArrayF& cf,
                                          E_Int nature, E_Int extrapOrder, E_Float constraint)=0;
    virtual
    short searchExtrapolationCellStruct(E_Int ni, E_Int nj, E_Int nk, 
                                        E_Float* xl, E_Float* yl, E_Float* zl,
                                        E_Float* cellNp,
                                        E_Float x, E_Float y, E_Float z,
                                        E_Int& ic, E_Int& jc, E_Int& kc,
                                        FldArrayF& cf,
                                        E_Int nature, E_Int extrapOrder, E_Float constraint)=0;

    virtual
    short searchInterpolationCellCartO2(E_Int ni, E_Int nj, E_Int nk,
                                        E_Float x, E_Float y, E_Float z,
                                        E_Int& ic, E_Int& jc, E_Int& kc,
                                        FldArrayF& cf)=0;

    virtual
    short searchInterpolationCellCartO3(E_Int ni, E_Int nj, E_Int nk,
                                        E_Float x, E_Float y, E_Float z,
                                        E_Int& icHO, E_Int& jcHO, E_Int& kcHO,
                                        FldArrayF& cf)=0;
    virtual
    short searchInterpolationCellCartO4(E_Int ni, E_Int nj, E_Int nk,
                                        E_Float x, E_Float y, E_Float z,
                                        E_Int& icHO, E_Int& jcHO, E_Int& kcHO,
                                        FldArrayF& cf)=0;
    virtual 
    short searchExtrapolationCellCart(E_Int ni, E_Int nj, E_Int nk, 
                                      E_Float* xl, E_Float* yl, E_Float* zl,
                                      E_Float* cellNp,
                                      E_Float x, E_Float y, E_Float z,
                                      E_Int& ic, E_Int& jc, E_Int& kc,
                                      FldArrayF& cf,
                                      E_Int nature, E_Int extrapOrder, E_Float constraint)=0;
    short getExtrapolationCoeffForCell(E_Float x, E_Float y, E_Float z,
                                       E_Int ic, E_Int jc, E_Int kc,
                                       FldArrayF& cf,
                                       E_Int ni, E_Int nj, E_Int nk,
                                       E_Float* xl, E_Float* yl, E_Float* zl,
                                       E_Float* cellNp, 
                                       E_Int& is, E_Int& js, E_Int& ks,
                                       E_Int nature, E_Float constraint, E_Float& diff_coeff);

    void coordHexa(E_Int ind, E_Int ni, E_Int nj, E_Int nk,
                   E_Float* xl, E_Float* yl, E_Float* zl,
                   E_Int& ic, E_Int& jc, E_Int& kc,
                   E_Float* xt, E_Float* yt, E_Float* zt);
    E_Boolean coeffInterpHexa(E_Float x, E_Float y, E_Float z,
                              E_Float* xt, E_Float* yt, E_Float* zt,
                              FldArrayF& cf);

    E_Boolean getCellJump(E_Float x, E_Float y, E_Float z,
                          E_Float* xt, E_Float* yt, E_Float* zt,
                          E_Int& isomm,
                          E_Float& xi, E_Float& yi, E_Float& zi);
    void coeffInterpTetra(E_Float x, E_Float y, E_Float z,
                          E_Float xp, E_Float yp, E_Float zp,
                          E_Float xq, E_Float yq, E_Float zq,
                          E_Float xr, E_Float yr, E_Float zr,
                          E_Float xs, E_Float ys, E_Float zs,
                          E_Float& xi, E_Float& yi, E_Float& zi);

    E_Boolean getCoeffInterpHexa(E_Float x, E_Float y, E_Float z,
                                 E_Int isomm,
                                 E_Float xi, E_Float yi, E_Float zi, 
                                 E_Float* xt, E_Float* yt, E_Float* zt,
                                 FldArrayF& cf);

  /* Recherche de la liste des cellules candidates. 
     Retourne la taille de listOfCandidateCells */
    virtual E_Int getListOfCandidateCells(E_Float x, E_Float y, E_Float z,
                                          std::list<E_Int>& listOfCandidateCells,
                                          E_Float alphaTol=0.)=0;


  public:
  E_Int _topology; //0 : cart, 1 : struct, 2 : non struct
  const E_Float _EPS_DET;    // EPS for determinant in tetrahedra inversion
  const E_Float _EPS_TETRA;  // EPS for tetrahedra belonging test
  const E_Float _EPS_GEOM;   // EPS for Geometric intersection
  E_Float _xmin, _ymin, _zmin;
  E_Float _xmax, _ymax, _zmax;
};
}
#endif
