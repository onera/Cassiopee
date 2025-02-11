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
#ifndef _INTERP_CART_H_
#define _INTERP_CART_H_

# include "Def/DefTypes.h"
# include "Interp/InterpData.h"
namespace K_INTERP
{
// ============================================================================
// Class pour interpoler d'une grille cartesienne reguliere 
//=============================================================================
class InterpCart  : public InterpData
{
public: 
  
  ///+ 1- Constructor / Destructor
  
  /* 1- Destructor. */
  virtual ~InterpCart();
  
  /* Build data for a donor cell search on an uniform per direction Cartesian grid*/
  InterpCart(E_Int ni, E_Int nj, E_Int nk,
             E_Float hi, E_Float hj, E_Float hk,
             E_Float x0, E_Float y0, E_Float z0,
             E_Int ioff=0, E_Int joff=0, E_Int koff=0);

  virtual
  short searchInterpolationCellCartO2(E_Int ni, E_Int nj, E_Int nk,
                                      E_Float x, E_Float y, E_Float z,
                                      E_Int& ic, E_Int& jc, E_Int& kc,
                                      FldArrayF& cf);
  virtual
  short searchInterpolationCellCartO3(E_Int ni, E_Int nj, E_Int nk,
                                      E_Float x, E_Float y, E_Float z,
                                      E_Int& icHO, E_Int& jcHO, E_Int& kcHO,
                                      FldArrayF& cf);
  virtual 
  short searchExtrapolationCellCart(E_Int ni, E_Int nj, E_Int nk, 
                                    E_Float* xl, E_Float* yl, E_Float* zl,
                                    E_Float* cellNp,
                                    E_Float x, E_Float y, E_Float z,
                                    E_Int& ic, E_Int& jc, E_Int& kc,
                                    FldArrayF& cf,
                                    E_Int nature, E_Int extrapOrder, E_Float constraint);

  virtual
  short searchInterpolationCellUnstruct(E_Float* xt, E_Float* yt, E_Float* zt,
                                        FldArrayI& connectEV,
                                        E_Float x, E_Float y, E_Float z,
                                        E_Int& noelt, FldArrayF& cf);
  virtual 
  short searchInterpolationCellStruct(E_Int ni, E_Int nj, E_Int nk, 
                                      E_Float* xl, E_Float* yl, E_Float* zl,
                                      E_Float x, E_Float y, E_Float z,
                                      E_Int& ic, E_Int& jc, E_Int& kc,
                                      FldArrayF& cf);
  
  virtual
  short searchExtrapolationCellUnstruct(E_Float* xt, E_Float* yt, E_Float* zt,
                                        E_Float* cellNp,
                                        FldArrayI& connectEV,
                                        E_Float x, E_Float y, E_Float z,
                                        E_Int& noet, FldArrayF& cf,
                                        E_Int nature, E_Int extrapOrder, E_Float constraint);
  
  virtual
  short searchExtrapolationCellStruct(E_Int ni, E_Int nj, E_Int nk, 
                                      E_Float* xl, E_Float* yl, E_Float* zl,
                                      E_Float* cellNp,
                                      E_Float x, E_Float y, E_Float z,
                                      E_Int& ic, E_Int& jc, E_Int& kc,
                                      FldArrayF& cf,
                                      E_Int nature, E_Int extrapOrder, E_Float constraint);
  
  virtual 
  E_Int getListOfCandidateCells(E_Float x, E_Float y, E_Float z,
                                std::list<E_Int>& getListOfCandidateCells, 
                                E_Float alphaTol=0.);
private:
  // E_Float _xmin, _ymin, _zmin;
  // E_Float _xmax, _ymax, _zmax;
  E_Int   _ni, _nj, _nk;
  E_Float _hi, _hj, _hk;
  E_Float _his2, _hjs2, _hks2; // for optimization
  E_Float _hii, _hji, _hki;
  E_Float _his2i, _hjs2i, _hks2i;
  E_Float _xminp, _xminm, _yminp, _yminm, _zminp, _zminm;
  E_Int _ioff, _joff, _koff;
};
}
#endif
