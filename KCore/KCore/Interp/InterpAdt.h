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
# ifndef _INTERP_INTERPADT_H
#define _INTERP_INTERPADT_H

# include "Interp/IntTreeNode.h"
# include "Interp/InterpData.h"
using namespace K_FLD;

namespace K_INTERP 
{
// ============================================================================
// Remarque : en non structure tetra, les noeuds de l adt correspondent
// aux elements et non aux sommets des elts 
//=============================================================================
class InterpAdt  : public InterpData
{
public: 

  ///+ 1- Constructor / Destructor
  
  /* 1- Destructor. */
  virtual ~InterpAdt();
  void destroy();

   InterpAdt();  
  /* Construit l'ADT a partir du maillage donné par xD,yD,zD 
   a1, a2, a3: ni,nj,nk en structure, connect,eltType,NULL en non structure
   OUT: built: 0 si echec, 1 si succes */
   InterpAdt(E_Int npts, E_Float* xD, E_Float* yD, E_Float* zD,
             void* a1, void* a2, void* a3, E_Int& built);

 private:
  /* Construit l'adt a partir d'un maillage structure 
     Retourne 0 si nk=1 mais maillage non plan, 1 dans les autres cas. */
  E_Int buildStructAdt(E_Int ni, E_Int nj, E_Int nk,
                       E_Float* x, E_Float* y, E_Float* z);
  /* construit l'adt a partir d'un maillage TETRA */
  E_Int buildUnstrAdt(E_Int npts, FldArrayI& cEV,
                      E_Float* x, E_Float* y, E_Float* z);
  
  void insert(E_Int ind, 
              E_Float xmin, E_Float ymin, E_Float zmin,
              E_Float xmax, E_Float ymax, E_Float zmax);

  public:
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
  /* Recherche de la liste des cellules candidates. 
     Retourne la taille de listOfCandidateCells */
    virtual E_Int getListOfCandidateCells(E_Float x, E_Float y, E_Float z,
                                          std::list<E_Int>& listOfCandidateCells,
                                          E_Float alphaTol=0.);


  // public:
  // E_Float _xmax;   // bounding box of mesh
  // E_Float _ymax;
  // E_Float _zmax;
  // E_Float _xmin;
  // E_Float _ymin;
  // E_Float _zmin;

  protected:
  IntTreeNode* _tree;
};
}
#endif
