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

#ifndef _KCORE_LOC_H
#define _KCORE_LOC_H
# include "Def/DefTypes.h"
# include "Fld/FldArray.h"
# include "Def/DefFunction.h"
# include <vector>
  
namespace K_LOC
{
  /* algo: traitement pour le cellN.
     algo=0: Si au moins une valeur de cellN=0 ==> 0
     Sinon 1
     algo=1: Si toutes les valeurs des cellN=0 ==> 0
     Sinon 1 */ 
  E_Int center2nodeStruct(K_FLD::FldArrayF& FCenter, 
                          E_Int ni, E_Int nj, E_Int nk,
                          E_Int cellN, E_Int mod,
                          E_Int posx, E_Int posy, E_Int posz,
                          K_FLD::FldArrayF& FNode,
                          E_Int& nin, E_Int& njn, E_Int& nkn,
                          E_Int algo=0);
  E_Int center2nodeStruct_OLD(K_FLD::FldArrayF& FCenter, 
                          E_Int ni, E_Int nj, E_Int nk,
                          E_Int cellN, E_Int mod,
                          E_Int posx, E_Int posy, E_Int posz,
                          K_FLD::FldArrayF& FNode,
                          E_Int& nin, E_Int& njn, E_Int& nkn,
                          E_Int algo=0);

  E_Int center2nodeUnstruct(K_FLD::FldArrayF& FCenter, 
                            K_FLD::FldArrayI& c,
                            E_Int cellN, E_Int mod,
                            E_Int posx, E_Int posy, E_Int posz,
                            K_FLD::FldArrayF& FNode,
                            E_Int algo=0);
  E_Int center2nodeUnstruct_OLD(K_FLD::FldArrayF& FCenter, 
                            K_FLD::FldArrayI& c,
                            E_Int cellN, E_Int mod,
                            E_Int posx, E_Int posy, E_Int posz,
                            K_FLD::FldArrayF& FNode,
                            E_Int algo=0);

  E_Int node2centerStruct(K_FLD::FldArrayF& FNode, 
                          E_Int ni, E_Int nj, E_Int nk,
                          E_Int cellN, E_Int mod, 
                          K_FLD::FldArrayF& FCenter);
  E_Int node2centerStruct_OLD(K_FLD::FldArrayF& FNode, 
                          E_Int ni, E_Int nj, E_Int nk,
                          E_Int cellN, E_Int mod, 
                          K_FLD::FldArrayF& FCenter);

  E_Int node2centerUnstruct(K_FLD::FldArrayF& FNode, 
                            K_FLD::FldArrayI& c,
                            E_Int cellN, E_Int mod, 
                            K_FLD::FldArrayF& FCenter);
  E_Int node2centerUnstruct_OLD(K_FLD::FldArrayF& FNode, 
                            K_FLD::FldArrayI& c,
                            E_Int cellN, E_Int mod, 
                            K_FLD::FldArrayF& FCenter);

  E_Int node2centerNGon(K_FLD::FldArrayF& FNode, K_FLD::FldArrayI& cNG,
                        K_FLD::FldArrayF& FCenter, E_Int sorted=0);
  E_Int node2centerNGon_OLD(K_FLD::FldArrayF& FNode, K_FLD::FldArrayI& cNG,
                        K_FLD::FldArrayF& FCenter, E_Int sorted=0);

  E_Int center2nodeNGon(K_FLD::FldArrayF& FCenter, K_FLD::FldArrayI& cNG,
                        std::vector< std::vector<E_Int> >& cEV, 
                        K_FLD::FldArrayF& FNode, E_Int cellN, E_Int mod, 
                        E_Int algo=0);
  E_Int center2nodeNGon_OLD(K_FLD::FldArrayF& FCenter, K_FLD::FldArrayI& cNG,
                        std::vector< std::vector<E_Int> >& cEV, 
                        K_FLD::FldArrayF& FNode, E_Int cellN, E_Int mod, 
                        E_Int algo=0);

  E_Int node2ExtCenterStruct(E_Int imo, E_Int jmo, E_Int kmo,
                             K_FLD::FldArrayF& FNode,
                             E_Int ime, E_Int jme, E_Int kme, 
                             K_FLD::FldArrayF& FExtCenter);
  E_Int center2ExtCenterStruct(E_Int im, E_Int jm, E_Int km,
                               K_FLD::FldArrayF& FCenter,
                               E_Int ime, E_Int jme, E_Int kme, 
                               K_FLD::FldArrayF& FExtCenter);
  E_Int extCenters2NodeStruct(E_Int ime, E_Int jme, E_Int kme,
                              K_FLD::FldArrayF& FextCenter,
                              E_Int im, E_Int jm, E_Int km,
                              K_FLD::FldArrayF& FNode);

  /* Retourne indTab sous la forme [i,i+1, ..., j,j+1] des centres 
     a partir de la cellule indExt en centres etendus
     Si une des cellules en centres est au bord extrapB = 1
     indTab est alloue ici */
  short fromExtCenters2StdCenters(E_Int ime, E_Int jme, E_Int kme, 
                                  E_Int indExt, E_Int type,  K_FLD::FldArrayI& indTab,
                                  E_Int& extrapB);

  /* Transformation repere cart -> repere cylindrique 
   if depth > 0, theta must be continuous between two adjacent points - work for structured zones only (ni,nj,nk must be provided)
  */
  E_Int cart2Cyl(E_Int npts, E_Float* xt, E_Float* yt, E_Float* zt,
                 E_Float X0, E_Float Y0, E_Float Z0,
                 E_Float ex, E_Float ey, E_Float ez,
                 E_Float* rt, E_Float* thetat, 
                 E_Int ni=0, E_Int nj=0, E_Int nk=0, E_Int depth=0, 
                 E_Float thetaShift=0.);
  /* Transformation repere cylindrique -> repere cartesien */
  E_Int cyl2Cart(E_Int npts, E_Float* rt, E_Float* thetat, 
                 E_Float X0, E_Float Y0, E_Float Z0,
                 E_Float ex, E_Float ey, E_Float ez,
                 E_Float* xt, E_Float* yt, E_Float* zt);
}  
#endif
