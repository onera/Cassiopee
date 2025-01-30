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

#ifndef _KCORE_METRIC_H
#define _KCORE_METRIC_H
# include "Def/DefTypes.h"
# include "Fld/FldArray.h"
# include "Def/DefFunction.h"

namespace K_METRIC
{
  /* Calcul les volumes pour des elements NGon
     IN: (xt, yt, zt): pointeurs sur les coordonnees du maillage
     IN: cn: connectivite NGon
     OUT: volp: pointeur sur le tableau des volumes calcules aux centres des elements
  */
  E_Int CompNGonVol(E_Float* xt, E_Float* yt, E_Float* zt, 
                    K_FLD::FldArrayI& cn, E_Float* volp); 
  
  E_Int CompNGonVolOfElement(
    E_Float* xt,E_Float* yt,E_Float* zt, 
    K_FLD::FldArrayI& cn, E_Int indE, std::vector<std::vector<E_Int> > cnEV, 
    K_FLD::FldArrayI& posElt, K_FLD::FldArrayI& posFace, 
    K_FLD::FldArrayI& dimElt, E_Float& vol);

  E_Int compNGonSurf(E_Float* xt, E_Float* yt, E_Float* zt, 
                     K_FLD::FldArrayI& cn, 
                     E_Float* sxp, E_Float* syp,  E_Float* szp); 
  
  /* Calcul des surfaces orientees des faces et la norme associee
     On suppose que le NGON est deja correctement oriente
     IN: (xt, yt, zt): pointeurs sur les coordonnees du maillage
     IN:  cnp: pointeur sur la connectivite NGon
     OUT: sxp, syp, szp, snp: surface orientee calculee pour les faces et norme associee
     Return 0 (OK), 1 (Failed)
  */
  E_Int compNGonFacesSurf(
    E_Float* xt, E_Float* yt, E_Float* zt, K_FLD::FldArrayI& cn,
    E_Float* sxp, E_Float* syp,  E_Float* szp, E_Float* snp, 
    K_FLD::FldArrayI* cFE=NULL);  

  /* Calcule l aire d'une cellule d un maillage surfacique nk=1. N'est pas n�cessairement dans le plan 
  On rentre soit l indice de la cellule indcell, soit indnode l indice du premier point
  d indices i et j min de la cellule. Si indnode est different de -1, c'est lui qui prime
  */
  E_Float compVolOfStructCell2D(E_Int ni, E_Int nj, E_Float* xt, E_Float* yt, E_Float* zt,
                                E_Int indcell=-1, E_Int indnode=-1);

  // Compute cell volumes for NGons
  void compute_face_center_and_area(E_Int id, E_Int stride, E_Int *pn,
    E_Float *x, E_Float *y, E_Float *z, E_Float *fc, E_Float *fa);
  
  E_Int compute_volumes_ngon(E_Float *x, E_Float *y, E_Float *z,
    K_FLD::FldArrayI &cn, E_Float *cellVols);

  void compute_cell_volume(E_Int, K_FLD::FldArrayI &, E_Float *, E_Float *,
    E_Float *, E_Float &, E_Int refIdx=0);
  
  void compute_face_centers_and_areas(K_FLD::FldArrayI &cn, E_Float *x,
    E_Float *y, E_Float *z, E_Float *fcenters, E_Float *fareas);
  
  void compute_cell_centers_and_vols(K_FLD::FldArrayI &cn, E_Float *x,
    E_Float *y, E_Float *z, E_Int *owner, E_Int *neigh, E_Float *fcenters,
    E_Float *fareas, E_Float *cx, E_Float *cy, E_Float *cz, E_Float *volumes);
  
  E_Int compute_gradients_ngon(K_FLD::FldArrayI &cn, E_Float *x, E_Float *y,
    E_Float *z, E_Int *owner, E_Int *neigh, E_Float *centers,
    const std::vector<E_Float *> &flds, std::vector<E_Float *> &Gs);
  
  void compute_cell_center_and_volume(E_Int id, E_Int stride,
    E_Int *pf, E_Float *x, E_Float *y, E_Float *z, E_Float *fc, E_Float *fa,
    E_Int *owner, E_Float &cx, E_Float &cy, E_Float &cz, E_Float &vol);
}
#endif
