/*    
    Copyright 2013-2026 ONERA.

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

#include"kcore.h"

# include "Connect/connect.h"
# include "metric.h"
# include <stdio.h>

using namespace K_FUNC;
using namespace K_FLD;

//=============================================================================
// IN: (xt, yt, zt): pointeurs sur les coordonnees du maillage
// IN: cnp: pointeur sur la connectivite NGon
// IN: cFE: optionnel. Si present, on oriente les surfaces a l'exterieur de indg
// OUT: sxp, syp, szp, snp: surface orientee calculee pour les faces et 
// norme associee (deja alloue)
// Return 0 (OK), 1 (Failed)
//=============================================================================
E_Int K_METRIC::compNGonFacesSurf(
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  FldArrayI& cn, E_Float* sxp, E_Float* syp, E_Float* szp, E_Float* snp,
  FldArrayI* cFE
)
{
  E_Int* ngon = cn.getNGon();
  E_Int* indPG = cn.getIndPG();
  E_Int nfaces = cn.getNFaces();
  E_Int dim = cn.getDim();

  if (dim < 2) return 1;  // only 2D and 3D cases are supported
  
  #pragma omp parallel
  {
    E_Float xbf, ybf, zbf; // coordonnees du barycentre d'une face
    E_Float surfnx, surfny, surfnz; // normale a la surface d un triangle 
    E_Float l1x, l1y, l1z, l2x, l2y, l2z; // delta de coordonnees de noeuds
    E_Float inv;
    E_Int nv, ind1, ind2, ind, kp;

    // parcours des faces
    if (dim == 2)
    {
      #pragma omp for
      for (E_Int f = 0; f < nfaces; f++)
      {
        E_Int* face = cn.getFace(f, nv, ngon, indPG);
        ind1 = face[0]-1; ind2 = face[1]-1;
        l1x = xt[ind2] - xt[ind1];
        l1y = yt[ind2] - yt[ind1];
        l1z = zt[ind2] - zt[ind1];
        sxp[f] = -l1y;
        syp[f] = l1x;
        szp[f] = 0;
        snp[f] = std::sqrt(l1x*l1x+l1y*l1y+l1z*l1z);
      }
    }
    else  // 3D: surface d'un NGON 2D
    {
      #pragma omp for
      for (E_Int f = 0; f < nfaces; f++)
      {
        E_Int* face = cn.getFace(f, nv, ngon, indPG);

        // calcul du barycentre de la face
        xbf = 0.; ybf = 0.; zbf = 0.;
        for (E_Int k = 0; k < nv; k++)
        {
          ind = face[k]-1;
          xbf += xt[ind]; ybf += yt[ind]; zbf += zt[ind];
        }
        inv = 1./nv;
        xbf *= inv; ybf *= inv; zbf *= inv;
      
        // parcours des noeuds de la face
        sxp[f] = 0.; syp[f] = 0.; szp[f] = 0.; snp[f] = 0.;
        for (E_Int k = 0; k < nv; k++)
        {
          if (k == nv-1) kp = 0;
          else kp = k+1;
          ind1 = face[k]-1; ind2 = face[kp]-1;
          // calcul de la normale au triangle (k, k+1, bf)
          l1x = xt[ind2] - xt[ind1];
          l1y = yt[ind2] - yt[ind1];
          l1z = zt[ind2] - zt[ind1];
          l2x = xt[ind2] - xbf; l2y = yt[ind2] - ybf; l2z = zt[ind2] - zbf;
          surfnx = l1y*l2z - l1z*l2y;
          surfny = l1z*l2x - l1x*l2z;
          surfnz = l1x*l2y - l1y*l2x;
          sxp[f] += surfnx; syp[f] += surfny; szp[f] += surfnz;
        }
        sxp[f] *= 0.5; syp[f] *= 0.5; szp[f] *= 0.5;
        // norm
        snp[f] = sqrt(sxp[f]*sxp[f] + syp[f]*syp[f] + szp[f]*szp[f]);
      }
    }
  }

  // si cFE est present, on oriente les normales a l'exterieur de indg
  if (dim == 3 && cFE != NULL)
  {
    E_Int nelts = cn.getNElts();
    std::vector<std::vector<E_Int> > cnEV(nelts);
    K_CONNECT::connectNG2EV(cn, cnEV);

    E_Int* cFE1 = cFE->begin(1);
    E_Int* cFE2 = cFE->begin(2);

    #pragma omp parallel
    {
      E_Int ind, indg, indd, nbNodes, nv;
      E_Float xbe, ybe, zbe, sign, inv, scal;
      E_Float xbf, ybf, zbf;  // coordonnees du barycentre d'une face

      #pragma omp for
      for (E_Int f = 0; f < nfaces; f++)
      {
        // elements gauche et droit de la face
        indg = cFE1[f];
        indd = cFE2[f];
        if (indg > 0) { ind = indg-1; sign = -1.; }
        else { ind = indd-1; sign = +1.; }

        // calcul du barycentre de l element
        const std::vector<E_Int>& vertices = cnEV[ind];
        nbNodes = vertices.size();
        xbe = 0.; ybe = 0.; zbe = 0.;
        for (E_Int k = 0; k < nbNodes; k++)
        {
          ind = vertices[k]-1;
          xbe += xt[ind]; ybe += yt[ind]; zbe += zt[ind];
        }
        inv = 1./nbNodes;
        xbe *= inv; ybe *= inv; zbe *= inv;

        // calcul du barycentre de la face
        xbf = 0.; ybf = 0.; zbf = 0.;
        E_Int* face = cn.getFace(f, nv, ngon, indPG);
        for (E_Int k = 0; k < nv; k++)
        {
          ind = face[k]-1;
          xbf += xt[ind]; ybf += yt[ind]; zbf += zt[ind];
        }
        inv = 1./nv;
        xbf *= inv; ybf *= inv; zbf *= inv;

        // Orientation de la face suivant le produit scalaire xbf-xbe
        scal = (xbe-xbf)*sxp[f] + (ybe-ybf)*syp[f] + (zbe-zbf)*szp[f];
        if (scal*sign < 0) 
        {
          sxp[f] = -sxp[f]; syp[f] = -syp[f]; szp[f] = -szp[f];
        }
      }
    }
  }

  return 0;
}
