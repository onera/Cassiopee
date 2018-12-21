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

# include "Connect/connect.h"
# include "metric.h"
# include <stdio.h>
using namespace K_FUNC;
using namespace K_FLD;
using namespace std;

//=============================================================================
// IN: (xt, yt, zt): pointeurs sur les coordonnees du maillage
// IN: cnp: pointeur sur la connectivite NGon
// IN: cFE: optionnel. Si present, on oriente les surfaces a l'exterieur de ig
// OUT: sxp, syp, szp, snp: surface orientee calculee pour les faces et 
// norme associee (deja alloue)
// Return 0 (OK), 1 (Failed)
//=============================================================================
E_Int K_METRIC::compNGonFacesSurf(
  E_Float* xt, E_Float* yt, E_Float* zt, FldArrayI& cn,
  E_Float* sxp, E_Float* syp, E_Float* szp, E_Float* snp,
  FldArrayI* cFE)
{
  // Donnees liees a la connectivite
  E_Int* cnp = cn.begin(); // pointeur sur la connectivite NGon
  E_Int nfaces = cnp[0]; // nombre total de faces
  
  FldArrayI posFace(nfaces); // tableau de position des faces dans la connectivite
  K_CONNECT::getPosFaces(cn, posFace);
  E_Float xbf, ybf, zbf; // coordonnees du barycentre d'une face
  E_Float surfnx, surfny, surfnz; // normale a la surface d un triangle 
  E_Float l1x, l1y, l1z, l2x, l2y, l2z; // delta de coordonnees de noeuds
  E_Float inv;

  // parcours des faces
  E_Int* ptrf = NULL;
  E_Int nvertices, ind1, ind2, ind, dim;
  for (E_Int fa = 0; fa < nfaces; fa++)
  {
    E_Int posf = posFace[fa];
    ptrf = cnp+posf;
    nvertices = ptrf[0];
    dim = K_FUNC::E_min(nvertices,3);
    if (dim == 1) 
    {
      return 1;
    }
    else if (dim == 2) // surface = longueur
    {
      ind1 = ptrf[1]-1; ind2 = ptrf[2]-1;
      l1x = xt[ind2]-xt[ind1];
      l1y = yt[ind2]-yt[ind1];
      l1z = zt[ind2]-zt[ind1];
      sxp[fa] = l1x;
      syp[fa] = l1y;
      szp[fa] = l1z;
      snp[fa] = sqrt(l1x*l1x+l1y*l1y+l1z*l1z);
    }
    else // 3D: surface d'un NGON 2D
    {
      sxp[fa] = 0.; syp[fa] = 0.; szp[fa] = 0.; snp[fa] = 0.;
      xbf = 0.; ybf = 0.; zbf = 0.;
      // calcul du barycentre de la face
      for (E_Int nv = 1; nv <= nvertices; nv++)
      {
        ind = ptrf[nv]-1;
        xbf += xt[ind]; ybf += yt[ind]; zbf += zt[ind];
      }
      inv = 1./nvertices;
      xbf = xbf*inv; ybf = ybf*inv; zbf = zbf*inv;
    
      // parcours des noeuds de la face
      for (E_Int nv = 1; nv < nvertices; nv++)
      {
        ind1 = ptrf[nv]-1; ind2 = ptrf[nv+1]-1;
        //printf("face %d : %d %d\n", fa,ind1,ind2);
        // calcul de la normale au triangle (nv, nv+1, bf)
        l1x = xt[ind2]-xt[ind1]; l1y = yt[ind2]-yt[ind1]; l1z = zt[ind2]-zt[ind1];
        l2x = xt[ind2]-xbf; l2y = yt[ind2]-ybf; l2z = zt[ind2]-zbf;
        surfnx = l1y*l2z-l1z*l2y;
        surfny = l1z*l2x-l1x*l2z;
        surfnz = l1x*l2y-l1y*l2x;
        sxp[fa] += surfnx; syp[fa] += surfny; szp[fa] += surfnz;       
      }
      // le dernier et le premier pour boucler
      ind1 = ptrf[nvertices]-1; ind2 = ptrf[1]-1;
      //printf("facef %d : %d %d\n", fa,ind1,ind2);
      l1x = xt[ind2]-xt[ind1]; l1y = yt[ind2]-yt[ind1]; l1z = zt[ind2]-zt[ind1];
      l2x = xt[ind2]-xbf; l2y = yt[ind2]-ybf; l2z = zt[ind2]-zbf;
      surfnx = l1y*l2z-l1z*l2y;
      surfny = l1z*l2x-l1x*l2z; 
      surfnz = l1x*l2y-l1y*l2x;
      sxp[fa] += surfnx; syp[fa] += surfny; szp[fa] += surfnz;
      sxp[fa] = 0.5*sxp[fa]; syp[fa] = 0.5*syp[fa]; szp[fa] = 0.5*szp[fa];  
      // norme
      snp[fa] = sqrt(sxp[fa]*sxp[fa]+syp[fa]*syp[fa]+szp[fa]*szp[fa]);
    }
  }

  if (cFE != NULL)
  {
    E_Int* cFE1 = cFE->begin(1);
    E_Int* cFE2 = cFE->begin(2);
    // si cFE est present, on oriente les normales a l'exterieur de ig
    E_Int nelts = cn[2+cn[1]];
    vector< vector<E_Int> > cnEV(nelts);
    K_CONNECT::connectNG2EV(cn, cnEV);

    ptrf = NULL;
    E_Int nbNodes;
    vector<E_Int> vertices;
    for (E_Int fa = 0; fa < nfaces; fa++)
    {
      E_Int posf = posFace[fa];
      ptrf = cnp+posf;
      nvertices = ptrf[0];
      E_Float xbe, ybe, zbe;
      dim = K_FUNC::E_min(nvertices,3);
      
      if (dim != 1 && dim != 2) // 3D
      {
        // elements gauche et droit de la face
        E_Int ig = cFE1[fa];
        E_Int id = cFE2[fa];
        E_Int icell; E_Float isign;
        if (ig != 0) { icell = ig-1; isign = -1.; }
        else { icell = id-1; isign = +1.; }
        // Barycenter of cell
        vertices = cnEV[icell];
        nbNodes = vertices.size();
        xbe = 0.; ybe = 0.; zbe = 0.;
        for (E_Int n = 0; n < nbNodes; n++)
        {
          ind = vertices[n]-1;
          xbe += xt[ind]; ybe += yt[ind]; zbe += zt[ind];
        }
        inv = 1./nbNodes;
        xbe = xbe*inv; ybe = ybe*inv; zbe = zbe*inv;
        xbf = 0.; ybf = 0.; zbf = 0.;
        // calcul du barycentre de la face
        for (E_Int nv = 1; nv <= nvertices; nv++)
        {
          ind = ptrf[nv]-1;
          xbf += xt[ind]; ybf += yt[ind]; zbf += zt[ind];
        }
        inv = 1./nvertices;
        xbf = xbf*inv; ybf = ybf*inv; zbf = zbf*inv;
        // Orientation de la face suivant le produit scalaire xbf-xbe
        E_Float scal = (xbe-xbf)*sxp[fa]+(ybe-ybf)*syp[fa]+(zbe-zbf)*szp[fa];
        if (scal*isign < 0) 
        { sxp[fa] = -sxp[fa]; syp[fa] = -syp[fa]; szp[fa] = -szp[fa]; }
        //printf("face %d (%f %f %f) -> %d %d\n", fa, xbf,ybf,zbf,ig,id);
        //printf("sxyz %f %f %f\n", sxp[fa],syp[fa],szp[fa]);
        //printf("snorm %f\n", snp[fa]);
      }
    }
  }
  
  return 0;
}
