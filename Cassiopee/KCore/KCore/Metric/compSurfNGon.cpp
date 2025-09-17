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

# include "Connect/connect.h"
# include "metric.h"
# include <vector>

using namespace K_FUNC;
using namespace K_CONST;
using namespace K_FLD;
using namespace std;

//=============================================================================
// IN: (xt, yt, zt): pointeurs sur les coordonnees du maillage NGON surfacique
// IN: cnp: pointeur sur la connectivite NGon
// OUT: sxp, syp, szp: surface orientee calculee aux centres des elements
// Return 0 (OK), 1 (Failed)
//=============================================================================
E_Int K_METRIC::compSurfNGon(
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  FldArrayI& cn,
  E_Float* sxp, E_Float* syp, E_Float* szp)
{
  // Acces non universel sur le ptrs
  E_Int* ngon = cn.getNGon();
  E_Int* nface = cn.getNFace();
  E_Int* indPG = cn.getIndPG();
  E_Int* indPH = cn.getIndPH();
  
  // Acces universel nbre d'elements
  E_Int nelts = cn.getNElts(); // nombre total d elements

  FldArrayI dimElt(nelts); // tableau de la dimension des elements
  K_CONNECT::getDimElts(cn, dimElt);

  E_Int ierr = 0; // error index

  #pragma omp parallel
  {
    // sommets associes a l'element
    vector<E_Int> vertices;
    E_Int nbNodes; // nombre de noeuds pour une face donnee
    E_Int dim; // dimension d'un element
    E_Float xbe, ybe, zbe; // coordonnees du barycentre d'un element
    E_Int ind, ind1, ind2; // indices de noeuds
    E_Float surfnx, surfny, surfnz; // normale a la surface d un triangle
    E_Float l1x, l1y, l1z, l2x, l2y, l2z; // delta de coordonnees de noeuds
    
    // parcours des elements
    #pragma omp for
    for (E_Int elt = 0; elt < nelts; elt++)
    {
      sxp[elt] = 0.; syp[elt] = 0.; szp[elt] = 0.;
      dim = dimElt[elt]; // dimension de l'element
      vertices.clear();
      switch (dim)
      {
        case 1: // NGon 1D
          ierr = 1;
        case 3:
          ierr = 1;

        case 2: // NGon 2D
        {
          // sommets associes a l'elt dans l'ordre
          K_CONNECT::getVertexIndices(cn, ngon, nface, indPG, indPH,
                                      elt, vertices);

          // calcul du barycentre be (xbe, ybe, zbe) de l'element
          nbNodes = vertices.size();
          xbe = 0.; ybe = 0.; zbe = 0.;
          for (E_Int n = 0; n < nbNodes; n++)
          {
            ind = vertices[n]-1;
            xbe += xt[ind]; ybe += yt[ind]; zbe += zt[ind];
          }
          xbe = xbe/nbNodes; ybe = ybe/nbNodes; zbe = zbe/nbNodes;

          // parcours des faces de l'element elt
          for (E_Int fa = 0; fa < nbNodes-1; fa++)
          {
            ind1 = vertices[fa]-1; ind2 = vertices[fa+1]-1;
            // calcul de la normale au triangle (n, n+1, be)
            l1x = xt[ind1]-xbe; l1y = yt[ind1]-ybe; l1z = zt[ind1]-zbe;
            l2x = xt[ind2]-xbe; l2y = yt[ind2]-ybe; l2z = zt[ind2]-zbe;
            surfnx = l1y*l2z-l1z*l2y;
            surfny = l1z*l2x-l1x*l2z;
            surfnz = l1x*l2y-l1y*l2x;
            sxp[elt] += surfnx; syp[elt] += surfny; szp[elt] += surfnz;
          }
          // dernier pour boucler
          ind1 = vertices[nbNodes-1]-1; ind2 = vertices[0]-1;
          // calcul de la normale au triangle (n, n+1, be)
          l1x = xt[ind1]-xbe; l1y = yt[ind1]-ybe; l1z = zt[ind1]-zbe;
          l2x = xt[ind2]-xbe; l2y = yt[ind2]-ybe; l2z = zt[ind2]-zbe;
          surfnx = l1y*l2z-l1z*l2y;
          surfny = l1z*l2x-l1x*l2z;
          surfnz = l1x*l2y-l1y*l2x;
          sxp[elt] += surfnx;  syp[elt] += surfny;  szp[elt] += surfnz;

          sxp[elt] = 0.5*sxp[elt]; syp[elt] = 0.5*syp[elt]; szp[elt] = 0.5*szp[elt];
        }
        break;

        default:
          ierr = 1;
      }
    }
  }
  return ierr;
}