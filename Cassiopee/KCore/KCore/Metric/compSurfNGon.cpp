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
  E_Int dim = cn.getDim();
  if (dim != 2) return 1;

  E_Int* ngon = cn.getNGon(); E_Int* nface = cn.getNFace();
  E_Int* indPG = cn.getIndPG(); E_Int* indPH = cn.getIndPH();
  E_Int nelts = cn.getNElts();

  #pragma omp parallel
  {
    // sommets associes a l'element
    vector<E_Int> vertices;
    E_Int nbNodes; // nombre de noeuds pour un element donne
    E_Float xbe, ybe, zbe; // coordonnees du barycentre d'un element
    E_Int ind, ind1, ind2; // indices de noeuds
    E_Float surfnx, surfny, surfnz; // normale a la surface d un triangle
    E_Float l1x, l1y, l1z, l2x, l2y, l2z; // delta de coordonnees de noeuds
    
    // parcours des elements
    #pragma omp for
    for (E_Int i = 0; i < nelts; i++)
    {
      sxp[i] = 0.; syp[i] = 0.; szp[i] = 0.;

      // sommets associes a l'element dans l'ordre
      K_CONNECT::getVertexIndices(cn, ngon, nface, indPG, indPH, i, vertices);

      // calcul du barycentre be (xbe, ybe, zbe) de l'element
      nbNodes = vertices.size();
      xbe = 0.; ybe = 0.; zbe = 0.;
      for (E_Int n = 0; n < nbNodes; n++)
      {
        ind = vertices[n]-1;
        xbe += xt[ind]; ybe += yt[ind]; zbe += zt[ind];
      }
      xbe = xbe/nbNodes; ybe = ybe/nbNodes; zbe = zbe/nbNodes;

      // parcours des faces de l'element i
      for (E_Int fa = 0; fa < nbNodes-1; fa++)
      {
        ind1 = vertices[fa] - 1; ind2 = vertices[fa+1] - 1;

        // calcul de la normale au triangle (n, n+1, be)
        l1x = xt[ind1] - xbe; l1y = yt[ind1] - ybe; l1z = zt[ind1] - zbe;
        l2x = xt[ind2] - xbe; l2y = yt[ind2] - ybe; l2z = zt[ind2] - zbe;
        surfnx = l1y * l2z - l1z * l2y;
        surfny = l1z * l2x - l1x * l2z;
        surfnz = l1x * l2y - l1y * l2x;
        sxp[i] += surfnx; syp[i] += surfny; szp[i] += surfnz;
      }

      // dernier pour boucler
      ind1 = vertices[nbNodes-1] - 1; ind2 = vertices[0] - 1;

      // calcul de la normale au triangle (n, n+1, be)
      l1x = xt[ind1] - xbe; l1y = yt[ind1] - ybe; l1z = zt[ind1] - zbe;
      l2x = xt[ind2] - xbe; l2y = yt[ind2] - ybe; l2z = zt[ind2] - zbe;

      surfnx = l1y * l2z - l1z * l2y;
      surfny = l1z * l2x - l1x * l2z;
      surfnz = l1x * l2y - l1y * l2x;

      sxp[i] += surfnx;  syp[i] += surfny;  szp[i] += surfnz;
      sxp[i] = 0.5*sxp[i]; syp[i] = 0.5*syp[i]; szp[i] = 0.5*szp[i];
    }
  }

  return 0;
}