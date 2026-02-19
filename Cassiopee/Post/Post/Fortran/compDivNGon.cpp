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
# include "post.h"

// ============================================================================
// Calcul de la divergence d'un champ defini aux noeuds d'une grille
// polyedrique
// Retourne la divergence definie aux centres des elements
// ============================================================================
E_Int K_POST::computeDivNGon(
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  const E_Float* fpx, const E_Float* fpy, const E_Float* fpz, FldArrayI& cn,
  E_Float* div
)
{
  E_Int  dim = cn.getDim();
  if (dim < 3)
  {
    printf("computeDiv: not valid for " SF_D_ "D NGONs\n", dim);
    return 1;
  }
  
  // Donnees liees a la connectivite
  E_Int nfaces = cn.getNFaces(); E_Int nelts = cn.getNElts();
  E_Int* ngon = cn.getNGon(); E_Int* indPG = cn.getIndPG();
  E_Int* nface = cn.getNFace(); E_Int* indPH = cn.getIndPH();

  // calcul de la metrique
  E_Float* sxp = new E_Float [3*nfaces];
  E_Float* syp = new E_Float [3*nfaces];
  E_Float* szp = new E_Float [3*nfaces];
  E_Float* snp = new E_Float [nfaces];

  FldArrayI* cFE = new FldArrayI();
  K_CONNECT::connectNG2FE(cn, *cFE);
  K_METRIC::compNGonFacesSurf(xt, yt, zt, cn, sxp, syp, szp, snp, cFE);
  delete cFE;
  E_Float* volp = new E_Float [nelts];
  K_METRIC::compVolNGon(xt, yt, zt, cn, volp);
  // Connectivite Element/Noeuds
  std::vector<std::vector<E_Int> > cnEV(nelts);
  K_CONNECT::connectNG2EV(cn, cnEV); //deja calculee dans NGONVol

  #pragma omp parallel
  {
    E_Int ind, noface, indnode, nbFaces, nbNodes, nbNodesPerFace;
    E_Float fpxmeanface, fpymeanface, fpzmeanface, invvol;
    E_Float xbe, ybe, zbe; // coordonnees du barycentre d un element
    E_Float xbf, ybf, zbf; // coordonnees du barycentre d une face
    E_Float sens, sx, sy, sz;
    std::vector<E_Int> vertices; // sommets associes a l'element

    // parcours des elements
    #pragma omp for
    for (E_Int et = 0; et < nelts; et++)
    {
      invvol = 1./volp[et];
      div[et] = 0.;

      // calcul du barycentre be (xbe, ybe, zbe) de l'element
      vertices = cnEV[et];
      nbNodes = vertices.size();
      xbe = 0.; ybe = 0.; zbe = 0.;
      for (E_Int n = 0; n < nbNodes; n++)
      {
        ind = vertices[n]-1;
        xbe += xt[ind]; ybe += yt[ind]; zbe += zt[ind];
      }
      xbe = xbe/nbNodes; ybe = ybe/nbNodes; zbe = zbe/nbNodes;

      // parcours des faces de l element et
      E_Int* elt = cn.getElt(et, nbFaces, nface, indPH);
      for (E_Int fa = 0; fa < nbFaces; fa++)
      {
        noface = elt[fa]-1;
        E_Int* face = cn.getFace(noface, nbNodesPerFace, ngon, indPG);
        // valeur moyenne de fpx,fpy,fpz pour la face
        fpxmeanface = 0.; fpymeanface = 0.; fpzmeanface = 0.;
        // calcul du barycentre bf (xbf, ybf, zbf) de la face
        xbf = 0.; ybf = 0.; zbf = 0.;
        for (E_Int n = 0; n < nbNodesPerFace; n++)
        {
          indnode = face[n]-1;
          xbf += xt[indnode]; ybf += yt[indnode]; zbf += zt[indnode];
          fpxmeanface += fpx[indnode];
          fpymeanface += fpy[indnode];
          fpzmeanface += fpz[indnode];
        }
        xbf = xbf/nbNodesPerFace; ybf = ybf/nbNodesPerFace; zbf = zbf/nbNodesPerFace;
        fpxmeanface /= nbNodesPerFace;
        fpymeanface /= nbNodesPerFace;
        fpzmeanface /= nbNodesPerFace;
        // bilan
        // verification du sens de la normale. Celle-ci doit etre exterieure
        sx = sxp[noface]; sy = syp[noface]; sz = szp[noface];
        sens = (xbe-xbf)*sx + (ybe-ybf)*sy + (zbe-zbf)*sz;
        if (sens > 0.) {sx=-sx; sy=-sy; sz=-sz;}
        div[et] += fpxmeanface*sx + fpymeanface*sy + fpzmeanface*sz;
      }
      div[et] *= invvol;
    }
  }

  delete [] volp;
  delete [] sxp;
  delete [] syp;
  delete [] szp;
  delete [] snp;
  return 0;
}