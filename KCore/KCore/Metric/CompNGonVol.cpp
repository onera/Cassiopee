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
# include <vector>

using namespace K_FUNC;
using namespace K_CONST;
using namespace K_FLD;
using namespace std;

//=============================================================================
// IN: (xt, yt, zt): pointeurs sur les coordonnees du maillage
// IN: cnp: pointeur sur la connectivite NGon
// OUT: volp: pointeur sur le tableau des volumes calcules aux centres 
// des elements (deja alloue)
// Return 0 (OK), 1 (Failed)
//=============================================================================
E_Int K_METRIC::CompNGonVol(E_Float* xt,E_Float* yt,E_Float* zt, 
                            FldArrayI& cn, E_Float* volp)
{
  // Donnees liees a la connectivite
  E_Int* cnp = cn.begin(); // pointeur sur la connectivite NGon
  E_Int nfaces = cnp[0]; // nombre total de faces
  E_Int sizeFN = cnp[1]; //  taille de la connectivite Face/Noeuds
  E_Int nelts = cnp[sizeFN+2];  // nombre total d elements
  E_Int* cEFp = cnp+4+sizeFN;// debut connectivite Elmt/Faces

  // Connectivite Element/Noeuds
  vector< vector<E_Int> > cnEV(nelts);
  K_CONNECT::connectNG2EV(cn, cnEV);

  // sommets associes a l'element
  vector<E_Int> vertices;
        
  E_Int nbFaces; // nombre de faces pour un element donne
  E_Int nbNodes; // nombre de noeuds pour une face donnee
  FldArrayI posFace(nfaces); // tableau de position des faces dans la connectivite
  K_CONNECT::getPosFaces(cn, posFace);
  E_Int* posFacep = posFace.begin(); // pointeur sur posFace
  FldArrayI dimElt(nelts); // tableau de la dimensions des elements
  K_CONNECT::getDimElts(cn, posFace, dimElt);
  E_Int dim; // dimension d un element
  E_Int pos; // position d une face donnee dans la connectivite
  E_Float xbe, ybe, zbe; // coordonnees du barycentre d un element
  E_Float xbf, ybf, zbf; // coordonnees du barycentre d une face
  E_Float xc, yc, zc; // coordonnees du barycentre d un triangle composant une face
  E_Int ind, ind1, ind2; // indices de noeuds
  E_Float surfnx, surfny, surfnz; // normale a la surface d un triangle 
  E_Float sens; // sens de la normale d une face (>0 : interieure, <0 : exterieure) 
  E_Float l1x, l1y, l1z, l2x, l2y, l2z; // delta de coordonnees de noeuds

  // parcours des elements
  for (E_Int elt = 0; elt < nelts; elt++)
  {
    volp[elt] = 0.;
    dim = dimElt[elt]; // dimension de l element
    switch(dim)
    {
      case 1: // NGon 1D
      {
        E_Float dx, dy, dz; // delta  de coordonnees de noeuds
        volp[elt] = 0.;
        nbFaces = cEFp[0];
        ind1 = cnp[posFacep[cEFp[1]-1]+1]-1; // indice du point de la premiere face
        for (E_Int fa = 1; fa < nbFaces; fa++) // parcours des faces de l'element elt
        {
          pos = posFacep[cEFp[fa+1]-1];
          nbNodes = cnp[pos];
          ind2 = cnp[pos+1]-1; // indice du point associe a la face
          dx = xt[ind2]-xt[ind1]; dy = yt[ind2]-yt[ind1]; dz = zt[ind2]-zt[ind1];
          volp[elt] += sqrt(dx*dx + dy*dy + dz*dz); //  "volume" = longueur du segment
          ind1 = ind2;
        }
        cEFp += nbFaces + 1;
     }
      break;
      case 2: // NGon 2D
      {
        // sommets associes a l'element
        vertices = cnEV[elt];

        // calcul du barycentre be (xbe, ybe, zbe) de l'element
        nbNodes = vertices.size();
        xbe = 0.; ybe = 0.; zbe = 0.;
        for (E_Int n = 0; n < nbNodes; n++)
        {
          ind = vertices[n]-1;
          xbe += xt[ind]; ybe += yt[ind]; zbe += zt[ind];
        }
        xbe = xbe/nbNodes; ybe = ybe/nbNodes; zbe = zbe/nbNodes;

        // parcours des faces de l element elt
        nbFaces = cEFp[0];
        for (E_Int fa = 0; fa < nbFaces; fa++)
        {
          pos = posFacep[cEFp[fa+1]-1];
          ind1 = cnp[pos+1]-1;
          ind2 = cnp[pos+2]-1;
          // calcul de la normal au triangle (n, n+1, be)
          l1x = xt[ind2]-xt[ind1]; l1y = yt[ind2]-yt[ind1]; l1z = zt[ind2]-zt[ind1];
          l2x = xt[ind2]-xbe; l2y = yt[ind2]-ybe; l2z = zt[ind2]-zbe;
          surfnx = l1y*l2z-l1z*l2y; 
          surfny = l1z*l2x-l1x*l2z; 
          surfnz = l1x*l2y-l1y*l2x;
          volp[elt] += ONE_HALF*sqrt(surfnx*surfnx+surfny*surfny+surfnz*surfnz);
        }
        cEFp += nbFaces + 1;
      }
      break;
      case 3: // NGon 3D
      {
        // sommets associes a l'element
        vertices = cnEV[elt];

        // calcul du barycentre be (xbe, ybe, zbe) de l'element
        nbNodes = vertices.size();
        xbe = 0.; ybe = 0.; zbe = 0.;
        for (E_Int n = 0; n < nbNodes; n++)
        {
          ind = vertices[n]-1;
          xbe += xt[ind]; ybe += yt[ind]; zbe += zt[ind];
        }
        xbe = xbe/nbNodes; ybe = ybe/nbNodes; zbe = zbe/nbNodes;
              
        // parcours des faces de l element elt
        nbFaces = cEFp[0]; 
        for (E_Int fa = 0; fa < nbFaces; fa++)
        {
          pos = posFacep[cEFp[fa+1]-1];
          nbNodes = cnp[pos]; pos++;
          // calcul du barycentre bf (xbf, ybf, zbf) de la face
          xbf = 0.; ybf = 0.; zbf = 0.;
          for (E_Int n = 0; n < nbNodes; n++) // parcours des noeuds de la face fa
          {
            ind = cnp[pos+n]-1;
            xbf += xt[ind]; ybf += yt[ind]; zbf += zt[ind];
          }
          xbf = xbf/nbNodes; ybf = ybf/nbNodes; zbf = zbf/nbNodes;            
          for (E_Int n = 0; n < nbNodes; n++ ) // parcours des noeuds de la face fa
          {
            ind1 = cnp[pos+n]-1; // indice du point n
            ind2 = cnp[pos+(n+1)%(nbNodes)]-1; // indice du point n+1
            // calcul de la normal au triangle (n, n+1, bf)
            l1x = xt[ind2]-xt[ind1]; l1y = yt[ind2]-yt[ind1]; l1z = zt[ind2]-zt[ind1];
            l2x = xt[ind2]-xbf; l2y = yt[ind2]-ybf; l2z = zt[ind2]-zbf;
            surfnx = ONE_HALF * (l1y*l2z-l1z*l2y); 
            surfny = ONE_HALF * (l1z*l2x-l1x*l2z); 
            surfnz = ONE_HALF * (l1x*l2y-l1y*l2x);

            // verification du sens de la normale. Celle-ci doit etre exterieure
            sens = (xbe - xbf)*surfnx + (ybe - ybf)*surfny + (zbe - zbf)*surfnz;
            if (sens > 0.) {surfnx = -surfnx; surfny = -surfny; surfnz = -surfnz;}

            // ajout de la contribution du triangle (n, n+1, bf) au volume de l element elt
            xc = (xt[ind1]+xt[ind2]+xbf)/3.; 
            yc = (yt[ind1]+yt[ind2]+ybf)/3.; // coordonnees du centre du triangle(n, n+1, bf) 
            zc = (zt[ind1]+zt[ind2]+zbf)/3.; 
            volp[elt] += xc * surfnx + yc * surfny + zc * surfnz;
          }
        }
        volp[elt] = volp[elt]/3.;
        cEFp += nbFaces + 1;
      }
      break;
      default:
        return 1;

    }
  }
  return 0;
}



//=============================================================================
// IN: (xt, yt, zt): pointeurs sur les coordonnees du maillage
// IN: cnp: pointeur sur la connectivite NGon
// IN: cnEV: Connectivite Element/Noeuds pour l element considere 
// IN: posFace: tableau de position des elements dans la connectivite
// IN: posFace: tableau de position des faces dans la connectivite
// IN: dim: dimension de l element
// OUT: vol: volume de l element calcule au centre
// Return 0 (OK), 1 (Failed)
//=============================================================================
E_Int K_METRIC::CompNGonVolOfElement(E_Float* xt,E_Float* yt,E_Float* zt, FldArrayI& cn,
                                     E_Int indE, vector< vector<E_Int> > cnEV, FldArrayI& posElt, 
                                     FldArrayI& posFace, FldArrayI& dimElt,
                                     E_Float& vol)
{

  // Donnees liees a la connectivite
  E_Int* cnp = cn.begin(); // pointeur sur la connectivite NGon
  E_Int sizeFN = cnp[1]; //  taille de la connectivite Face/Noeuds
  E_Int* cEFp = cnp+4+sizeFN;// debut connectivite Elmt/Faces
  
  // sommets associes a l element
  vector<E_Int> vertices;
        
  E_Int nbFaces; // nombre de faces pour un element donne
  E_Int nbNodes; // nombre de noeuds pour une face donnee
  E_Int* posFacep = posFace.begin(); // pointeur sur posFace
  E_Int* posEltp = posElt.begin(); // pointeur sur posElt

  E_Int pos; // position d une face donnee dans la connectivite
  E_Float xbe, ybe, zbe; // coordonnees du barycentre de l element
  E_Float xbf, ybf, zbf; // coordonnees du barycentre d une face
  E_Float xc, yc, zc; // coordonnees du barycentre d un triangle composant une face
  E_Int ind, ind1, ind2; // indices de noeuds
  E_Float surfnx, surfny, surfnz; // normale a la surface d un triangle 
  E_Float sens; // sens de la normale d une face (>0 : interieure, <0 : exterieure) 
  E_Float l1x, l1y, l1z, l2x, l2y, l2z; // delta de coordonnees de noeuds
  
  vol = 0.;
  E_Int dim = dimElt[indE];// dimension de l element
  switch(dim)
  {
    case 1: // NGon 1D
    {
      E_Float dx, dy, dz; // delta  de coordonnees de noeuds
      vol = 0.;
      cEFp += posEltp[indE]-posEltp[0];
      nbFaces = cEFp[0];
      ind1 = cnp[posFacep[cEFp[1]-1]+1]-1; // indice du point de la premiere face
      for (E_Int fa = 1; fa < nbFaces; fa++) // parcours des faces de l element elt
      {
        pos = posFacep[cEFp[fa+1]-1];
        nbNodes = cnp[pos];
        ind2 = cnp[pos+1]-1; // indice du point associe a la face
        dx = xt[ind2]-xt[ind1]; dy = yt[ind2]-yt[ind1]; dz = zt[ind2]-zt[ind1];
        vol += sqrt(dx*dx + dy*dy + dz*dz); //  "volume" = longueur du segment
        ind1 = ind2;
      }
    }
    break;
    case 2: // NGon 2D
    {
      // sommets associes a l element
      vertices = cnEV[indE];
      
      // calcul du barycentre be (xbe, ybe, zbe) de l element
      nbNodes = vertices.size();
      xbe = 0.; ybe = 0.; zbe = 0.;
      for (E_Int n=0; n<nbNodes; n++)
      {
        ind = vertices[n]-1;
        xbe += xt[ind]; ybe += yt[ind]; zbe += zt[ind];
      }
      xbe = xbe/nbNodes; ybe = ybe/nbNodes; zbe = zbe/nbNodes;
      
      // parcours des faces de l element elt
      cEFp += posEltp[indE]-posEltp[0];
      nbFaces = cEFp[0];
      for (E_Int fa = 0; fa < nbFaces; fa++)
      {
        pos = posFacep[cEFp[fa+1]-1];
        ind1 = cnp[pos+1]-1;
        ind2 = cnp[pos+2]-1;
        // calcul de la normal au triangle (n, n+1, be)
        l1x = xt[ind2]-xt[ind1]; l1y = yt[ind2]-yt[ind1]; l1z = zt[ind2]-zt[ind1];
        l2x = xt[ind2]-xbe; l2y = yt[ind2]-ybe; l2z = zt[ind2]-zbe;
        surfnx = l1y*l2z-l1z*l2y; 
        surfny = l1z*l2x-l1x*l2z; 
        surfnz = l1x*l2y-l1y*l2x;
        vol += ONE_HALF * sqrt(surfnx*surfnx+surfny*surfny+surfnz*surfnz);
      }
    }
    break;
    case 3: // NGon 3D
    {
      // sommets associes a l element
      vertices = cnEV[indE];
      
      // calcul du barycentre be (xbe, ybe, zbe) de l element
      nbNodes = vertices.size();
      xbe = 0.; ybe = 0.; zbe = 0.;
      for (E_Int n=0; n<nbNodes; n++)
      {
        ind = vertices[n]-1;
        xbe += xt[ind]; ybe += yt[ind]; zbe += zt[ind];
      }
      xbe = xbe/nbNodes; ybe = ybe/nbNodes; zbe = zbe/nbNodes;
              
      // parcours des faces de l element elt
      cEFp += posEltp[indE]-posEltp[0];  // position de l element
      nbFaces = cEFp[0];
      for (E_Int fa = 0; fa < nbFaces; fa++)
      {
        pos = posFacep[cEFp[fa+1]-1];
        nbNodes = cnp[pos]; pos++;
        // calcul du barycentre bf (xbf, ybf, zbf) de la face
        xbf = 0.; ybf = 0.; zbf = 0.;
        for (E_Int n=0; n <nbNodes; n++ ) // parcours des noeuds de la face fa
        {
          ind = cnp[pos+n]-1;
          xbf += xt[ind]; ybf += yt[ind]; zbf += zt[ind];
        }
        xbf = xbf/nbNodes; ybf = ybf/nbNodes; zbf = zbf/nbNodes;                
        for (E_Int n=0; n <nbNodes; n++ ) // parcours des noeuds de la face fa
        {
          ind1 = cnp[pos+n]-1; // indice du point n
          ind2 = cnp[pos+(n+1)%(nbNodes)]-1; // indice du point n+1
          // calcul de la normal au triangle (n, n+1, bf)
          l1x = xt[ind2]-xt[ind1]; l1y = yt[ind2]-yt[ind1]; l1z = zt[ind2]-zt[ind1];
          l2x = xt[ind2]-xbf; l2y = yt[ind2]-ybf; l2z = zt[ind2]-zbf;
          surfnx = ONE_HALF * (l1y*l2z-l1z*l2y); 
          surfny = ONE_HALF * (l1z*l2x-l1x*l2z); 
          surfnz = ONE_HALF * (l1x*l2y-l1y*l2x);

          // verification du sens de la normale. Celle-ci doit etre exterieure
          sens = (xbe - xbf)*surfnx + (ybe - ybf)*surfny + (zbe - zbf)*surfnz;
          if (sens > 0.) {surfnx = -surfnx; surfny = -surfny; surfnz = -surfnz;}

          // ajout de la contribution du triangle (n, n+1, bf) au volume de l element elt
          xc = (xt[ind1]+xt[ind2]+xbf)/3.; 
          yc = (yt[ind1]+yt[ind2]+ybf)/3.; // coordonnees du centre du triangle(n, n+1, bf) 
          zc = (zt[ind1]+zt[ind2]+zbf)/3.; 
          vol += xc * surfnx + yc * surfny + zc * surfnz;
        }
      }
      vol = vol/3.;
    }
    break;
    default:
      return 1;
      
  }
  return 0;
}
