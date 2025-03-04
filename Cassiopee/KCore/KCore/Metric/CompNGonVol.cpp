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
# include <iostream>

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
  // Donnees liees a la connectivite - Acces non universel sur le ptrs
  E_Int* ngon = cn.getNGon();
  E_Int* nface = cn.getNFace();
  E_Int* indPG = cn.getIndPG();
  E_Int* indPH = cn.getIndPH();
  // Acces universel nbre d'elements
  E_Int nelts = cn.getNElts(); // nombre total d elements

  // Connectivite Element/Noeuds
  vector< vector<E_Int> > cnEV(nelts);
  K_CONNECT::connectNG2EV(cn, cnEV);
  // Tableau de la dimensions des elements
  FldArrayI dimElt(nelts); 
  K_CONNECT::getDimElts(cn, dimElt);

  E_Int ierr = 0; // error index

  #pragma omp parallel
  {
    E_Int nbFaces; // nombre de faces pour un element donne
    E_Int nbNodes; // nombre de noeuds pour une face donnee
    E_Int dim; // dimension d un element
    E_Float xbe, ybe, zbe; // coordonnees du barycentre d un element
    E_Float xbf, ybf, zbf; // coordonnees du barycentre d une face
    E_Float xc, yc, zc; // coordonnees du barycentre d un triangle composant une face
    E_Int ind, ind1, ind2; // indices de noeuds
    E_Float surfnx, surfny, surfnz; // normale a la surface d un triangle 
    E_Float sens; // sens de la normale d une face (>0 : interieure, <0 : exterieure) 
    E_Float l1x, l1y, l1z, l2x, l2y, l2z; // delta de coordonnees de noeuds
    E_Int dummy; // dummy variable
    
    // parcours des elements
    #pragma omp for
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
          // Acces universel element elt
          E_Int* elem = cn.getElt(elt, nbFaces, nface, indPH);
          // Acces universel face elem[0]-1
          E_Int* face = cn.getFace(elem[0]-1, dummy, ngon, indPG);
          ind1 = face[0]-1; // indice du point de la premiere face
          for (E_Int fa = 1; fa < nbFaces; fa++) // parcours des faces de l'element elt
          {
            // Acces universel face elem[fa]-1
            E_Int* face = cn.getFace(elem[fa]-1, dummy, ngon, indPG);
            ind2 = face[0]-1; // indice du point associe a la face
            dx = xt[ind2]-xt[ind1]; dy = yt[ind2]-yt[ind1]; dz = zt[ind2]-zt[ind1];
            volp[elt] += sqrt(dx*dx + dy*dy + dz*dz); //  "volume" = longueur du segment
            ind1 = ind2;
          }
        }
        break;
        case 2: // NGon 2D
        {
          // sommets associes a l'element
          const vector<E_Int>& vertices = cnEV[elt];

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
          E_Int* elem = cn.getElt(elt, nbFaces, nface, indPH);
          for (E_Int fa = 0; fa < nbFaces; fa++)
          {
            E_Int* face = cn.getFace(elem[fa]-1, dummy, ngon, indPG);
            ind1 = face[0]-1;
            ind2 = face[1]-1;
            // calcul de la normal au triangle (n, n+1, be)
            l1x = xt[ind2]-xt[ind1]; l1y = yt[ind2]-yt[ind1]; l1z = zt[ind2]-zt[ind1];
            l2x = xt[ind2]-xbe; l2y = yt[ind2]-ybe; l2z = zt[ind2]-zbe;
            surfnx = l1y*l2z-l1z*l2y; 
            surfny = l1z*l2x-l1x*l2z; 
            surfnz = l1x*l2y-l1y*l2x;
            volp[elt] += ONE_HALF*sqrt(surfnx*surfnx+surfny*surfny+surfnz*surfnz);
          }
        }
        break;
        case 3: // NGon 3D
        {
          // sommets associes a l'element
          const vector<E_Int>& vertices = cnEV[elt];

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
          E_Int* elem = cn.getElt(elt, nbFaces, nface, indPH);
          for (E_Int fa = 0; fa < nbFaces; fa++)
          {
            E_Int* face = cn.getFace(elem[fa]-1, nbNodes, ngon, indPG);
            // calcul du barycentre bf (xbf, ybf, zbf) de la face
            xbf = 0.; ybf = 0.; zbf = 0.;
            for (E_Int n = 0; n < nbNodes; n++) // parcours des noeuds de la face fa
            {
              ind = face[n]-1;
              xbf += xt[ind]; ybf += yt[ind]; zbf += zt[ind];
            }
            xbf = xbf/nbNodes; ybf = ybf/nbNodes; zbf = zbf/nbNodes;            
            for (E_Int n = 0; n < nbNodes; n++ ) // parcours des noeuds de la face fa
            {
              ind1 = face[n]-1; // indice du point n
              ind2 = face[(n+1)%(nbNodes)]-1; // indice du point n+1
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
        }
        break;
        default:
          ierr = 1;
      }
    }
  }
  return ierr;
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

  // Donnees liees a la connectivite - Acces non universel sur le ptrs
  E_Int* ngon = cn.getNGon();
  E_Int* nface = cn.getNFace();
  E_Int* indPG = cn.getIndPG();
  E_Int* indPH = cn.getIndPH();

  // sommets associes a l element
  vector<E_Int> vertices;
      
  E_Int nbFaces; // nombre de faces pour un element donne
  E_Int nbNodes; // nombre de noeuds pour une face donnee

  E_Float xbe, ybe, zbe; // coordonnees du barycentre de l element
  E_Float xbf, ybf, zbf; // coordonnees du barycentre d une face
  E_Float xc, yc, zc; // coordonnees du barycentre d un triangle composant une face
  E_Int ind, ind1, ind2; // indices de noeuds
  E_Float surfnx, surfny, surfnz; // normale a la surface d un triangle 
  E_Float sens; // sens de la normale d une face (>0 : interieure, <0 : exterieure) 
  E_Float l1x, l1y, l1z, l2x, l2y, l2z; // delta de coordonnees de noeuds
  E_Int dummy; // dummy variable
  
  vol = 0.;
  E_Int dim = dimElt[indE];// dimension de l element
  switch(dim)
  {
    case 1: // NGon 1D
    {
      E_Float dx, dy, dz; // delta  de coordonnees de noeuds
      vol = 0.;
      // Acces universel element indE
      E_Int* elem = cn.getElt(indE, nbFaces, nface, indPH);
      // Acces universel à la première face, elem[0]-1
      E_Int* face = cn.getFace(elem[0]-1, dummy, ngon, indPG);
      ind1 = face[0]-1; // indice du point de la premiere face
      for (E_Int fa = 1; fa < nbFaces; fa++) // parcours des faces de l element indE
      {
        // Acces universel face elem[fa]-1
        E_Int* face = cn.getFace(elem[fa]-1, dummy, ngon, indPG);
        ind2 = face[0]-1; // indice du point associe a la face
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
      
      // parcours des faces de l element indE
      E_Int* elem = cn.getElt(indE, nbFaces, nface, indPH);
      for (E_Int fa = 0; fa < nbFaces; fa++)
      {
        E_Int* face = cn.getFace(elem[fa]-1, dummy, ngon, indPG);
        ind1 = face[0]-1;
        ind2 = face[1]-1;
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
              
      // parcours des faces de l element indE
      E_Int* elem = cn.getElt(indE, nbFaces, nface, indPH);
      for (E_Int fa = 0; fa < nbFaces; fa++)
      {
        E_Int* face = cn.getFace(elem[fa]-1, nbNodes, ngon, indPG);
        // calcul du barycentre bf (xbf, ybf, zbf) de la face
        xbf = 0.; ybf = 0.; zbf = 0.;
        for (E_Int n=0; n <nbNodes; n++ ) // parcours des noeuds de la face fa
        {
          ind = face[n]-1;
          xbf += xt[ind]; ybf += yt[ind]; zbf += zt[ind];
        }
        xbf = xbf/nbNodes; ybf = ybf/nbNodes; zbf = zbf/nbNodes;                
        for (E_Int n=0; n <nbNodes; n++ ) // parcours des noeuds de la face fa
        {
          ind1 = face[n]-1; // indice du point n
          ind2 = face[(n+1)%(nbNodes)]-1; // indice du point n+1
          // calcul de la normal au triangle (n, n+1, bf)
          l1x = xt[ind2]-xt[ind1]; l1y = yt[ind2]-yt[ind1]; l1z = zt[ind2]-zt[ind1];
          l2x = xt[ind2]-xbf; l2y = yt[ind2]-ybf; l2z = zt[ind2]-zbf;
          surfnx = ONE_HALF * (l1y*l2z-l1z*l2y); 
          surfny = ONE_HALF * (l1z*l2x-l1x*l2z); 
          surfnz = ONE_HALF * (l1x*l2y-l1y*l2x);

          // verification du sens de la normale. Celle-ci doit etre exterieure
          sens = (xbe - xbf)*surfnx + (ybe - ybf)*surfny + (zbe - zbf)*surfnz;
          if (sens > 0.) {surfnx = -surfnx; surfny = -surfny; surfnz = -surfnz;}

          // ajout de la contribution du triangle (n, n+1, bf) au volume de l element indE
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