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

# include "CompGeom/compGeom.h"
# include <stdio.h>

using namespace K_FLD;
using namespace std;

//=============================================================================
/* IN: npts: nv de pts sur la surface a projeter (fx,fy,fz)
   IN: nelts2: nb d elements sur la surface de projection
   IN: fx2, fy2, fz2: coordonnees de la surface de projection
   IN: cn2: connectivite elements/noeuds de la surface de projection
   IN: Px, Py, Pz: pt de depart des rayons
   IN/OUT: fx, fy, fz: coordonnees des pts a projeter */
//=============================================================================
void K_COMPGEOM::projectRay(E_Int npts, E_Int nelts2,
                            E_Float Px, E_Float Py, E_Float Pz,
                            E_Float* fx2, E_Float* fy2, E_Float* fz2,
                            K_FLD::FldArrayI& cn2, 
                            E_Float* fx, E_Float* fy, E_Float* fz)
{
  E_Int* cn2p1 = cn2.begin(1);
  E_Int* cn2p2 = cn2.begin(2);
  E_Int* cn2p3 = cn2.begin(3);

#pragma omp parallel
  {
    E_Float p[3]; E_Float pr1[3]; E_Float pr2[3]; E_Float pi[3];
    E_Float p0[3]; E_Float p1[3]; E_Float p2[3];
    E_Float dist; E_Float distc;
    E_Int ret; E_Int ind1, ind2, ind3;

#pragma omp for
    for (E_Int ind = 0; ind < npts; ind++)
    {
      p[0] = fx[ind]; p[1] = fy[ind]; p[2] = fz[ind];
      pr1[0] = p[0]; pr1[1] = p[1]; pr1[2] = p[2];
      pr2[0] = Px; pr2[1] = Py; pr2[2] = Pz;
      distc = 1e6;
      for (E_Int e = 0; e < nelts2; e++)
      {
        ind1 = cn2p1[e]-1; ind2 = cn2p2[e]-1; ind3 = cn2p3[e]-1;
        p0[0] = fx2[ind1]; p0[1] = fy2[ind1]; p0[2] = fz2[ind1];
        p1[0] = fx2[ind2]; p1[1] = fy2[ind2]; p1[2] = fz2[ind2];
        p2[0] = fx2[ind3]; p2[1] = fy2[ind3]; p2[2] = fz2[ind3];        
     
        ret = K_COMPGEOM::intersectRayTriangle(p0, p1, p2,
                                               pr1, pr2,
                                               pi);
        if (ret == 1)
        {
          dist = (pi[0]-p[0])*(pi[0]-p[0]) + (pi[1]-p[1])*(pi[1]-p[1]) +
          (pi[2]-p[2])*(pi[2]-p[2]);
          if (dist < distc)
          {fx[ind] = pi[0]; fy[ind] = pi[1]; fz[ind] = pi[2]; distc = dist;}
        }
      }
    }
  }
}
//=============================================================================
/* IN: nelts2: nb d'elements sur la surface de projection
   IN: fx2, fy2, fz2: coordonnees de la surface de projection
   IN: cn2: connectivite elements/noeuds de la surface de projection
   IN: Px, Py, Pz: pt de depart des rayons
   IN/OUT: fx, fy, fz: coordonnees des pts a projeter */
//=============================================================================
void K_COMPGEOM::projectRay(E_Int nelts2,
                            E_Float Px, E_Float Py, E_Float Pz,
                            E_Float* fx2, E_Float* fy2, E_Float* fz2,
                            K_FLD::FldArrayI& cn2, 
                            vector<E_Int>& sizet, vector<E_Float*>& fxt, vector<E_Float*>& fyt, vector<E_Float*>& fzt)
{
  
  E_Int* cn2p1 = cn2.begin(1);
  E_Int* cn2p2 = cn2.begin(2);
  E_Int* cn2p3 = cn2.begin(3);
  E_Int nzones = sizet.size();

#pragma omp parallel
  {
    E_Float p[3]; E_Float pr1[3]; E_Float pr2[3]; E_Float pi[3];
    E_Float p0[3]; E_Float p1[3]; E_Float p2[3];
    E_Float dist; E_Float distc; 
    E_Int ret; E_Int ind1, ind2, ind3;

    for (E_Int v = 0; v < nzones; v++)
    {
      E_Int npts = sizet[v];
      E_Float* fx = fxt[v];
      E_Float* fy = fyt[v];
      E_Float* fz = fzt[v];
#pragma omp for
      for (E_Int ind = 0; ind < npts; ind++)
      {
        p[0] = fx[ind]; p[1] = fy[ind]; p[2] = fz[ind];
        pr1[0] = p[0]; pr1[1] = p[1]; pr1[2] = p[2];
        pr2[0] = Px; pr2[1] = Py; pr2[2] = Pz;
        distc = 1e6;
        for (E_Int e = 0; e < nelts2; e++) // no procond!
        {
          ind1 = cn2p1[e]-1; ind2 = cn2p2[e]-1; ind3 = cn2p3[e]-1;
          p0[0] = fx2[ind1]; p0[1] = fy2[ind1]; p0[2] = fz2[ind1];
          p1[0] = fx2[ind2]; p1[1] = fy2[ind2]; p1[2] = fz2[ind2];
          p2[0] = fx2[ind3]; p2[1] = fy2[ind3]; p2[2] = fz2[ind3];        

          ret = K_COMPGEOM::intersectRayTriangle(p0, p1, p2,
                                                 pr1, pr2,
                                                 pi);
          if (ret == 1)
          {
            dist = (pi[0]-p[0])*(pi[0]-p[0]) + (pi[1]-p[1])*(pi[1]-p[1]) +
            (pi[2]-p[2])*(pi[2]-p[2]);
            if (dist < distc) 
            {fx[ind] = pi[0]; fy[ind] = pi[1]; fz[ind] = pi[2]; distc = dist;}
          }
        }
      }
    }
  }
}


