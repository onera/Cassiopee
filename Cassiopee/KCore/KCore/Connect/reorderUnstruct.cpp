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

#include "Connect/connect.h"
#include <stdlib.h>

using namespace K_FLD;
using namespace std;

#define SWAPEDGESQUAD {                         \
  temp = cn(ie, edge20);                        \
  cn(ie, edge20) = cn(ie, edge21);              \
  cn(ie, edge21) = temp;                        \
  if ((edge20 == 1 && edge21 == 2)||            \
      (edge20 == 2 && edge21 == 1))             \
  {                                             \
    temp = cn3[ie];                             \
    cn3[ie] = cn4[ie];                          \
    cn4[ie] = temp;                             \
  }                                             \
  else if ((edge20 == 1 && edge21 == 3)||       \
           (edge20 == 3 && edge21 == 1))        \
  {                                             \
    temp = cn2[ie];                             \
    cn2[ie] = cn4[ie];                          \
    cn4[ie] = temp;                             \
  }                                             \
  else if ((edge20 == 1 && edge21 == 4)||       \
           (edge20 == 4 && edge21 == 1))        \
  {                                             \
    temp = cn2[ie];                             \
    cn2[ie] = cn3[ie];                          \
    cn3[ie] = temp;                             \
  }                                             \
  else if ((edge20 == 2 && edge21 == 3)||       \
           (edge20 == 3 && edge21 == 2))        \
  {                                             \
    temp = cn1[ie];                             \
    cn1[ie] = cn4[ie];                          \
    cn4[ie] = temp;                             \
  }                                             \
  else if ((edge20 == 2 && edge21 == 4)||       \
           (edge20 == 4 && edge21 == 2))        \
  {                                             \
    temp = cn1[ie];                             \
    cn1[ie] = cn3[ie];                          \
    cn3[ie] = temp;                             \
  }                                             \
  else if ((edge20 == 3 && edge21 == 4)||       \
           (edge20 == 4 && edge21 == 3))        \
  {                                             \
    temp = cn1[ie];                             \
    cn1[ie] = cn2[ie];                          \
    cn2[ie] = temp;                             \
  }                                             \
  }                                             
         
#define SWAPTWOEDGES(i, j) {                                      \
    t1 = edge1[i]; t2 = edge2[i];                               \
    edge1[i] = edge1[j]; edge2[i] = edge2[j];                   \
    edge1[j] = t1; edge2[j] = t2; }

// ============================================================================
/* 
   Reorder the numerotation of unstructured array.
   Marche pour les TRI-array et les QUAD-array.
   Sa connectivite doit etre propre (cleanConnectivity).
   Retourne 1 si OK. 
 */
// ============================================================================
E_Int K_CONNECT::reorderUnstruct2D(FldArrayF& f, FldArrayI& cn, E_Int dir)
{
  vector< vector<E_Int> > cEEN(cn.getSize());
  if (cn.getNfld() == 3)
    K_CONNECT::connectEV2EENbrs("TRI", f.getSize(), cn, cEEN);
  else if (cn.getNfld() == 4)
    K_CONNECT::connectEV2EENbrs("QUAD", f.getSize(), cn, cEEN);
  else return 0;

  E_Int nt = cn.getNfld();
  E_Int ne = cn.getSize(); // nbre d'elements
  E_Int nev = 0; // nbre d'elements deja visites
  char* isVisited = (char*)calloc(ne, sizeof(char)); // elt deja visite?
  E_Int* mustBeVisited = (E_Int*)malloc(ne * sizeof(E_Int));
  E_Int mbv, p, iv, i, ie, elt, i1, i2, sign1, sign2, temp, ind1, ind2;
  E_Int size;
  E_Int edge1[8], edge2[8]; // a cause des cas degeneres
  E_Int edge10=0, edge11=0, edge20=0, edge21=0;
  E_Int t1, t2, found;

  mbv = 0;
  E_Int* cn1 = NULL; E_Int* cn2 = NULL;
  E_Int* cn3 = NULL; E_Int* cn4 = NULL;
  cn1 = cn.begin(1); cn2 = cn.begin(2);
  if (nt == 4) 
  {
    cn3 = cn.begin(3); cn4 = cn.begin(4);
  }

  while (nev < ne)
  {
    // Recherche le premier elt pas encore visite
    for (p = 0; (isVisited[p] != 0); p++);
    
    // C'est un nouveau composant connexe
    if (dir == -1)
    {
      // Reorder de l'element de reference
      if (nt == 3) // TRI
      {
        temp = cn1[p];
        cn1[p] = cn2[p];
        cn2[p] = temp;
      } 
      else // QUAD
      {
        temp = cn1[p];
        cn1[p] = cn2[p];
        cn2[p] = temp;
        temp = cn3[p];
        cn3[p] = cn4[p];
        cn4[p] = temp;
      }
    }

    mustBeVisited[mbv] = p;
    mbv++; nev++;
    isVisited[p] = 1;
    
    while (mbv > 0)
    {
      mbv--;
      elt = mustBeVisited[mbv];
      size = cEEN[elt].size();
      //printf("elt %d visiting size %d\n", elt, size);

      for (iv = 0; iv < size; iv++)
      {
        ie = cEEN[elt][iv];
        //printf("elt %d visiting elt %d\n", elt, ie);

        if (isVisited[ie] == 0)
        {
          // Trouve l'arrete commune entre les deux elements
          i = 0;
          for (i1 = 1; i1 <= nt; i1++)
          {
            ind1 = cn(elt, i1);
            for (i2 = 1; i2 <= nt; i2++)
            {
              ind2 = cn(ie, i2);
              //printf("inds %d %d\n", ind1, ind2);
              if (ind1 == ind2)
              {
                edge1[i] = i1; edge2[i] = i2; i++;
              }
            }
          }

          // On a trouve un matching edge
          if (nt == 3) // TRIANGLES
          {
            if (i == 2)
            {
              //printf("Cas 2\n");
              edge10 = edge1[0]; edge11 = edge1[1];
              edge20 = edge2[0]; edge21 = edge2[1];
              ind1 = cn(elt, edge10);
              ind2 = cn(elt, edge11);
              if (ind1 == ind2) goto next; // edge degenere -> skipped
              sign1 = edge11 - edge10;
              sign2 = edge21 - edge20;
              if (sign1 == nt-1) sign1 = -1;
              else if (sign1 == -nt+1) sign1 = +1;
              else sign1 = K_FUNC::E_sign(sign1);
              if (sign2 == nt-1) sign2 = -1;
              else if (sign2 == -nt+1) sign2 = +1;
              else sign2 = K_FUNC::E_sign(sign2);
            
              if (sign1 == sign2) // reorder
              {
                temp = cn(ie, edge20);
                cn(ie, edge20) = cn(ie, edge21);
                cn(ie, edge21) = temp;
              }
            }
            //else // elts identiques
            //{ 
            //  cn(ie, 1) = cn(elt, 1);
            //  cn(ie, 2) = cn(elt, 2);
            //  cn(ie, 3) = cn(elt, 3);
            //}
          }
          else // QUADRANGLES
          {
            if (i == 2)
            {
              //printf("Cas 2\n");
              edge10 = edge1[0]; edge11 = edge1[1];
              edge20 = edge2[0]; edge21 = edge2[1];
              ind1 = cn(elt, edge10);
              ind2 = cn(elt, edge11);
              if (ind1 == ind2) goto next; // edge degenere -> skipped
              sign1 = edge11 - edge10;
              sign2 = edge21 - edge20;
              if (sign1 == nt-1) sign1 = -1;
              else if (sign1 == -nt+1) sign1 = +1;
              else sign1 = K_FUNC::E_sign(sign1);
              if (sign2 == nt-1) sign2 = -1;
              else if (sign2 == -nt+1) sign2 = +1;
              else sign2 = K_FUNC::E_sign(sign2);
            
              if (sign1 == sign2) // reorder
              {
                SWAPEDGESQUAD;
              }
            }
            else if (i == 3)
            { // degenerescence QUAD/TRI
              //printf("Cas 3\n");
              //printf("edge1 %d %d %d\n", edge1[0], edge1[1], edge1[2]);
              //printf("edge2 %d %d %d\n", edge2[0], edge2[1], edge2[2]);
              if (edge1[0] == edge1[1])
              {
                edge20 = edge2[0]; edge21 = edge2[2];
                edge10 = edge1[0]; edge11 = edge1[2];
                if (abs(edge20-edge21) == 2)
                {
                  edge20 = edge2[1]; edge21 = edge2[2];
                  edge10 = edge1[1]; edge11 = edge1[2];
                }
              }
              else if (edge1[0] == edge1[2])
              {
                edge20 = edge2[0]; edge21 = edge2[1];
                edge10 = edge1[0]; edge11 = edge1[1];
                if (abs(edge20-edge21) == 2)
                {
                  edge20 = edge2[1]; edge21 = edge2[2];
                  edge10 = edge1[1]; edge11 = edge1[2];
                }
              }
              else if (edge1[1] == edge1[2])
              {
                edge20 = edge2[0]; edge21 = edge2[1];
                edge10 = edge1[0]; edge11 = edge1[1];
                if (abs(edge20-edge21) == 2)
                {
                  edge20 = edge2[0]; edge21 = edge2[2];
                  edge10 = edge1[0]; edge11 = edge1[2];
                }
              }
              else if (edge2[0] == edge2[1])
              {
                edge10 = edge1[0]; edge11 = edge1[2];
                edge20 = edge2[0]; edge21 = edge2[2];
                if (abs(edge10-edge11) == 2)
                {
                  edge10 = edge1[1]; edge11 = edge1[2];
                  edge20 = edge2[1]; edge21 = edge2[2];
                }
              }
              else if (edge2[0] == edge2[2])
              {
                edge10 = edge1[0]; edge11 = edge1[1];
                edge20 = edge2[0]; edge21 = edge2[1];
                if (abs(edge10-edge11) == 2)
                {
                  edge10 = edge1[1]; edge11 = edge1[2];
                  edge20 = edge2[1]; edge21 = edge2[2];
                }
              }
              else if (edge2[1] == edge2[2])
              {
                edge10 = edge1[0]; edge11 = edge1[1];
                edge20 = edge2[0]; edge21 = edge2[1];
                if (abs(edge10-edge11) == 2)
                {
                  edge10 = edge1[0]; edge11 = edge1[2];
                  edge20 = edge2[0]; edge21 = edge2[2];
                }
              }
              else // check only first edge
              {
                edge10 = edge1[0]; edge11 = edge1[1];
                edge20 = edge2[0]; edge21 = edge2[1];
              }

              //printf("final edge1 %d %d\n", edge10, edge11);
              //printf("final edge2 %d %d\n", edge20, edge21);
              sign1 = edge11 - edge10;
              sign2 = edge21 - edge20;
              if (sign1 == nt-1) sign1 = -1;
              else if (sign1 == -nt+1) sign1 = +1;
              else sign1 = K_FUNC::E_sign(sign1);
              if (sign2 == nt-1) sign2 = -1;
              else if (sign2 == -nt+1) sign2 = +1;
              else sign2 = K_FUNC::E_sign(sign2);
            
              if (sign1 == sign2) SWAPEDGESQUAD;
            }
            else if (i == 4)
            { // degenerescence TRI/TRI sur 2 noeuds
              //printf("cas 4\n");
              //printf("edge1 %d %d %d %d\n", edge1[0], edge1[1], edge1[2], edge1[3]);
              //printf("edge2 %d %d %d %d\n", edge2[0], edge2[1], edge2[2], edge2[3]);
              found = 0;
              for (i = 1; i < 4; i++)
              { if (edge1[i] == edge1[0]) found = 1; break; }
              if (found == 0)
              {
                for (i = 2; i < 4; i++)
                { if (edge1[i] == edge1[1]) found = 1; break; }
                if (found == 1) { SWAPTWOEDGES(0, 1); }
                else { SWAPTWOEDGES(0, 2); }
              }
            
              if (edge1[0] == edge1[2]) { SWAPTWOEDGES(1, 2); }
              else if (edge1[0] == edge1[3]) { SWAPTWOEDGES(1, 3); }

              if (edge1[2] == edge1[3]) goto next;

              //printf("edgetrie1 %d %d %d %d\n", edge1[0], edge1[1], edge1[2], edge1[3]);
              //printf("edgetrie2 %d %d %d %d\n", edge2[0], edge2[1], edge2[2], edge2[3]);
            
              edge20 = edge2[0]; edge21 = edge2[2];
              edge10 = edge1[0]; edge11 = edge1[2];
              if (abs(edge20-edge21) == 2)
              {
                edge20 = edge2[1]; edge21 = edge2[2];
                edge10 = edge1[1]; edge11 = edge1[2];
              }
              if (abs(edge11-edge10) == 2)
              {
                edge21 = edge2[3];
                edge11 = edge1[3];
              }

              //printf("final edge1 %d %d\n", edge10, edge11);
              //printf("final edge2 %d %d\n", edge20, edge21);
              sign1 = edge11 - edge10;
              sign2 = edge21 - edge20;
              if (sign1 == nt-1) sign1 = -1;
              else if (sign1 == -nt+1) sign1 = +1;
              else sign1 = K_FUNC::E_sign(sign1);
              if (sign2 == nt-1) sign2 = -1;
              else if (sign2 == -nt+1) sign2 = +1;
              else sign2 = K_FUNC::E_sign(sign2);
            
              if (sign1 == sign2) SWAPEDGESQUAD;
            }
            else if (i == 5)
            { // degenerescence TRI/TRI sur 1 noeuds
              //printf("cas 5\n");
              //printf("edge1 %d %d %d %d %d\n", edge1[0], edge1[1], edge1[2], edge1[3], edge1[4]);
              //printf("edge2 %d %d %d %d %d\n", edge2[0], edge2[1], edge2[2], edge2[3], edge2[4]);
              found = 0;
              for (i = 1; i < 5; i++)
              { if (edge1[i] == edge1[0]) found = 1; break; }
              if (found == 0) { SWAPTWOEDGES(0, 1); }

              if (edge1[0] == edge1[2]) { SWAPTWOEDGES(1, 2); }
             
              else if (edge1[0] == edge1[3]) { SWAPTWOEDGES(1, 3); }
            
              else if (edge1[0] == edge1[4]) { SWAPTWOEDGES(1, 4); }
            
              found = 0;
              for (i = 3; i < 5; i++)
              { if (edge1[i] == edge1[2]) found = 1; break; }
              if (found == 0) { SWAPTWOEDGES(2, 4); }

              if (edge1[2] == edge1[4]) { SWAPTWOEDGES(3, 4); }
              
              if (edge2[3] == edge2[0]) { SWAPTWOEDGES(2, 3); }
            
              //printf("edgetrie1 %d %d %d %d %d\n", edge1[0], edge1[1], edge1[2], edge1[3], edge1[4]);
              //printf("edgetrie2 %d %d %d %d %d\n", edge2[0], edge2[1], edge2[2], edge2[3], edge2[4]);
              edge20 = edge2[0]; edge21 = edge2[4];
              edge10 = edge1[0]; edge11 = edge1[4];
              if (abs(edge20-edge21) == 2 || abs(edge11-edge10) == 2)
              {
                edge20 = edge2[1]; edge21 = edge2[4];
                edge10 = edge1[1]; edge11 = edge1[4];
                if (abs(edge20-edge21) == 2 || abs(edge11-edge10) == 2)
                {
                  edge20 = edge2[2]; edge21 = edge2[4];
                  edge10 = edge1[2]; edge11 = edge1[4];
                  if (abs(edge20-edge21) == 2 || abs(edge11-edge10) == 2)
                  {
                    edge20 = edge2[3]; edge21 = edge2[4];
                    edge10 = edge1[3]; edge11 = edge1[4];
                  }
                }
              }
          
              sign1 = edge11 - edge10;
              sign2 = edge21 - edge20;
              if (sign1 == nt-1) sign1 = -1;
              else if (sign1 == -nt+1) sign1 = +1;
              else sign1 = K_FUNC::E_sign(sign1);
              if (sign2 == nt-1) sign2 = -1;
              else if (sign2 == -nt+1) sign2 = +1;
              else sign2 = K_FUNC::E_sign(sign2);
            
              if (sign1 == sign2) SWAPEDGESQUAD;
            }
            //if (mbv > ne) printf("danger\n");
            //if (ie > ne) printf("danger2\n");
          }
          mustBeVisited[mbv] = ie;
          mbv++; nev++;
          isVisited[ie] = 1;
          next: ;
        } // isvisited
      } // iv
    } // mbv
  }
  free(isVisited);
  free(mustBeVisited);

  return 1;
}

// ============================================================================
// Reorder the vertex indices of a ME connectivity such that each facet normal
// is pointing outward.
// The input connectivity must be such that the image points of the 'base' facet
// are correct. Only 3D element types are supported.
// Return
//    -1 if the element type is not supported.
//     0 if no vertex reordering took place (incl. when the coordinates are
//                                           missing in varString).
//     1 if vertex reordering took place.
// ============================================================================
E_Int K_CONNECT::reorderUnstruct3D(const char* varString, FldArrayF& f, 
                                   FldArrayI& cn, const char* eltType)
{
  // Get dimensionality - 3D elements only
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);
  E_Int dim = K_CONNECT::getDimME(eltTypes);
  if(dim < 3) return -1;

  // Check that coordinates are present
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_WarnEx(PyExc_Warning,
                 "reorderUnstruct: coords must be present in array.", 1);
    return 0;
  }
  posx++; posy++; posz++;

  E_Float* xp = f.begin(posx);
  E_Float* yp = f.begin(posy);
  E_Float* zp = f.begin(posz);

  // List of vertex indices forming each of the element's facets, for each
  // connectivity
  std::vector<std::vector<std::vector<E_Int> > > facetspc(nc);

  // List the pairs of vertex indices to swap for each connectivity if a face
  // normal is found to be pointing inward (i.e., in the element's volume)
  std::vector<std::vector<std::pair<E_Int, E_Int> > > ind2Swap(nc);
  for (E_Int ic = 0; ic < nc; ic++)
  {
    if (K_STRING::cmp(eltTypes[ic], 5, "TETRA") == 0)
      ind2Swap[ic].push_back({2, 3});
    else if (K_STRING::cmp(eltTypes[ic], 4, "PYRA") == 0)
      ind2Swap[ic].push_back({2, 4});
    else if (K_STRING::cmp(eltTypes[ic], 5, "PENTA") == 0)
    {
      ind2Swap[ic].push_back({2, 3});
      ind2Swap[ic].push_back({5, 6});
    }
    else if (K_STRING::cmp(eltTypes[ic], 4, "HEXA") == 0)
    {
      ind2Swap[ic].push_back({2, 4});
      ind2Swap[ic].push_back({6, 8});
    }

    // Fill the list of facets for this BE
    E_Int ierr = K_CONNECT::getEVFacets(facetspc[ic], eltTypes[ic]);
    if (ierr == 1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "reorderUnstruct: element type not supported.");
      for (E_Int ic = 0; ic < nc; ic++) delete [] eltTypes[ic];
      return -1;
    }
  }
  
  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];

  E_Int reoriented = 0;

  #pragma omp parallel reduction(+: reoriented)
  {
    E_Int nelts, nvpe;
    E_Int ind, ind1, ind2, ind3;
    E_Float xb, yb, zb;
    E_Float x1, y1, z1, x2, y2, z2, x3, y3, z3;
    E_Float x12, y12, z12, x13, y13, z13;
    E_Float xf, yf, zf, nx, ny, nz, xbf, ybf, zbf;
    E_Float res;
    
    // Loop over each element of all connectivities
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cn.getConnect(ic));
      nelts = cm.getSize();
      nvpe = cm.getNfld();
      const std::vector<E_Int>& facet0 = facetspc[ic][0];

      #pragma omp for
      for (E_Int i = 0; i < nelts; i++)
      {
        // Compute the position of the cell barycenter (xb, yb, zb)
        xb = 0.; yb = 0.; zb = 0.;
        for (E_Int j = 1; j <= nvpe; j++)
        {
          ind = cm(i, j)-1;
          xb += xp[ind]; yb += yp[ind]; zb += zp[ind];
        }
        xb /= nvpe; yb /= nvpe; zb /= nvpe;

        // Compute the facet normal (nx, ny, nz) using the first 3 points of the
        // first facet, facet0
        ind1 = cm(i, facet0[0]) - 1;
        ind2 = cm(i, facet0[1]) - 1;
        ind3 = cm(i, facet0[2]) - 1;
        x1 = xp[ind1]; y1 = yp[ind1]; z1 = zp[ind1];
        x2 = xp[ind2]; y2 = yp[ind2]; z2 = zp[ind2];
        x3 = xp[ind3]; y3 = yp[ind3]; z3 = zp[ind3];
        x12 = x2 - x1; x13 = x3 - x1;
        y12 = y2 - y1; y13 = y3 - y1;
        z12 = z2 - z1; z13 = z3 - z1;
        K_MATH::cross(x12, y12, z12, x13, y13, z13, nx, ny, nz);

        // Compute the face center and the vector from element barycenter to the
        // face center (xbf, ybf, zbf)
        xf = K_MATH::ONE_THIRD * (x1 + x2 + x3);
        yf = K_MATH::ONE_THIRD * (y1 + y2 + y3);
        zf = K_MATH::ONE_THIRD * (z1 + z2 + z3);
        xbf = xf - xb; ybf = yf - yb; zbf = zf - zb;

        // If dot(n, bf) is positive, the facet normal is pointing outward
        // which is the correct orientation
        K_MATH::dot(nx, ny, nz, xbf, ybf, zbf, res);

        if (res < 0.)
        {
          reoriented = true;
          for (size_t j = 0; j < ind2Swap[ic].size(); j++)
          {
            std::swap(
              cm(i, ind2Swap[ic][j].first),
              cm(i, ind2Swap[ic][j].second)
            );
          }
        }
      }
    }
  }

  if (reoriented > 0) reoriented = 1;
  return reoriented;
}
