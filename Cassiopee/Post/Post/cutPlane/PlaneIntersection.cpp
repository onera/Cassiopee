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
# include <stdio.h>
# include <stdlib.h>
# include <vector>

# include "PlaneIntersection.h"

using namespace std;
using namespace K_FLD;
using namespace K_FUNC;

extern "C"
{
  void k6compvoloftetracell_(const E_Int& nelts, 
                             const E_Int& ind1, const E_Int& ind2, 
                             const E_Int& ind3, const E_Int& ind4,
                             const E_Float* xt, const E_Float* yt, 
                             const E_Float* zt, E_Float& vol);
} 

//=============================================================================
/* Given the points defined by coord, for each segment [i,i+1] or [j,j+1] ...
   search for the intersection with the plane. If the point is found, 
   store it in plane Field */ 
//=============================================================================
short K_POST::computeStructIntersectionWithPlane( 
  K_INTERP::InterpData* interpData,
  K_INTERP::InterpData::InterpolationType interpType,
  E_Float coefa, E_Float coefb, 
  E_Float coefc, E_Float coefd, 
  E_Int ni, E_Int nj, E_Int nk,
  E_Int posx, E_Int posy, E_Int posz, E_Int posc,
  FldArrayF& field, FldArrayI& tagC,
  FldArrayF& intersectPts,
  FldArrayF& volOfIntersectPts)
{ 
  E_Int dim;
  if (nk == 1) dim = 2;
  else dim = 3;
  E_Int ninj = ni*nj;
  E_Int nfld = field.getNfld();
  E_Int npts = field.getSize();
  E_Int sizeMax = 3*npts;
  intersectPts.malloc(sizeMax, nfld);
  volOfIntersectPts.malloc(sizeMax);

  E_Int cnt = 0;
  E_Int nic = K_FUNC::E_max(1,ni-1);
  E_Int njc = K_FUNC::E_max(1,nj-1);
  E_Int nkc = K_FUNC::E_max(1,nk-1);
  E_Int nicnjc = nic*njc;

  switch (dim)
  {
    case 2:
      for (E_Int j = 0; j < njc; j++)
        for (E_Int i = 0; i < nic; i++)
        {
          E_Int indv  = i+j*nic;
          if (tagC[indv] == 1) 
          {
            E_Int indv1, indv2, indv3, indv4;
            indv1 = i + j*ni; //(i,j)
            indv2 = indv1 + 1;         //(i+1,j)   
            indv3 = indv2 + ni;        //(i+1,j+1)
            indv4 = indv3 - 1;         //(i,j+1)
            // (i,j) et (i+1,j)
            searchStructIntersectForSegment(interpData, interpType,
                                            coefa, coefb, coefc, coefd,
                                            ni, nj, nk, indv1, 
                                            posx, posy, posz,
                                            posc, cnt, indv2, field, 
                                            intersectPts, volOfIntersectPts);
            // (i,j) et (i,j+1)
            searchStructIntersectForSegment(interpData, interpType,
                                            coefa, coefb, coefc, coefd,
                                            ni, nj, nk, indv1, 
                                            posx, posy, posz,
                                            posc, cnt, indv4, field, 
                                            intersectPts, volOfIntersectPts);
            // (i+1,j) et (i+1,j+1)
            searchStructIntersectForSegment(interpData, interpType,
                                            coefa, coefb, coefc, coefd,
                                            ni, nj, nk, indv2, 
                                            posx, posy, posz,
                                            posc, cnt, indv3, field, 
                                            intersectPts, volOfIntersectPts);
            // (i+1,j+1) et (i,j+1)
            searchStructIntersectForSegment(interpData, interpType,
                                            coefa, coefb, coefc, coefd,
                                            ni, nj, nk, indv3, 
                                            posx, posy, posz,
                                            posc, cnt, indv4, field, 
                                            intersectPts, volOfIntersectPts);
          }
        }
    
    break;
    
    case 3:

//#pragma omp parallel default(shared) if (nkc > __MIN_SIZE_MEAN__)
  {    
      E_Int indv, indv1, indv2, indv3, indv4, indv5, indv6, indv7, indv8;
//#pragma omp for
      for (E_Int k = 0; k < nkc; k++)
        for (E_Int j = 0; j < njc; j++)
          for (E_Int i = 0; i < nic; i++)
          {
            indv  = i+j*nic+k*nicnjc;
            if (tagC[indv] == 1) 
            {
              indv1 = i + j*ni + k*ninj; 
              indv2 = indv1 + 1;            
              indv3 = indv2 + ni;           
              indv4 = indv3 - 1;
              indv5 = indv1 + ninj;  
              indv6 = indv2 + ninj;   
              indv7 = indv3 + ninj;  
              indv8 = indv4 + ninj;  
              
              //Warning: cnt is a counter that is modified !!!
              // OpenMP algorithm must be adapted
              // Plan k 
              searchStructIntersectForSegment(interpData, interpType,
                                              coefa, coefb, coefc, coefd,
                                              ni, nj, nk, indv1,
                                              posx, posy, posz, 
                                              posc, cnt, indv2, field, 
                                              intersectPts, volOfIntersectPts);
              searchStructIntersectForSegment(interpData, interpType,
                                              coefa, coefb, coefc, coefd,
                                              ni, nj, nk, indv2,
                                              posx, posy, posz, 
                                              posc, cnt, indv3, field, 
                                              intersectPts, volOfIntersectPts);
              searchStructIntersectForSegment(interpData, interpType,
                                              coefa, coefb, coefc, coefd,
                                              ni, nj, nk, indv3,
                                              posx, posy, posz, 
                                              posc, cnt, indv4, field, 
                                              intersectPts, volOfIntersectPts);
              searchStructIntersectForSegment(interpData, interpType,
                                              coefa, coefb, coefc, coefd,
                                              ni, nj, nk, indv4,
                                              posx, posy, posz, 
                                              posc, cnt, indv1, field, 
                                              intersectPts, volOfIntersectPts);
              // Plan k+1
              searchStructIntersectForSegment(interpData, interpType,
                                              coefa, coefb, coefc, coefd,
                                              ni, nj, nk, indv5,
                                              posx, posy, posz, 
                                              posc, cnt, indv6, field, 
                                              intersectPts, volOfIntersectPts);
              searchStructIntersectForSegment(interpData, interpType,
                                              coefa, coefb, coefc, coefd,
                                              ni, nj, nk, indv6,
                                              posx, posy, posz, 
                                              posc, cnt, indv7, field, 
                                              intersectPts, volOfIntersectPts);
              searchStructIntersectForSegment(interpData, interpType,
                                              coefa, coefb, coefc, coefd,
                                              ni, nj, nk, indv7,
                                              posx, posy, posz, 
                                              posc, cnt, indv8, field, 
                                              intersectPts, volOfIntersectPts);
              searchStructIntersectForSegment(interpData, interpType,
                                              coefa, coefb, coefc, coefd,
                                              ni, nj, nk, indv8,
                                              posx, posy, posz, 
                                              posc, cnt, indv5, field, 
                                              intersectPts, volOfIntersectPts);
              // aretes laterales
              searchStructIntersectForSegment(interpData, interpType,
                                              coefa, coefb, coefc, coefd,
                                              ni, nj, nk, indv1,
                                              posx, posy, posz, 
                                              posc, cnt, indv5, field, 
                                              intersectPts, volOfIntersectPts);
              searchStructIntersectForSegment(interpData, interpType,
                                              coefa, coefb, coefc, coefd,
                                              ni, nj, nk, indv2,
                                              posx, posy, posz, 
                                              posc, cnt, indv6, field, 
                                              intersectPts, volOfIntersectPts);
              searchStructIntersectForSegment(interpData, interpType,
                                              coefa, coefb, coefc, coefd,
                                              ni, nj, nk, indv3,
                                              posx, posy, posz, 
                                              posc, cnt, indv7, field, 
                                              intersectPts, volOfIntersectPts);
              searchStructIntersectForSegment(interpData, interpType,
                                              coefa, coefb, coefc, coefd,
                                              ni, nj, nk, indv4,
                                              posx, posy, posz, 
                                              posc, cnt, indv8, field, 
                                              intersectPts, volOfIntersectPts);
            }
          }      
  }
  break;

    default:
      printf("Warning: extractPlane: compIntersection: not a valid value for dim.\n");
      return 0;
  }

  // Resize the intersectPts array
  if (cnt > 0)
  {
    intersectPts.reAllocMat(cnt, nfld);
    volOfIntersectPts.resize(cnt);
  }
  else
  {
    intersectPts.malloc(0);
    volOfIntersectPts.malloc(0);
  }
  return 1;
}
//=============================================================================
/* Calcul de l'intersection des grilles non structurees TETRA avec le plan
   Retourne la liste des pts d intersection et le volume de la cellule 
   d'interpolation correspondante */
//=============================================================================
void K_POST::computeUnstrIntersectionWithPlane(
  E_Float coefa, E_Float coefb, E_Float coefc, E_Float coefd, 
  K_INTERP::InterpData* interpData, 
  K_INTERP::InterpData::InterpolationType interpType,
  FldArrayI& connect,
  E_Int posx, E_Int posy, E_Int posz, E_Int posc, FldArrayF& field, 
  FldArrayI& tagC,
  FldArrayF& intersectPts, FldArrayF& volOfIntersectPts)
{
  E_Int nfld = field.getNfld();
  E_Int nnodes = field.getSize();
  E_Int nelts = connect.getSize();
  E_Int sizeMax = 6*nelts;

  intersectPts.malloc(sizeMax, nfld);
  volOfIntersectPts.malloc(sizeMax);
  E_Int cnt = 0;
  
  E_Int* indt1 = connect.begin(1);// 1ers noeuds des tetras
  E_Int* indt2 = connect.begin(2);
  E_Int* indt3 = connect.begin(3);
  E_Int* indt4 = connect.begin(4);
  E_Float cellVol = 0.;

  for (E_Int et = 0; et < connect.getSize(); et++)
  {
    if (tagC[et] == 1) 
    {
      E_Int ind1 = indt1[et]-1;
      E_Int ind2 = indt2[et]-1;
      E_Int ind3 = indt3[et]-1;
      E_Int ind4 = indt4[et]-1;
          
      k6compvoloftetracell_(nnodes, ind1, ind2, ind3, ind4, 
                            field.begin(posx), field.begin(posy),
                            field.begin(posz), cellVol);
      // segment 12
      searchUnstrIntersectForSegment(coefa, coefb, coefc, coefd, 
                                     ind1, ind2, posx, posy, posz, posc,
                                     cellVol, connect, field, interpData, interpType,
                                     cnt, intersectPts, volOfIntersectPts);
      // segment 23
      searchUnstrIntersectForSegment(coefa, coefb, coefc, coefd, 
                                     ind2, ind3, posx, posy, posz, posc,
                                     cellVol, connect, field, interpData,  interpType,
                                     cnt, intersectPts, volOfIntersectPts);
      // segment 31
      searchUnstrIntersectForSegment(coefa, coefb, coefc, coefd,
                                     ind3, ind1, posx, posy, posz, posc,
                                     cellVol, connect, field, interpData,  interpType,
                                     cnt, intersectPts, volOfIntersectPts);
      // segment 41
      searchUnstrIntersectForSegment(coefa, coefb, coefc, coefd,
                                     ind4, ind1, posx, posy, posz, posc,
                                     cellVol, connect, field, interpData,  interpType,
                                     cnt, intersectPts, volOfIntersectPts);
      // segment 42
      searchUnstrIntersectForSegment(coefa, coefb, coefc, coefd, 
                                     ind4, ind2, posx, posy, posz, posc,
                                     cellVol, connect, field, interpData,  interpType,
                                     cnt, intersectPts, volOfIntersectPts);
      // segment 43
      searchUnstrIntersectForSegment(coefa, coefb, coefc, coefd, 
                                     ind4, ind3, posx, posy, posz, posc,
                                     cellVol, connect, field, interpData,  interpType,
                                     cnt, intersectPts, volOfIntersectPts);
    }
  }
  // Resize the intersectPts array
  if (cnt > 0)
  {
    intersectPts.reAllocMat(cnt, nfld);
    volOfIntersectPts.resize(cnt);
  }
  else
  {
    intersectPts.malloc(0);
    volOfIntersectPts.malloc(0);
  }
}
// //=============================================================================
// /* Given the points defined by coord, for each segment [i,i+1] or [j,j+1] ...
//    search for the intersection with the plane. If the point is found, 
//    store it in plane Field */ 
// //=============================================================================
// short K_POST::computeStructIntersectionWithPlane( 
//   K_INTERP::BlkInterpData* interpData,
//   K_INTERP::BlkInterpData::InterpolationType interpType,
//   E_Float coefa, E_Float coefb, 
//   E_Float coefc, E_Float coefd, 
//   E_Int ni, E_Int nj, E_Int nk,
//   E_Int posx, E_Int posy, E_Int posz, E_Int posc,
//   FldArrayF& field, FldArrayI& tagC,
//   FldArrayF& intersectPts,
//   FldArrayF& volOfIntersectPts)
// { 
//   E_Int dim;
//   if (nk == 1) dim = 2;
//   else dim = 3;
//   E_Int ind;
//   E_Int indp; // voisin de ind en i+1, j+1, k+1
//   E_Int ninj = ni * nj;

//   E_Int nfld = field.getNfld();
//   E_Int npts = field.getSize();
//   E_Int sizeMax = 3*npts;
//   intersectPts.malloc(sizeMax, nfld);
//   volOfIntersectPts.malloc(sizeMax);

//   E_Int cnt = 0;
  
//   switch (dim)
//   {
//     case 2:
//       for (E_Int j = 0; j < nj; j++)
//         for (E_Int i = 0; i < ni-1; i++)
//         {
//           ind = i + j * ni ;  // (i,j)
//           indp = ind + 1;  // (i+1,j)
//           searchStructIntersectForSegment(interpData, interpType,
//                                           coefa, coefb, coefc, coefd,
//                                           ni, nj, nk, ind, 
//                                           posx, posy, posz,
//                                           posc, cnt, indp, field, 
//                                           intersectPts, volOfIntersectPts);
//         }
//       for (E_Int j = 0; j < nj-1; j++)
//         for (E_Int i = 0; i < ni; i++)
//         {
//           ind = i + j * ni ;  // (i,j)
//           indp = ind + ni;  // (i,j+1)
//           searchStructIntersectForSegment(interpData, interpType,
//                                           coefa, coefb, coefc, coefd,
//                                           ni, nj, nk, ind, 
//                                           posx, posy, posz,
//                                           posc, cnt, indp, field, 
//                                           intersectPts, volOfIntersectPts);
//         }
//       break;
//     case 3:
//       for (E_Int k = 0; k < nk; k++)
//         for (E_Int j = 0; j < nj; j++)
//           for (E_Int i = 0; i < ni-1; i++)
//           {
//             ind = i + j * ni + k * ninj; // (i,j,k)
//             indp = ind + 1; // (i+1,j,k)
//             searchStructIntersectForSegment(interpData, interpType,
//                                             coefa, coefb, coefc, coefd,
//                                             ni, nj, nk, ind,
//                                             posx, posy, posz, 
//                                             posc, cnt, indp, field, 
//                                             intersectPts, volOfIntersectPts);
//           }

//       for (E_Int k = 0; k < nk; k++)
//         for (E_Int j = 0; j < nj-1; j++)
//           for (E_Int i = 0; i < ni; i++)
//           {
//             ind = i + j * ni + k * ninj; // (i,j,k)
//             indp = ind + ni; // (i,j+1,k)
//             searchStructIntersectForSegment(interpData, interpType,
//                                             coefa, coefb, coefc, coefd,
//                                             ni, nj, nk, ind,
//                                             posx, posy, posz, 
//                                             posc, cnt, indp, field, 
//                                             intersectPts, volOfIntersectPts);
//           }

//       for (E_Int k = 0; k < nk-1; k++)
//         for (E_Int j = 0; j < nj; j++)
//           for (E_Int i = 0; i < ni; i++)
//           {
//             ind = i + j * ni + k * ninj; // (i,j,k)
//             indp = ind + ninj; // (i,j,k+1)
//             searchStructIntersectForSegment(interpData, interpType,
//                                             coefa, coefb, coefc, coefd,
//                                             ni, nj, nk, ind,
//                                             posx, posy, posz, 
//                                             posc, cnt, indp, field, 
//                                             intersectPts, volOfIntersectPts);
//           }
//       break;
//     default:
//       printf("Warning: extractPlane: compIntersection: not a valid value for dim.\n");
//       return 0;
//   }
  
//   // Resize the intersectPts array
//   if (cnt > 0)
//   {
//     intersectPts.reAllocMat(cnt, nfld);
//     volOfIntersectPts.resize(cnt);
//   }
//   else
//   {
//     intersectPts.malloc(0);
//     volOfIntersectPts.malloc(0);
//   }
//   return 1;
// }
//===========================================================================
/* Etant donnes 2 pts indA et indB formant un segment de la grille non struct
   calcule l'intersection de ce segment avec le plan. Si intersection,
   cnt est incremente et le pt d'intersection est insere dans intersectPts*/
//===========================================================================
void K_POST::searchUnstrIntersectForSegment(
  E_Float coefa, E_Float coefb, E_Float coefc, E_Float coefd,
  E_Int indA, E_Int indB, E_Int posx, E_Int posy, E_Int posz, E_Int posc,
  E_Float cellVol, FldArrayI& connect, FldArrayF& field, 
  K_INTERP::InterpData* interpData, 
  K_INTERP::InterpData::InterpolationType interpType,
  E_Int& cnt, FldArrayF& intersectPts, FldArrayF& volOfIntersectPts)
{
  E_Float eps = 1.e-12;

  E_Int nfld = field.getNfld();
  E_Float* xp = field.begin(posx);
  E_Float* yp = field.begin(posy);
  E_Float* zp = field.begin(posz);

  // Test si A est dans le plan
  E_Float xA = xp[indA];
  E_Float yA = yp[indA];
  E_Float zA = zp[indA];

  // interpData non structure
  E_Int nindi = 5;
  E_Int ncf = 4;
  FldArrayI indi(nindi); FldArrayF cf(ncf);
  FldArrayI tmpIndi(nindi); FldArrayF tmpCf(ncf);
  
  E_Float res = coefa * xA + coefb * yA + coefc * zA + coefd;  
  E_Float xAB = xp[indB] - xA;
  E_Float yAB = yp[indB] - yA;
  E_Float zAB = zp[indB] - zA;
  E_Float det = coefa * xAB + coefb * yAB + coefc * zAB;

  if (fEqualZero(det, 1.e-12) == true) 
  {
    // pts A et B sont soit dans le plan ou soit dans un plan parallele
    // il suffit donc de le savoir pour le point A
    if (fEqualZero(res, 1.e-12) == true)// point A dans le plan 
    { 
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        intersectPts(cnt,   eq) = field(indA, eq);
        intersectPts(cnt+1, eq) = field(indB, eq);
      }
     
      volOfIntersectPts[cnt] = cellVol;
      volOfIntersectPts[cnt+1] = cellVol;
      cnt = cnt+2;
    }
  }
  else 
  {
    E_Float invdet = (-1.)/det;
    E_Float k = coefd + coefa * xA + coefb * yA + coefc * zA;
    k = k * invdet;
    // Point H intersection
    if (k >= -eps && k <= 1.+ eps)
    {
      if (k <= eps)
      {
        for (E_Int v = 1; v <= nfld; v++)
          intersectPts(cnt, v) = field(indA, v);
        
        // valeur du celln si celln existe dans field:
        //if (posc > 0)
        //  intersectPts(cnt, posc) = field(indA, posc);

        volOfIntersectPts[cnt] = cellVol;
        cnt++;
      }
      else if ( k >= 1.-eps)
      {
        for (E_Int v = 1; v <= nfld; v++)
          intersectPts(cnt, v) = field(indB, v);
        
        //  valeur du celln si celln existe dans field:
        //if (posc > 0)
        //  intersectPts(cnt, posc) = field(indB, posc);

        volOfIntersectPts[cnt] = cellVol;
        cnt++;
      }
      else 
      {
        E_Float k1 = 1.-k;
        E_Float xH = k1*xp[indA] + k*xp[indB];
        E_Float yH = k1*yp[indA] + k*yp[indB];
        E_Float zH = k1*zp[indA] + k*zp[indB];
        E_Int type = 0;
        E_Float voli;
        E_Int noblk = 0; 
        short found0 = K_INTERP::getInterpolationCell(xH, yH, zH, interpData,
                                                      &field, &connect, NULL, NULL, NULL,
                                                      posx, posy, posz, posc,
                                                      voli, indi, cf, tmpIndi, tmpCf, type, noblk, interpType);
        if (found0 > 0) 
        {
          K_INTERP::compInterpolatedValues(indi.begin(), cf, field, &connect, NULL, NULL,
                                           cnt, type, intersectPts);
                                        
          // coordonnees
          intersectPts(cnt, posx) = xH;
          intersectPts(cnt, posy) = yH;
          intersectPts(cnt, posz) = zH;    
          //  valeur du celln si celln existe dans field :
          if (posc > 0)
          {
            E_Float cellNA = field(indA, posc);
            E_Float cellNB = field(indB, posc);
            if (cellNA == 0. || cellNB == 0.) intersectPts(cnt, posc) = 0.;
            else intersectPts(cnt, posc) = k1*cellNA + k*cellNB;
//               intersectPts(cnt, posc) = E_max(cellNA, cellNB);
          }

          volOfIntersectPts[cnt] = cellVol;
          cnt++;
        }
        else 
        {
          for (E_Int v = 1; v <= nfld; v++)
            intersectPts(cnt,v) = k1*field(indA,v)+k*field(indB,v);

          volOfIntersectPts[cnt] = cellVol;
          cnt++;
        }
      }
    }
  }
}

//=============================================================================
/* Given 1- 2 vertices A(field(indA,.)) and Bi(field(indp[i],.)) 
         2- the plane equation coefficients coefa, coefb...
   Compute the intersection pt H coordinates of (AB) and the plane 
   If (AB) is in the plane insert A and B, 
   else insert H if and only if k in [0,1], where k is such that  AH = k.AB */
//=============================================================================
void K_POST::searchStructIntersectForSegment( 
  K_INTERP::InterpData* interpData,
  K_INTERP::InterpData::InterpolationType interpType,
  E_Float coefa, E_Float coefb, 
  E_Float coefc, E_Float coefd,
  E_Int ni, E_Int nj, E_Int nk,
  E_Int indA, 
  E_Int posx, E_Int posy, E_Int posz,
  E_Int posc, 
  E_Int& cnt, 
  E_Int indp,
  FldArrayF& field,
  FldArrayF& intersectPts,
  FldArrayF& volOfIntersectPts)
{
  E_Int inddummy = -1;  // doit rester a -1 pour K_METRIC::compVolOfStructCell3D
  FldArrayI connect(0);
  E_Float eps = 1.e-12;
  E_Int nfld = field.getNfld();

  E_Float* xp = field.begin(posx);
  E_Float* yp = field.begin(posy);
  E_Float* zp = field.begin(posz);

  // Test si A est dans le plan
  E_Float xA = xp[indA];
  E_Float yA = yp[indA];
  E_Float zA = zp[indA];

  // interpolation data
  E_Int nindi, ncf;
  switch (interpType)
  {
    case K_INTERP::InterpData::O2CF:
      ncf = 8;
      nindi = 1;
      break; 
    case K_INTERP::InterpData::O3ABC: 
      ncf = 9;
      nindi = 1;
      break;
    case K_INTERP::InterpData::O5ABC: 
      ncf = 15;
      nindi = 1;
      break;
    default:
       ncf = 8;
       nindi = 1;
       interpType = K_INTERP::InterpData::O2CF;
  }
  FldArrayI indi(nindi); FldArrayF cf(ncf);
  FldArrayI tmpIndi(nindi); FldArrayF tmpCf(ncf);

  E_Float res = coefa * xA + coefb * yA + coefc * zA + coefd;  

  /* Recherche sur un des 3 voisins en i+1, j+1 et k+1 */
  E_Int indB = indp;// voisin de indA    
  E_Float xAB = xp[indB] - xA;
  E_Float yAB = yp[indB] - yA;
  E_Float zAB = zp[indB] - zA;

  E_Float det = coefa * xAB + coefb * yAB + coefc * zAB;

  if (fEqualZero(det, 1.e-12)) 
  {
    // pts A et B sont soit dans le plan ou soit dans un plan parallele
    // il suffit donc de le savoir pour le point A
    if (fEqualZero(res, 1.e-12))// point A dans le plan 
    { 
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        intersectPts(cnt,   eq) = field(indA, eq);
        intersectPts(cnt+1, eq) = field(indB, eq);
      }
      K_METRIC::compVolOfStructCell3D(
        ni, nj, nk, inddummy, indA, 
        xp, yp, zp, volOfIntersectPts[cnt]);
      K_METRIC::compVolOfStructCell3D(
        ni, nj, nk, inddummy, indB, 
        xp, yp, zp, volOfIntersectPts[cnt+1]);
      cnt = cnt+2;
    }
  }
  else 
  {
    E_Float invdet = (-1.)/det;
    E_Float k = coefd + coefa * xA + coefb * yA + coefc * zA;
    k = k * invdet;
    // Point H intersection
    if (k >= -eps && k <= 1.+ eps)
    {
      if (k <= eps)
      {
        for (E_Int v = 1; v <= nfld; v++)
          intersectPts(cnt,v) = field(indA,v);
        
        // valeur du celln si celln existe dans field :
        //if (posc != 0)
        //  intersectPts(cnt, posc) = field(indA, posc);
        
        K_METRIC::compVolOfStructCell3D(
          ni, nj, nk, inddummy, indA, 
          xp, yp, zp, volOfIntersectPts[cnt]);
        cnt++;
      }
      else if (k >= 1.-eps)
      {
        for (E_Int v = 1; v <= nfld; v++)
          intersectPts(cnt,v) = field(indB,v);
        
        //  valeur du celln si celln existe dans field :
        //if (posc != 0)
        //  intersectPts(cnt, posc) = field(indB, posc);
        
        K_METRIC::compVolOfStructCell3D(
          ni, nj, nk, inddummy, indB, 
          xp, yp, zp, volOfIntersectPts[cnt]);
        cnt++;
      }
      else 
      {
        E_Float k1 = 1.-k;
        E_Float xH = k1*xp[indA] + k*xp[indB];
        E_Float yH = k1*yp[indA] + k*yp[indB];
        E_Float zH = k1*zp[indA] + k*zp[indB];
        E_Float voli;
        E_Int type = 0, noblk = 0;
        short found0 = K_INTERP::getInterpolationCell(
          xH, yH, zH, interpData,
          &field, &ni, &nj, &nk, NULL,
          posx, posy, posz, posc,
          voli, indi, cf, tmpIndi, tmpCf, type, noblk, interpType);
        if (found0 > 0) 
        {          
          K_INTERP::compInterpolatedValues(
            indi.begin(), cf, field, &ni, &nj, &nk, cnt, type,intersectPts);
          
          // coordonnees
          intersectPts(cnt,posx) = xH;
          intersectPts(cnt,posy) = yH;
          intersectPts(cnt,posz) = zH;
          
          //  valeur du celln si celln existe dans field :
          if (posc != 0)
          {
            E_Float cellNA = field(indA, posc);
            E_Float cellNB = field(indB, posc);
            if (cellNA == 0. || cellNB == 0.) intersectPts(cnt, posc) = 0.;
            else intersectPts(cnt, posc) = k1*cellNA + k*cellNB;
//               intersectPts(cnt, posc) = E_max(cellNA, cellNB);
          }
          
          K_METRIC::compVolOfStructCell3D(
            ni, nj, nk, inddummy, indA,
            xp, yp, zp, volOfIntersectPts[cnt]);
          cnt++;
        }
        else 
        {
          for (E_Int v = 1; v <= nfld; v++)
            intersectPts(cnt,v) = k1*field(indA,v)+k*field(indB,v);
          K_METRIC::compVolOfStructCell3D(
            ni, nj, nk, inddummy, indA, 
            xp, yp, zp, volOfIntersectPts[cnt]);
          cnt++;  
        }
      }
    }
  }
}

//=================== Post/cutPlane/PlaneIntersection.cpp ====================
