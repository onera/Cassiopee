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

// Routines principales: celles qui sont appelees de extractPlane.cpp

# include "cutPlane.h"

using namespace std;
using namespace K_FLD;
using namespace K_FUNC;

//===========================================================================
/* Calcul des intersections des grilles  avec le plan  */
//===========================================================================
void K_POST::compIntersectionWithPlane(
  E_Float coefa, E_Float coefb, E_Float coefc, E_Float coefd,
  vector<K_INTERP::InterpData*>& structInterpDatas,
  vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks,
  vector<E_Int>& posxs, vector<E_Int>& posys, vector<E_Int>& poszs,
  vector<E_Int>& poscs,
  vector<FldArrayF*>& structFields, vector<FldArrayI*>& tagS, 
  vector<K_INTERP::InterpData*>& unstrInterpDatas,
  vector<FldArrayI*>& connectu,
  vector<E_Int>& posxu, vector<E_Int>& posyu, vector<E_Int>& poszu,
  vector<E_Int>& poscu,
  vector<FldArrayF*>& unstrFields, vector<FldArrayI*>& tagU, 
  vector<FldArrayF*>& vectOfIntersectPts,
  vector<FldArrayF*>& vectOfInterpCellVol,
  K_INTERP::InterpData::InterpolationType interpType)
{
  E_Int nzones = structInterpDatas.size();
  
  for (E_Int zone = 0; zone < nzones; zone++)
  {
    FldArrayF* intersectPts = new FldArrayF();
    FldArrayF* volOfIntersectPts = new FldArrayF();
    short ok = computeStructIntersectionWithPlane(
      structInterpDatas[zone], interpType, coefa, coefb, coefc, coefd,
      nis[zone], njs[zone], nks[zone], posxs[zone], posys[zone], 
      poszs[zone], poscs[zone], *structFields[zone], *tagS[zone],
      *intersectPts, *volOfIntersectPts);
    if (ok == 0)
    {
      printf("Warning: extractPlane: zone skipped.\n");
    }
    else 
    {
      vectOfIntersectPts.push_back(intersectPts);
      vectOfInterpCellVol.push_back(volOfIntersectPts);  
    }
  }

  nzones = unstrInterpDatas.size();
  for (E_Int zone = 0; zone < nzones; zone++)
  {
    FldArrayF* intersectPts = new FldArrayF();
    FldArrayF* volOfIntersectPts = new FldArrayF();
    computeUnstrIntersectionWithPlane(
      coefa, coefb, coefc, coefd, unstrInterpDatas[zone], interpType, *connectu[zone],
      posxu[zone], posyu[zone], poszu[zone], poscu[zone], 
      *unstrFields[zone], *tagU[zone], *intersectPts, *volOfIntersectPts);
    vectOfIntersectPts.push_back(intersectPts);
    vectOfInterpCellVol.push_back(volOfIntersectPts);
  }
}
//=============================================================================
/* Select points from overlapping zones: selected one is the one whose cell 
   volume is the smallest one
   Ici l ordre de la liste vectOfIntersectPts correspond a l ordre de la liste 
   des structFields+unstrFields
 */
//=============================================================================
void K_POST::selectPointsInOverlappingZones( 
  vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks,
  vector<E_Int>& posxs, vector<E_Int>& posys, 
  vector<E_Int>& poszs, vector<E_Int>& poscs,
  vector<K_INTERP::InterpData*>& structInterpDatas,
  vector<FldArrayF*>& structFields,
  vector<FldArrayI*>& connectu,
  vector<E_Int>& posxu, vector<E_Int>& posyu, 
  vector<E_Int>& poszu, vector<E_Int>& poscu,
  vector<K_INTERP::InterpData*>& unstrInterpDatas,
  vector<FldArrayF*>& unstrFields,
  vector<FldArrayF*>& vectOfIntersectPts,
  vector<FldArrayF*>& volOfIntersectPts,
  FldArrayF& selectedPts,
  K_INTERP::InterpData::InterpolationType interpType)
{
  E_Int nzones = vectOfIntersectPts.size();

  E_Int nfld = 0;
  E_Int vectSize = vectOfIntersectPts.size();
  for (E_Int v = 0; v < vectSize; v++)
  {
    if (vectOfIntersectPts[v]->getSize() != 0) 
    {
      nfld = E_max(nfld, vectOfIntersectPts[v]->getNfld());
    }
  }
  if (nzones == 0 || nfld == 0)
  {
    selectedPts.malloc(0);
    return;
  }
  // interpolation data
  E_Int nindi, ncf;
  switch (interpType)
  {
    case K_INTERP::InterpData::O2CF:
      ncf = 8;
      nindi = 1;
      //order = 2;
      break; 
    case K_INTERP::InterpData::O3ABC: 
      ncf = 9;
      nindi = 1;
      //order = 3;
      break;
    case K_INTERP::InterpData::O5ABC: 
      ncf = 15;
      nindi = 1;
      //order = 5;
      break;
    default:
      ncf = 8;
      nindi = 1;
      //order = 2;
      interpType = K_INTERP::InterpData::O2CF;
      break;
  }

  if (structFields.size()== 0) //purement non structure 
  { nindi = 1; ncf = 4; }
 
  FldArrayI indi(nindi);
  FldArrayF cf(ncf);
  FldArrayI tmpIndi(nindi); FldArrayF tmpCf(ncf);
  
  /* Compute sizeMax of selectedPts and allocate it */
  E_Float vol1, vol2;
  E_Int sizeMax = 0;
  for (E_Int zone1 = 0; zone1 < nzones; zone1++)
    sizeMax += vectOfIntersectPts[zone1]->getSize();  
  selectedPts.malloc(sizeMax, nfld);
  E_Float* xt = selectedPts.begin(1);
  E_Float* yt = selectedPts.begin(2);
  E_Float* zt = selectedPts.begin(3);
  E_Float* cellnt = selectedPts.begin(nfld);

  E_Int cnt = 0;
  /* Select best points in overlapping zones and keep them elsewhere*/
  E_Int nzoneStruct = structFields.size();
  E_Int nzoneUnstr = unstrFields.size();
  // par hypothese (construction ds getFromArrays): posx, posy, idem 
  // pour tous les blocs
  E_Int posx1, posy1, posz1, posc1;
  if (nzoneStruct > 0)
  {posx1 = posxs[0]; posy1 = posys[0]; posz1 = poszs[0]; posc1 = poscs[0];}
  else if (nzoneUnstr > 0)
  {posx1 = posxu[0]; posy1 = posyu[0]; posz1 = poszu[0]; posc1 = poscu[0];} 
  else {selectedPts.malloc(0); return;}

  vector<K_INTERP::InterpData*> allInterpDatas;
  vector<FldArrayF*> allFields;
  vector<void*> allA1;
  vector<void*> allA2;
  vector<void*> allA3;
  vector<void*> allA4;
  E_Int nzonesTot = posxs.size()+posxu.size();
  vector<E_Int> posxt; vector<E_Int> posyt; vector<E_Int> poszt; vector<E_Int> posct;
  posxt.reserve(nzonesTot);  // preallocate memory
  posxt.insert(posxt.end(), posxs.begin(), posxs.end()); 
  posxt.insert(posxt.end(), posxu.begin(), posxu.end()); 
  posyt.reserve(nzonesTot);  // preallocate memory
  posyt.insert(posyt.end(), posys.begin(), posys.end()); 
  posyt.insert(posyt.end(), posyu.begin(), posyu.end()); 
  poszt.reserve(nzonesTot);  // preallocate memory
  poszt.insert(poszt.end(), poszs.begin(), poszs.end()); 
  poszt.insert(poszt.end(), poszu.begin(), poszu.end());
  posct.reserve(nzonesTot);  // preallocate memory
  posct.insert(posct.end(), poscs.begin(), poscs.end()); 
  posct.insert(posct.end(), poscu.begin(), poscu.end());

  allFields.reserve(nzonesTot);  // preallocate memory
  allFields.insert(allFields.end(), structFields.begin(), structFields.end()); 
  allFields.insert(allFields.end(), unstrFields.begin(), unstrFields.end());

  allInterpDatas.reserve(nzonesTot);  // preallocate memory
  allInterpDatas.insert(allInterpDatas.end(), structInterpDatas.begin(), structInterpDatas.end()); 
  allInterpDatas.insert(allInterpDatas.end(), unstrInterpDatas.begin(), unstrInterpDatas.end());

  for (unsigned int noz = 0; noz < structFields.size(); noz++)
  {
    allA1.push_back(&nis[noz]); 
    allA2.push_back(&njs[noz]); 
    allA3.push_back(&nks[noz]); 
    allA4.push_back(NULL);
  }
  for (unsigned int noz = 0; noz < unstrFields.size(); noz++)
  {
    allA1.push_back(connectu[noz]); 
    allA2.push_back(NULL); 
    allA3.push_back(NULL); 
    allA4.push_back(NULL);
  }
  for (E_Int zone = 0; zone < nzones; zone++)
  {
    E_Int ninterpPts = vectOfIntersectPts[zone]->getSize();
        
    if (ninterpPts != 0)
    {
      FldArrayF& field1 = *(vectOfIntersectPts[zone]);
      E_Float* volArray = volOfIntersectPts[zone]->begin();
      E_Float* field1x = field1.begin(posx1);
      E_Float* field1y = field1.begin(posy1);
      E_Float* field1z = field1.begin(posz1);
      E_Float* field1c = NULL;
      E_Float x, y, z, celln1;
      E_Int noblk = 0;
      short found;
      if (posc1 > 0) field1c = field1.begin(posc1);
      for (E_Int ind1 = 0; ind1 < field1.getSize(); ind1++)
      {
        celln1 = 1.;
        if (posc1 != 0) celln1 = field1c[ind1];
        
        if (celln1 != 0.) // cas: pas de celln + si celln pas masque
        {
          x = field1x[ind1]; y = field1y[ind1]; z = field1z[ind1];
          E_Int type = 0;
          found = K_INTERP::getInterpolationCell(x, y, z, allInterpDatas,
                                                 allFields, allA1, allA2, allA3, allA4,
                                                 posxt, posyt, poszt, posct, 
                                                 vol2, indi, cf, tmpIndi, tmpCf, type, noblk, interpType);
          if (found > 0) // interpolable 
          {          
            noblk = noblk-1;//no du bloc d interpolation             
            if (noblk != zone)// interpolable depuis un autre bloc
            {
              vol1 = volArray[ind1]; // volume of current point's cell
              if (vol1 <= vol2 + 1.e-7)
              {
                E_Int eqn = 4;
                for (E_Int eq = 1; eq <= nfld; eq++)
                { 
                  if (eq != posx1 && eq != posy1 && 
                      eq != posz1 && eq != posc1)
                  {
                    selectedPts(cnt, eqn) = field1(ind1, eq);
                    eqn++;
                  }
                }
                xt[cnt] = field1x[ind1];
                yt[cnt] = field1y[ind1];
                zt[cnt] = field1z[ind1];
                if (posc1 != 0) cellnt[cnt] = field1c[ind1];//last: celln
                
                cnt++;
              }
            }
            else // interpolable depuis son propre bloc: ok
            {
              E_Int eqn = 4;
              for (E_Int eq = 1; eq <= nfld; eq++)
              {
                if (eq != posx1 && eq != posy1 && 
                    eq != posz1 && eq != posc1)
                {
                  selectedPts(cnt, eqn) = field1(ind1, eq);
                  eqn++;
                }
              }
              
              xt[cnt] = field1x[ind1];
              yt[cnt] = field1y[ind1];
              zt[cnt] = field1z[ind1];
              if (posc1 != 0) cellnt[cnt] = field1c[ind1];//last: celln
              cnt++;
            } 
          }
        }
      } //end of ind1
    }
  } //end of zone 
  /* Reallocate selectedPts array */
  if (cnt > 0) selectedPts.reAllocMat(cnt, nfld);
  else  selectedPts.malloc(0);
  allInterpDatas.clear(); allFields.clear(); allA1.clear(); allA2.clear(); allA3.clear(); allA4.clear();
  posxt.clear(); posyt.clear(); poszt.clear(); posct.clear();
}
//===========================================================================
/*
  IN: coefa, coefbm coefc, coefd : plane definition.
  IN: field: the field containing the coordinate of points to be 
                triangulated with the numerical solution.
  OUT: connect: the triangle connectivity. */
//===========================================================================
void K_POST::makeTriangulation(
  E_Float coefa, E_Float coefb,E_Float coefc, E_Float coefd,
  vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks,
  vector<E_Int>& posxs, vector<E_Int>& posys, vector<E_Int>& poszs, 
  vector<E_Int>& poscs, vector<FldArrayF*>& structF,
  vector<K_INTERP::InterpData*>& structInterpDatas,
  vector<FldArrayI*>& connectu, 
  vector<E_Int>& posxu, vector<E_Int>& posyu, vector<E_Int>& poszu, 
  vector<E_Int>& poscu, vector<FldArrayF*>& unstrF,
  std::vector<K_INTERP::InterpData*>& unstrInterpDatas,
  FldArrayF& field, FldArrayI& connect,
  K_INTERP::InterpData::InterpolationType interpType)
{
  K_COMPGEOM::delaunay(coefa, coefb, coefc, coefd, field,connect, 0);

  /* Remove triangles with one vertex with celln = 0 */
  //  if ( cellN != -1 )
  //    removeTrianglesWithBlankedVertices(field, triangles);
  removeTrianglesWithBlankedCenters(
    nis, njs, nks, posxs, posys, poszs, poscs, structF, structInterpDatas, 
    connectu, posxu, posyu, poszu, poscu, unstrF, unstrInterpDatas,
    field, connect, interpType);
}

//========================================================================
/* test if the barycenter of the triangle is interpolated or not 
   if no interpolation point with cellN > 0 is found the triangle is
   blanked */
//========================================================================
void K_POST::removeTrianglesWithBlankedCenters(
  vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks,
  vector<E_Int>& posxs, vector<E_Int>& posys, vector<E_Int>& poszs, 
  vector<E_Int>& poscs, vector<FldArrayF*>& structF,
  vector<K_INTERP::InterpData*>& structInterpDatas,
  vector<FldArrayI*>& connectu,
  vector<E_Int>& posxu, vector<E_Int>& posyu, vector<E_Int>& poszu,  
  vector<E_Int>& poscu, vector<FldArrayF*>& unstrF,
  vector<K_INTERP::InterpData*>& unstrInterpDatas,
  FldArrayF& field, FldArrayI& connect,
  K_INTERP::InterpData::InterpolationType interpType)
{
  // interpolation data
  E_Int nindi, ncf;
  switch (interpType)
  {
    case K_INTERP::InterpData::O2CF:
      ncf = 8; nindi = 1; //order = 2;
      break; 
    case K_INTERP::InterpData::O3ABC: 
      ncf = 9; nindi = 1; //order = 3;
      break;
    case K_INTERP::InterpData::O5ABC: 
      ncf = 15; nindi = 1; //order = 5;
      break;
    default:
        ncf = 8; nindi = 1; //order = 2;
       interpType = K_INTERP::InterpData::O2CF;
  }
  if (structInterpDatas.size() == 0) // purement non structure
    {ncf = 4; nindi = 1; } //order = 2;
  FldArrayI indi(nindi);
  FldArrayF cf(ncf);
  FldArrayI tmpIndi(nindi); FldArrayF tmpCf(ncf);

  E_Float inv = 1./3;
  E_Float* x = field.begin(1);
  E_Float* y = field.begin(2);
  E_Float* z = field.begin(3);

  E_Int ind;
  E_Float cellN;
  
  E_Int nelts = connect.getSize();
  E_Int* cn1 = connect.begin(1);
  E_Int* cn2 = connect.begin(2);
  E_Int* cn3 = connect.begin(3);
  FldArrayI connect2(nelts, 3);
  FldArrayIS tag(nelts); tag.setAllValuesAtNull();
  short* tagp = tag.begin();
  E_Int ind1, ind2, ind3, noblk = 0;
  E_Float vol2;
 

  vector<K_INTERP::InterpData*> allInterpDatas;
  vector<FldArrayF*> allFields;
  vector<void*> allA1;
  vector<void*> allA2;
  vector<void*> allA3;
  vector<void*> allA4;
  E_Int nzonesTot = posxs.size()+posxu.size();
  vector<E_Int> posxt; vector<E_Int> posyt; vector<E_Int> poszt; vector<E_Int> posct;
  posxt.reserve(nzonesTot);  // preallocate memory
  posxt.insert(posxt.end(), posxs.begin(), posxs.end()); 
  posxt.insert(posxt.end(), posxu.begin(), posxu.end()); 
  posyt.reserve(nzonesTot);  // preallocate memory
  posyt.insert(posyt.end(), posys.begin(), posys.end()); 
  posyt.insert(posyt.end(), posyu.begin(), posyu.end()); 
  poszt.reserve(nzonesTot);  // preallocate memory
  poszt.insert(poszt.end(), poszs.begin(), poszs.end()); 
  poszt.insert(poszt.end(), poszu.begin(), poszu.end());
  posct.reserve(nzonesTot);  // preallocate memory
  posct.insert(posct.end(), poscs.begin(), poscs.end()); 
  posct.insert(posct.end(), poscu.begin(), poscu.end());

  allFields.reserve(nzonesTot);  // preallocate memory
  allFields.insert(allFields.end(), structF.begin(), structF.end()); 
  allFields.insert(allFields.end(), unstrF.begin(), unstrF.end());

  allInterpDatas.reserve(nzonesTot);  // preallocate memory
  allInterpDatas.insert(allInterpDatas.end(), structInterpDatas.begin(), structInterpDatas.end()); 
  allInterpDatas.insert(allInterpDatas.end(), unstrInterpDatas.begin(), unstrInterpDatas.end());

  for (unsigned int noz = 0; noz < structF.size(); noz++)
  {
    allA1.push_back(&nis[noz]); 
    allA2.push_back(&njs[noz]); 
    allA3.push_back(&nks[noz]); 
    allA4.push_back(NULL);
  }
  for (unsigned int noz = 0; noz < unstrF.size(); noz++)
  {
    allA1.push_back(connectu[noz]); 
    allA2.push_back(NULL); 
    allA3.push_back(NULL); 
    allA4.push_back(NULL);
  }
  for (E_Int et = 0; et < nelts; et++)
  {
    ind1 = cn1[et]-1;
    ind2 = cn2[et]-1;
    ind3 = cn3[et]-1;

    E_Float xc = inv * (x[ind1] + x[ind2] + x[ind3]);
    E_Float yc = inv * (y[ind1] + y[ind2] + y[ind3]);
    E_Float zc = inv * (z[ind1] + z[ind2] + z[ind3]);

    E_Bool isInterp = false;
    E_Int type = 0;
    short found = K_INTERP::getInterpolationCell(xc, yc, zc, allInterpDatas,
                                                 allFields, allA1, allA2, allA3, allA4,
                                                 posxt, posyt, poszt, posct, 
                                                 vol2, indi, cf, tmpIndi, tmpCf, type, noblk, interpType);
//     short found = K_INTERP::getInterpolationCell(
//       xc, yc, zc, structInterpDatas, structF, nis, njs, nks, 
//       posxs, posys, poszs, poscs, unstrInterpDatas, unstrF,
//       connectu, posxu, posyu, poszu, poscu, indi, cf, tmpIndi, tmpCf,
//       interpMeshType, interpType);

    cellN = 1.;
    if (found > 0)
    { 
      noblk = noblk-1;
      E_Int size = structInterpDatas.size();
      if (noblk < size) // structure
      {
        E_Int posc = poscs[noblk];
        if (posc > 0) 
        {
          E_Float* cellN2 = structF[noblk]->begin(posc);
          E_Int ni = nis[noblk]; E_Int nj = njs[noblk]; E_Int ninj = ni*nj;
          E_Int ind0 = indi[0];
          E_Int k = ind0/(ninj);  E_Int j = (ind0-k*ninj)/ni; E_Int i = (ind0-j*ni-k*ninj);
          switch (type)
          {
            case 2:              
              for (E_Int k0 = 0; k0 < 2; k0++)
                for (E_Int j0 = 0; j0 < 2; j0++)
                  for (E_Int i0 = 0; i0 < 2; i0++)
                  {
                    ind = (i+i0) + (j+j0)*ni + (k+k0)*ninj;
                    cellN *= cellN2[ind];
                  }
              break;
            
            case 3:
            case 5:              
              for (E_Int k0 = 0; k0 < type; k0++)
                for (E_Int j0 = 0; j0 < type; j0++)
                  for (E_Int i0 = 0; i0 < type; i0++)
                  {
                    ind = (i+i0)+(j+j0)*ni+(k+k0)*ninj;
                    cellN *= cellN2[ind];  
                  }
              break;
  
            default:
              printf("cutPlane:removeTrianglesWithBlankedCenters : not a valid interpolation type\n");
              exit(0);
                
          }
        }
      }
      else // donneur non structure
      {
        E_Int posc = poscu[noblk];
        if (posc > 0) 
        {
          E_Float* cellN2 = unstrF[noblk]->begin(posc);
          FldArrayI& cEV = *connectu[noblk];
          switch (type)
          {
            case 4:   
              E_Int noet = indi[0];
              for (E_Int nov = 1; nov <= cEV.getNfld(); nov++)
              {
                ind = cEV(noet,nov)-1;
                cellN *= cellN2[ind];
              }             
          }
        }
      }      
      if (cellN != 0.) isInterp = true;
    }    
    if (isInterp == false) tagp[et] = 1;
  }
  
  E_Int c = 0;
  E_Int* c2p1 = connect2.begin(1);
  E_Int* c2p2 = connect2.begin(2);
  E_Int* c2p3 = connect2.begin(3);
  for (E_Int et = 0; et < nelts; et++)
  {
    if (tagp[et] != 1)
    {
      c2p1[c] = cn1[et];
      c2p2[c] = cn2[et];
      c2p3[c] = cn3[et];
      c++;
    }
  }
  connect2.reAllocMat(c,3);
  connect = connect2;
  allInterpDatas.clear(); allFields.clear(); allA1.clear(); allA2.clear(); allA3.clear(); allA4.clear();
  posxt.clear(); posyt.clear(); poszt.clear(); posct.clear();
}
// //=============================================================================
// /* Remove triangles with one or more blanked (masked) vertex  */
// //=============================================================================
// void K_POST::removeTrianglesWithBlankedVertices(FldArrayF& field, 
//                                         list<Triangle*>& triangles)
// {
//   E_Int nfld = field.getNfld();
//   E_Int ind1, ind2, ind3;
//   E_Float cellN1, cellN2, cellN3;
//   list<Triangle*>::iterator itr = triangles.begin();
//   list<Triangle*>::iterator itr2;

//   E_Float* cellNF = field.begin(nfld);

//   while ( itr != triangles.end())
//   {
//     (*itr)->getVertexIndices(ind1, ind2, ind3);

//     // Get cell nature of each vertex
//     cellN1 = cellNF[ind1];
//     cellN2 = cellNF[ind2];
//     cellN3 = cellNF[ind3];

//     if ( cellN1 * cellN2 * cellN3 == 0 )
//     {
//       itr2 = itr;
//       itr++;
//       delete *itr2;
//       triangles.erase(itr2);
//     }
//     else 
//       itr++;
//   }
// }
//=================== Post/cutPlane/cutPlane.cpp ==========================
