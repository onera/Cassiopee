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

// routines communes a streamLine et streamRibbon

# include <string.h>
# include <stdio.h>
# include <math.h>
# include "stream.h"

using namespace K_FLD;
using namespace std;

extern "C"
{
  void k6compmeanlengthofstructcell_(const E_Int& ni, const E_Int& nj, 
                                     const E_Int& nk, const E_Int& indA,
                                     const E_Float* xt, const E_Float* yt,
                                     const E_Float* zt, E_Float& meanl);
  
  void k6compmeanlengthoftetracell_(E_Int& npts, E_Int& indA, E_Int& indB,
                                   E_Int& indC, E_Int& indD, 
                                   E_Float* xt, E_Float* yt, E_Float* zt, 
                                   E_Float& meanl);  
}

//=============================================================================
/* Calcul du pas "de temps" initial pour le calcul de la ligne de courant
   noblk demarre a 1 */
//=============================================================================
void K_POST::compInitialStep(
  E_Int noblk, E_Int type, FldArrayI& indi, 
  vector<E_Int>& nis, vector<E_Int>& njs, 
  vector<E_Int>& nks, vector<E_Int>& posxs, 
  vector<E_Int>& posys, vector<E_Int>& poszs,
  vector<FldArrayF*>& structFields, 
  vector<FldArrayF*>& structVelocities,
  vector<E_Int>& posxu, vector<E_Int>& posyu, 
  vector<E_Int>& poszu,  
  vector<FldArrayF*>& unstrFields,
  vector<FldArrayI*>& connectu,
  vector<FldArrayF*>& unstrVelocities,
  E_Float& dt0)
{
  E_Int ns = structFields.size();
  E_Int noblk0 = noblk-1;
  
  E_Float l0 = 0.;
  E_Float u0 = K_CONST::E_MAX_FLOAT;
  if (noblk0 < ns)// structure
  {
    if (type != 2 && type != 3 && type != 5 )
    {
      printf("Error: stream: compInitialStep: not a valid interpolation type: " SF_D_ ".\n", type);
      exit(0);
    }
    FldArrayF* field = structFields[noblk0];
    FldArrayF* velo = structVelocities[noblk0];
    // Calcul de la vitesse de la cellule
    E_Float* vx = velo->begin(1);
    E_Float* vy = velo->begin(2);
    E_Float* vz = velo->begin(3);
    E_Int ind0 = indi[0];
    E_Int ni = nis[noblk0]; 
    E_Int nj = njs[noblk0];
    E_Int nk = nks[noblk0];
    E_Int posx = posxs[noblk0];
    E_Int posy = posys[noblk0];
    E_Int posz = poszs[noblk0];
    // Longueur de reference de la cellule d'interpolation
    k6compmeanlengthofstructcell_(ni, nj, nk, ind0, field->begin(posx),
                                  field->begin(posy), field->begin(posz), 
                                  l0);
    u0 = vx[ind0]*vx[ind0] + vy[ind0]*vy[ind0] + vz[ind0]*vz[ind0]; 
  }
  else
  {
    if (type != 4) 
    {
      printf("Error: stream: compInitialStep: not a valid interp type: " SF_D_ ".\n", type);
      exit(0);
    }
    noblk0 = noblk0-ns;// numero du bloc reel dans la liste des blocs non structures
    FldArrayI& cnEV = *(connectu[noblk0]);
    if (cnEV.getNfld() != 4)
    {
      printf("Error: stream: compInitialStep: unstructured zone must be TETRA.\n");
      exit(0);
    }
    FldArrayF* field = unstrFields[noblk0];
    FldArrayF* velo = unstrVelocities[noblk0];
    // Calcul de la vitesse de la cellule
    E_Float* vx = velo->begin(1);
    E_Float* vy = velo->begin(2);
    E_Float* vz = velo->begin(3);
    E_Int posx = posxu[noblk0];
    E_Int posy = posyu[noblk0];
    E_Int posz = poszu[noblk0];

    E_Int npts = field->getSize();
    E_Int noet = indi[0];
    E_Int indA = cnEV(noet,1)-1;
    E_Int indB = cnEV(noet,2)-1;
    E_Int indC = cnEV(noet,3)-1;
    E_Int indD = cnEV(noet,4)-1;

    k6compmeanlengthoftetracell_( npts, indA, indB, indC, indD, 
                                 field->begin(posx), field->begin(posy), field->begin(posz), l0);
    u0 = 0.;
    for (E_Int v = 1; v <= 4; v++)
    {
      E_Int ind = cnEV(noet,v)-1;
      E_Float u0x = vx[ind];
      E_Float u0y = vy[ind];
      E_Float u0z = vz[ind];
      u0 = u0 + u0x*u0x+ u0y*u0y + u0z*u0z;
    }
    u0 = u0 * 0.25;
  }

  u0 = sqrt(u0);
  u0 = K_FUNC::E_max(u0, 1.e-13);
  dt0 = l0/u0;
}
//=========================================================================
/* Calcul des champs du point (x0,y0,z0) par interpolation sur le bloc
   noblk (structure ou non), et des donnees indi, cf
   Mise a jour dans streamPt */
//=========================================================================
void K_POST::compStreamPtFields(
  E_Int nopt, E_Float x0, E_Float y0, E_Float z0, 
  E_Int noblk, E_Int type, FldArrayI& indi, FldArrayF& cf,
  vector<E_Int>& nis, vector<E_Int>& njs, 
  vector<E_Int>& nks, vector<E_Int>& posxs, 
  vector<E_Int>& posys, vector<E_Int>& poszs,
  vector<FldArrayF*>& structFields,
  vector<E_Int>& posxu, vector<E_Int>& posyu, 
  vector<E_Int>& poszu,  
  vector<FldArrayF*>& unstrFields, vector<FldArrayI*>& connectu,
  FldArrayF& streamPt,
  K_INTERP::InterpData::InterpolationType interpType)
{
  E_Int ns = structFields.size();
  //E_Int nu = unstrFields.size();
  E_Int noblk0 = noblk-1;
  if (noblk0 < ns) // structure
  {
    E_Int ni = nis[noblk0];
    E_Int nj = njs[noblk0];
    E_Int nk = nks[noblk0];
    K_INTERP::compInterpolatedValues(indi.begin(), cf, *structFields[noblk0],
                                     &ni, &nj, &nk, nopt, type, streamPt);
  }
  else
  {
    noblk0 = noblk0-ns;// numero du bloc reel dans la liste des blocs non structures
    FldArrayI& cnEV = *connectu[noblk0];
    if (cnEV.getNfld() != 4)
    {
      printf("Error: stream: compInitialStep: unstructured zone must be TETRA.\n");
      exit(0);
    }

    K_INTERP::compInterpolatedValues(indi.begin(), cf, *unstrFields[noblk0], 
                                     connectu[noblk0], NULL, NULL,
                                     nopt, type, streamPt);
  }

  //mise a jour des coordonnees (non interpolees...)
  streamPt(nopt,1) = x0; streamPt(nopt,2) = y0; streamPt(nopt,3) = z0;
}
//=============================================================================
// Calcul des coefficients de Runge-Kutta RK4 xn = xp + dt/6 * (k1+k2+k3+k4)
// retourne 0 si un sous-pas ne s est pas passe correctement
//=============================================================================
short K_POST::compRungeKutta4(
  E_Float xp, E_Float yp, E_Float zp,
  E_Float up, E_Float vp, E_Float wp, 
  E_Float& dt, E_Float& xn, E_Float& yn, E_Float& zn,
  vector<K_INTERP::InterpData*>& listOfStructInterpData,
  vector<FldArrayF*>& listOfStructFields,
  vector<FldArrayF*>& listOfStructVelocities,
  vector<E_Int>& nis,vector<E_Int>& njs,
  vector<E_Int>& nks, vector<E_Int>& posxs, 
  vector<E_Int>& posys, vector<E_Int>& poszs, 
  vector<E_Int>& poscs, 
  vector<K_INTERP::InterpData*>& listOfUnstrInterpData,
  vector<FldArrayF*>& listOfUnstrFields,
  vector<FldArrayF*>& listOfUnstrVelocities,
  vector<FldArrayI*>& connectu,
  vector<E_Int>& posxu, vector<E_Int>& posyu, 
  vector<E_Int>& poszu, vector<E_Int>& poscu, 
  FldArrayI& connectSurf, 
  E_Float* xSurf, E_Float* ySurf, E_Float* zSurf, E_Int sizeSurf, 
  K_INTERP::InterpData::InterpolationType interpType)
{
  E_Int nu = listOfUnstrInterpData.size();
  E_Int ns = listOfStructInterpData.size();
  vector<K_INTERP::InterpData*> allInterpDatas;
  vector<FldArrayF*> allFields;
  vector<void*> allA1;
  vector<void*> allA2;
  vector<void*> allA3;
  vector<void*> allA4;
  E_Int nzonesTot = ns+nu;
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
  allFields.insert(allFields.end(), listOfStructFields.begin(), listOfStructFields.end()); 
  allFields.insert(allFields.end(), listOfUnstrFields.begin(), listOfUnstrFields.end());
  
  allInterpDatas.reserve(nzonesTot);  // preallocate memory
  allInterpDatas.insert(allInterpDatas.end(), listOfStructInterpData.begin(), listOfStructInterpData.end()); 
  allInterpDatas.insert(allInterpDatas.end(), listOfUnstrInterpData.begin(), listOfUnstrInterpData.end());

  for (E_Int noz = 0; noz < ns; noz++)
  {
    allA1.push_back(&nis[noz]); 
    allA2.push_back(&njs[noz]);
    allA3.push_back(&nks[noz]);
    allA4.push_back(NULL);
  }
  for (E_Int noz = 0; noz < nu; noz++)
  {
    allA1.push_back(connectu[noz]); 
    allA2.push_back(NULL);
    allA3.push_back(NULL); 
    allA4.push_back(NULL);
  }

  E_Float dts6 = dt/6.;
  E_Float dts2 = dt/2.;
  E_Int noblk0 = 0;

  // 1- determination de k1 : k1 = u(xp)
  E_Float k1x = up; //coefficient "k1"
  E_Float k1y = vp;
  E_Float k1z = wp;
  
  //2- determination de k2 : k2 = u(xp2), xp2 = xp + dts2 * k1
  E_Float xp2 = xp + dts2 * k1x;
  E_Float yp2 = yp + dts2 * k1y;
  E_Float zp2 = zp + dts2 * k1z;
  
  // Cellule d interpolation pour calculer Up2
  FldArrayI indi;
  FldArrayF cf;
  if ( listOfStructVelocities.size() == 0 ) 
  {cf.malloc(4); indi.malloc(1);}
  else //ordre 2 structure
  {cf.malloc(8); indi.malloc(1);}
  FldArrayI tmpIndi(indi.getSize()); FldArrayF tmpCf(cf.getSize());
  FldArrayF up2(3);
  E_Float p0[3]; E_Float p1[3]; E_Float p2[3]; E_Float p[3];
  E_Float vol2 = 0.;
  // Si la surface n'est pas nulle, on projette le point sur cette surface
  if (sizeSurf != 0)
  {
    K_COMPGEOM::projectOrtho(xp2, yp2, zp2,
                             xSurf, ySurf, zSurf,
                             connectSurf, 
                             xp2, yp2, zp2,
                             p0, p1, p2, p);
  }
  E_Int type = 0; E_Int noblk = 0;
  short found = K_INTERP::getInterpolationCell(xp2, yp2, zp2, allInterpDatas,
                                               allFields, allA1, allA2, allA3, allA4,
                                               posxt, posyt, poszt, posct, 
                                               vol2, indi, cf, tmpIndi, tmpCf, type, noblk, interpType);
  
  if ( found == 0 ) return 0;// pas de pt d'interpolation trouve
  //Calcul de Up2
  noblk0 = noblk-1;
  if ( noblk0 < ns )
  {
    K_INTERP::compInterpolatedField(indi.begin(), cf, *listOfStructVelocities[noblk0],
                                    &nis[noblk0],&njs[noblk0],&nks[noblk0], type, up2);
  }
  else // non structure
  {
    E_Int noblku = noblk0-ns;
    K_INTERP::compInterpolatedField(indi.begin(), cf, *listOfUnstrVelocities[noblku],
                                    connectu[noblku], NULL, NULL, type, up2);
  }
  E_Float k2x = up2[0]; // "k2"
  E_Float k2y = up2[1];
  E_Float k2z = up2[2];
  
  //3- determination de k3 : k3 = u(xp3), xp3 = xp + dts2 * k2  
  xp2 = xp + dts2 * k2x;
  yp2 = yp + dts2 * k2y;
  zp2 = zp + dts2 * k2z;
  
  //Cellule d interpolation pour calculer Up3
  
  // Si la surface n'est pas nulle, on projette le point sur cette surface
  if (sizeSurf != 0)
  {
    K_COMPGEOM::projectOrtho(xp2, yp2, zp2,
                             xSurf, ySurf, zSurf,
                             connectSurf, 
                             xp2, yp2, zp2,
                             p0, p1, p2, p);
  }
  type = 0; noblk = 0;
  found = K_INTERP::getInterpolationCell( xp2, yp2, zp2, allInterpDatas,
                                          allFields, allA1, allA2, allA3, allA4,
                                          posxt, posyt, poszt, posct, 
                                          vol2, indi, cf, tmpIndi, tmpCf, type, noblk, interpType);


  if ( found == 0 ) return 0;// pas de pt d'interpolation trouve
    
  //Calcul de Up3
  noblk0 = noblk-1;
  if ( noblk0 < ns )
  {
    K_INTERP::compInterpolatedField(indi.begin(), cf, *listOfStructVelocities[noblk0],
                                    &nis[noblk0],&njs[noblk0],&nks[noblk0],
                                    type, up2);      
  }
  else // non structure
  {
    E_Int noblku = noblk0-ns;
    K_INTERP::compInterpolatedField(indi.begin(), cf, *listOfUnstrVelocities[noblku],
                                    connectu[noblku], NULL, NULL, type, up2);
  }
  E_Float k3x = up2[0]; // "k3"
  E_Float k3y = up2[1];
  E_Float k3z = up2[2];
  
  //4- determination de k4 : k4 = u(xp4), xp4 = xp + dt * k3  
  xp2 = xp + dt * k3x;
  yp2 = yp + dt * k3y;
  zp2 = zp + dt * k3z;

  //Cellule d interpolation pour calculer Up4

  // Si la surface n'est pas nulle, on projette le point sur cette surface
  if (sizeSurf != 0)
  {
    K_COMPGEOM::projectOrtho(xp2, yp2, zp2,
                             xSurf, ySurf, zSurf,
                             connectSurf, 
                             xp2, yp2, zp2,
                             p0, p1, p2, p);
  }
  type = 0; noblk = 0;
  found = K_INTERP::getInterpolationCell( xp2, yp2, zp2, allInterpDatas,
                                          allFields, allA1, allA2, allA3, allA4,
                                          posxt, posyt, poszt, posct, 
                                          vol2, indi, cf, tmpIndi, tmpCf, type, noblk, interpType);

  if ( noblk == 0 ) return 0;// pas de pt d interpolation trouve

  //Calcul de Up4
  noblk0 = noblk-1; 
  if ( noblk0 < ns )
  {
    K_INTERP::compInterpolatedField(indi.begin(), cf, *listOfStructVelocities[noblk0],
                                    &nis[noblk0],&njs[noblk0],&nks[noblk0], 
                                    type, up2);

  }
  else // non structure
  {
    E_Int noblku = noblk0-ns;
    K_INTERP::compInterpolatedField(indi.begin(), cf, *listOfUnstrVelocities[noblku],
                                    connectu[noblku], NULL, NULL, type, up2);

  }

  E_Float k4x = up2[0]; // "k4"
  E_Float k4y = up2[1];
  E_Float k4z = up2[2];

  // 5- Mise a jour de x(n+1)
  xn = xp + dts6 * (k1x+k2x+k3x+k4x);
  yn = yp + dts6 * (k1y+k2y+k3y+k4y);
  zn = zp + dts6 * (k1z+k2z+k3z+k4z); 

  allInterpDatas.clear(); allFields.clear(); allA1.clear(); allA2.clear(); allA3.clear(); allA4.clear();
  posxt.clear(); posyt.clear(); poszt.clear(); posct.clear();
  return 1;
}

//===========================================================================
/* Initialisation de la streamline pour des grilles structurees et 
   non structurees. Retourne le numero du bloc d interpolation dans la liste 
   des blocs structures si bloc d interpolation structure
   Sinon, si bloc d interpolation non structure, retourne le numero du bloc 
   dans liste non structuree
   Retourne 0 si pas trouve */
//===========================================================================
short K_POST::initStreamLine(
  E_Float xp, E_Float yp, E_Float zp,
  vector<K_INTERP::InterpData*>& listOfStructInterpData, 
  vector<FldArrayF*>& listOfStructFields,
  vector<FldArrayF*>& listOfStructVelocities,
  vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks,
  vector<E_Int>& posxs, vector<E_Int>& posys, 
  vector<E_Int>& poszs, vector<E_Int>& poscs,
  vector<K_INTERP::InterpData*>& listOfUnstrInterpData, 
  vector<FldArrayF*>& listOfUnstrFields,
  vector<FldArrayF*>& listOfUnstrVelocities,
  vector<FldArrayI*>& connectu,
  vector<E_Int>& posxu, vector<E_Int>& posyu, 
  vector<E_Int>& poszu, vector<E_Int>& poscu,
  FldArrayI& connectSurf, 
  E_Float* xSurf, E_Float* ySurf, E_Float* zSurf, E_Int sizeSurf, 
  E_Float& up, E_Float& vp, E_Float& wp, E_Float& dt,
  FldArrayI& indi, FldArrayF& cf, FldArrayF& streamPts, 
  K_INTERP::InterpData::InterpolationType interpType)
{
  FldArrayI tmpIndi(indi.getSize()); FldArrayF tmpCf(cf.getSize());
  E_Int noblkp0 = -1;
  E_Float p0[3]; E_Float p1[3]; E_Float p2[3]; E_Float p[3]; 
  E_Int ns = listOfStructInterpData.size();
  E_Int nu = listOfUnstrInterpData.size();
  vector<K_INTERP::InterpData*> allInterpDatas;
  vector<FldArrayF*> allFields;
  vector<void*> allA1;
  vector<void*> allA2;
  vector<void*> allA3;
  vector<void*> allA4;
  E_Int nzonesTot = ns+nu;
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
  allFields.insert(allFields.end(), listOfStructFields.begin(), listOfStructFields.end()); 
  allFields.insert(allFields.end(), listOfUnstrFields.begin(), listOfUnstrFields.end());
  
  allInterpDatas.reserve(nzonesTot);  // preallocate memory
  allInterpDatas.insert(allInterpDatas.end(), listOfStructInterpData.begin(), listOfStructInterpData.end()); 
  allInterpDatas.insert(allInterpDatas.end(), listOfUnstrInterpData.begin(), listOfUnstrInterpData.end());

  for (E_Int noz = 0; noz < ns; noz++)
  {
    allA1.push_back(&nis[noz]);
    allA2.push_back(&njs[noz]); 
    allA3.push_back(&nks[noz]); 
    allA4.push_back(NULL);
  }
  for (E_Int noz = 0; noz < nu; noz++)
  {
    allA1.push_back(connectu[noz]); 
    allA2.push_back(NULL); 
    allA3.push_back(NULL); 
    allA4.push_back(NULL);
  }
  // Si la surface n'est pas nulle, on projette le point sur cette surface
  if (sizeSurf != 0)
  {
    K_COMPGEOM::projectOrtho(xp, yp, zp,
                             xSurf, ySurf, zSurf,
                             connectSurf, 
                             xp, yp, zp,
                             p0, p1, p2, p);
  }
  // pt 0 : no du blk d interpolation dans interpDatas : demarre a 1 
  E_Int type = 0; E_Int noblkp = 0; E_Float voli = 0.;
  short found = K_INTERP::getInterpolationCell( xp, yp, zp, allInterpDatas,
                                                allFields, allA1, allA2, allA3, allA4,
                                                posxt, posyt, poszt, posct, 
                                                voli, indi, cf, tmpIndi, tmpCf, type, noblkp, interpType);

  if ( found < 1 ) //pas de pt trouve
  {
    streamPts.malloc(0);
    printf("Warning: streamLine: starting point not interpolable.\n");
    return 0;
  }

  /* Determination de dt */
  compInitialStep(noblkp, type, indi, nis, njs, nks, posxs, posys, poszs,
                  listOfStructFields, listOfStructVelocities,
                  posxu, posyu, poszu, 
                  listOfUnstrFields, connectu, listOfUnstrVelocities, dt);

  /* Calcul des champs du pt X(0)*/
  compStreamPtFields(0, xp, yp, zp, noblkp, type, indi, cf, 
                     nis, njs, nks, posxs, posys, poszs,
                     listOfStructFields,
                     posxu, posyu, poszu, listOfUnstrFields,connectu,
                     streamPts, interpType);
  FldArrayF u(3);
  noblkp0 = noblkp-1;
  if ( noblkp0 < ns )
  {      
    K_INTERP::compInterpolatedField(indi.begin(), cf, *listOfStructVelocities[noblkp0],
                                    &nis[noblkp0],&njs[noblkp0],&nks[noblkp0],
                                    type, u);
  }
  else // non structure
  {
    E_Int noblku = noblkp0-ns;
    K_INTERP::compInterpolatedField(indi.begin(), cf, *listOfUnstrVelocities[noblku],
                                    connectu[noblku], NULL, NULL, type, u);
  }

  up = u[0];
  vp = u[1];
  wp = u[2];
  allInterpDatas.clear(); allFields.clear(); allA1.clear(); allA2.clear(); allA3.clear(); allA4.clear();
  posxt.clear(); posyt.clear(); poszt.clear(); posct.clear();
  return noblkp;
}
//===========================================================================
/* Determine la liste des arrays structures qui ont la vitesse comme 
   information. Retourne la vitesse sur ces arrays */
//===========================================================================
E_Int K_POST::extractVectorFromStructArrays(
  E_Float signe, 
  vector<E_Int>& niIn, vector<E_Int>& njIn, vector<E_Int>& nkIn,
  vector<E_Int>& posxIn, vector<E_Int>& posyIn, vector<E_Int>& poszIn,
  vector<E_Int>& poscIn, vector<char*>& varStringIn, 
  vector<FldArrayF*>& fieldsIn, 
  vector<K_INTERP::InterpData*>& interpDataIn,
  vector<E_Int>& niOut, vector<E_Int>& njOut, vector<E_Int>& nkOut,
  vector<E_Int>& posxOut, vector<E_Int>& posyOut, vector<E_Int>& poszOut,
  vector<E_Int>& poscOut, vector<char*>& varStringOut, 
  vector<FldArrayF*>& fieldsOut, 
  vector<K_INTERP::InterpData*>& interpDataOut,
  vector<FldArrayF*>& vect, vector<char*>& vnames)
{
  E_Int size = fieldsIn.size();
  E_Int posv1, posv2, posv3;  
  for (E_Int v = 0; v < size; v++)
  {
    FldArrayF& f = *fieldsIn[v];
    posv1 = K_ARRAY::isNamePresent(vnames[0], varStringIn[v]);
    posv2 = K_ARRAY::isNamePresent(vnames[1], varStringIn[v]);
    posv3 = K_ARRAY::isNamePresent(vnames[2], varStringIn[v]);
    // variables pas presentes dans l'array v
    if ( posv1== -1 || posv2 == -1 || posv3 == -1) return -1;
    posv1++; posv2++; posv3++;

    E_Int ni = niIn[v];
    E_Int nj = njIn[v];
    E_Int nk = nkIn[v];
    if (ni == 1 || nj == 1 || nk == 1) return -2;
    else
    {
      E_Int npts = f.getSize();
      niOut.push_back(ni);
      njOut.push_back(nj);
      nkOut.push_back(nk);
      posxOut.push_back(posxIn[v]);
      posyOut.push_back(posyIn[v]);
      poszOut.push_back(poszIn[v]);
      poscOut.push_back(poscIn[v]);
      interpDataOut.push_back(interpDataIn[v]);
      fieldsOut.push_back(fieldsIn[v]);
      varStringOut.push_back(varStringIn[v]);
      
      FldArrayF* v = new FldArrayF(npts, 3);
      E_Float* v1 = v->begin(1);
      E_Float* v2 = v->begin(2);
      E_Float* v3 = v->begin(3);
      
      E_Float* fv1 = f.begin(posv1);
      E_Float* fv2 = f.begin(posv2);
      E_Float* fv3 = f.begin(posv3);
      for (E_Int i = 0; i < npts; i++)
      {
        v1[i] = signe * fv1[i];
        v2[i] = signe * fv2[i];
        v3[i] = signe * fv3[i];
      }
      
      vect.push_back(v);
    }
  } //3d mesh
  return 1;//variables du vecteur presentes dans variables
} 
//===========================================================================
/* Determine la liste des arrays non structures qui ont la vitesse comme 
   information. Retourne la vitesse sur ces arrays */
//===========================================================================
E_Int K_POST::extractVectorFromUnstrArrays(
  E_Float signe, vector<E_Int>& posxIn, vector<E_Int>& posyIn, 
  vector<E_Int>& poszIn, vector<E_Int>& poscIn, 
  vector<char*>& varStringIn, vector<FldArrayF*>& fieldsIn, 
  vector<FldArrayI*>& cntIn, vector<char*>& eltTypeIn,
  vector<K_INTERP::InterpData*>& interpDataIn,
  vector<E_Int>& posxOut, vector<E_Int>& posyOut, 
  vector<E_Int>& poszOut, vector<E_Int>& poscOut, 
  vector<char*>& varStringOut, vector<FldArrayF*>& fieldsOut, 
  vector<FldArrayI*>& cntOut, vector<char*>& eltTypeOut,
  vector<K_INTERP::InterpData*>& interpDataOut,
  vector<FldArrayF*>& vect, vector<char*>& vnames)
{
  E_Int size = fieldsIn.size();
  E_Int posv1, posv2, posv3;
  
  for (E_Int v = 0; v < size; v++)
  {
    FldArrayF& f = *fieldsIn[v];
    posv1 = K_ARRAY::isNamePresent(vnames[0], varStringIn[v]);
    posv2 = K_ARRAY::isNamePresent(vnames[1], varStringIn[v]); 
    posv3 = K_ARRAY::isNamePresent(vnames[2], varStringIn[v]); 

    // variables pas presentes dans l'array v
    if ( posv1== -1 || posv2 == -1 || posv3 == -1) return -1;
    posv1++; posv2++; posv3++;

    if ( strcmp(eltTypeIn[v],"TETRA") != 0 ) return -2;
    else
    {
      E_Int npts = f.getSize();
      posxOut.push_back(posxIn[v]);
      posyOut.push_back(posyIn[v]);
      poszOut.push_back(poszIn[v]);
      poscOut.push_back(poscIn[v]);
      varStringOut.push_back(varStringIn[v]);
      fieldsOut.push_back(fieldsIn[v]);
      cntOut.push_back(cntIn[v]);
      eltTypeOut.push_back(eltTypeIn[v]);
      interpDataOut.push_back(interpDataIn[v]);
      
      FldArrayF* v = new FldArrayF(npts, 3);
      E_Float* v1 = v->begin(1);
      E_Float* v2 = v->begin(2);
      E_Float* v3 = v->begin(3);
      
      E_Float* fv1 = f.begin(posv1);
      E_Float* fv2 = f.begin(posv2);
      E_Float* fv3 = f.begin(posv3);
            
      for (E_Int i = 0; i < npts; i++)
      {
        v1[i] = signe * fv1[i];
        v2[i] = signe * fv2[i];
        v3[i] = signe * fv3[i];
      }
      
      vect.push_back(v);
    } //3d mesh
  } 
  return 1;
}
//===========================================================================
/* Initialisation de la streamsurf pour des grilles structurees et 
   non structurees. Retourne le numero des blocs d interpolation dans la liste 
   des blocs structures si bloc d interpolation structure
   Sinon, si bloc d interpolation non structure, retourne le numero des blocs 
   dans liste non structuree
   Retourne 0 si pas trouve
*/
//===========================================================================
void K_POST::initStreamSurf(
  E_Float xp, E_Float yp, E_Float zp,
  vector<K_INTERP::InterpData*>& listOfStructInterpData, 
  vector<FldArrayF*>& listOfStructFields,
  vector<FldArrayF*>& listOfStructVelocities,
  vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks,
  vector<E_Int>& posxs, vector<E_Int>& posys, 
  vector<E_Int>& poszs, vector<E_Int>& poscs,
  vector<K_INTERP::InterpData*>& listOfUnstrInterpData, 
  vector<FldArrayF*>& listOfUnstrFields,
  vector<FldArrayF*>& listOfUnstrVelocities,
  vector<FldArrayI*>& connectu,
  vector<E_Int>& posxu, vector<E_Int>& posyu, 
  vector<E_Int>& poszu, vector<E_Int>& poscu,
  E_Float& up, E_Float& vp, E_Float& wp, E_Float& dt,
  FldArrayI& indi, FldArrayF& cf, FldArrayF& streamPts, 
  K_INTERP::InterpData::InterpolationType interpType)
{
  FldArrayI tmpIndi(indi.getSize()); FldArrayF tmpCf(cf.getSize());  
  E_Int ns = listOfStructFields.size();
  E_Int nu = listOfUnstrFields.size();
  vector<K_INTERP::InterpData*> allInterpDatas;
  vector<FldArrayF*> allFields;
  vector<void*> allA1;
  vector<void*> allA2;
  vector<void*> allA3;
  vector<void*> allA4;
  E_Int nzonesTot = ns+nu;
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
  allFields.insert(allFields.end(), listOfStructFields.begin(), listOfStructFields.end()); 
  allFields.insert(allFields.end(), listOfUnstrFields.begin(), listOfUnstrFields.end());
    
  allInterpDatas.reserve(nzonesTot);  // preallocate memory
  allInterpDatas.insert(allInterpDatas.end(), listOfStructInterpData.begin(), listOfStructInterpData.end()); 
  allInterpDatas.insert(allInterpDatas.end(), listOfUnstrInterpData.begin(), listOfUnstrInterpData.end());
    
  for (E_Int noz = 0; noz < ns; noz++)
  {
    allA1.push_back(&nis[noz]); 
    allA2.push_back(&njs[noz]); 
    allA3.push_back(&nks[noz]); 
    allA4.push_back(NULL);
  }
  for (E_Int noz = 0; noz < nu; noz++)
  {
    allA1.push_back(connectu[noz]); 
    allA2.push_back(NULL); 
    allA3.push_back(NULL); 
    allA4.push_back(NULL);
  }
  E_Int sizeBAR = streamPts.getSize();
  vector<E_Int> noblkp(sizeBAR);
  E_Int noblkp0;
  FldArrayI cnDum(0);

  for (E_Int p = 0; p < sizeBAR; p++)
  {
    // pt 0 : no du blk d interpolation dans interpDatas : demarre a 1 
    E_Int type = 0; noblkp[p] = 0; E_Float voli = 0.;
    short found = K_INTERP::getInterpolationCell(xp, yp, zp, allInterpDatas,
                                                 allFields, allA1, allA2, allA3, allA4,
                                                 posxt, posyt, poszt, posct, 
                                                 voli, indi, cf, tmpIndi, tmpCf, type, noblkp[p], interpType);
    if (found<1 ) //pas de pt trouve
    {
      streamPts.malloc(0);
      printf("Warning: streamSurf: starting BAR-array has not interpolable points.\n");
      return;
    }

    /* Determination de dt */
    compInitialStep(noblkp[p], type, indi, nis, njs, nks, posxs, posys, poszs,
                    listOfStructFields, listOfStructVelocities,
                    posxu, posyu, poszu, 
                    listOfUnstrFields, connectu, listOfUnstrVelocities, dt);
    /* Calcul des champs du pt X(0)*/
    compStreamPtFields(0, xp, yp, zp, noblkp[p], type, indi, cf, 
                       nis, njs, nks, posxs, posys, poszs,
                       listOfStructFields,
                       posxu, posyu, poszu, listOfUnstrFields, connectu,
                       streamPts, interpType);
    FldArrayF u(3);
    noblkp0 = noblkp[p]-1;
    if ( noblkp0 < ns ) 
    {
      K_INTERP::compInterpolatedField(indi.begin(), cf, *listOfStructVelocities[noblkp0],
                                      &nis[noblkp0],&njs[noblkp0],&nks[noblkp0],
                                      type, u);
    }
    else// non structure 
    {
      E_Int noblku = noblkp0-ns;
      K_INTERP::compInterpolatedField(indi.begin(), cf, *listOfUnstrVelocities[noblku],
                                      connectu[noblku], NULL, NULL, type, u);
    }
    up = u[0]; vp = u[1]; wp = u[2];
  }
  return;
}
