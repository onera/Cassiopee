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

# include "BlkInterp.h"
# include "Def/DefCplusPlusConst.h"

using namespace K_FUNC;
using namespace std;
using namespace K_FLD;

//=============================================================================
/* Find the interpolation coefficients in a given cell. Cut a cell
   in 24 tetrahedras made of 2 edges of the cell, the center of the cell and
   one center of a face. Uses a "jump" technique to find quickly the most 
   probable tetrahedra.
   Taken from FLU3M. */
//=============================================================================
E_Boolean
K_KINTERP::BlkInterpData::coeffInterpHexa(
  E_Float x, E_Float y, E_Float z,
  E_Float* xt, E_Float* yt, E_Float* zt,
  FldArrayF& cf)
{
  E_Float xi, yi, zi;
  E_Float xp, yp, zp;
  E_Float xr, yr, zr;
  E_Float xs, ys, zs;
  E_Float xq, yq, zq;
  E_Int isom, its, inds, indq, indr, itq, ibsom, ind, i;
  E_Float cf0;
  E_Int isomm, ixp, iyp, izp;
  E_Int is1, is2, is3, itri;
  const E_Float EPS = _EPS_TETRA;
  const E_Float EPS2 = _EPS_GEOM;
  
  /* For each edge of the cell, the index of the neigbour edges (3 x 8) */ 
  static E_Int neighbour[24] = {1,2,4,0,3,5,3,0,6,2,1,7,5,0,6,4,7,1,7,4,
                                2,6,5,3};
  
  /* For each edge and its neighbour, the proper center of face for 
     tetrahedra */
  static E_Int center[32] = {13,8,9,13,13,8,11,13,12,8,9,12,12,8,11,12,10,13,9,
                             10,13,10,11,13,12,10,9,12,12,10,11,12};
  /* (4 x 8) */
  static E_Int indss[32] = {0,3,5,6,1,2,4,7,2,4,7,1,3,5,6,0,4,7,1,2,
                            5,6,0,3,6,0,3,5,7,1,2,4};
  
  /* For each face of the cell, the number of edges (4 x 6) */
  static E_Int indfc[24] = {0,1,3,2,4,0,2,6,4,5,7,6,5,1,3,7,6,7,3,2,4,5,1,0};
    
  /* Index of edge of cell following interp. coeff. in the tetrahedra 
     (C,7,6,5,3) */
  static E_Int indsom[8] = {0,1,2,3,4,5,6,7};
  
  /* Index of most probable triangle */
  static E_Int indtr[64] = {5,13,5,13,7,15,7,15,
                            5,13,5,13,7,15,7,15,
                            23,21,19,17,23,21,19,17,
                            4,12,6,14,4,12,6,14,
                            0,0,2,2,8,8,10,10,
                            3,1,3,1,11,9,11,9,
                            22,22,18,18,20,20,16,16,
                            22,22,18,18,20,20,16,16};
  
  /* Index of points in tetrahedra from triangle number */
  static E_Int indtd[96] = {0,1,8,14,
                            1,3,8,14,
                            2,3,8,14,
                            0,2,8,14,
                            0,4,9,14,
                            0,2,9,14,
                            2,6,9,14,
                            4,6,9,14,
                            4,5,10,14,
                            5,7,10,14,
                            6,7,10,14,
                            4,6,10,14,
                            1,5,11,14,
                            1,3,11,14,
                            3,7,11,14,
                            5,7,11,14,
                            6,7,12,14,
                            3,7,12,14,
                            2,3,12,14,
                            2,6,12,14,
                            4,5,13,14,
                            1,5,13,14,
                            0,1,13,14,
                            0,4,13,14};

/*--------------------------------------------------------------*/
/* Get the most probable tetrahedra                             */
/* test a first tetrahedra, following the xi,yi,zi, deduce the  */
/* most probable tetrahedra for interpolation                   */
/*--------------------------------------------------------------*/
  
  xp = xt[14];
  yp = yt[14];
  zp = zt[14];

  xq = xt[11];
  yq = yt[11];
  zq = zt[11];

  xr = xt[12];
  yr = yt[12];
  zr = zt[12];

  xs = xt[10];
  ys = yt[10];
  zs = zt[10];

  coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                   xr, yr, zr, xs, ys, zs, xi, yi, zi);

  /* Compute index to choose the most probable edge of cell */
  ixp = 0;
  iyp = 0;
  izp = 0;
    
  if ((E_abs(xi) <= K_CONST::ONE+EPS2) &&
      (E_abs(yi) <= K_CONST::ONE+EPS2) &&
      (E_abs(zi) <= K_CONST::ONE+EPS2))
  {
    if (xi >= 0) ixp = 1;
    if (yi >= 0) iyp = 1;
    if (zi >= 0) izp = 1;
  }
    
  /* Most probable edge of the cell */
  isomm = indsom[ixp+2*iyp+4*izp];
    
  switch (isomm)
  {
    case 0:
      xi = -xi;
      yi = -yi;
      zi = -zi;
      break;
      
    case 1:
      yi = -yi;
      zi = -zi;
      break;
      
    case 2:
      xi = -xi;
      zi = -zi;
      break;
      
    case 3:
      zi = -zi;
      break;
      
    case 4:
      xi = -xi;
      yi = -yi;
      break;
      
    case 5:
      yi = -yi;
      break;
      
    case 6:
      xi = -xi;
      break;
  }
    
  /* Compute index to find the most probable triangle */
  is1 = E_Int((K_CONST::ONE + E_sign(xi-yi)) * K_CONST::ONE_HALF);
  is2 = E_Int((K_CONST::ONE + E_sign(yi-zi)) * K_CONST::ONE_HALF);
  is3 = E_Int((K_CONST::ONE + E_sign(zi-xi)) * K_CONST::ONE_HALF);
  
  itri = indtr[isomm+8*is1+16*is2+32*is3];

  indr = indtd[0+itri*4];
  inds = indtd[1+itri*4];
  indq = indtd[2+itri*4];

  /* Most probable tetrahedra */
  xq = xt[indq];
  yq = yt[indq];
  zq = zt[indq];
    
  xr = xt[indr];
  yr = yt[indr];
  zr = zt[indr];
  
  xs = xt[inds];
  ys = yt[inds];
  zs = zt[inds];

  coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                   xr, yr, zr, xs, ys, zs, xi, yi, zi);
  
  if ((xi > -EPS) && (yi > -EPS) &&
      (zi > -EPS) && (xi+yi+zi < K_CONST::ONE+3*EPS))
  {
    /* Transforming tetrahedra interpolation coefficient into */
    /* hexahedra interpolation coefficients */
    cf0 = K_CONST::ONE_EIGHT*(K_CONST::ONE-xi-yi-zi);
    
    cf[0] = cf0;
    cf[1] = cf0;
    cf[2] = cf0;
    cf[3] = cf0;
    cf[4] = cf0;
    cf[5] = cf0;
    cf[6] = cf0;
    cf[7] = cf0;

    E_Int indl = (indq-8)*4;
    E_Float tl = xi*K_CONST::ONE_FOURTH;
    ind = indfc[indl];
    cf[ind] = cf[ind] + tl;
    indl++;
    ind = indfc[indl];
    cf[ind] = cf[ind] + tl;
    indl++;
    ind = indfc[indl];
    cf[ind] = cf[ind] + tl;
    indl++;
    ind = indfc[indl];
    cf[ind] = cf[ind] + tl;

    cf[indr] = cf[indr] + yi;
    cf[inds] = cf[inds] + zi;
    
    return true;
  }
    
/*--------------------------------------------------------*/
/* If we cannot interpolate from the previous tetrahedra, */
/* we test all tetrahedra (24) in a most probable order   */
/*--------------------------------------------------------*/
    
  for (ibsom = 0; ibsom < 4; ibsom++) /* edges */
  {
    isom = indss[ibsom+4*isomm];
    indr = isom;
    xr = xt[indr];
    yr = yt[indr];
    zr = zt[indr];
    
    for (its = 0; its < 3; its++) /* neighbour of edge */
    {
      inds = neighbour[its+isom*3];
      xs = xt[inds];
      ys = yt[inds];
      zs = zt[inds];
      
      for (itq = 0; itq < 2; itq++) /* center */
      {
        indq = center[itq+its+isom*4];
        xq = xt[indq];
        yq = yt[indq];
        zq = zt[indq];
        
        coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                         xr, yr, zr, xs, ys, zs, xi, yi, zi);

        if ((xi > -EPS)&&(yi > -EPS)&&
            (zi > -EPS)&&(xi+yi+zi < K_CONST::ONE+3*EPS))
        {
          /* Transforming to hexahedra coefficients */ 
          cf0 = K_CONST::ONE_EIGHT*(1.-xi-yi-zi);
          for (i = 0; i < 8; i++) cf[i] = cf0;
          
          for (i = 0; i < 4; i++)
          {
            ind = indfc[i+(indq-8)*4];
            cf[ind] = cf[ind]+xi*K_CONST::ONE_FOURTH;
          }
          
          cf[indr] = cf[indr]+yi;
          cf[inds] = cf[inds]+zi;
          
          return true;
        }           
      }
    }
  }
  return false; /* Cannot interpolate from this cell */
}
//=============================================================================
/* Find the interpolation coefficients in a given cell. Cuts a cell
   in 24 tetrahedras made of 2 edges of the cell, the center of the cell and
   one center of a face. Uses a "jump" technique to find quickly the most 
   probable tetrahedra.
   For use with vectorized version */
//=============================================================================
E_Boolean K_KINTERP::BlkInterpData::
coeffInterpHexav(E_Int index, E_Int index2,
                 E_Float x, E_Float y, E_Float z,
                 FldArrayF& xtv, FldArrayF& ytv, FldArrayF& ztv,
                 FldArrayF& cfv)
{
  E_Float xi, yi, zi;
  E_Float xp, yp, zp;
  E_Float xr, yr, zr;
  E_Float xs, ys, zs;
  E_Float xq, yq, zq;
  E_Int isom, its, inds, indq, indr, itq, ibsom, ind, i;
  E_Float cf0;
  E_Int isomm, ixp, iyp, izp;
  E_Int is1, is2, is3, itri;
  E_Int isom3, isom4;
  E_Float xis4;
  E_Float* xtvp = xtv.begin() + index;
  E_Float* ytvp = ytv.begin() + index;
  E_Float* ztvp = ztv.begin() + index;

  const E_Float EPS = _EPS_TETRA;
  const E_Float EPS2 = _EPS_GEOM;
  E_Int off1 = xtv.getSize();
  E_Int off2 = 2*off1;
  E_Int off10 = 10*off1;
  E_Int itri4, indqm, isomm4;
  E_Float *ptx, *pty, *ptz;

  /* For each edge of the cell, the index of the neigbour edges (3 x 8) */ 
  static E_Int neighbour[24] = {1,2,4,0,3,5,3,0,6,2,1,7,5,0,6,4,7,1,7,4,
                                2,6,5,3};
  
  /* For each edge and its neighbour, the proper center of face for 
     tetrahedra */
  static E_Int center[32] = {13,8,9,13,13,8,11,13,12,8,9,12,12,8,11,12,10,13,9,
                             10,13,10,11,13,12,10,9,12,12,10,11,12};
  /* (4 x 8) */
  static E_Int indss[32] = {0,3,5,6,1,2,4,7,2,4,7,1,3,5,6,0,4,7,1,2,
                            5,6,0,3,6,0,3,5,7,1,2,4};
  
  /* For each face of the cell, the number of edges (4 x 6) */
  static E_Int indfc[24] = {0,1,3,2,4,0,2,6,4,5,7,6,5,1,3,7,6,7,3,2,4,5,1,0};
  
  /* Index of most probable triangle */
  static E_Int indtr[64] = {5,13,5,13,7,15,7,15,
                            5,13,5,13,7,15,7,15,
                            23,21,19,17,23,21,19,17,
                            4,12,6,14,4,12,6,14,
                            0,0,2,2,8,8,10,10,
                            3,1,3,1,11,9,11,9,
                            22,22,18,18,20,20,16,16,
                            22,22,18,18,20,20,16,16};
  
  /* Index of points in tetrahedra from triangle number */
  static E_Int indtd[96] = {0,1,8,14,
                            1,3,8,14,
                            2,3,8,14,
                            0,2,8,14,
                            0,4,9,14,
                            0,2,9,14,
                            2,6,9,14,
                            4,6,9,14,
                            4,5,10,14,
                            5,7,10,14,
                            6,7,10,14,
                            4,6,10,14,
                            1,5,11,14,
                            1,3,11,14,
                            3,7,11,14,
                            5,7,11,14,
                            6,7,12,14,
                            3,7,12,14,
                            2,3,12,14,
                            2,6,12,14,
                            4,5,13,14,
                            1,5,13,14,
                            0,1,13,14,
                            0,4,13,14};

/*--------------------------------------------------------------*/
/* Get the most probable tetrahedra                             */
/* test a first tetrahedra, following the xi,yi,zi, deduce the  */
/* most probable tetrahedra for interpolation                   */
/*--------------------------------------------------------------*/

  ptx = xtvp + off10;
  pty = ytvp + off10;
  ptz = ztvp + off10;
  
  xs = *ptx; 
  ys = *pty;
  zs = *ptz; 

  ptx = ptx + off1;
  pty = pty + off1;
  ptz = ptz + off1;
  xq = *ptx;
  yq = *pty; 
  zq = *ptz; 

  ptx = ptx + off1;
  pty = pty + off1;
  ptz = ptz + off1;

  xr = *ptx; 
  yr = *pty;
  zr = *ptz; 
  
  ptx = ptx + off2;
  pty = pty + off2;
  ptz = ptz + off2;
  xp = *ptx; 
  yp = *pty;
  zp = *ptz; 

  coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                   xr, yr, zr, xs, ys, zs, xi, yi, zi);

  /* Compute index to choose the most probable edge of cell */
  ixp = 0;
  iyp = 0;
  izp = 0;
    
  if ((E_abs(xi) <= K_CONST::ONE+EPS2)&&
      (E_abs(yi) <= K_CONST::ONE+EPS2)&&
      (E_abs(zi) <= K_CONST::ONE+EPS2))
  {
    if (xi >= 0) ixp = 1;
    if (yi >= 0) iyp = 2;
    if (zi >= 0) izp = 4;
  }
    
  /* Most probable edge of the cell */
  isomm = ixp + iyp + izp;
    
  switch (isomm)
  {
    case 0:
      xi = -xi;
      yi = -yi;
      zi = -zi;
      break;
      
    case 1:
      yi = -yi;
      zi = -zi;
      break;
      
    case 2:
      xi = -xi;
      zi = -zi;
      break;
      
    case 3:
      zi = -zi;
      break;
      
    case 4:
      xi = -xi;
      yi = -yi;
      break;
      
    case 5:
      yi = -yi;
      break;
      
    case 6:
      xi = -xi;
      break;
  }

  /* Compute index to find the most probable triangle */
  is1 = E_Int((K_CONST::ONE + E_sign(xi-yi)) * K_CONST::ONE_HALF);
  is2 = E_Int((K_CONST::ONE + E_sign(yi-zi)) * K_CONST::ONE_HALF);
  is3 = E_Int((K_CONST::ONE + E_sign(zi-xi)) * K_CONST::ONE_HALF);
  
  itri = indtr[isomm+8*is1+16*is2+32*is3];
  itri4 = 4*itri;

  indr = indtd[0+itri4];
  inds = indtd[1+itri4];
  indq = indtd[2+itri4];

  /* Most probable tetrahedra */
  ptx = xtvp + indq*off1;
  pty = ytvp + indq*off1;
  ptz = ztvp + indq*off1;

  xq = *ptx;
  yq = *pty; 
  zq = *ptz;
    
  ptx = xtvp + indr*off1;
  pty = ytvp + indr*off1;
  ptz = ztvp + indr*off1;
  xr = *ptx;
  yr = *pty; 
  zr = *ptz; 
  
  ptx = xtvp + inds*off1;
  pty = ytvp + inds*off1;
  ptz = ztvp + inds*off1;
  xs = *ptx;
  ys = *pty;
  zs = *ptz;

  coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                   xr, yr, zr, xs, ys, zs, xi, yi, zi);
    
  if ((xi>-EPS)&&(yi>-EPS)&&
      (zi>-EPS)&&(xi+yi+zi<K_CONST::ONE+3*EPS))
  {
    /* transforming tetrahedra interpolation coefficient into */
    /* hexahedra interpolation coefficients */
    xis4 = xi*K_CONST::ONE_FOURTH;
    indqm = (indq-8)*4;
    cf0 = K_CONST::ONE_EIGHT*(K_CONST::ONE-xi-yi-zi);
    for (i = 0; i < 8; i++) cfv(index2,i+1) = cf0;
    
    for (i = 0; i < 4; i++)
    {
      ind = indfc[i+indqm];
      cfv(index2, ind+1) = cfv(index2, ind+1)+ xis4;
    }
    
    cfv(index2, indr+1) = cfv(index2, indr+1)+yi;
    cfv(index2, inds+1) = cfv(index2, inds+1)+zi;
    
    return true;
  }
    
/*--------------------------------------------------------*/
/* If we cannot interpolate from the previous tetrahedra, */
/* we test all tetrahedra (24) in a most probable order   */
/*--------------------------------------------------------*/

  isomm4 = 4*isomm;
  for (ibsom = 0; ibsom < 4; ibsom++) /* edges */
  {
    isom = indss[ibsom+isomm4];
    indr = isom;
    
    ptx = xtvp + indr*off1;
    pty = ytvp + indr*off1;
    ptz = ztvp + indr*off1;

    xr = *ptx; 
    yr = *pty; 
    zr = *ptz; 
    
    isom3 = 3*isom;
    for (its = 0; its < 3; its++) /* neighbour of edge */
    {
      inds = neighbour[its+isom3];
      ptx = xtvp + inds*off1;
      pty = ytvp + inds*off1;
      ptz = ztvp + inds*off1;
      xs = *ptx; 
      ys = *pty; 
      zs = *ptz; 
      
      isom4 = its + 4*isom;
      for (itq = 0; itq < 2; itq++) /* center */
      {
        indq = center[itq+isom4];
        ptx = xtvp + indq*off1;
        pty = ytvp + indq*off1;
        ptz = ztvp + indq*off1;
        xq = *ptx; 
        yq = *pty; 
        zq = *ptz; 
        
        coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                         xr, yr, zr, xs, ys, zs, xi, yi, zi);
        
        if ((xi>-EPS)&&(yi>-EPS)&&
            (zi>-EPS)&&(xi+yi+zi<K_CONST::ONE+3*EPS))
        {
          indqm = (indq-8)*4;
          xis4 = xi*K_CONST::ONE_FOURTH;
          /* transforming to hexahedra coefficients */ 
          cf0 = K_CONST::ONE_EIGHT*(1.-xi-yi-zi);
          for (i=0; i<8; i++) cfv(index2,i+1) = cf0;
          
          for (i=0; i<4; i++)
          {
            ind = indfc[i+indqm];
            cfv(index2, ind+1) = cfv(index2, ind+1)+xis4;
          }
          
          cfv(index2, indr+1) = cfv(index2, indr+1)+yi;
          cfv(index2, inds+1) = cfv(index2, inds+1)+zi;
          
          return true;
        }           
      }
    }
  }
  return false; /* Cannot interpolate from this cell */
}  
//=============================================================================
/* Test if the cell contains the point to interpolate */
//=============================================================================
E_Boolean K_KINTERP::BlkInterpData::
getCellJump(E_Float x, E_Float y, E_Float z,
            E_Float* xt, E_Float* yt, E_Float* zt,
            E_Int& isomm,
            E_Float& xi, E_Float& yi, E_Float& zi)
{
  E_Int ixp, iyp, izp;
  E_Float xp, yp, zp;
  E_Float xr, yr, zr;
  E_Float xs, ys, zs;
  E_Float xq, yq, zq;

  const E_Float EPS = _EPS_GEOM;

  /* Index of edge of cell following interp. coeff. in the tetrahedra (C,7,6,5,3) */
  static E_Int indsom[8] = {0,1,2,3,4,5,6,7};

  /*-------------------------------------------------------*/
  /* Test a tetraedra, following the xi,yi,zi, deduce the  */
  /* probable  cell by technique of jump for interpolation */
  /*-------------------------------------------------------*/

  xp = xt[14];
  yp = yt[14];
  zp = zt[14];

  xq = xt[11];
  yq = yt[11];
  zq = zt[11];

  xr = xt[12];
  yr = yt[12];
  zr = zt[12];

  xs = xt[10];
  ys = yt[10];
  zs = zt[10];

  coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                   xr, yr, zr, xs, ys, zs, xi, yi, zi);

  /* Compute index to choose the most probable edge of cell */
  ixp = 0;
  iyp = 0;
  izp = 0;
  if ((E_abs(xi) <= K_CONST::ONE+EPS)&&
      (E_abs(yi) <= K_CONST::ONE+EPS)&&
      (E_abs(zi) <= K_CONST::ONE+EPS))
  {
    if (xi >= 0) ixp = 1;
    if (yi >= 0) iyp = 1;
    if (zi >= 0) izp = 1;
    isomm = indsom[ixp+2*iyp+4*izp];
    return true;
  }
  else
    return false;
}
//=============================================================================
/*
  Find the interpolation coefficients in a given cell. Cuts a cell
  in 24 tetrahedras made of 2 edges of the cell, the center of the cell and
  one center of a face. Uses a "jump" technique to find quickly the most 
  probable tetrahedra.
  Taken from FLU3M.
*/
//=============================================================================
E_Boolean K_KINTERP::BlkInterpData::
getCoeffInterpHexa(E_Float x, E_Float y, E_Float z,
                   E_Int isomm,
                   E_Float xi, E_Float yi, E_Float zi, 
                   E_Float* xt, E_Float* yt, E_Float* zt,
                   FldArrayF& cf)
{
  E_Int isom, ibsom;
  E_Int ind, i;

  E_Float xp, yp, zp;
  E_Float xr, yr, zr;
  E_Float xs, ys, zs;
  E_Float xq, yq, zq;
  E_Float cf0;
  const E_Float EPS = _EPS_TETRA;
  
  /* for each edge of the cell, the index of the neigbour edges (3 x 8) */ 
  static E_Int neighbour[24] = {1,2,4,0,3,5,3,0,6,2,1,7,5,0,6,4,7,1,7,4,2,6,5,3};
  
  /* for each edge and its neighbour, the proper center of face for tetrahedra */
  static E_Int center[32] = {13,8,9,13,13,8,11,13,12,8,9,12,12,8,11,12,10,13,9,
                             10,13,10,11,13,12,10,9,12,12,10,11,12};
  /* (4 x 8) */
  static E_Int indss[32] = {0,3,5,6,1,2,4,7,2,4,7,1,3,5,6,0,4,7,1,2,
                            5,6,0,3,6,0,3,5,7,1,2,4};
  
  /* for each face of the cell, the number of edges (4 x 6) */
  static E_Int indfc[24] = {0,1,3,2,4,0,2,6,4,5,7,6,5,1,3,7,6,7,3,2,4,5,1,0};
      
  /* index of most probable triangle */
  static E_Int indtr[64] = {5,13,5,13,7,15,7,15,
                            5,13,5,13,7,15,7,15,
                            23,21,19,17,23,21,19,17,
                            4,12,6,14,4,12,6,14,
                            0,0,2,2,8,8,10,10,
                            3,1,3,1,11,9,11,9,
                            22,22,18,18,20,20,16,16,
                            22,22,18,18,20,20,16,16};
  
  /* index of points in tetrahedra from triangle number */
  static E_Int indtd[96] = {0,1,8,14,
                            1,3,8,14,
                            2,3,8,14,
                            0,2,8,14,
                            0,4,9,14,
                            0,2,9,14,
                            2,6,9,14,
                            4,6,9,14,
                            4,5,10,14,
                            5,7,10,14,
                            6,7,10,14,
                            4,6,10,14,
                            1,5,11,14,
                            1,3,11,14,
                            3,7,11,14,
                            5,7,11,14,
                            6,7,12,14,
                            3,7,12,14,
                            2,3,12,14,
                            2,6,12,14,
                            4,5,13,14,
                            1,5,13,14,
                            0,1,13,14,
                            0,4,13,14};

  /* most probable edge of the cell */    
  switch (isomm)
  {
    case 0:
      xi = -xi;
      yi = -yi;
      zi = -zi;
      break;
      
    case 1:
      yi = -yi;
      zi = -zi;
      break;
      
    case 2:
      xi = -xi;
      zi = -zi;
      break;
      
    case 3:
      zi = -zi;
      break;
      
    case 4:
      xi = -xi;
      yi = -yi;
      break;
      
    case 5:
      yi = -yi;
      break;
      
    case 6:
      xi = -xi;
      break;
  }
    
  /* Compute index to find the most probable triangle */
  E_Int is1 = static_cast<E_Int>((K_CONST::ONE + E_sign(xi-yi)) * K_CONST::ONE_HALF);
  E_Int is2 = static_cast<E_Int>((K_CONST::ONE + E_sign(yi-zi)) * K_CONST::ONE_HALF);
  E_Int is3 = static_cast<E_Int>((K_CONST::ONE + E_sign(zi-xi)) * K_CONST::ONE_HALF);
  
  E_Int itri = indtr[isomm+8*is1+16*is2+32*is3];

  E_Int indr = indtd[0+itri*4];
  E_Int inds = indtd[1+itri*4];
  E_Int indq = indtd[2+itri*4];
  E_Int indp = indtd[3+itri*4];

  /* Most probable tetrahedra */
  xp = xt[indp];
  yp = yt[indp];
  zp = zt[indp];

  xq = xt[indq];
  yq = yt[indq];
  zq = zt[indq];
    
  xr = xt[indr];
  yr = yt[indr];
  zr = zt[indr];
  
  xs = xt[inds];
  ys = yt[inds];
  zs = zt[inds];

  coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                   xr, yr, zr, xs, ys, zs, xi, yi, zi);

  if ((xi > -EPS)&&(yi > -EPS)&&
      (zi > -EPS)&&(xi+yi+zi < K_CONST::ONE+3*EPS))
  {
    /* transforming tetrahedra interpolation coefficient into */
    /* hexahedra interpolation coefficients */
      
    cf0 = K_CONST::ONE_EIGHT*(K_CONST::ONE-xi-yi-zi);
    for (i = 0; i < 8; i++) cf[i] = cf0;
    
    for (i = 0; i < 4; i++)
    {
      ind = indfc[i+(indq-8)*4];
      cf[ind] = cf[ind]+xi*K_CONST::ONE_FOURTH;
    }
      
    cf[indr] = cf[indr]+yi;
    cf[inds] = cf[inds]+zi;
      
    return true;
  }
  
  /*--------------------------------------------------------*/
  /* if we cannot interpolate from the previous tetrahedra, */
  /* we test all tetrahedra (24) in a most probable order   */
  /*--------------------------------------------------------*/

  /* first point of the tetrahedra is the center of cell (supp) */
   
  for (ibsom = 0; ibsom < 4; ibsom++) /* edges */
  {
    isom = indss[ibsom+4*isomm];
    indr = isom;
    xr = xt[indr];
    yr = yt[indr];
    zr = zt[indr];
    
    for (E_Int its = 0; its < 3; its++) /* neighbour of edge */
    {
      inds = neighbour[its+isom*3];
      xs = xt[inds];
      ys = yt[inds];
      zs = zt[inds];
      
      for (E_Int itq = 0; itq < 2; itq++) /* center */
      {
        indq = center[itq+its+isom*4];
        xq = xt[indq];
        yq = yt[indq];
        zq = zt[indq];
        
        coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                         xr, yr, zr, xs, ys, zs, xi, yi, zi);


        if ((xi > -EPS)&&(yi > -EPS) &&
            (zi > -EPS)&&(xi+yi+zi < K_CONST::ONE+3*EPS))
        {
          /* transforming to hexahedra coefficients */ 
          cf0 = K_CONST::ONE_EIGHT*(1.-xi-yi-zi);
          for (i = 0; i < 8; i++) cf[i] = cf0;

          for (i = 0; i < 4; i++)
          {
            ind = indfc[i+(indq-8)*4];
            cf[ind] = cf[ind]+xi*K_CONST::ONE_FOURTH;
          }
          
          cf[indr] = cf[indr]+yi;
          cf[inds] = cf[inds]+zi;

          return true;
        }           
      }
    }
  }
  return false; /* Cannot interpolate from this cell */
}
//=============================================================================
/* Find the interpolation coefficients in all cells. Cuts each cell
   in 24 tetrahedras made of 2 edges of the cell, the center of the cell and
   one center of a face. Uses a "jump" technique to find quickly the most 
   probable tetrahedra. */
//=============================================================================
void K_KINTERP::BlkInterpData::
coeffInterpHexav(const E_Int istart, const E_Int nbI,
                 const FldArrayF& coord, 
                 FldArrayF& xtv, FldArrayF& ytv, FldArrayF& ztv,
                 FldArrayF& cfv, FldArrayIS& found, 
                 FldArrayI& indv, E_Int& c)
{
  E_Int index, index2;
  E_Float x, y, z;
  E_Boolean flg;
  E_Float xi, yi, zi;
  E_Float xp, yp, zp;
  E_Float xr, yr, zr;
  E_Float xs, ys, zs;
  E_Float xq, yq, zq;
  E_Int isom, its, inds, indq, indr, itq, ibsom, ind, i;
  E_Float cf0;
  E_Int isomm, ixp, iyp, izp;
  E_Int is1, is2, is3, itri;
  E_Int isom3, isom4;
  E_Float xis4;
  const E_Float EPS = _EPS_TETRA;
  const E_Float EPS2 = _EPS_GEOM;
  E_Int itri4, indqm, isomm4;
  E_Float *ptx, *pty, *ptz;

  /* For each edge of the cell, the index of the neigbour edges (3 x 8) */ 
  static E_Int neighbour[24] = {1,2,4,0,3,5,3,0,6,2,1,7,5,0,6,4,7,1,7,4,
                                2,6,5,3};
  
  /* For each edge and its neighbour, the proper center of face for 
     tetrahedra */
  static E_Int center[32] = {13,8,9,13,13,8,11,13,12,8,9,12,12,8,11,12,10,13,9,
                             10,13,10,11,13,12,10,9,12,12,10,11,12};
  /* (4 x 8) */
  static E_Int indss[32] = {0,3,5,6,1,2,4,7,2,4,7,1,3,5,6,0,4,7,1,2,
                            5,6,0,3,6,0,3,5,7,1,2,4};
  
  /* For each face of the cell, the number of edges (4 x 6) */
  static E_Int indfc[24] = {0,1,3,2,4,0,2,6,4,5,7,6,5,1,3,7,6,7,3,2,4,5,1,0};
  
  /* Index of most probable triangle */
  static E_Int indtr[64] = {5,13,5,13,7,15,7,15,
                            5,13,5,13,7,15,7,15,
                            23,21,19,17,23,21,19,17,
                            4,12,6,14,4,12,6,14,
                            0,0,2,2,8,8,10,10,
                            3,1,3,1,11,9,11,9,
                            22,22,18,18,20,20,16,16,
                            22,22,18,18,20,20,16,16};
  
  /* Index of points in tetrahedra from triangle number */
  static E_Int indtd[96] = {0,1,8,14,
                            1,3,8,14,
                            2,3,8,14,
                            0,2,8,14,
                            0,4,9,14,
                            0,2,9,14,
                            2,6,9,14,
                            4,6,9,14,
                            4,5,10,14,
                            5,7,10,14,
                            6,7,10,14,
                            4,6,10,14,
                            1,5,11,14,
                            1,3,11,14,
                            3,7,11,14,
                            5,7,11,14,
                            6,7,12,14,
                            3,7,12,14,
                            2,3,12,14,
                            2,6,12,14,
                            4,5,13,14,
                            1,5,13,14,
                            0,1,13,14,
                            0,4,13,14};
 
  E_Int off1 = xtv.getSize();
  E_Int off2 = 2*off1;
  E_Int off10 = 10*off1;
  
  for (E_Int j = 0; j < nbI; j++)
  {
    if (indv[j] != -1)
    {
      flg = false;
      index = j;
      index2 = j+istart;

      E_Float* xtvp = xtv.begin() + index;
      E_Float* ytvp = ytv.begin() + index;
      E_Float* ztvp = ztv.begin() + index;

      x = coord(index2,1);
      y = coord(index2,2);
      z = coord(index2,3);
      
      /*--------------------------------------------------------------*/
      /* Get the most probable tetrahedra                             */
      /* test a first tetrahedra, following the xi,yi,zi, deduce the  */
      /* most probable tetrahedra for interpolation                   */
      /*--------------------------------------------------------------*/
      
      ptx = xtvp + off10;
      pty = ytvp + off10;
      ptz = ztvp + off10;
      
      xs = *ptx; 
      ys = *pty; 
      zs = *ptz;
      
      ptx = ptx + off1;
      pty = pty + off1;
      ptz = ptz + off1;
      xq = *ptx; 
      yq = *pty; 
      zq = *ptz; 
      
      ptx = ptx + off1;
      pty = pty + off1;
      ptz = ptz + off1;
      
      xr = *ptx; 
      yr = *pty; 
      zr = *ptz; 
      
      ptx = ptx + off2;
      pty = pty + off2;
      ptz = ptz + off2;
      xp = *ptx; 
      yp = *pty; 
      zp = *ptz; 
      
      coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                   xr, yr, zr, xs, ys, zs, xi, yi, zi);
      
      /* Compute index to choose the most probable edge of cell */
      ixp = 0;
      iyp = 0;
      izp = 0;
      
      if ((E_abs(xi) <= K_CONST::ONE+EPS2)&&
          (E_abs(yi) <= K_CONST::ONE+EPS2)&&
          (E_abs(zi) <= K_CONST::ONE+EPS2))
      {
        if (xi >= 0) ixp = 1;
        if (yi >= 0) iyp = 2;
        if (zi >= 0) izp = 4;
      }
      
      /* Most probable edge of the cell */
      isomm = ixp + iyp + izp;
    
      switch (isomm)
      {
        case 0:
          xi = -xi;
          yi = -yi;
          zi = -zi;
          break;
          
        case 1:
          yi = -yi;
          zi = -zi;
          break;
          
        case 2:
          xi = -xi;
          zi = -zi;
          break;
          
        case 3:
          zi = -zi;
          break;
          
        case 4:
          xi = -xi;
          yi = -yi;
          break;
          
        case 5:
          yi = -yi;
          break;
          
        case 6:
          xi = -xi;
          break;
      }
      
      /* Compute index to find the most probable triangle */
      is1 = E_Int((K_CONST::ONE + E_sign(xi-yi)) * K_CONST::ONE_HALF);
      is2 = E_Int((K_CONST::ONE + E_sign(yi-zi)) * K_CONST::ONE_HALF);
      is3 = E_Int((K_CONST::ONE + E_sign(zi-xi)) * K_CONST::ONE_HALF);
      
      itri = indtr[isomm+8*is1+16*is2+32*is3];
      itri4 = 4*itri;
      
      indr = indtd[0+itri4];
      inds = indtd[1+itri4];
      indq = indtd[2+itri4];
      
      /* Most probable tetrahedra */
      ptx = xtvp + indq*off1;
      pty = ytvp + indq*off1;
      ptz = ztvp + indq*off1;
      
      xq = *ptx;
      yq = *pty;
      zq = *ptz;
      
      ptx = xtvp + indr*off1;
      pty = ytvp + indr*off1;
      ptz = ztvp + indr*off1;
      xr = *ptx; 
      yr = *pty; 
      zr = *ptz; 
  
      ptx = xtvp + inds*off1;
      pty = ytvp + inds*off1;
      ptz = ztvp + inds*off1;
      xs = *ptx;
      ys = *pty;
      zs = *ptz;
      
 
      coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                       xr, yr, zr, xs, ys, zs, xi, yi, zi);
      
      if ((xi>-EPS)&&(yi>-EPS)&&
          (zi>-EPS)&&(xi+yi+zi<K_CONST::ONE+3*EPS))
      {
        /* transforming tetrahedra interpolation coefficient into */
        /* hexahedra interpolation coefficients */
        xis4 = xi*K_CONST::ONE_FOURTH;
        indqm = (indq-8)*4;
        cf0 = K_CONST::ONE_EIGHT*(K_CONST::ONE-xi-yi-zi);
        for (i = 0; i < 8; i++) cfv(index2,i+1) = cf0;
        
        for (i = 0; i < 4; i++)
        {
          ind = indfc[i+indqm];
          cfv(index2, ind+1) = cfv(index2, ind+1)+ xis4;
        }
        
        cfv(index2, indr+1) = cfv(index2, indr+1)+yi;
        cfv(index2, inds+1) = cfv(index2, inds+1)+zi;
        
        flg = true;
        goto end;
      }
      
      /*--------------------------------------------------------*/
      /* If we cannot interpolate from the previous tetrahedra, */
      /* we test all tetrahedra (24) in a most probable order   */
      /*--------------------------------------------------------*/

      isomm4 = 4*isomm;
      for (ibsom = 0; ibsom < 4; ibsom++) /* edges */
      {
        isom = indss[ibsom+isomm4];
        indr = isom;
        
        ptx = xtvp + indr*off1;
        pty = ytvp + indr*off1;
        ptz = ztvp + indr*off1;
        
        xr = *ptx;
        yr = *pty; 
        zr = *ptz; 
        
        isom3 = 3*isom;
        for (its = 0; its < 3; its++) /* neighbour of edge */
        {
          inds = neighbour[its+isom3];
          ptx = xtvp + inds*off1;
          pty = ytvp + inds*off1;
          ptz = ztvp + inds*off1;
          xs = *ptx; 
          ys = *pty; 
          zs = *ptz; 
          
          isom4 = its + 4*isom;
          for (itq = 0; itq < 2; itq++) /* center */
          {
            indq = center[itq+isom4];
            ptx = xtvp + indq*off1;
            pty = ytvp + indq*off1;
            ptz = ztvp + indq*off1;
            xq = *ptx; 
            yq = *pty; 
            zq = *ptz; 
            
            coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                             xr, yr, zr, xs, ys, zs, xi, yi, zi);
           
            if ((xi>-EPS)&&(yi>-EPS)&&
                (zi>-EPS)&&(xi+yi+zi<K_CONST::ONE+3*EPS))
            {
              indqm = (indq-8)*4;
              xis4 = xi*K_CONST::ONE_FOURTH;
              /* transforming to hexahedra coefficients */ 
              cf0 = K_CONST::ONE_EIGHT*(1.-xi-yi-zi);
              for (i=0; i<8; i++) cfv(index2,i+1) = cf0;
              
              for (i=0; i<4; i++)
              {
                ind = indfc[i+indqm];
                cfv(index2, ind+1) = cfv(index2, ind+1)+xis4;
              }
              
              cfv(index2, indr+1) = cfv(index2, indr+1)+yi;
              cfv(index2, inds+1) = cfv(index2, inds+1)+zi;
              
              flg = true;
              goto end;
            }           
          }
        }
      }
      // flg = false; /* Cannot interpolate from this cell */
     
      end:;
      if ( flg == true) // really found in candidate cells list
      { 
        found[j+istart] = 1;
        indv[j] = -1;
        c++;
      }
    } //test indv[j] != -1 
  }//boucle j 
}  
//========================================================================
/* compLagrangeCoefs : for an array of points
   corr[i] vaut 1 si correction de type O2CF pour les interpType OiABC
*/
//========================================================================
void K_KINTERP::BlkInterpData::
compLagrangeCoefsv(FldArrayF& coord, E_Int istart, E_Int iend,
                   FldArrayI& icv, FldArrayI& jcv, FldArrayI& kcv,
                   FldArrayIS& found, FldArrayF& cf,
                   FldArrayIS& corr,
                   InterpolationType interpType, 
                   InterpMeshType interpMeshType)
{
  corr.setAllValuesAtNull();

  FldArrayF cfloc(cf.getNfld());
  
  for (E_Int i = istart; i < iend; i++)
  {
    if ( found[i] > 0 )
    {
      corr[i] = compLagrangeCoefs(coord(i,1), coord(i,2), coord(i,3),
                                  icv[i], jcv[i], kcv[i], cfloc, interpType,
                                  interpMeshType);
      for (E_Int v = 1; v <= cfloc.getSize(); v++)
        cf(i,v) = cfloc[v-1];
    }
  }
}
//=============================================================================
/* Correction des coefficients d interpolation a l ordre 2 de type O2CF
   pour les points dont l'approximation n a pas marche 
   Maj de indi a 2 localement
*/
//=============================================================================
void K_KINTERP::BlkInterpData::correctInterpToO2CFv(
  FldArrayF& coord, E_Int istart, E_Int iend,
  FldArrayI& indiv,
  FldArrayIS& foundv, FldArrayF& cfv,
  FldArrayIS& corrv,
  InterpMeshType interpMeshType)
{
  E_Int ic, jc, kc;
  E_Float x, y, z;
  FldArrayF cf(8);
  E_Int nfld = cfv.getNfld();
  FldArrayI tmpIndi(7);

  //correction pour les points en O2CF
  for (E_Int i = istart; i < iend; i++)
  {
    if ( corrv[i] == 1 ) 
    {
      x = coord(i,1);
      y = coord(i,2);
      z = coord(i,3);
      
      foundv[i] = searchInterpolationCellO2CF( x, y, z, ic, jc, kc, cf );
      if ( foundv[i] < 1 ) goto next;
      
      //stockage en premier dans le tableau cfv
      for (E_Int j = 0; j < 8; j++)
        cfv(i,j+1) = cf[j];
      
      //autres coefs nuls
      for (E_Int j = 8; j < nfld; j++)
        cfv(i,j+1) = 0.;

      // maj de indi
      if (interpMeshType == EXT_CENTERS)
        fromExtendedToStandardCenters(ic, jc, kc, tmpIndi, O2CF);
      else 
        compStandardIndices(ic, jc, kc, tmpIndi, O2CF);
    
      for (E_Int j = 1; j <= 7; j++) indiv(i,j) = tmpIndi[j-1];
    }
    next:;
  }
}

//=============================================================================
/* Corrige a l ordre 2 O2CF les coefs d interpolation et les indi
   Attention : les tableaux ne sont pas redimensionnes
   ex : ordre 5 : cf de taille 15 toujours, mais les 8 premiers coefs seult
   sont valables, idem pour indi.
*/
//=============================================================================
short K_KINTERP::BlkInterpData::correctInterpToO2CF(
  E_Float x, E_Float y, E_Float z, 
  FldArrayI& indi, FldArrayF& cf,
  InterpMeshType interpMeshType)
{
  E_Int ic, jc, kc;    
  FldArrayF tmpcf(8);

  short found = searchInterpolationCellO2CF( x, y, z, ic, jc, kc, tmpcf );
  if ( found < 1) return found;  
  // 1- modif des cf
  for (E_Int i = 0; i < 8; i++)
    cf[i] = tmpcf[i];

  // 2- modif de indi
  if (interpMeshType == EXT_CENTERS)
    fromExtendedToStandardCenters(ic, jc, kc, indi, O2CF);
  else compStandardIndices(ic, jc, kc, indi, O2CF);
  
  return found;
}
// ============= Interp/BlkInterpDataStruct.cpp ================================
