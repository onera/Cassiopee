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
# include "Interp/InterpData.h"

//=============================================================================
/* Destructor */
//=============================================================================
K_INTERP::InterpData::~InterpData()
{
}
// ============================================================================
/* Constructor */
// ============================================================================
K_INTERP::InterpData::InterpData() :
  _EPS_DET(1.e-16), _EPS_TETRA(1.e-4), _EPS_GEOM(K_CONST::E_GEOM_CUTOFF)
{
  _topology = -1;
  _xmin = -K_CONST::E_MAX_FLOAT;
  _ymin = -K_CONST::E_MAX_FLOAT;
  _zmin = -K_CONST::E_MAX_FLOAT;
  _xmax =  K_CONST::E_MAX_FLOAT;
  _ymax =  K_CONST::E_MAX_FLOAT;
  _zmax = K_CONST::E_MAX_FLOAT;
}
K_INTERP::InterpData::InterpData(E_Int topology, E_Float xmin, E_Float ymin, E_Float zmin) :
  _EPS_DET(1.e-16), _EPS_TETRA(1.e-4), _EPS_GEOM(K_CONST::E_GEOM_CUTOFF)
{
  _topology = topology;
  _xmin = xmin;
  _ymin = ymin;
  _zmin = zmin;
  _xmax = K_CONST::E_MAX_FLOAT;
  _ymax = K_CONST::E_MAX_FLOAT;
  _zmax = K_CONST::E_MAX_FLOAT;
}
// ============================================================================
/* Get the hexahedra coordinates+coord of centers of faces+center
   This is the reference for the cell edges numerotation. */
// ============================================================================
void K_INTERP::InterpData::coordHexa(E_Int ind, E_Int ni, E_Int nj, E_Int nk,
                                     E_Float* xl, E_Float* yl, E_Float* zl,
                                     E_Int& ic, E_Int& jc, E_Int& kc,
                                     E_Float* xt, E_Float* yt, E_Float* zt)
{
  E_Int nij = ni*nj;
  kc = ind/nij+1;
  E_Int kcnij = kc*nij;
  E_Int kcmnij = (kc-1)*nij;
  
  jc = (ind-kcmnij)/ni+1;
  E_Int jcni = jc*ni;
  E_Int jcmni = (jc-1)*ni;

  ic = ind-jcmni-kcmnij+1;
  xt[0] = xl[ind]; yt[0] = yl[ind]; zt[0] = zl[ind];

  ind = ic+jcmni+kcmnij;
  xt[1] = xl[ind]; yt[1] = yl[ind]; zt[1] = zl[ind];

  ind = ic-1+jcni+kcmnij;
  xt[2] = xl[ind]; yt[2] = yl[ind]; zt[2] = zl[ind];

  ind = ic+jcni+kcmnij;
  xt[3] = xl[ind]; yt[3] = yl[ind]; zt[3] = zl[ind];

  if ( nk > 1 ) 
  {
    ind = ic-1+jcmni+kcnij;
    xt[4] = xl[ind]; yt[4] = yl[ind]; zt[4] = zl[ind];

    ind = ic+jcmni+kcnij;
    xt[5] = xl[ind]; yt[5] = yl[ind]; zt[5] = zl[ind];

    ind = ic-1+jcni+kcnij;
    xt[6] = xl[ind]; yt[6] = yl[ind]; zt[6] = zl[ind];

    ind = ic+jcni+kcnij;
    xt[7] = xl[ind]; yt[7] = yl[ind]; zt[7] = zl[ind];
  }
  else 
  {
    xt[4] = xt[0]; yt[4] = yt[0]; zt[4] = zt[0]+1.;
    xt[5] = xt[1]; yt[5] = yt[1]; zt[5] = zt[1]+1.;
    xt[6] = xt[2]; yt[6] = yt[2]; zt[6] = zt[2]+1.;
    xt[7] = xt[3]; yt[7] = yt[3]; zt[7] = zt[3]+1.;
  }

          
  /* Compute the center of faces */
  xt[8] = K_CONST::ONE_FOURTH*(xt[0]+xt[1]+xt[2]+xt[3]);
  yt[8] = K_CONST::ONE_FOURTH*(yt[0]+yt[1]+yt[2]+yt[3]);
  zt[8] = K_CONST::ONE_FOURTH*(zt[0]+zt[1]+zt[2]+zt[3]);
  
  xt[9] = K_CONST::ONE_FOURTH*(xt[0]+xt[2]+xt[4]+xt[6]);
  yt[9] = K_CONST::ONE_FOURTH*(yt[0]+yt[2]+yt[4]+yt[6]);
  zt[9] = K_CONST::ONE_FOURTH*(zt[0]+zt[2]+zt[4]+zt[6]);
  
  xt[10] = K_CONST::ONE_FOURTH*(xt[4]+xt[5]+xt[6]+xt[7]);
  yt[10] = K_CONST::ONE_FOURTH*(yt[4]+yt[5]+yt[6]+yt[7]);
  zt[10] = K_CONST::ONE_FOURTH*(zt[4]+zt[5]+zt[6]+zt[7]);
  
  xt[11] = K_CONST::ONE_FOURTH*(xt[1]+xt[3]+xt[5]+xt[7]);
  yt[11] = K_CONST::ONE_FOURTH*(yt[1]+yt[3]+yt[5]+yt[7]);
  zt[11] = K_CONST::ONE_FOURTH*(zt[1]+zt[3]+zt[5]+zt[7]);
  
  xt[12] = K_CONST::ONE_FOURTH*(xt[2]+xt[3]+xt[6]+xt[7]);
  yt[12] = K_CONST::ONE_FOURTH*(yt[2]+yt[3]+yt[6]+yt[7]);
  zt[12] = K_CONST::ONE_FOURTH*(zt[2]+zt[3]+zt[6]+zt[7]);
  
  xt[13] = K_CONST::ONE_FOURTH*(xt[0]+xt[1]+xt[4]+xt[5]);
  yt[13] = K_CONST::ONE_FOURTH*(yt[0]+yt[1]+yt[4]+yt[5]);
  zt[13] = K_CONST::ONE_FOURTH*(zt[0]+zt[1]+zt[4]+zt[5]);
    
  /* Compute the center of the cell */
  xt[14] = K_CONST::ONE_HALF*(xt[8]+xt[10]);
  yt[14] = K_CONST::ONE_HALF*(yt[8]+yt[10]);
  zt[14] = K_CONST::ONE_HALF*(zt[8]+zt[10]);
}
//=============================================================================
/* Find the interpolation coefficients in a given cell. Cut a cell
   in 24 tetrahedras made of 2 edges of the cell, the center of the cell and
   one center of a face. Uses a "jump" technique to find quickly the most 
   probable tetrahedron */
//=============================================================================
E_Bool K_INTERP::InterpData::coeffInterpHexa(E_Float x, E_Float y, E_Float z,
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
  E_Float* cfp = cf.begin();
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
  
  xp = xt[14]; yp = yt[14]; zp = zt[14];
  xq = xt[11]; yq = yt[11]; zq = zt[11];
  xr = xt[12]; yr = yt[12]; zr = zt[12];
  xs = xt[10]; ys = yt[10]; zs = zt[10];

  coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                   xr, yr, zr, xs, ys, zs, xi, yi, zi);

  /* Compute index to choose the most probable edge of cell */
  ixp = 0; iyp = 0; izp = 0;
    
  if ((K_FUNC::E_abs(xi) <= K_CONST::ONE+EPS2) &&
      (K_FUNC::E_abs(yi) <= K_CONST::ONE+EPS2) &&
      (K_FUNC::E_abs(zi) <= K_CONST::ONE+EPS2))
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
#ifdef E_ADOLC
  E_Float vx = (K_CONST::ONE + K_FUNC::E_sign(xi-yi)) * K_CONST::ONE_HALF;
  E_Float vy = (K_CONST::ONE + K_FUNC::E_sign(yi-zi)) * K_CONST::ONE_HALF; 
  E_Float vz = (K_CONST::ONE + K_FUNC::E_sign(zi-xi)) * K_CONST::ONE_HALF; 
  is1 = E_Int(vx.value());
  is2 = E_Int(vy.value());
  is3 = E_Int(vz.value());
#else
  is1 = E_Int((K_CONST::ONE + K_FUNC::E_sign(xi-yi)) * K_CONST::ONE_HALF);
  is2 = E_Int((K_CONST::ONE + K_FUNC::E_sign(yi-zi)) * K_CONST::ONE_HALF);
  is3 = E_Int((K_CONST::ONE + K_FUNC::E_sign(zi-xi)) * K_CONST::ONE_HALF);
#endif

  itri = indtr[isomm+8*is1+16*is2+32*is3];

  indr = indtd[0+itri*4];
  inds = indtd[1+itri*4];
  indq = indtd[2+itri*4];

  /* Most probable tetrahedra */
  xq = xt[indq]; yq = yt[indq]; zq = zt[indq];
  xr = xt[indr]; yr = yt[indr]; zr = zt[indr];
  xs = xt[inds]; ys = yt[inds]; zs = zt[inds];

  coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                   xr, yr, zr, xs, ys, zs, xi, yi, zi);
  
  if ((xi > -EPS) && (yi > -EPS) &&
      (zi > -EPS) && (xi+yi+zi < K_CONST::ONE+3*EPS))
  {
    /* Transforming tetrahedra interpolation coefficient into */
    /* hexahedra interpolation coefficients */
    cf0 = K_CONST::ONE_EIGHTH*(K_CONST::ONE-xi-yi-zi);
    
    cfp[0] = cf0;
    cfp[1] = cf0;
    cfp[2] = cf0;
    cfp[3] = cf0;
    cfp[4] = cf0;
    cfp[5] = cf0;
    cfp[6] = cf0;
    cfp[7] = cf0;

    E_Int indl = (indq-8)*4;
    E_Float tl = xi*K_CONST::ONE_FOURTH;
    ind = indfc[indl];
    cfp[ind] += tl;
    indl++;
    ind = indfc[indl];
    cfp[ind] += tl;
    indl++;
    ind = indfc[indl];
    cfp[ind] += tl;
    indl++;
    ind = indfc[indl];
    cfp[ind] += tl;

    cfp[indr] += yi;
    cfp[inds] += zi;
    
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
    xr = xt[indr]; yr = yt[indr]; zr = zt[indr];
    
    for (its = 0; its < 3; its++) /* neighbour of edge */
    {
      inds = neighbour[its+isom*3];
      xs = xt[inds]; ys = yt[inds]; zs = zt[inds];
      
      for (itq = 0; itq < 2; itq++) /* center */
      {
        indq = center[itq+its+isom*4];
        xq = xt[indq]; yq = yt[indq]; zq = zt[indq];
        
        coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                         xr, yr, zr, xs, ys, zs, xi, yi, zi);

        if ((xi > -EPS)&&(yi > -EPS)&&
            (zi > -EPS)&&(xi+yi+zi < K_CONST::ONE+3*EPS))
        {
          /* Transforming to hexahedra coefficients */ 
          cf0 = K_CONST::ONE_EIGHTH*(1.-xi-yi-zi);
          for (i = 0; i < 8; i++) cfp[i] = cf0;
          
          for (i = 0; i < 4; i++)
          {
            ind = indfc[i+(indq-8)*4];
            cfp[ind] += xi*K_CONST::ONE_FOURTH;
          }
          
          cfp[indr] += yi;
          cfp[inds] += zi;
          
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
E_Bool K_INTERP::InterpData::getCellJump(E_Float x, E_Float y, E_Float z,
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

  xp = xt[14]; yp = yt[14]; zp = zt[14];
  xq = xt[11]; yq = yt[11]; zq = zt[11];
  xr = xt[12]; yr = yt[12]; zr = zt[12];
  xs = xt[10]; ys = yt[10]; zs = zt[10];

  coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                   xr, yr, zr, xs, ys, zs, xi, yi, zi);

  /* Compute index to choose the most probable edge of cell */
  ixp = 0; iyp = 0; izp = 0;
  if ((K_FUNC::E_abs(xi) <= K_CONST::ONE+EPS)&&
      (K_FUNC::E_abs(yi) <= K_CONST::ONE+EPS)&&
      (K_FUNC::E_abs(zi) <= K_CONST::ONE+EPS))
  {
    if (xi >= 0) ixp = 1;
    if (yi >= 0) iyp = 1;
    if (zi >= 0) izp = 1;
    isomm = indsom[ixp+2*iyp+4*izp];
    return true;
  }
  else return false;
}
//=============================================================================
/* Find the interp. coeff of point (x,y,z) in the given tetrahedron */
//=============================================================================
void K_INTERP::InterpData::coeffInterpTetra(E_Float x, E_Float y, E_Float z,
                                            E_Float xp, E_Float yp, E_Float zp,
                                            E_Float xq, E_Float yq, E_Float zq,
                                            E_Float xr, E_Float yr, E_Float zr,
                                            E_Float xs, E_Float ys, E_Float zs,
                                            E_Float& xi, E_Float& yi, E_Float& zi)
{
  E_Float a11, a12, a13, a21, a22, a23, a31, a32, a33;
  E_Float c11, c12, c13, c21, c22, c23, c31, c32, c33;
  E_Float det;
  E_Float xpm, ypm, zpm;
  const E_Float EPS = _EPS_DET;
  
  /* Computation of the coefficient of transfer matrix */
  a11 = xq-xp;
  a12 = xr-xp;
  a13 = xs-xp;
  a21 = yq-yp;
  a22 = yr-yp;
  a23 = ys-yp;
  a31 = zq-zp;
  a32 = zr-zp;
  a33 = zs-zp;
  
  /* Computation of the coefficient of the comatrix */
  c11 =   a22*a33-a32*a23;
  c12 = -(a21*a33-a31*a23);
  c13 =   a21*a32-a31*a22;
  c21 = -(a12*a33-a32*a13);
  c22 =   a11*a33-a31*a13;
  c23 = -(a11*a32-a31*a12);
  c31 =   a12*a23-a22*a13; 
  c32 = -(a11*a23-a21*a13);
  c33 =   a11*a22-a21*a12;
    
  det = a11*c11+a12*c12+a13*c13;

  /* When det is null, the routine should declare */
  /* this tetrahedra as not candidate for interpolation */
  if (K_FUNC::E_abs(det) < EPS)
  {
    xi = -10.;
    yi = -10.;
    zi = -10.;
  }
  else
  {
    xpm = x-xp;
    ypm = y-yp;
    zpm = z-zp;
  
    det = K_CONST::ONE/det;
    xi = (c11*xpm+c21*ypm+c31*zpm)*det;
    yi = (c12*xpm+c22*ypm+c32*zpm)*det;
    zi = (c13*xpm+c23*ypm+c33*zpm)*det;
  }
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
E_Bool K_INTERP::InterpData::getCoeffInterpHexa(E_Float x, E_Float y, E_Float z,
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
  E_Float* cfp = cf.begin();
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
      xi = -xi; yi = -yi; zi = -zi;
      break;
      
    case 1:
      yi = -yi; zi = -zi;
      break;
      
    case 2:
      xi = -xi; zi = -zi;
      break;
      
    case 3:
      zi = -zi;
      break;
      
    case 4:
      xi = -xi; yi = -yi;
      break;
      
    case 5:
      yi = -yi;
      break;
      
    case 6:
      xi = -xi;
      break;
  }
    
  /* Compute index to find the most probable triangle */
#ifdef E_ADOLC
  E_Float vx = (K_CONST::ONE + K_FUNC::E_sign(xi-yi)) * K_CONST::ONE_HALF;
  E_Float vy = (K_CONST::ONE + K_FUNC::E_sign(yi-zi)) * K_CONST::ONE_HALF;
  E_Float vz = (K_CONST::ONE + K_FUNC::E_sign(zi-xi)) * K_CONST::ONE_HALF;
  E_Int is1 = static_cast<E_Int>(vx.value());
  E_Int is2 = static_cast<E_Int>(vy.value());
  E_Int is3 = static_cast<E_Int>(vz.value());
#else
  E_Int is1 = static_cast<E_Int>((K_CONST::ONE + K_FUNC::E_sign(xi-yi)) * K_CONST::ONE_HALF);
  E_Int is2 = static_cast<E_Int>((K_CONST::ONE + K_FUNC::E_sign(yi-zi)) * K_CONST::ONE_HALF);
  E_Int is3 = static_cast<E_Int>((K_CONST::ONE + K_FUNC::E_sign(zi-xi)) * K_CONST::ONE_HALF);
#endif

  E_Int itri = indtr[isomm+8*is1+16*is2+32*is3];

  E_Int indr = indtd[0+itri*4];
  E_Int inds = indtd[1+itri*4];
  E_Int indq = indtd[2+itri*4];
  E_Int indp = indtd[3+itri*4];

  /* Most probable tetrahedra */
  xp = xt[indp]; yp = yt[indp]; zp = zt[indp];
  xq = xt[indq]; yq = yt[indq]; zq = zt[indq];
  xr = xt[indr]; yr = yt[indr]; zr = zt[indr];
  xs = xt[inds]; ys = yt[inds]; zs = zt[inds];

  coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                   xr, yr, zr, xs, ys, zs, xi, yi, zi);

  if ((xi > -EPS)&&(yi > -EPS)&&
      (zi > -EPS)&&(xi+yi+zi < K_CONST::ONE+3*EPS))
  {
    /* transforming tetrahedra interpolation coefficient into */
    /* hexahedra interpolation coefficients */
      
    cf0 = K_CONST::ONE_EIGHTH*(K_CONST::ONE-xi-yi-zi);
    for (i = 0; i < 8; i++) cfp[i] = cf0;
    
    for (i = 0; i < 4; i++)
    {
      ind = indfc[i+(indq-8)*4];
      cfp[ind] += xi*K_CONST::ONE_FOURTH;
    }
      
    cfp[indr] += yi;
    cfp[inds] += zi;
      
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
    xr = xt[indr]; yr = yt[indr]; zr = zt[indr];
    
    for (E_Int its = 0; its < 3; its++) /* neighbour of edge */
    {
      inds = neighbour[its+isom*3];
      xs = xt[inds]; ys = yt[inds]; zs = zt[inds];
      
      for (E_Int itq = 0; itq < 2; itq++) /* center */
      {
        indq = center[itq+its+isom*4];
        xq = xt[indq]; yq = yt[indq]; zq = zt[indq];
        
        coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                         xr, yr, zr, xs, ys, zs, xi, yi, zi);


        if ((xi > -EPS)&&(yi > -EPS) &&
            (zi > -EPS)&&(xi+yi+zi < K_CONST::ONE+3*EPS))
        {
          /* transforming to hexahedra coefficients */ 
          cf0 = K_CONST::ONE_EIGHTH*(1.-xi-yi-zi);
          for (i = 0; i < 8; i++) cfp[i] = cf0;

          for (i = 0; i < 4; i++)
          {
            ind = indfc[i+(indq-8)*4];
            cfp[ind] += xi*K_CONST::ONE_FOURTH;
          }
          
          cfp[indr] += yi;
          cfp[inds] += zi;

          return true;
        }           
      }
    }
  }
  return false; /* Cannot interpolate from this cell */
}

//=============================================================================
/* Compute the extrapolation coefficients (storage OiABC) for point (x,y,z) 
   inside cell (ic,jc,kc) for structured meshes.
   It is based on the coefficients in the best valid tetrahedra. 
   In this version (default), the mesh in centers is supposed to be known.
   IN: (x,y,z): coordonnees du point a extrapoler
   IN: (ic,jc,kc): indices de la cellule a tester (debut a 1)
   IN: cellNatureField: champ de la nature Chimere des cellules sur la grille d interpolation
   IN: nature: type de cellN a tester
   IN: constraint: contrainte sue la val abs de la somme des coeff
   OUT: cf: coefficient d'extrapolation */
//=============================================================================
short K_INTERP::InterpData::getExtrapolationCoeffForCell(
  E_Float x, E_Float y, E_Float z,
  E_Int ic, E_Int jc, E_Int kc,
  FldArrayF& cf,
  E_Int ni, E_Int nj, E_Int nk,
  E_Float* xl, E_Float* yl, E_Float* zl,
  E_Float* cellNp, 
  E_Int& is, E_Int& js, E_Int& ks,
  E_Int nature, E_Float constraint, E_Float& diff_coeff)
{
  E_Float min_coeff = K_CONST::E_MAX_FLOAT;
  E_Float max_coeff =-K_CONST::E_MAX_FLOAT;

  E_Int dim = 3; 
  E_Int ncells = 8;
  if (nk == 1) { dim = 2; ncells = 4;}

  // vecteur pour stocker les coordonnees du centre, des sommets et des centres des faces du hexa d'interpolation
  E_Float xt[15]; E_Float yt[15]; E_Float zt[15];
  E_Float* cfp = cf.begin();
  // variables pour le calcul des coefficients tetraedriques
  E_Float xi, yi, zi;
  E_Float xp, yp, zp;
  E_Float xr, yr, zr;
  E_Float xs, ys, zs;
  E_Float xq, yq, zq;

  // Coord of interpolation cell
  ic = ic-1; jc = jc-1; kc = K_FUNC::E_max(0,kc-1);
  E_Int ind = ic + jc*ni + kc*ni*nj;
  E_Int icdummy,jcdummy,kcdummy;
  coordHexa(ind,ni,nj,nk,xl,yl,zl,icdummy,jcdummy,kcdummy,xt,yt,zt);

  // Try to find a correct tetrahedra
  // that is a tetrahedra where cellNatureField is equal to 1.
  E_Int isom, its, inds, indq, indr, itq, ibsom, isomm;
  E_Float cf0;
  
  /* for each edge of the cell, the index of the neigbour edges (3 x 8) */ 
  static E_Int neighbour[24] = {1,2,4,0,3,5,3,0,6,2,1,7,5,0,6,4,7,1,7,4,
                                2,6,5,3};
  
  /* for each edge and its neighbour, the proper center of face for tetrahedra */
  static E_Int center[32] = {13,8,9,13,13,8,11,13,12,8,9,12,12,8,11,12,10,13,9,
                             10,13,10,11,13,12,10,9,12,12,10,11,12};
  /* (4 x 8) */
  static E_Int indss[32] = {0,3,5,6,1,2,4,7,2,4,7,1,3,5,6,0,4,7,1,2,
                            5,6,0,3,6,0,3,5,7,1,2,4};
  
  /* for each face of the cell, the number of edges (4 x 6) */
  static E_Int indfc[24] = {0,1,3,2,4,0,2,6,4,5,7,6,5,1,3,7,6,7,3,2,4,5,1,0};

  E_Float cfSumMin = K_CONST::E_MAX_FLOAT;
  E_Float cfSum;
  E_Float cfSumMin2 = K_CONST::E_MAX_FLOAT;
  E_Int ncf = 8;
  FldArrayF cfSav(ncf);
  cfSav.setAllValuesAtNull();
  FldArrayF cfSav2(ncf);
  cfSav2.setAllValuesAtNull();
  E_Float sum, sum2;
  E_Float report, cellN0;
  E_Float val = 1.;
  E_Bool valid = false;
  xp = xt[14]; yp = yt[14]; zp = zt[14];

  // test les 24 tetraedres a la recherche d'un tetraedre
  // valide pour l'interpolation (ne possedant pas de points
  // interpoles) et du tetraedre presentant le min de la somme des 
  // val. abs des coeff.
  isomm = 0;
  for (ibsom = 0; ibsom < 4; ibsom++) /* edges */
  {
    isom = indss[ibsom+4*isomm];
    indr = isom;
    xr = xt[indr]; yr = yt[indr]; zr = zt[indr];
    
    for (its = 0; its < 3; its++) /* neighbour of edge */
    {
      inds = neighbour[its+isom*3];
      xs = xt[inds]; ys = yt[inds]; zs = zt[inds];
      
      for (itq = 0; itq < 2; itq++) /* center */
      {
        indq = center[itq+its+isom*4];
        xq = xt[indq]; yq = yt[indq]; zq = zt[indq];
        coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                         xr, yr, zr, xs, ys, zs, xi, yi, zi);       
        /* transforming to hexahedra coefficients */ 
        cf0 = K_CONST::ONE_EIGHTH*(K_CONST::ONE-xi-yi-zi);
        for (E_Int i = 0; i < 8; i++) cfp[i] = cf0;
        
        for (E_Int i = 0; i < 4; i++)
        {
          ind = indfc[i+(indq-8)*4];
          cfp[ind] += xi*K_CONST::ONE_FOURTH;
        }
        cfp[indr] += yi;
        cfp[inds] += zi;
        cfSum = K_FUNC::E_abs(cfp[0])+K_FUNC::E_abs(cfp[1])+K_FUNC::E_abs(cfp[2])+K_FUNC::E_abs(cfp[3])+
        K_FUNC::E_abs(cfp[4])+K_FUNC::E_abs(cfp[5])+K_FUNC::E_abs(cfp[6])+K_FUNC::E_abs(cfp[7]);
        
        // keep the tetrahedra with smallest sum of abs. coeff.
        if (cfSum < cfSumMin2)
        {
          cfSumMin2 = cfSum;
          cfSav2 = cf; // copy the best extrap. coeff (best tetra)
        }

        if ((xi+yi+zi <= K_CONST::ONE)&&(xi >= 0.)&&(yi >= 0.)&&(zi >= 0.))
        {
          // indeed the best (it contains P)
          cfSumMin2 = 0.;
          cfSav2 = cf;
        }

        val = 1.;
        if (cellNp != NULL)
        {
          if ( dim == 3)
          {
          // criterion: all cells are cellN=1
          for (E_Int kk = 0; kk < 2; kk++)
            for (E_Int jj = 0; jj < 2; jj++)
              for (E_Int ii = 0; ii < 2; ii++)
              {
                ind = (ic+ii) + (jc+jj)*ni + (kc+kk)*ni*nj;
                cellN0 = cellNp[ind];
                if (nature == 0) val *= cellN0; // pas de pt masque ds la molecule donneuse
                else val *= cellN0*(cellN0-2.); // pas de pt masque ou interpole dans la cellule donneuse
              }
          }
          else
          {
            // criterion: all cells are cellN=1
            for (E_Int jj = 0; jj < 2; jj++)
              for (E_Int ii = 0; ii < 2; ii++)
              {
                ind = (ic+ii) + (jc+jj)*ni;
                cellN0 = cellNp[ind];
                if (nature == 0) val *= cellN0; // pas de pt masque ds la molecule donneuse
                else val *= cellN0*(cellN0-2.); // pas de pt masque ou interpole dans la cellule donneuse
              }
          }
        }
        if (K_FUNC::fEqualZero(val, K_CONST::E_GEOM_CUTOFF) == false) // valid 
        {
          // keep the best valid tetra
          if (cfSum < cfSumMin)
          {
            cfSumMin = cfSum;
            cfSav = cf;
            valid = true;
          }
        }
      }       
    }
  }

  if (valid == true && cfSumMin < constraint) // OK, We found a good valid tetra
  {
    cf = cfSav;
    min_coeff = K_CONST::E_MAX_FLOAT;
    max_coeff =-K_CONST::E_MAX_FLOAT;
    for (E_Int nocf = 0; nocf < 8; nocf++)
    {
      if (K_FUNC::E_abs(cf[nocf])>1.e-10)
      {
        min_coeff = K_FUNC::E_min(min_coeff,cf[nocf]);
        max_coeff = K_FUNC::E_max(max_coeff,cf[nocf]);
      }
    }
    diff_coeff = max_coeff-min_coeff;
    return 1;
  }
  else if (cfSumMin2 < constraint)
  {
    // Une cellule est trouvee, mais elle contient des cellN invalides
    cf = cfSav2; 
    min_coeff = K_CONST::E_MAX_FLOAT;
    max_coeff =-K_CONST::E_MAX_FLOAT;
    for (E_Int nocf = 0; nocf < 8; nocf++)
    {
      if (K_FUNC::E_abs(cf[nocf])>1.e-10)
      {
        min_coeff = K_FUNC::E_min(min_coeff,cf[nocf]);
        max_coeff = K_FUNC::E_max(max_coeff,cf[nocf]);
      }
    }
    diff_coeff = max_coeff-min_coeff;
    // On essaie de reporter les coeffs sur les sommets valides
    E_Int icell = 0;
    E_Float locCellN[8]; // 1: on le garde, 0: on l'enleve
    for (E_Int i = 0; i < 8; i++) locCellN[i] = 1.;
      
    if (cellNp != NULL)
    {
      sum = 0.;
      if (dim == 3)
      {
        for (E_Int kk = 0; kk < 2; kk++)
          for (E_Int jj = 0; jj < 2; jj++)
            for (E_Int ii = 0; ii < 2; ii++)
            {
              E_Int ind = (ic + ii) + (jc + jj)*ni + (kc+kk)*ni*nj;
              if (nature == 1) locCellN[icell] = 1.-K_FUNC::E_abs(cellNp[ind]-1.);
              else if (nature == 0) locCellN[icell] = K_FUNC::E_min(1.,cellNp[ind]);
              sum += locCellN[icell];
              icell++;
            }
      }
      else
      {
        for (E_Int jj = 0; jj < 2; jj++)
          for (E_Int ii = 0; ii < 2; ii++)
          {
            E_Int ind = (ic + ii) + (jc + jj)*ni;
            if (nature == 1) locCellN[icell] = 1.-K_FUNC::E_abs(cellNp[ind]-1.);
            else if (nature == 0) locCellN[icell] = K_FUNC::E_min(1.,cellNp[ind]);
            sum += locCellN[icell];
            icell++;
          }       
      }

    }
    else sum = 8.;

    if (sum < K_CONST::E_CUTOFF) goto next; // sum=0 (que des pts invalides)
    if ( nk == 1)//report des coefs en 2D
    {
      for (E_Int nocf = 0; nocf < 4; nocf++)
      {
        cfp[nocf] += cfp[nocf+4];
        cfp[nocf+4] = 0.;
      }
    }
    sum2 = 0.;
    for (icell = 0; icell < ncells; icell++)
      sum2 += (1.-locCellN[icell])*cfp[icell];

    report = sum2 / sum;
    for (icell = 0; icell < ncells; icell++)
      cfp[icell] = (cfp[icell]+report)*locCellN[icell];
    return 1;
  }

  next:
  return 0; /* Cannot extrapolate from this cell */
}

