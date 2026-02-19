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
# include "Interp/InterpCart.h"
# include "Def/DefFunction.h"

using namespace std;
//using namespace K_FLD;

//=============================================================================
/* Destructor */
//=============================================================================
K_INTERP::InterpCart::~InterpCart()
{
}

//=============================================================================
/* Constructor 
  IN: ni,nj,nk: nbre de pts de la grille cartesienne reguliere
  IN: hi,hj,hk: le pas dans chaque direction
  IN: x0,y0,z0: les coords du premier point
  IN: ioff,joff,koff: par defaut 0
  Si la grille est une grille cartesienne reguliere avec 2 ghost cells de chaque
  cote, il faut mettre ioff=2, ni au nbre de pts-2, x0 a x[2] et hi a x[3]-x[2].
  */
//=============================================================================
K_INTERP::InterpCart::InterpCart(E_Int ni, E_Int nj, E_Int nk,
                                 E_Float hi, E_Float hj, E_Float hk,
                                 E_Float x0, E_Float y0, E_Float z0,
                                 E_Int ioff, E_Int joff, E_Int koff):
   InterpData(0, x0, y0, z0), _ni(ni), _nj(nj), _nk(nk)
{
  _hi = hi;
  _hj = hj;
  _hk = hk;
  _hii = 1.0/_hi;
  _hji = 1.0/_hj;
  _hki = 1.0/_hk;
  _his2 = 0.5*hi;
  _hjs2 = 0.5*hj;
  _hks2 = 0.5*hk;
  _xmax = _xmin+(_ni-1-ioff)*_hi;
  _ymax = _ymin+(_nj-1-joff)*_hj;
  _zmax = _zmin+(_nk-1-koff)*_hk;
  _ioff = ioff+1; // offset si decalage de la grille (ghostcells)
  _joff = joff+1;
  _koff = koff+1;
}

//=============================================================================
/* 2nd order interpolation - (ic,jc,kc) demarre a 1 - type O2CF */
//=============================================================================
short K_INTERP::InterpCart::searchInterpolationCellCartO2(E_Int ni, E_Int nj, E_Int nk,
                                                          E_Float x, E_Float y, E_Float z,
                                                          E_Int& ic, E_Int& jc, E_Int& kc,
                                                          FldArrayF& cf)
{ 
  const E_Float EPS = K_CONST::E_GEOM_CUTOFF;
  if (x < _xmin-EPS) return 0;
  if (x > _xmax+EPS) return 0;
  if (y < _ymin-EPS) return 0;
  if (y > _ymax+EPS) return 0;
  if (z < _zmin-EPS) return 0;
  if (z > _zmax+EPS) return 0;

#ifdef E_ADOLC
  E_Float vx = (x-_xmin)*_hii;
  E_Float vy = (y-_ymin)*_hji;
  E_Float vz = (z-_zmin)*_hki;
  ic = E_Int(vx.value())+_ioff;
  jc = E_Int(vy.value())+_joff;
  kc = E_Int(vz.value())+_koff;
#else
  ic = E_Int((x-_xmin)*_hii)+_ioff;
  jc = E_Int((y-_ymin)*_hji)+_joff;
  kc = E_Int((z-_zmin)*_hki)+_koff;
#endif

  ic = std::min(ic,ni-1); ic = std::max(ic,E_Int(1)); 
  jc = std::min(jc,nj-1); jc = std::max(jc,E_Int(1)); 
  kc = std::min(kc,nk-1); kc = std::max(kc,E_Int(1));
        
  E_Float* cfp = cf.begin();
  
  E_Float x0 = _xmin+(ic-1)*_hi;
  E_Float y0 = _ymin+(jc-1)*_hj;
  E_Float z0 = _zmin+(kc-1)*_hk;

  E_Float x1 = (x-x0)*_hii;
  E_Float x2 = (y-y0)*_hji;
  E_Float x3 = (z-z0)*_hki;

  E_Float t1 = (1.-x1);
  E_Float t2 = (1.-x2);
  E_Float t3 = (1.-x3);

  E_Float t1t2 = t1*t2;
  E_Float x2t3 = x2*t3;
  E_Float x1t2 = x1*t2;
  E_Float x2x3 = x2*x3;

  cfp[0] = t1t2*t3;
  cfp[1] = x1t2*t3;
  cfp[2] = t1*x2t3;
  cfp[3] = x1*x2t3;

  if (_nk <=2 && K_FUNC::E_abs(_zmax-_zmin)<EPS)
  {
    cfp[4]=0.; cf[5]=0.; cfp[6]=0.; cfp[7]=0.;
  }
  else 
  {
    cfp[4] = t1t2*x3;
    cfp[5] = x1t2*x3;
    cfp[6] = t1*x2x3;
    cfp[7] = x1*x2x3;
  }
  return 1;
}

// ============================================================================
/* 3rd order interpolation on Cartesian grids 
  Returns coefs stored by direction (O3ABC) 
  e.g. : coef(P_0) = cf[0]*cf[3]*cf[6] */
// ============================================================================
short K_INTERP::InterpCart::searchInterpolationCellCartO3(E_Int ni, E_Int nj, E_Int nk,
                                                          E_Float x, E_Float y, E_Float z,
                                                          E_Int& ic, E_Int& jc, E_Int& kc,
                                                          FldArrayF& cf)
{ 
  E_Int icHO, jcHO, kcHO;
  const E_Float EPS = K_CONST::E_GEOM_CUTOFF;
  if (x < _xmin-EPS) return 0;
  if (x > _xmax+EPS) return 0;
  if (y < _ymin-EPS) return 0;
  if (y > _ymax+EPS) return 0;
  if (z < _zmin-EPS) return 0;
  if (z > _zmax+EPS) return 0;
  
#ifdef E_ADOLC
  E_Float vx = (x-_xmin)*_hii;
  E_Float vy = (y-_ymin)*_hji;
  E_Float vz = (z-_zmin)*_hki;  
  ic = E_Int(vx.value())+_ioff;
  jc = E_Int(vy.value())+_joff;
  kc = E_Int(vz.value())+_koff;
#else
  ic = E_Int((x-_xmin)*_hii)+_ioff;
  jc = E_Int((y-_ymin)*_hji)+_joff;
  kc = E_Int((z-_zmin)*_hki)+_koff;
#endif

  ic = std::max(ic,E_Int(1)); ic = std::min(ic,ni-1);
  jc = std::max(jc,E_Int(1)); jc = std::min(jc,nj-1);
  kc = std::max(kc,E_Int(1)); kc = std::min(kc,nk-1);

  //E_Float* cfp = cf.begin();
  
  E_Float x_i = _xmin+(ic-1)*_hi; E_Float x_ip = x_i +_hi;
  E_Float y_j = _ymin+(jc-1)*_hj; E_Float y_jp = y_j +_hj;
  E_Float z_k = _zmin+(kc-1)*_hk; E_Float z_kp = z_k +_hk;
  E_Float x_iHO, x_ipHO, x_ippHO;
  E_Float y_jHO, y_jpHO, y_jppHO;
  E_Float z_kHO, z_kpHO, z_kppHO;

  // Interpolation cell
  if (x-x_i >= _his2) // cellule d'interpolation à droite
  {
    if (ic >= _ni-1)
    {
      icHO = ic-1;
      x_iHO = x_i-_hi;
      x_ipHO = x_i;
      x_ippHO = x_ip;
    }
    else
    {
      icHO = ic;
      x_iHO = x_i;
      x_ipHO = x_ip;
      x_ippHO = x_ipHO+_hi;
    }
  }
  else
  {
    if (ic < 3)
    {
      icHO = ic;
      x_iHO = x_i;
      x_ipHO = x_ip;
      x_ippHO = x_ipHO + _hi;
    }
    else
    {
      icHO = ic-1;
      x_iHO = x_i - _hi;
      x_ipHO = x_i;
      x_ippHO = x_ip;
    } 
  }
  
  if ( y-y_j >= _hjs2)
  {
    if (jc >= _nj-1)
    {
      jcHO = jc-1;
      y_jHO = y_j - _hj;
      y_jpHO = y_j;
      y_jppHO = y_jp;
    }
    else
    {
      jcHO = jc;
      y_jHO = y_j;
      y_jpHO = y_jp;
      y_jppHO = y_jpHO + _hj;
    }
  }
  else
  {
    if (jc < 3)
    {
      jcHO = jc;
      y_jHO = y_j;
      y_jpHO = y_jp;
      y_jppHO = y_jpHO + _hj;
    }
    else
    {
      jcHO = jc-1;
      y_jHO = y_j - _hj;
      y_jpHO = y_j;
      y_jppHO = y_jp;
    } 
  }

  if (z-z_k >= _hks2)
  {
    if (kc >= _nk-1)
    {
      kcHO = kc-1;
      z_kHO = z_k - _hk;
      z_kpHO = z_k;
      z_kppHO = z_kp;
    }
    else
    {
      kcHO = kc;
      z_kHO = z_k;
      z_kpHO = z_kp;
      z_kppHO = z_kpHO + _hk;
    }
  }
  else
  {
    if (kc < 3)
    {
      kcHO = kc;
      z_kHO = z_k;
      z_kpHO = z_kp;
      z_kppHO = z_kpHO + _hk;
    }
    else
    {
      kcHO = kc-1;
      z_kHO = z_k - _hk;
      z_kpHO = z_k;
      z_kppHO = z_kp;
    } 
  }

  if (_nk == 2) // dimension 2
  {
    kcHO = 1;
    if (kc == 1)
    {
      z_kHO = z_k;
      z_kpHO = z_kp;
      z_kppHO = z_kp + _hks2;
    }
    else
    {
      z_kHO = z_k - _hks2;
      z_kpHO = z_k;
      z_kppHO = z_kp;
    }    
  }
  
  E_Float hciHO = x_ipHO-x_iHO;
  E_Float hcjHO = y_jpHO-y_jHO;
  E_Float hckHO = z_kpHO-z_kHO;
 
  E_Float hcipHO = x_ippHO-x_ipHO;
  E_Float hcjpHO = y_jppHO-y_jpHO;
  E_Float hckpHO = z_kppHO-z_kpHO;

  E_Float hcippHO = x_ippHO-x_iHO;
  E_Float hcjppHO = y_jppHO-y_jHO;
  E_Float hckppHO = z_kppHO-z_kHO;
  
  E_Float inv_hciHO = 1./hciHO;
  E_Float inv_hcipHO = 1./hcipHO;
  E_Float inv_hcippHO = 1./hcippHO;

  E_Float inv_hcjHO = 1./hcjHO;
  E_Float inv_hcjpHO = 1./hcjpHO;
  E_Float inv_hcjppHO = 1./hcjppHO;

  E_Float inv_hckHO = 1./hckHO;
  E_Float inv_hckpHO = 1./hckpHO;
  E_Float inv_hckppHO = 1./hckppHO;

  E_Float lx0 = (x-x_ipHO)*(x-x_ippHO)*inv_hciHO*inv_hcippHO;
  E_Float lx1 =-(x-x_iHO) *(x-x_ippHO)*inv_hciHO*inv_hcipHO;
  E_Float lx2 = (x-x_iHO) *(x-x_ipHO) *inv_hcipHO*inv_hcippHO;

  E_Float ly0 = (y-y_jpHO)*(y-y_jppHO)*inv_hcjHO*inv_hcjppHO;
  E_Float ly1 =-(y-y_jHO) *(y-y_jppHO)*inv_hcjHO*inv_hcjpHO;
  E_Float ly2 = (y-y_jHO) *(y-y_jpHO) *inv_hcjpHO*inv_hcjppHO;

  E_Float lz0 = (z-z_kpHO)*(z-z_kppHO)*inv_hckHO*inv_hckppHO;
  E_Float lz1 =-(z-z_kHO) *(z-z_kppHO)*inv_hckHO*inv_hckpHO;
  E_Float lz2 = (z-z_kHO) *(z-z_kpHO) *inv_hckpHO*inv_hckppHO;

  cf[0] = lx0; cf[1] = lx1; cf[2] = lx2;
  cf[3] = ly0; cf[4] = ly1; cf[5] = ly2;
  cf[6] = lz0; cf[7] = lz1; cf[8] = lz2;
  ic = icHO; jc = jcHO; kc = kcHO;
  return 1;
}

// ============================================================================
/* 4rd order interpolation on Cartesian grids 
  Returns coefs stored by direction (O4ABC) 
  e.g. : coef(P_0) = cf[0]*cf[3]*cf[6] */
// ============================================================================
short K_INTERP::InterpCart::searchInterpolationCellCartO4(E_Int ni, E_Int nj, E_Int nk,
                                                          E_Float x, E_Float y, E_Float z,
                                                          E_Int& ic, E_Int& jc, E_Int& kc,
                                                          FldArrayF& cf)
{
  const E_Float EPS = K_CONST::E_GEOM_CUTOFF;
  if (x < _xmin-EPS) return 0;
  if (x > _xmax+EPS) return 0;
  if (y < _ymin-EPS) return 0;
  if (y > _ymax+EPS) return 0;
  if (z < _zmin-EPS) return 0;
  if (z > _zmax+EPS) return 0;

#ifdef E_ADOLC
  E_Float vx = (x-_xmin)*_hii;
  E_Float vy = (y-_ymin)*_hji;
  E_Float vz = (z-_zmin)*_hki;
  ic = E_Int(vx.value())+_ioff;
  jc = E_Int(vy.value())+_joff;
  kc = E_Int(vz.value())+_koff;
#else
  ic = E_Int((x-_xmin)*_hii)+_ioff;
  jc = E_Int((y-_ymin)*_hji)+_joff;
  kc = E_Int((z-_zmin)*_hki)+_koff;
#endif

  ic = std::max(ic,E_Int(1)); ic = std::min(ic,ni-1);
  jc = std::max(jc,E_Int(1)); jc = std::min(jc,nj-1);
  kc = std::max(kc,E_Int(1)); kc = std::min(kc,nk-1);

  //Dir X(I)
  if (ic >1 ) ic -=1; //on decale le donneur de 1
  if (ic+3 >  _ni) ic -= ic+3-_ni;

  E_Float x_i = _xmin+(ic-1)*_hi; E_Float x_ip = x_i +_hi;  E_Float x_ipp = x_i + 2*_hi; E_Float x_ippp = x_i + 3*_hi;

  //xi
  E_Float norm= -6.*_hi*_hi*_hi;
  E_Float lx0 = (x-x_ip)*(x-x_ipp)*(x-x_ippp)/norm;
  //xip
  norm=  2.*_hi*_hi*_hi;
  E_Float lx1 =( x-x_i) *(x-x_ipp)*(x-x_ippp)/norm;
  //xipp
  norm= -2.*_hi*_hi*_hi;
  E_Float lx2 = (x-x_i)*(x-x_ip)*(x-x_ippp)/norm;
  //xippp
  norm=  6.*_hi*_hi*_hi;
  E_Float lx3 = (x-x_i)*(x-x_ip)*(x-x_ipp)/norm;

  //Dir Y(J)
  //if (jc >= _nj-2) jc -=2;
  if (jc >1 ) jc -=1; //on decale le donneur de 1
  if (jc+3 >  _nj) jc -= jc+3-_nj;

  E_Float y_j = _ymin+(jc-1)*_hj; E_Float y_jp = y_j +_hj;  E_Float y_jpp = y_j + 2*_hj; E_Float y_jppp = y_j + 3*_hj;
  //yj
  norm=-6.*_hj*_hj*_hj;
  E_Float ly0 = (y-y_jp)*(y-y_jpp)*(y-y_jppp)/norm;
  //yjp
  norm= 2.*_hj*_hj*_hj;
  E_Float ly1 = (y-y_j) *(y-y_jpp)*(y-y_jppp)/norm;
  //yjpp
  norm=-2.*_hj*_hj*_hj;
  E_Float ly2 =(y-y_j)*(y-y_jp)*(y-y_jppp)/norm;
  //yjppp
  norm= 6.*_hj*_hj*_hj;
  E_Float ly3 = (y-y_j)*(y-y_jp)*(y-y_jpp)/norm;

  //Dir Z(K)
  if (kc >1 ) kc -=1; //on decale le donneur de 1
  if (kc+3 >  _nk) kc -= kc+3-_nk;

  E_Float z_k = _zmin+(kc-1)*_hk; E_Float z_kp = z_k +_hk;  E_Float z_kpp = z_k + 2*_hk; E_Float z_kppp = z_k + 3*_hk;
  //zk
  norm=-6.*_hk*_hk*_hk;
  E_Float lz0 = (z-z_kp)*(z-z_kpp)*(z-z_kppp)/norm;
  //zkp
  norm= 2.*_hk*_hk*_hk;
  E_Float lz1 = (z-z_k) *(z-z_kpp)*(z-z_kppp)/norm;
  //zkpp
  norm=-2.*_hk*_hk*_hk;
  E_Float lz2 = (z-z_k)*(z-z_kp)*(z-z_kppp)/norm;
  //zkppp
  norm= 6.*_hk*_hk*_hk;
  E_Float lz3 = (z-z_k)*(z-z_kp)*(z-z_kpp)/norm;

  cf[0] = lx0; cf[1] = lx1; cf[2 ] = lx2;  cf[3] = lx3;
  cf[4] = ly0; cf[5] = ly1; cf[6 ] = ly2;  cf[7] = ly3;
  cf[8] = lz0; cf[9] = lz1; cf[10] = lz2;  cf[11] =lz3;
 
  return 1;
}

// =====================================================================================
/* Retourne la liste des cellules contenant (x,y,z) a la tolerance alphaTol pres */
// =====================================================================================
E_Int K_INTERP::InterpCart::getListOfCandidateCells(E_Float x, E_Float y, E_Float z,
                                                    list<E_Int>& listIndices, 
                                                    E_Float alphaTol)
{
  //E_Float dx = 0.1*alphaTol*(_xmax-_xmin);
  //E_Float dy = 0.1*alphaTol*(_ymax-_ymin);
  //E_Float dz = 0.1*alphaTol*(_zmax-_zmin);
  E_Float dx = 0.1*alphaTol*(_hi);
  E_Float dy = 0.1*alphaTol*(_hj);
  E_Float dz = 0.1*alphaTol*(_hk);
  if (x < _xmin - dx) return 0;
  if (y < _ymin - dy) return 0;
  if (_nk > 1 && z < _zmin - dz) return 0;
  if (x > _xmax + dx) return 0;
  if (y > _ymax + dy) return 0;
  if (_nk > 1 && z > _zmax + dz) return 0;

  E_Float x0 = x-dx; E_Float x1 = x+dx;
  E_Float y0 = y-dy; E_Float y1 = y+dy;
  E_Float z0 = z-dz; E_Float z1 = z+dz;

#ifdef E_ADOLC
  E_Float vx = (x0-_xmin)*_hii;
  E_Float vy = (y0-_ymin)*_hji;
  E_Float vz = (z0-_zmin)*_hki;
  E_Int icmin = E_Int(vx.value())+_ioff;
  E_Int jcmin = E_Int(vy.value())+_joff;
  E_Int kcmin = E_Int(vz.value())+_koff;
#else
  E_Int icmin = E_Int((x0-_xmin)*_hii)+_ioff;
  E_Int jcmin = E_Int((y0-_ymin)*_hji)+_joff;
  E_Int kcmin = E_Int((z0-_zmin)*_hki)+_koff;
#endif

#ifdef E_ADOLC
  vx = (x1-_xmin)*_hii;
  vy = (y1-_ymin)*_hji;
  vz = (z1-_zmin)*_hki;
  E_Int icmax = E_Int(vx.value())+_ioff;
  E_Int jcmax = E_Int(vy.value())+_joff;
  E_Int kcmax = E_Int(vz.value())+_koff;
#else
  E_Int icmax = E_Int((x1-_xmin)*_hii)+_ioff;
  E_Int jcmax = E_Int((y1-_ymin)*_hji)+_joff;
  E_Int kcmax = E_Int((z1-_zmin)*_hki)+_koff;
#endif

  if (icmin < 1 && icmax < 1) return 0;
  if (jcmin < 1 && jcmax < 1) return 0;
  if (kcmin < 1 && kcmax < 1) return 0;
  if (icmin > _ni && icmax > _ni) return 0;
  if (jcmin > _nj && jcmax > _nj) return 0;
  if (kcmin > _nk && kcmax > _nk) return 0;

  icmin = K_FUNC::E_min(_ni-1,icmin); icmin = K_FUNC::E_max(1,icmin); 
  icmax = K_FUNC::E_min(_ni-1,icmax); icmax = K_FUNC::E_max(1,icmax); 
  jcmin = K_FUNC::E_min(_nj-1,jcmin); jcmin = K_FUNC::E_max(1,jcmin);  
  jcmax = K_FUNC::E_min(_nj-1,jcmax); jcmax = K_FUNC::E_max(1,jcmax); 
  kcmin = K_FUNC::E_min(_nk-1,kcmin); kcmin = K_FUNC::E_max(1,kcmin);  
  kcmax = K_FUNC::E_min(_nk-1,kcmax); kcmax = K_FUNC::E_max(1,kcmax);
  for (E_Int k = kcmin; k <= kcmax; k++)
    for (E_Int j = jcmin; j <= jcmax; j++)
      for (E_Int i = icmin; i <= icmax; i++)
      {
        E_Int ind = (i-1) + (j-1)*_ni + (k-1)*_ni*_nj;
        listIndices.push_back(ind);
      }
  // max indice
  //printf("size=%d - tol=%f %f %f - alphaTol=%f, xmin=%f, xmax=%f\n", listIndices.size(),dx,dy,dz,alphaTol,_xmin,_xmax); 
  return listIndices.size();
}
// ============================================================================
/* Find the interpolation cell for point (x,y,z) in case of a structured donor
   zone */
// ============================================================================
short K_INTERP::InterpCart::searchExtrapolationCellCart(
  E_Int ni, E_Int nj, E_Int nk, 
  E_Float* xl, E_Float* yl, E_Float* zl,
  E_Float* cellNp,
  E_Float x, E_Float y, E_Float z,
  E_Int& ic, E_Int& jc, E_Int& kc,
  FldArrayF& cf,
  E_Int nature, E_Int extrapOrder, E_Float constraint)
{
  //E_Int dim = 3;
  //if (nk == 1) dim = 2;
  E_Int i,j,k;
  E_Int icsav = 0;
  E_Int jcsav = 0;
  E_Int kcsav = 0;
  E_Int foundInDomain = 0;
  E_Int ncf = 8;
  FldArrayF cf_sav(ncf);
  E_Float* cfp = cf.begin();
  cf.setAllValuesAtNull();

  // search list of candidate cells
  list<E_Int> listOfCandidateCells;
  E_Float alphaTol = K_FUNC::E_max( (2*constraint-5.)*0.1, 0.);
  E_Int found = getListOfCandidateCells(x,y,z,listOfCandidateCells, alphaTol);
  if (found == 0) return 0;

  list<E_Int>::iterator itr;
  E_Float sum_coef;
  E_Float saved_sum_coef = K_CONST::E_MAX_FLOAT;
  E_Float saved_max_diff_coef = K_CONST::E_MAX_FLOAT;
  E_Float diff_coeff = K_CONST::E_MAX_FLOAT;
  E_Int ind;
  // parameters for extrapolation routine
  E_Int is, js, ks, ret;
  E_Int nij = _ni*_nj;

  // Loop on all candidate cells: computation of coefficients cf
  for (itr = listOfCandidateCells.begin(); itr != listOfCandidateCells.end(); itr++)
  {
    // 1D-index of interpolation cell
    ind = *itr;
 
    // (i,j,k)-indices of interpolation cell
    k = ind/nij; 
    j = (ind-k*nij)/_ni;
    i = ind-j*_ni-k*nij;
    k++; j++; i++;

    // compute interpolation coefficients for hexahedra
    // (is, js, ks): neighbour cell (not used here)
    ret = getExtrapolationCoeffForCell(x, y, z, i, j, k, cf, ni, nj, nk, 
                                       xl, yl, zl, cellNp, is, js, ks, 
                                       nature, constraint, diff_coeff);

    // Keep extrapolation cell with minimum distance (~ minimum sum of absolute coefficients)
    // from the interpolated point
    sum_coef = K_FUNC::E_abs(cfp[0])+K_FUNC::E_abs(cfp[1])+K_FUNC::E_abs(cfp[2])+K_FUNC::E_abs(cfp[3])
    +K_FUNC::E_abs(cfp[4])+K_FUNC::E_abs(cfp[5])+K_FUNC::E_abs(cfp[6])+K_FUNC::E_abs(cfp[7]);

    if (ret == 1 && sum_coef < saved_sum_coef+1.e-6 && diff_coeff< saved_max_diff_coef)
    {
      saved_max_diff_coef = diff_coeff;
      foundInDomain = 1;
      saved_sum_coef = sum_coef;
      icsav = i; jcsav = j; kcsav = k;
      cf_sav = cf;
    }
  }
  ic = icsav; jc = jcsav; kc = kcsav;
  cf = cf_sav;

  if (extrapOrder == 0) // passage a l'ordre 0
  {
    // On reconstruit un xi,eta,zeta comme si les coeff avaient
    // ete obtenus en tri-lineaire
    E_Float xi, eta, zeta;
    // eval, xi, eta, zeta
    xi = cfp[1]+cfp[3]+cfp[5]+cfp[7];
    eta = cfp[2]+cfp[3]+cfp[6]+cfp[7];
    zeta = cfp[4]+cfp[5]+cfp[6]+cfp[7];

    //printf("xi: %f %f %f\n", xi,eta,zeta);
    if (xi < 0) xi = 0.;
    if (xi > 1) xi = 1.;
    if (eta < 0) eta = 0.;
    if (eta > 1) eta = 1.;
    if (zeta < 0) zeta = 0.;
    if (zeta > 1) zeta = 1.;

    cfp[0] = (1-xi)*(1-eta)*(1-zeta);
    cfp[1] = xi*(1-eta)*(1-zeta);
    cfp[2] = (1-xi)*eta*(1-zeta);
    cfp[3] = xi*eta*(1-zeta);
    cfp[4] = (1-xi)*(1-eta)*zeta;
    cfp[5] = xi*(1-eta)*zeta;
    cfp[6] = (1-xi)*eta*zeta;
    cfp[7] = xi*eta*zeta;
  }

  return foundInDomain;
}
//===================================================================
/* Fonctions virtuelles de InterpAdt */
//===================================================================
short K_INTERP::InterpCart::searchInterpolationCellUnstruct(E_Float* xt, E_Float* yt, E_Float* zt,
                                                            FldArrayI& connectEV,
                                                            E_Float x, E_Float y, E_Float z,
                                                            E_Int& noelt, FldArrayF& cf)
{
  return -1;
}

short K_INTERP::InterpCart::searchInterpolationCellStruct(E_Int ni, E_Int nj, E_Int nk, 
                                                          E_Float* xl, E_Float* yl, E_Float* zl,
                                                          E_Float x, E_Float y, E_Float z,
                                                          E_Int& ic, E_Int& jc, E_Int& kc,
                                                          FldArrayF& cf)
{
  return -1;
}

short K_INTERP::InterpCart::searchExtrapolationCellUnstruct(E_Float* xt, E_Float* yt, E_Float* zt,
                                                            E_Float* cellNp,
                                                            FldArrayI& connectEV,
                                                            E_Float x, E_Float y, E_Float z,
                                                            E_Int& noet, FldArrayF& cf,
                                                            E_Int nature, E_Int extrapOrder, E_Float constraint)
{
  return -1;
}

short K_INTERP::InterpCart::searchExtrapolationCellStruct(E_Int ni, E_Int nj, E_Int nk, 
                                                          E_Float* xl, E_Float* yl, E_Float* zl,
                                                          E_Float* cellNp,
                                                          E_Float x, E_Float y, E_Float z,
                                                          E_Int& ic, E_Int& jc, E_Int& kc,
                                                          FldArrayF& cf,
                                                          E_Int nature, E_Int extrapOrder, E_Float constraint)
{
  return -1;
}
