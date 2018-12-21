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
# include <cstdio>
# include <cstdlib>

# include "KMesh.h"
# include "Def/DefCplusPlusConst.h"
# include "Def/DefFunction.h"

using namespace K_FUNC;
using namespace std;
using namespace K_FLD;

extern "C"
{
  void k6rotatemesh_(const E_Int& dim, const E_Float* center,
                     const E_Float* axis, const E_Float& teta,
                     E_Float* x, E_Float* y, E_Float* z);
  
  void k6boundbox_(const E_Int& im, const E_Int& jm, const E_Int& km,
                   const E_Float* x, const E_Float* y, const E_Float* z,
                   E_Float& xmax, E_Float& ymax, E_Float& zmax, 
                   E_Float& xmin, E_Float& ymin, E_Float& zmin );
  
  void k6boundbox2_(const E_Int& im, const E_Int& jm, const E_Int& km,
                    const E_Float* x, const E_Float* y, const E_Float* z,
                    const E_Float* m,const E_Float* r0,const E_Float* xc0, 
                    E_Float& xmax, E_Float& ymax, E_Float& zmax, 
                    E_Float& xmin, E_Float& ymin, E_Float& zmin );

  void k6boundboxunstr_(const E_Int& npts, 
                        const E_Float* x, const E_Float* y, const E_Float* z, 
                        E_Float& xmax, E_Float& ymax, E_Float& zmax, 
                        E_Float& xmin, E_Float& ymin, E_Float& zmin);

  void k6boundboxunstr2_(const E_Int& npts, 
                         const E_Float* x, const E_Float* y, const E_Float* z, 
                         const E_Float* m, const E_Float* r0, const E_Float* xc0,
                         E_Float& xmax, E_Float& ymax, E_Float& zmax, 
                         E_Float& xmin, E_Float& ymin, E_Float& zmin);

  void k6compcartelembox_(const E_Int& is1, const E_Int& is2,
                          const E_Int& js1, const E_Int& js2,
                          const E_Int& ks1, const E_Int& ks2,
                          const E_Int& i2, const E_Int& j2,
                          const E_Int& k2,
                          const E_Int& im, const E_Int& jm,
                          const E_Int& km,
                          const E_Float* x, const E_Float* y, const E_Float* z,
                          const E_Float* rotMat,
                          const E_Float& du, const E_Float& dv,
                          const E_Int& nbElts1, const E_Int& nbElts2,
                          const E_Int& dir, const E_Int& area,
                          const E_Float* nodeMin,
                          const E_Float* nodeMax,
                          const E_Int& szCartElt,
                          E_Float* cartEltMin, E_Float* cartEltMax );

  void k6compstructmetric_(
    const E_Int& ni, const E_Int& nj, const E_Int& nk, 
    const E_Int& nbCells, const E_Int& nbInt,
    const E_Int& nbInti, const E_Int& nbIntj, 
    const E_Int& nbIntk, 
    const E_Float* x, const E_Float* y, const E_Float* z,  
    E_Float* vol,
    E_Float* surfx, E_Float* surfy, E_Float* surfz,
    E_Float* snorm, E_Float* cix, E_Float* ciy, E_Float* ciz); 

  void k6compunstrmetric_(const E_Int& npts, const E_Int& nelts, 
                          const E_Int& nedges, const E_Int& nnodes, 
                          const E_Int* cn, const E_Float* xt, 
                          const E_Float* yt, const E_Float* zt, 
                          const E_Float* xint, const E_Float* yint, const E_Float* zint, 
                          const E_Float* snx, const E_Float* sny, 
                          const E_Float* snz, const E_Float* surf, 
                          const E_Float* vol);

  void k6compstructcellcenter_(E_Int& im, E_Int& jm, E_Int& km, 
                               E_Int& nbNode, E_Int& nbcell, 
                               const E_Float* xt, const E_Float* yt, 
                               const E_Float* zt, E_Float* bary);

  void k6comptetracellcenter_(const E_Int& npts, const E_Int& nelts, 
                              const E_Int* cn, const E_Float* xt,
                              const E_Float* yt, const E_Float* zt,
                              E_Float* bary);
}

//=============================================================================
K_KINTERP::KMesh::KMesh()
  : _im(0), _jm(0), _km(0), _imjm(0), 
    _coord(0,3), _npts(0), _cellVol(0),_bary(0), _cn(0)
{
  //printf("Warning : structured mesh is set by default.\n");
  _isStruct = true;
}

//=============================================================================
/*  Maillage structure */ 
//=============================================================================
K_KINTERP::KMesh::KMesh(E_Int im, E_Int jm, E_Int km)
  : _im(im), _jm(jm), _km(km), _imjm(im*jm), 
    _coord(im*jm*km, 3), _npts(im*jm*km), _cellVol(0),_bary(0), _cn(0)
{
  _isStruct = true;
}
//=============================================================================
/*  Maillage structure */ 
//=============================================================================
K_KINTERP::KMesh::KMesh(E_Int im, E_Int jm, E_Int km,
             const FldArrayF& coord)
  : _im(im), _jm(jm), _km(km), _imjm(im*jm),
    _coord(coord), _npts(im*jm*km), _cellVol(0),_bary(0), _cn(0)
{
  _isStruct = true;
}
//=============================================================================
/*  Maillage non structure */ 
//=============================================================================
K_KINTERP::KMesh::KMesh(FldArrayF& coord, FldArrayI& cn):
  _coord(coord), _cellVol(0),_bary(0), _cn(cn)
{
  _npts = _coord.getSize();
  _isStruct = false;
}
//=============================================================================
K_KINTERP::KMesh::~KMesh()
{
}

//=============================================================================
/* Retourne la connectivite dans le cas non structure */
//=============================================================================
FldArrayI& K_KINTERP::KMesh::getConnectivity()
{
  if ( _isStruct == true ) 
  {
    printf("Warning: getConnectivity: not valid for a structured KMesh.\n");
  }
  return _cn;
}
// ============================================================================
/* Compute the coordinates of the cell barycenters */
// ============================================================================
FldArrayF& K_KINTERP::KMesh::getCellCenter()
{
  if ( _npts == 0 ) 
  {
    printf("Error: getCellCenter: KMesh is not initialised.\n");
    exit(0);
  }
  
  if ( _bary.getSize() != 0 ) return _bary;

  if ( _isStruct == true )
  {
    E_Int imc, jmc, kmc;
 
    if ( _im == 1) 
      imc = 1;
    else 
      imc = _im-1;
    if ( _jm == 1) 
      jmc = 1;
    else jmc = _jm -1;
    if ( _km == 1) 
      kmc = 1;
    else kmc = _km-1;
    
    E_Int nbCell = imc*jmc*kmc;
    E_Int nbNode = _npts;
    
    _bary.malloc(nbCell,3);
    k6compstructcellcenter_( _im, _jm, _km, nbNode, nbCell, 
                             _coord.begin(1), _coord.begin(2), 
                             _coord.begin(3), _bary.begin());
  } 
  else 
  {
    E_Int nelts = _cn.getSize();
    _bary.malloc(nelts, 3);
    
    k6comptetracellcenter_( _npts, nelts, _cn.begin(), _coord.begin(1),
                            _coord.begin(2), _coord.begin(3), 
                            _bary.begin());
  }
  return _bary;
}
//=============================================================================
/* Get the cell volume. If not already computed, compute it. */
//=============================================================================
FldArrayF& K_KINTERP::KMesh::getCellVol()
{
  if ( _cellVol.getSize() != 0) return _cellVol;

  if ( _isStruct == true ) //structure
  {
    E_Int im1 = _im-1;
    E_Int jm1 = _jm-1;
    E_Int km1 = _km-1;
    E_Int nbCells = im1 * jm1 * km1;
    E_Int nbInti = _im*jm1*km1;
    E_Int nbIntj = im1*_jm*km1;
    E_Int nbIntk = im1*jm1*_km;
    E_Int nbInt = nbInti + nbIntj + nbIntk;
 
    _cellVol.malloc(nbCells);
    FldArrayF surf(nbInt,3);
    FldArrayF snorm(nbInt);
    FldArrayF centerInt(nbInt,3);
    k6compstructmetric_(
      _im, _jm, _km, nbCells, nbInt,
      nbInti, nbIntj, nbIntk,
      _coord.begin(1), _coord.begin(2), _coord.begin(3), 
      _cellVol.begin(),
      surf.begin(1), surf.begin(2), surf.begin(3), 
      snorm.begin(),
      centerInt.begin(1), centerInt.begin(1), centerInt.begin(3) );
  }
  else 
  {
    E_Int nelts = _cn.getSize();
    E_Int nedges = 4;  E_Int nnodes = 4;

    _cellVol.malloc(nelts);
    FldArrayF surf(nelts, nedges);
    FldArrayF snx(nelts, nedges);
    FldArrayF sny(nelts, nedges);
    FldArrayF snz(nelts, nedges);
    //tableau local au fortran
    FldArrayF xint(nelts,nedges);
    FldArrayF yint(nelts,nedges);
    FldArrayF zint(nelts,nedges);
    //
    k6compunstrmetric_(_npts, nelts, nedges, nnodes, _cn.begin(), 
                       _coord.begin(1), _coord.begin(2), _coord.begin(3), 
                       xint.begin(), yint.begin(), zint.begin(),
                       snx.begin(), sny.begin(), snz.begin(), 
                       surf.begin(), _cellVol.begin());
  }
  return _cellVol;
}
//=============================================================================
/* Create a new mesh which is made of the centers of the original mesh 
   and which contains points in the middle of boundary interfaces */
//=============================================================================
void K_KINTERP::KMesh::createExtendedCenterMesh(const KMesh& origMesh)
{  
  if ( _isStruct == false )
  {
    printf("KMesh: createExtendedCenterMesh: not valid for an unstructured kmesh. Nothing done.\n");
    return;
  }

  E_Int i, j, k;
  E_Int pos, pos2, pos2a, pos2b, pos2c, pos2d, pos2e, pos2f, pos2g;

  E_Int imo = origMesh.getIm();
  E_Int jmo = origMesh.getJm();
  E_Int kmo = origMesh.getKm();
  
  _im = imo+1;
  _jm = jmo+1;
  _km = kmo+1;
  _imjm = _im*_jm;
  _npts = _imjm*_km;
  E_Int imojmo = imo*jmo;

  _coord.malloc(_npts, 3);
  
  const E_Float* xo = origMesh.getXVector();
  const E_Float* yo = origMesh.getYVector();
  const E_Float* zo = origMesh.getZVector();
  
  E_Float* xt = _coord.begin(1);
  E_Float* yt = _coord.begin(2);
  E_Float* zt = _coord.begin(3);

  E_Int im1 = _im-1; E_Int jm1 = _jm-1; E_Int km1 = _km-1;
  E_Int imo1 = imo-1; E_Int jmo1 = jmo-1; E_Int kmo1 = kmo-1;
  
  /* corners */
  pos  = 0;
  pos2 = 0;
  xt[pos] = xo[pos2];
  yt[pos] = yo[pos2];
  zt[pos] = zo[pos2];

  pos  = im1; 
  pos2 = imo1; 
  xt[pos] = xo[pos2];
  yt[pos] = yo[pos2];
  zt[pos] = zo[pos2];
  
  pos  = jm1*_im; 
  pos2 = jmo1*imo;
  
  xt[pos] = xo[pos2];
  yt[pos] = yo[pos2];
  zt[pos] = zo[pos2];

  pos  = km1*_imjm; 
  pos2 = kmo1*imojmo;
  xt[pos] = xo[pos2];
  yt[pos] = yo[pos2];
  zt[pos] = zo[pos2];

  pos  = im1 + jm1*_im; 
  pos2 = imo1 + jmo1*imo; 
  xt[pos] = xo[pos2];
  yt[pos] = yo[pos2];
  zt[pos] = zo[pos2]; 

  pos  = im1 + km1*_imjm; 
  pos2 = imo1 + kmo1*imojmo;
  xt[pos] = xo[pos2];
  yt[pos] = yo[pos2];
  zt[pos] = zo[pos2];

  pos  = jm1*_im + km1*_imjm; 
  pos2 = jmo1*imo + kmo1*imojmo;
  xt[pos] = xo[pos2];
  yt[pos] = yo[pos2];
  zt[pos] = zo[pos2];
  
  pos  = im1 + jm1*_im + km1*_imjm;
  pos2 = imo1 + jmo1*imo + kmo1*imojmo;
  xt[pos] = xo[pos2];
  yt[pos] = yo[pos2];
  zt[pos] = zo[pos2];
  
  /* border lines */
  for (i = 2; i < _im; i++)
  {
    pos   = getPos(i,1,1);
    pos2  = origMesh.getPos(i,1,1);
    pos2a = origMesh.getPos(i-1,1,1);
    xt[pos] = K_CONST::ONE_HALF*(xo[pos2] + xo[pos2a]);
    yt[pos] = K_CONST::ONE_HALF*(yo[pos2] + yo[pos2a]);
    zt[pos] = K_CONST::ONE_HALF*(zo[pos2] + zo[pos2a]);

    pos   = getPos(i,_jm,1);
    pos2  = origMesh.getPos(i,jmo,1);
    pos2a = origMesh.getPos(i-1,jmo,1);
    xt[pos] = K_CONST::ONE_HALF*(xo[pos2] + xo[pos2a]);
    yt[pos] = K_CONST::ONE_HALF*(yo[pos2] + yo[pos2a]);
    zt[pos] = K_CONST::ONE_HALF*(zo[pos2] + zo[pos2a]);
  
    pos   = getPos(i,1,_km);
    pos2  = origMesh.getPos(i,1,kmo);
    pos2a = origMesh.getPos(i-1,1,kmo);
    xt[pos] = K_CONST::ONE_HALF*(xo[pos2] + xo[pos2a]);
    yt[pos] = K_CONST::ONE_HALF*(yo[pos2] + yo[pos2a]);
    zt[pos] = K_CONST::ONE_HALF*(zo[pos2] + zo[pos2a]);

    pos   = getPos(i,_jm,_km);
    pos2  = origMesh.getPos(i,jmo,kmo);
    pos2a = origMesh.getPos(i-1,jmo,kmo);
    xt[pos] = K_CONST::ONE_HALF*(xo[pos2] + xo[pos2a]);
    yt[pos] = K_CONST::ONE_HALF*(yo[pos2] + yo[pos2a]);
    zt[pos] = K_CONST::ONE_HALF*(zo[pos2] + zo[pos2a]);
  }
  
  for (j = 2; j < _jm; j++)
  {
    pos   = getPos(1,j,1);
    pos2  = origMesh.getPos(1,j,1);
    pos2a = origMesh.getPos(1,j-1,1);
    xt[pos] = K_CONST::ONE_HALF*(xo[pos2] + xo[pos2a]);
    yt[pos] = K_CONST::ONE_HALF*(yo[pos2] + yo[pos2a]);
    zt[pos] = K_CONST::ONE_HALF*(zo[pos2] + zo[pos2a]);

    pos   = getPos(_im,j,1);
    pos2  = origMesh.getPos(imo,j,1);
    pos2a = origMesh.getPos(imo,j-1,1);
    xt[pos] = K_CONST::ONE_HALF*(xo[pos2] + xo[pos2a]);
    yt[pos] = K_CONST::ONE_HALF*(yo[pos2] + yo[pos2a]);
    zt[pos] = K_CONST::ONE_HALF*(zo[pos2] + zo[pos2a]);
    
    pos   = getPos(1,j,_km);
    pos2  = origMesh.getPos(1,j,kmo);
    pos2a = origMesh.getPos(1,j-1,kmo);
    xt[pos] = K_CONST::ONE_HALF*(xo[pos2] + xo[pos2a]);
    yt[pos] = K_CONST::ONE_HALF*(yo[pos2] + yo[pos2a]);
    zt[pos] = K_CONST::ONE_HALF*(zo[pos2] + zo[pos2a]);

    pos   = getPos(_im,j,_km);
    pos2  = origMesh.getPos(imo,j,kmo);
    pos2a = origMesh.getPos(imo,j-1,kmo);
    xt[pos] = K_CONST::ONE_HALF*(xo[pos2] + xo[pos2a]);
    yt[pos] = K_CONST::ONE_HALF*(yo[pos2] + yo[pos2a]);
    zt[pos] = K_CONST::ONE_HALF*(zo[pos2] + zo[pos2a]);
  }

  for (k = 2; k < _km; k++)
  {
    pos   = getPos(1,1,k);
    pos2  = origMesh.getPos(1,1,k);
    pos2a = origMesh.getPos(1,1,k-1);
    xt[pos] = K_CONST::ONE_HALF*(xo[pos2] + xo[pos2a]);
    yt[pos] = K_CONST::ONE_HALF*(yo[pos2] + yo[pos2a]);
    zt[pos] = K_CONST::ONE_HALF*(zo[pos2] + zo[pos2a]);
  
    pos   = getPos(_im,1,k);
    pos2  = origMesh.getPos(imo,1,k);
    pos2a = origMesh.getPos(imo,1,k-1);
    xt[pos] = K_CONST::ONE_HALF*(xo[pos2] + xo[pos2a]);
    yt[pos] = K_CONST::ONE_HALF*(yo[pos2] + yo[pos2a]);
    zt[pos] = K_CONST::ONE_HALF*(zo[pos2] + zo[pos2a]);

    pos   = getPos(1,_jm,k);
    pos2  = origMesh.getPos(1,jmo,k);
    pos2a = origMesh.getPos(1,jmo,k-1);
    xt[pos] = K_CONST::ONE_HALF*(xo[pos2] + xo[pos2a]);
    yt[pos] = K_CONST::ONE_HALF*(yo[pos2] + yo[pos2a]);
    zt[pos] = K_CONST::ONE_HALF*(zo[pos2] + zo[pos2a]);
    
    pos   = getPos(_im,_jm,k);
    pos2  = origMesh.getPos(imo,jmo,k);
    pos2a = origMesh.getPos(imo,jmo,k-1);
    xt[pos] = K_CONST::ONE_HALF*(xo[pos2] + xo[pos2a]);
    yt[pos] = K_CONST::ONE_HALF*(yo[pos2] + yo[pos2a]);
    zt[pos] = K_CONST::ONE_HALF*(zo[pos2] + zo[pos2a]);
  }
  
  /* inside plan */
  E_Int inc = km1 * _imjm;
  E_Int inco = kmo1 * imojmo;

  for (j = 1; j < jm1; j++)
    for (i = 1; i < im1; i++)
    {
      pos   = i + j * _im;
      pos2  = i + j * imo;
      pos2a = pos2-1;
      pos2b = pos2-imo;
      pos2c = pos2b-1;
      xt[pos] =
        K_CONST::ONE_FOURTH*( xo[pos2]  + xo[pos2a] + xo[pos2b] + xo[pos2c] );
      yt[pos] =
        K_CONST::ONE_FOURTH*( yo[pos2]  + yo[pos2a] + yo[pos2b] + yo[pos2c] );
      zt[pos] =
        K_CONST::ONE_FOURTH*( zo[pos2]  + zo[pos2a] + zo[pos2b] + zo[pos2c] );

      pos   = pos + inc;
      pos2  = pos2 + inco;
      pos2a = pos2a + inco;
      pos2b = pos2b + inco;
      pos2c = pos2c + inco;
      xt[pos] =
        K_CONST::ONE_FOURTH*(  xo[pos2]  + xo[pos2a] + xo[pos2b] + xo[pos2c] );
      yt[pos] =
        K_CONST::ONE_FOURTH*(  yo[pos2]  + yo[pos2a] + yo[pos2b] + yo[pos2c] );
      zt[pos] =
        K_CONST::ONE_FOURTH*(  zo[pos2]  + zo[pos2a] + zo[pos2b] + zo[pos2c] );
    }
  
  inc = jm1 * _im;
  inco = jmo1 * imo;
  for (k = 1; k < km1; k++)
    for (i = 1; i < im1; i++)
    {
      pos   = i + k * _imjm; 
      pos2  = i + k * imojmo; 
      pos2a = pos2 - 1;
      pos2b = pos2 - imojmo;
      pos2c = pos2b - 1;
      xt[pos] =
        K_CONST::ONE_FOURTH*(  xo[pos2]  + xo[pos2a] + xo[pos2b] + xo[pos2c] );
      yt[pos] =
        K_CONST::ONE_FOURTH*(  yo[pos2]  + yo[pos2a] + yo[pos2b] + yo[pos2c] );
      zt[pos] =
        K_CONST::ONE_FOURTH*(  zo[pos2]  + zo[pos2a] + zo[pos2b] + zo[pos2c] );

      pos   = pos + inc;
      pos2  = pos2 + inco;
      pos2a = pos2a + inco;
      pos2b = pos2b + inco;
      pos2c = pos2c + inco;
      xt[pos] =
        K_CONST::ONE_FOURTH*(  xo[pos2]  + xo[pos2a] + xo[pos2b] + xo[pos2c] );
      yt[pos] =
        K_CONST::ONE_FOURTH*(  yo[pos2]  + yo[pos2a] + yo[pos2b] + yo[pos2c] );
      zt[pos] =
        K_CONST::ONE_FOURTH*(  zo[pos2]  + zo[pos2a] + zo[pos2b] + zo[pos2c] );
    }

  inc = im1;
  inco = imo1;
  for (k = 1; k < km1; k++)
    for (j = 1; j < jm1; j++)
    {
      pos   = j*_im + k*_imjm;
      pos2  = j*imo + k*imojmo;
      pos2a = pos2 - imojmo;
      pos2b = pos2 - imo;
      pos2c = pos2a - imo;
      xt[pos] =
        K_CONST::ONE_FOURTH*(  xo[pos2]  + xo[pos2a] + xo[pos2b] + xo[pos2c] );
      yt[pos] =
        K_CONST::ONE_FOURTH*(  yo[pos2]  + yo[pos2a] + yo[pos2b] + yo[pos2c] );
      zt[pos] =
        K_CONST::ONE_FOURTH*(  zo[pos2]  + zo[pos2a] + zo[pos2b] + zo[pos2c] );
      pos   = pos + inc;
      pos2  = pos2 + inco;
      pos2a = pos2a + inco;
      pos2b = pos2b + inco;
      pos2c = pos2c + inco;
     
      xt[pos] =
        K_CONST::ONE_FOURTH*(  xo[pos2]  + xo[pos2a] + xo[pos2b] + xo[pos2c] );
      yt[pos] =
        K_CONST::ONE_FOURTH*(  yo[pos2]  + yo[pos2a] + yo[pos2b] + yo[pos2c] );
      zt[pos] =
        K_CONST::ONE_FOURTH*(  zo[pos2]  + zo[pos2a] + zo[pos2b] + zo[pos2c] );
    }

  /* current points */
  for (k = 1; k < km1; k++)
    for (j = 1; j < jm1; j++)
      for (i = 1; i < im1; i++)
      {
        pos   = i + j*_im + k*_imjm;
        pos2  = i + j*imo + k*imojmo;
        pos2a = pos2-1;
        pos2b = pos2 - imo;
        pos2c = pos2 - imojmo;
        pos2d = pos2a - imo;
        pos2e = pos2a - imojmo;
        pos2f = pos2d - imojmo;
        pos2g = pos2b - imojmo;
        xt[pos] =
          K_CONST::ONE_EIGHT*(xo[pos2]  + xo[pos2a] + xo[pos2b] + xo[pos2c] +
                              xo[pos2d]  + xo[pos2e] + xo[pos2f] + xo[pos2g] );
        yt[pos] =
          K_CONST::ONE_EIGHT*(yo[pos2]  + yo[pos2a] + yo[pos2b] + yo[pos2c] +
                              yo[pos2d] + yo[pos2e] + yo[pos2f] + yo[pos2g] );
        zt[pos] =
          K_CONST::ONE_EIGHT*(zo[pos2]  + zo[pos2a] + zo[pos2b] + zo[pos2c] +
                              zo[pos2d] + zo[pos2e] + zo[pos2f] + zo[pos2g] );
      }
}
//=============================================================================
/*  Creation of an extended mesh  from of the original mesh using a rotation
    given the direction dir, duplicate the original mesh into one (or two)
    mesh(es) on the left or/and right side(s) of the origMesh  */
//=============================================================================
void K_KINTERP::KMesh::createDuplicatedExtendedPeriodMesh(const KMesh& origMesh,
                                                         FldArrayF& axisVct,
                                                         FldArrayF& axisPnt,
                                                         E_Float theta)
{
  if ( _isStruct == false )
  {
    printf("KMesh: createDuplicatedExtendedCenterMesh: not valid for an unstructured kmesh. Nothing done.\n");
    return;
  }
  E_Int pos;
  E_Int imo = origMesh.getIm();
  E_Int jmo = origMesh.getJm();
  E_Int kmo = origMesh.getKm();
  E_Int dim1 = imo * jmo * kmo ;
  _im = imo;
  _jm = jmo;
  _km = kmo;
  _imjm = _im*_jm;
  _npts = _imjm*_km;

  FldArrayF coord(dim1, 3);
  _coord.malloc(_npts, 3);  
  E_Float* xt = coord.begin(1);
  E_Float* yt = coord.begin(2);
  E_Float* zt = coord.begin(3);
  const E_Float* xo = origMesh.getXVector();
  const E_Float* yo = origMesh.getYVector();
  const E_Float* zo = origMesh.getZVector();

  E_Int imojmo = imo * jmo;
  /* Points of the original mesh */
  for ( E_Int k = 0; k < kmo; k++)
    for ( E_Int j = 0; j < jmo; j++)
      for ( E_Int i = 0; i < imo; i++)
      {
        pos  = i + j * imo + k*imojmo;
        xt[pos] = xo[pos];
        yt[pos] = yo[pos];
        zt[pos] = zo[pos];
      }
  
  k6rotatemesh_(dim1, axisPnt.begin(), axisVct.begin(), theta, 
                coord.begin(1), coord.begin(2), coord.begin(3));
  
  E_Float* xn = _coord.begin(1);
  E_Float* yn = _coord.begin(2);
  E_Float* zn = _coord.begin(3);
  for (E_Int i = 0; i < dim1; i++)
  {
    xn[i] = xt[i]; 
    yn[i] = yt[i]; 
    zn[i] = zt[i];   
  } 
}

//=============================================================================
/* Find the bounding box of a mesh */
//=============================================================================
void
K_KINTERP::KMesh::boundingBox(E_Float& xmax, E_Float& ymax, E_Float& zmax, 
                             E_Float& xmin, E_Float& ymin, E_Float& zmin) const
{
  if ( _isStruct == true)
    k6boundbox_(_im, _jm, _km, _coord.begin(1), _coord.begin(2), 
                _coord.begin(3), xmax, ymax, zmax, xmin, ymin, zmin);
  else 
    k6boundboxunstr_(_npts, 
                     _coord.begin(1), _coord.begin(2), _coord.begin(3), 
                     xmax, ymax, zmax, xmin, ymin, zmin);
}

//=============================================================================
/* Find the bounding box of a mesh in the absolute frame */
//=============================================================================
void
K_KINTERP::KMesh::boundingBox(E_Float& xmax, E_Float& ymax, E_Float& zmax, 
                             E_Float& xmin, E_Float& ymin, E_Float& zmin,
                             const FldArrayF& m, 
                             const FldArrayF& r0,
                             const FldArrayF& xc0) const
{
  if ( _isStruct == true)
    k6boundbox2_(_im, _jm, _km, 
                 _coord.begin(1), _coord.begin(2), _coord.begin(3), 
                 m.begin(), r0.begin(), xc0.begin(),
                 xmax, ymax, zmax, xmin, ymin, zmin);
  else 
    k6boundboxunstr2_(_npts, 
                      _coord.begin(1), _coord.begin(2), _coord.begin(3), 
                      m.begin(), r0.begin(), xc0.begin(),
                      xmax, ymax, zmax, xmin, ymin, zmin);
}
//=============================================================================
/* Find the cartesian elements bounding box (CEBB) in the Z direction */
//=============================================================================
void K_KINTERP::KMesh::computeZCEBB(E_Int ni, E_Int nj,
                                   E_Float& xmin, E_Float& ymin, E_Float& zmin,
                                   E_Float& xmax, E_Float& ymax, E_Float& zmax,
                                   E_Float& deltax, E_Float& deltay,
                                   FldArrayF& cartZmin, FldArrayF& cartZmax)
{
  if ( _isStruct == false )
  {
    printf("KMesh: computeZCEBB: not valid for an unstructured kmesh.\n");
    exit(0);
  }
  // Compute the standard bounding box of mesh
  boundingBox(xmax, ymax, zmax, xmin, ymin, zmin);

  // Compute discretization step
  deltax = (xmax-xmin)/ni;
  deltay = (ymax-ymin)/nj;

  // Allocate the field and init
  cartZmin.malloc(ni*nj);
  cartZmax.malloc(ni*nj);
  cartZmin.setAllValuesAt(+K_CONST::E_MAX_FLOAT);
  cartZmax.setAllValuesAt(-K_CONST::E_MAX_FLOAT);

  FldArrayF rotMat(3,3);
  rotMat.setAllValuesAtNull();
  for (E_Int i = 0; i < 3; i++) rotMat(i,i+1) = 1.;

  FldArrayF nodeMin(3);
  FldArrayF nodeMax(3);
  nodeMin[0] = xmin;
  nodeMin[1] = ymin;
  nodeMin[2] = zmin;
  nodeMax[0] = xmax;
  nodeMax[1] = ymax;
  nodeMax[2] = zmax;

  // Compute the heigh of each cartesian elements
  k6compcartelembox_( 1, _im, 1, _jm, 1, _km,
                      _im, _jm, _km,
                      _im, _jm, _km,
                      _coord.begin(1), _coord.begin(2), _coord.begin(3),
                      rotMat.begin(),
                      deltax, deltay, ni, nj,
                      3, 0,
                      nodeMin.begin(), nodeMax.begin(),
                      ni*nj,
                      cartZmin.begin(), 
                      cartZmax.begin() );
}

//=============================================================================
/* Find the cartesian elements bounding box (CEBB) in the Y direction */
//=============================================================================
void K_KINTERP::KMesh::computeYCEBB(E_Int ni, E_Int nk,
                                   E_Float& xmin, E_Float& ymin, E_Float& zmin,
                                   E_Float& xmax, E_Float& ymax, E_Float& zmax,
                                   E_Float& deltax, E_Float& deltaz,
                                   FldArrayF& cartYmin, FldArrayF& cartYmax)
{
  if ( _isStruct == false )
  {
    printf("KMesh: computeYCEBB: not valid for an unstructured kmesh.\n");
    exit(0);
  }
  // Compute the standard bounding box of mesh
  boundingBox(xmax, ymax, zmax, xmin, ymin, zmin);

  // Compute discretization step
  deltax = (xmax-xmin)/ni;
  deltaz = (zmax-zmin)/nk;

  // Allocate the field and init
  cartYmin.malloc(ni*nk);
  cartYmax.malloc(ni*nk);
  cartYmin.setAllValuesAt(+K_CONST::E_MAX_FLOAT);
  cartYmax.setAllValuesAt(-K_CONST::E_MAX_FLOAT);

  FldArrayF rotMat(3,3);
  rotMat.setAllValuesAtNull();
  for (E_Int i = 0; i < 3; i++)
    rotMat(i,i+1) = 1.;

  FldArrayF nodeMin(3);
  FldArrayF nodeMax(3);
  nodeMin[0] = xmin;
  nodeMin[1] = ymin;
  nodeMin[2] = zmin;
  nodeMax[0] = xmax;
  nodeMax[1] = ymax;
  nodeMax[2] = zmax;
  
  // Compute the heigh of each cartesian elements
  k6compcartelembox_( 1, _im, 1, _jm, 1, _km,
                      _im, _jm, _km,
                      _im, _jm, _km,
                      _coord.begin(1), _coord.begin(2), _coord.begin(3),
                      rotMat.begin(),
                      deltaz, deltax, nk, ni,
                      2, 0,
                      nodeMin.begin(), nodeMax.begin(),
                      nk*ni,
                      cartYmin.begin(), 
                      cartYmax.begin() );
}

//=============================================================================
/* Find the cartesian elements bounding box (CEBB) in the X direction */
//=============================================================================
void K_KINTERP::KMesh::computeXCEBB(E_Int nj, E_Int nk,
                                   E_Float& xmin, E_Float& ymin, E_Float& zmin,
                                   E_Float& xmax, E_Float& ymax, E_Float& zmax,
                                   E_Float& deltay, E_Float& deltaz,
                                   FldArrayF& cartXmin, FldArrayF& cartXmax)
{
  if ( _isStruct == false )
  {
    printf("KMesh: computeXCEBB: valid for an unstructured kmesh.\n");
    exit(0);
  }
  // Compute the standard bounding box of mesh
  boundingBox(xmax, ymax, zmax, xmin, ymin, zmin);

  // Compute discretization step
  deltay = (ymax-ymin)/nj;
  deltaz = (zmax-zmin)/nk;

  // Allocate the field and init
  cartXmin.malloc(nj*nk);
  cartXmax.malloc(nj*nk);
  cartXmin.setAllValuesAt(+K_CONST::E_MAX_FLOAT);
  cartXmax.setAllValuesAt(-K_CONST::E_MAX_FLOAT);

  FldArrayF rotMat(3,3);
  rotMat.setAllValuesAtNull();
  for (E_Int i = 0; i < 3; i++)
    rotMat(i,i+1) = 1.;

  FldArrayF nodeMin(3);
  FldArrayF nodeMax(3);
  nodeMin[0] = xmin;
  nodeMin[1] = ymin;
  nodeMin[2] = zmin;
  nodeMax[0] = xmax;
  nodeMax[1] = ymax;
  nodeMax[2] = zmax;

  // Compute the heigh of each cartesian elements
  k6compcartelembox_( 1, _im, 1, _jm, 1, _km,
                      _im, _jm, _km,
                      _im, _jm, _km,
                      _coord.begin(1), _coord.begin(2), _coord.begin(3),
                      rotMat.begin(),
                      deltay, deltaz, nj, nk,
                      1, 0,
                      nodeMin.begin(), nodeMax.begin(),
                      nj*nk,
                      cartXmin.begin(), 
                      cartXmax.begin() );
}
//============================================================================
/*

Design
compute the greatest length in each physical direction x, y, z.
returns these lengths.

Synopsis
KMesh mesh;
E_Float dxmax, dymax, dzmax;
mesh.computeCellGreatLength(dxmax,dymax,dzmax);

Description
Compute dimensions of the greatest virtual orthogonal cell of the mesh.
*/
// ----------------------------------------------------------------------------
void K_KINTERP::KMesh::computeCellGreatLength(E_Float& dxmax, E_Float& dymax, 
                                             E_Float& dzmax)
{
  E_Float dxm, dym, dzm;
  dxmax = 0.;
  dymax = 0.;
  dzmax = 0.;

  E_Int   ijk1, ijk2;
  E_Float x1, y1, z1;
  E_Float x2, y2, z2;

  if ( _isStruct == false )
  {
    printf("KMesh: computeCellGreatLength: not valid for an unstructured kmesh.\n");
    exit(0);
  }
  for (E_Int k = 1; k <= _km; k++)
    for (E_Int j = 1; j <= _jm; j++)
      for (E_Int i = 1; i < _im; i++)
      {
        ijk1 = getPos(i,j,k);
        x1 = getX(ijk1);
        y1 = getY(ijk1);
        z1 = getZ(ijk1);

        ijk2 = getPos(i+1,j,k);
        x2 = getX(ijk2);
        y2 = getY(ijk2);
        z2 = getZ(ijk2);

        dxm = E_abs(x2-x1);
        dym = E_abs(y2-y1);
        dzm = E_abs(z2-z1);

        dxmax = E_max(dxmax,dxm);
        dymax = E_max(dymax,dym);
        dzmax = E_max(dzmax,dzm);
      }

  for (E_Int i = 1; i <= _im; i++)
    for (E_Int k = 1; k <= _km; k++)
      for (E_Int j = 1; j < _jm; j++)
      {
        ijk1 = getPos(i,j,k);
        x1 = getX(ijk1);
        y1 = getY(ijk1);
        z1 = getZ(ijk1);

        ijk2 = getPos(i,j+1,k);
        x2 = getX(ijk2);
        y2 = getY(ijk2);
        z2 = getZ(ijk2);

        dxm = E_abs(x2-x1);
        dym = E_abs(y2-y1);
        dzm = E_abs(z2-z1);

        dxmax = E_max(dxmax,dxm);
        dymax = E_max(dymax,dym);
        dzmax = E_max(dzmax,dzm);
      }

  for (E_Int j = 1; j <= _jm; j++)
    for (E_Int i = 1; i <= _im; i++)
      for (E_Int k = 1; k < _km; k++)
      {
        ijk1 = getPos(i, j, k);
        x1 = getX(ijk1);
        y1 = getY(ijk1);
        z1 = getZ(ijk1);

        ijk2 = getPos(i, j, k+1);
        x2 = getX(ijk2);
        y2 = getY(ijk2);
        z2 = getZ(ijk2);

        dxm = E_abs(x2-x1);
        dym = E_abs(y2-y1);
        dzm = E_abs(z2-z1);
 
        dxmax = E_max(dxmax,dxm);
        dymax = E_max(dymax,dym);
        dzmax = E_max(dzmax,dzm);
      }
}
//=============================================================================
/* determine la position d un noeud (i,j,k) demarrant a 1*/
//=============================================================================
E_Int K_KINTERP::KMesh::getPos(E_Int i, E_Int j, E_Int k) const 
{ 

  if ( _isStruct == false )
  {
    printf("KMesh: getPos: not valid for an unstructured KMesh.\n");
    exit(0);
  }

  assert(i >= 1);
  assert(i <= _im);
  assert(j >= 1);
  assert(j <= _jm);
  assert(k >= 1);
  assert(k <= _km);

  return ((i-1) + (j-1)*_im + (k-1)*_imjm);     
}

//=========================Interp/KMesh.cpp================================
