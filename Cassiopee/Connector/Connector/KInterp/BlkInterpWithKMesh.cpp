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
 
# include "BlkInterp.h"
# include "Def/DefCplusPlusConst.h"
# include "String/kstring.h"
# include <stdio.h>
# include <stdlib.h>

using namespace std;
using namespace K_FLD;
using namespace K_FUNC;

extern "C" 
{
  void compinterpolatedptinrefelt_(
    const E_Float* xt, const E_Float* yt, const E_Float* zt, 
    const E_Int& npts, const E_Int& ic, const E_Int& jc, const E_Int& kc, 
    const E_Int& ni, const E_Int& nj, 
    const E_Float& x0, const E_Float& y0, const E_Float& z0, 
    const E_Int& npts_interp_1D, const E_Int& npts_interp_3D,   
    E_Float* x_interp, E_Float* y_interp, E_Float* z_interp, 
    E_Float* base, E_Float* A1, E_Float* B1, 
    E_Int* indxc, E_Int* indxr, E_Int* ipiv,
    E_Float& xx, E_Float& yy, E_Float& zz, E_Int& err);
}

// ============================================================================
/* Constructor */
// ============================================================================
K_KINTERP::BlkInterpWithKMesh::
BlkInterpWithKMesh(K_KINTERP::KMesh& mesh) :
  BlkInterpData(),
  _mesh(mesh)
{
}

// ============================================================================
/* Destructor */
// ============================================================================
K_KINTERP::BlkInterpWithKMesh::~BlkInterpWithKMesh()
{
}

// ============================================================================
/* Search and return interpolation cells with coefficients. 
   Cas  non structure*/
// ============================================================================
short K_KINTERP::BlkInterpWithKMesh::
searchInterpolationCellUnstruct( E_Float x, E_Float y, E_Float z,
                                 E_Int& noelt,
                                 FldArrayI& indi, FldArrayF& cf )   
{
  short found = searchInterpolationCell(x, y, z, noelt, cf);
  if ( found == 0 ) return found;

  indi[0] = 4; // marqueur du type de stockage
  indi[1] = noelt;
  return found;
}

//=============================================================================
/* Find (if possible) the cell containing point(x,y,z). 
   2nd order interpolations. Version : O2CF */
//=============================================================================
short K_KINTERP::BlkInterpWithKMesh::
searchInterpolationCellO2CF( E_Float x, E_Float y, E_Float z,
                             E_Int& ic, E_Int& jc, E_Int& kc,
                             FldArrayF& cf)
{
  return searchInterpolationCell(x, y, z, ic, jc, kc, cf);
}
//=============================================================================
/* O3CF not relevant for extended centers interpolations  */
//=============================================================================
short K_KINTERP::BlkInterpWithKMesh::
searchInterpolationCellO3CF(E_Float x, E_Float y, E_Float z,
                            E_Int& ic, E_Int& jc, E_Int& kc,
                            FldArrayF& cf)
{
  printf("BlkInterpWithKMesh: searchInterpolationCellO3CF not relevant.\n");
  return -1;
}
//=============================================================================
/* Find (if possible) the cell containing point(x,y,z). 
   2nd order interpolations.  */
//=============================================================================
short K_KINTERP::BlkInterpWithKMesh::
searchInterpolationCellO2ABC(E_Float x, E_Float y, E_Float z,
                             E_Int& ic, E_Int& jc, E_Int& kc)
{
  FldArrayF cfloc(8);
  return searchInterpolationCell(x, y, z, ic, jc, kc, cfloc);
}

//=============================================================================
/* O3ABC OK  */
//=============================================================================
short K_KINTERP::BlkInterpWithKMesh::
searchInterpolationCellO3ABC(E_Float x, E_Float y, E_Float z,
                             E_Int& ic, E_Int& jc, E_Int& kc)
{
  FldArrayF cf(8);  
  E_Int ni = _mesh.getIm();
  E_Int nj = _mesh.getJm();
  E_Int nk = _mesh.getKm();
  if ( ni < 3 || nj < 3 || nk < 3)
  { 
    printf("Error: 3rd order interpolation requires at least 3 points per direction.\n");
    return -1;
  }

  short ret = searchInterpolationCell(x, y, z, ic, jc, kc, cf);

  ic = ic-1;
  jc = jc-1;
  kc = kc-1;
  //pour les cas 2d
  if ( ic < 1) ic = 1;
  if ( jc < 1) jc = 1;
  if ( kc < 1) kc = 1;
  return ret;
}
//=============================================================================
/* 5th order interpolation cell search : O5ABC type */
//=============================================================================
short K_KINTERP::BlkInterpWithKMesh::
searchInterpolationCellO5ABC(E_Float x, E_Float y, E_Float z, 
                             E_Int& ic, E_Int& jc, E_Int& kc)
{

  FldArrayF cf(8);
  E_Int ni = _mesh.getIm();
  E_Int nj = _mesh.getJm();
  E_Int nk = _mesh.getKm();
  if ( ni < 5 || nj < 5 || nk < 5)
  { 
    printf("Error: 5th order interpolation requires at least 5 points per direction.\n");
    return -1;
  }
  short ret = searchInterpolationCell(x, y, z, ic, jc, kc, cf);

  //decalage pour frontiere max
  if (ic >= ni-1)
    ic = ni-4;
  else ic = ic-2;

  if (jc >= nj-1)
    jc = nj-4;
  else jc = jc-2;

  if (kc >= nk-1)
    kc = nk-4;
  else kc = kc-2;

  //pour les cas 2d
  if ( ic < 1) ic = 1;
  if ( jc < 1) jc = 1;
  if ( kc < 1) kc = 1;
  return ret;
}

//=============================================================================
/* O2CF type : vectorized */
//=============================================================================
void K_KINTERP::BlkInterpWithKMesh::
searchInterpolationCellO2CFv(FldArrayF& coord,
                             E_Int istart, E_Int iend,
                             FldArrayI& ic, FldArrayI& jc, FldArrayI& kc,
                             FldArrayF& cf, FldArrayIS& found)
{
  searchInterpolationCellv(coord, istart, iend, ic, jc, kc, cf, found);
}
//=============================================================================
/* O3CF type : not relevant for InterpWithExtendedCenters */
//=============================================================================
void K_KINTERP::BlkInterpWithKMesh::
searchInterpolationCellO3CFv(FldArrayF& coord,
                             E_Int istart, E_Int iend,
                             FldArrayI& ic, FldArrayI& jc, FldArrayI& kc,
                             FldArrayF& cf, FldArrayIS& found)
{
  printf("BlkInterpWithKMesh: searchInterpolationCellO3CFv not relevant.\n");
  exit(0);
}
//=============================================================================
/* O2ABC type : vectorized */
//=============================================================================
void K_KINTERP::BlkInterpWithKMesh::
searchInterpolationCellO2ABCv(FldArrayF& coord,
                              E_Int istart, E_Int iend,
                              FldArrayI& ic, FldArrayI& jc, FldArrayI& kc,
                              FldArrayIS& found)
{
  E_Int ic0, jc0, kc0;
  for (E_Int i = istart; i < iend; i++)
  {
    found[i] = searchInterpolationCellO2ABC(
      coord(i,1), coord(i,2), coord(i,3), ic0, jc0, kc0);
    ic[i] = ic0;
    jc[i] = jc0;
    kc[i] = kc0;
  }
}
//=============================================================================
/* O3ABC type : vectorized */
//=============================================================================
void K_KINTERP::BlkInterpWithKMesh::
searchInterpolationCellO3ABCv(FldArrayF& coord,
                              E_Int istart, E_Int iend,
                              FldArrayI& ic, FldArrayI& jc, FldArrayI& kc,
                              FldArrayIS& found)
{
  E_Int ni = _mesh.getIm();
  E_Int nj = _mesh.getJm();
  E_Int nk = _mesh.getKm();

  if ( ni < 3 || nj < 3 || nk < 3)
  { 
    printf("Error: 3rd order interpolation requires at least 3 points per direction.\n");
    exit(0);
  }

  E_Int ic0, jc0, kc0;
  for (E_Int i = istart; i < iend; i++)
  {
    found[i] = searchInterpolationCellO3ABC(
      coord(i,1), coord(i,2), coord(i,3), ic0, jc0, kc0);
    ic[i] = ic0;
    jc[i] = jc0;
    kc[i] = kc0;
  }
}
//=============================================================================
/* 5th order interpolation. cf are not computed here*/
//=============================================================================
void K_KINTERP::BlkInterpWithKMesh::
searchInterpolationCellO5ABCv(FldArrayF& coord,
                              E_Int istart, E_Int iend,
                              FldArrayI& ic, FldArrayI& jc, FldArrayI& kc,
                              FldArrayIS& found)
{
  E_Int ni = _mesh.getIm();
  E_Int nj = _mesh.getJm();
  E_Int nk = _mesh.getKm();
  if ( ni < 5 || nj < 5 || nk < 5)
  { 
    printf("Error: 5th order interpolation requires at least 5 points per direction.\n");
    exit(0);
  }
  E_Int ic0, jc0, kc0;
  for (E_Int i = istart; i < iend; i++)
  {
    found[i] = searchInterpolationCellO5ABC(
      coord(i,1), coord(i,2), coord(i,3), ic0, jc0, kc0);
    
    ic[i] = ic0;
    jc[i] = jc0;
    kc[i] = kc0;
  }
}
// ============================================================================
/*
  Get the hexahedra coordinates+coord of centers of faces+center
  This is the reference for the cell edges numerotation.
*/
// ============================================================================
void K_KINTERP::BlkInterpWithKMesh::
coordHexa(E_Int ind, E_Int ni, E_Int nj,
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
    
  xt[0] = xl[ind];
  yt[0] = yl[ind];
  zt[0] = zl[ind];

  ind = ic+jcmni+kcmnij;
  xt[1] = xl[ind];
  yt[1] = yl[ind];
  zt[1] = zl[ind];

  ind = ic-1+jcni+kcmnij;
  xt[2] = xl[ind];
  yt[2] = yl[ind];
  zt[2] = zl[ind];

  ind = ic+jcni+kcmnij;
  xt[3] = xl[ind];
  yt[3] = yl[ind];
  zt[3] = zl[ind];

  ind = ic-1+jcmni+kcnij;
  xt[4] = xl[ind];
  yt[4] = yl[ind];
  zt[4] = zl[ind];

  ind = ic+jcmni+kcnij;
  xt[5] = xl[ind];
  yt[5] = yl[ind];
  zt[5] = zl[ind];

  ind = ic-1+jcni+kcnij;
  xt[6] = xl[ind];
  yt[6] = yl[ind];
  zt[6] = zl[ind];

  ind = ic+jcni+kcnij;
  xt[7] = xl[ind];
  yt[7] = yl[ind];
  zt[7] = zl[ind];
          
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

// ============================================================================
/*
  Get the hexahedra coordinates
  This is the reference for the cell edges numerotation.
*/
// ============================================================================
void K_KINTERP::BlkInterpWithKMesh::coordHexal(
  E_Int ind, E_Int ni, E_Int nj,
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
    
  xt[0] = xl[ind];
  yt[0] = yl[ind];
  zt[0] = zl[ind];

  ind = ic+jcmni+kcmnij;
  xt[1] = xl[ind];
  yt[1] = yl[ind];
  zt[1] = zl[ind];

  ind = ic-1+jcni+kcmnij;
  xt[2] = xl[ind];
  yt[2] = yl[ind];
  zt[2] = zl[ind];

  ind = ic+jcni+kcmnij;
  xt[3] = xl[ind];
  yt[3] = yl[ind];
  zt[3] = zl[ind];

  ind = ic-1+jcmni+kcnij;
  xt[4] = xl[ind];
  yt[4] = yl[ind];
  zt[4] = zl[ind];

  ind = ic+jcmni+kcnij;
  xt[5] = xl[ind];
  yt[5] = yl[ind];
  zt[5] = zl[ind];

  ind = ic-1+jcni+kcnij;
  xt[6] = xl[ind];
  yt[6] = yl[ind];
  zt[6] = zl[ind];

  ind = ic+jcni+kcnij;
  xt[7] = xl[ind];
  yt[7] = yl[ind];
  zt[7] = zl[ind];
}

//=============================================================================
/* Get interpolation cell by jump from initialized ic, jc, kc */
//=============================================================================
// ic, jc, kc contiennent les anciens indices de la cellule d'interpolation
// et cf les coefficients correspondant.
short K_KINTERP::BlkInterpWithKMesh::
searchInterpolationCellByJump(E_Float x, E_Float y, E_Float z,
                              E_Int& ic, E_Int& jc, E_Int& kc,
                              FldArrayF& cf )
{
//   KMesh& mesh = const_cast<KMesh&>(_mesh);

  K_KINTERP::KMesh& mesh = _mesh;
  /* Get informations on the KMesh on which is built the InterpCartGrid */
  E_Float* xl = mesh.getXVector();
  E_Float* yl = mesh.getYVector();
  E_Float* zl = mesh.getZVector();
  E_Float xt[15], yt[15], zt[15];
  E_Int ni = mesh.getIm();
  E_Int nj = mesh.getJm();
  E_Int nk = mesh.getKm();
  E_Int nij = ni*nj;
  
  E_Int isomm;
  E_Int is, js, ks;  
  E_Float xi, yi, zi;
    
  E_Bool jump = false;
  if ( ic != -1 )
  {
    E_Int ind = (ic-1) + (jc-1)*ni + (kc-1)*nij;
      
    const E_Int maxItJmp = 2;
    for ( E_Int i = 0; i < maxItJmp; i++ )
    {
      coordHexa(ind, ni, nj, xl, yl, zl,
                ic, jc, kc, xt, yt, zt);
      if (getCellJump(x, y, z,
                      xt, yt, zt,
                      isomm,
                      xi, yi, zi) == false)
      {
        if (i == maxItJmp-1)
        {
          jump = false;
          break;
        }
        /* Apply a technique of jump */
        is = E_min(8, E_Int((xi+E_sign(xi))*K_CONST::ONE_HALF));
        is = E_max(-8, is);
        js = E_min(8, E_Int((yi+E_sign(yi))*K_CONST::ONE_HALF));
        js = E_max(-8, js);
        ks = E_min(8, E_Int((zi+E_sign(zi))*K_CONST::ONE_HALF));
        ks = E_max(-8, ks);
          
        ind = ind+is+js*ni+ks*nij;
        kc = ind/nij;
        E_Int kcnij = kc*nij;
        jc = (ind-kcnij)/ni;
        E_Int jcni = jc*ni;
        ic = ind-jcni-kcnij;
      
        if ( ic<0 || ic>ni-2 || jc<0 || jc>nj-2 || kc<0 || kc>nk-2)
        {
          jump = false;
          break;
        }
      }
      else
      {
        jump = true;
        break;
      }
    }
  }
    
  if ( jump == true )
  {
    if ( getCoeffInterpHexa(x, y, z,
                            isomm,
                            xi, yi, zi, 
                            xt, yt, zt,
                            cf) == true )
    {
      return 1;
    }
    else
    {
      jump = false;
    }      
  }
  if ( jump == false )
  {
    return searchInterpolationCell( x, y, z, ic, jc, kc, cf );
  }
  return 0;
}
 
//=============================================================================
/* Get interpolation cell by jump from initialized ic,jc,kc */
//=============================================================================
void K_KINTERP::BlkInterpWithKMesh::
searchInterpolationCellByJumpv(FldArrayF& coord,E_Int istart, E_Int iend,
                               FldArrayI& ic, FldArrayI& jc, FldArrayI& kc,
                               FldArrayF& cf, FldArrayIS& found )
{
  E_Int nocell = 0;
  FldArrayF coordStore( coord.getSize(), coord.getNfld() );
  FldArrayI icStore   ( ic.getSize(), ic.getNfld() );
  FldArrayI jcStore   ( jc.getSize(), jc.getNfld() );
  FldArrayI kcStore   ( kc.getSize(), kc.getNfld() );
  FldArrayF cfStore   ( cf.getSize(), cf.getNfld() );
  FldArrayIS foundStore(found.getSize());
  foundStore.setAllValuesAtNull();
  FldArrayI indirect  (coord.getSize());
  // ic, jc, kc contiennent les anciens indices de la cellule d'interpolation
  // et cf les coefficients correspondant.
  K_KINTERP::KMesh& mesh = _mesh;//const_cast<KMesh&>(_mesh);
  
  /* get informations on the KMesh on which is built the InterpCartGrid */
  E_Float* xl = mesh.getXVector();
  E_Float* yl = mesh.getYVector();
  E_Float* zl = mesh.getZVector();
  E_Float xt[15], yt[15], zt[15];
  E_Int ni = mesh.getIm();
  E_Int nj = mesh.getJm();
  E_Int nk = mesh.getKm();
  E_Int nij = ni*nj;
  
  E_Int storeEnd = 0;
  
  for ( E_Int ii = istart; ii < iend; ii++ )
  {
    E_Int icl = ic[ii];
    E_Int jcl = jc[ii];
    E_Int kcl = kc[ii];
    E_Bool jump = false;
    
    if ( icl != -1 )
    {
      E_Float x = coord(ii,1);
      E_Float y = coord(ii,2);
      E_Float z = coord(ii,3);
      E_Int ind = (icl-1) + (jcl-1)*ni + (kcl-1)*nij;
        
      E_Int isomm;
      E_Int is, js, ks;  
      E_Float xi, yi, zi;
      
      const E_Int maxItJmp = 1;
      for ( E_Int i = 0; i < maxItJmp; i++ )
      {
        coordHexa(ind, ni, nj,
                  xl, yl, zl,
                  icl, jcl, kcl,
                  xt, yt, zt);
        if (getCellJump(x, y, z,
                        xt, yt, zt,
                        isomm,
                        xi, yi, zi) == false)
        {
          if (i == maxItJmp-1)
          {
            jump = false;
            break;
          }
          /* Apply a technique of jump */
          is = E_min(8, E_Int((xi+E_sign(xi))*K_CONST::ONE_HALF));
          is = E_max(-8, is);
          js = E_min(8, E_Int((yi+E_sign(yi))*K_CONST::ONE_HALF));
          js = E_max(-8, js);
          ks = E_min(8, E_Int((zi+E_sign(zi))*K_CONST::ONE_HALF));
          ks = E_max(-8, ks);
      
          ind = ind+is+js*ni+ks*nij;
          kcl = ind/nij;
          E_Int kcnij = kcl*nij;
          jcl = (ind-kcnij)/ni;
          E_Int jcni = jcl*ni;
          icl = ind-jcni-kcnij;
      
          if ( icl<0 || icl>ni-2 || jcl<0 || jcl>nj-2 || kcl<0 || kcl>nk-2)
          {
            jump = false;
            break;
          }
        }
        else
        {
          // the point is always in the cell
          // printf("Point is always in the cell (jump = " SF_D_ ").\n", i);
          jump = true;
          break;
        }     
      }
      
      if ( jump == true )
      {
        FldArrayF cfLoc(8);
        if ( getCoeffInterpHexa(x, y, z,
                                isomm,
                                xi, yi, zi, 
                                xt, yt, zt,
                                cfLoc) == true )
        {
          ic[ii]    = icl;
          jc[ii]    = jcl;
          kc[ii]    = kcl;
          cf(ii,1)  = cfLoc[0];
          cf(ii,2)  = cfLoc[1];
          cf(ii,3)  = cfLoc[2];
          cf(ii,4)  = cfLoc[3];
          cf(ii,5)  = cfLoc[4];
          cf(ii,6)  = cfLoc[5];
          cf(ii,7)  = cfLoc[6];
          cf(ii,8)  = cfLoc[7];
          found[ii] = 1;
        }
        else
        {
          jump = false;
        }
      }
    }
    else
      nocell++;
    
    // nbJump += ( jump ? 1 : 0 );
    if ( false == jump )
    {
      indirect[storeEnd]   = ii;
      coordStore(storeEnd,1) = coord(ii,1);
      coordStore(storeEnd,2) = coord(ii,2);
      coordStore(storeEnd,3) = coord(ii,3);
      storeEnd++;
    }
  }
  //cerr << "Nb Jump valides : " << nbJump << endl;

  printf("Nbre de points sans cellule previous: " SF_D_ " \n", nocell);
  printf("Nbre de points pour recherche classique: " SF_D_ "\n",storeEnd);
  
  if ( storeEnd != 0 )
  {
    searchInterpolationCellv(
      coordStore,0,storeEnd,icStore,jcStore,kcStore,cfStore,foundStore);

    for ( E_Int ii = 0; ii < storeEnd; ii++ )
    {
      ic[indirect[ii]]      = icStore[ii];
      jc[indirect[ii]]      = jcStore[ii];
      kc[indirect[ii]]      = kcStore[ii];
      cf(indirect[ii],1)    = cfStore(ii,1);
      cf(indirect[ii],2)    = cfStore(ii,2);
      cf(indirect[ii],3)    = cfStore(ii,3);
      cf(indirect[ii],4)    = cfStore(ii,4);
      cf(indirect[ii],5)    = cfStore(ii,5);
      cf(indirect[ii],6)    = cfStore(ii,6);
      cf(indirect[ii],7)    = cfStore(ii,7);
      cf(indirect[ii],8)    = cfStore(ii,8);
      found[indirect[ii]]   = foundStore[ii];
    }
  }
}

//=============================================================================
/* Find the best extrapolation cell for point (x,y,z) on structured meshes
   Suppose that the mesh in centers exists 
IN: (x,y,z) : coordonnees du point interpole
IN: (ic,jc,kc) : indices de la cellule d interpolation
OUT: cf : coefficient d extrapolation
IN: cellNatureField : champ de la nature Chimere des cellules sur la grille d interpolation
IN: testNature : si testNature = 0, interpolation aux centres des cellules
                 si testNature = 1, interpolation aux points EX
IN: interpType : type d interpolation (ordre 2 par tetraedre : O2CF, ordre 2,3,5, lagrangienne :  O2ABC, O3ABC, O5ABC)
IN: interpMeshType: type du maillage d interpolation (aux noeuds ou en centres etendus) */
//=============================================================================
short K_KINTERP::BlkInterpWithKMesh::
getExtrapolationCell(E_Float x, E_Float y, E_Float z,
                     E_Int& ic, E_Int& jc, E_Int& kc,
                     FldArrayF& cf,
                     const FldArrayI& cellNatureField,
                     E_Int testNature,
                     E_Float& test,
                     K_KINTERP::BlkInterpData::InterpolationType interpType,
                     K_KINTERP::BlkInterpData::InterpMeshType interpMeshType,
                     E_Float cfMax)
{
  K_KINTERP::KMesh& mesh = _mesh;//const_cast<KMesh&>(_mesh);
  // Find most probable surface
  E_Int nic = mesh.getIm();
  E_Int njc = mesh.getJm();
  E_Int nkc = mesh.getKm();

  E_Float x1, y1, z1;
  E_Int ind1, c;
  E_Int isav = 0;
  E_Int jsav = 0;
  E_Int ksav = 0;
  E_Int is,js,ks;
  E_Float distance[6];
  E_Float distance1,distance2;
  E_Int dir[6];
  E_Int bnd[6];
  short found;
  E_Int i1, i2, j1, j2, k1, k2;
  E_Int dir1, bnd1;

  // Init cf
  cf.setAllValuesAtNull();
  
  // 1st point
  i1 = nic>>1;
  j1 = 1;
  k1 = nkc>>1;

  ind1 = mesh.getPos(i1,j1,k1);
  x1 = mesh.getX(ind1);
  y1 = mesh.getY(ind1);
  z1 = mesh.getZ(ind1);

  distance[0] = (x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1);
  dir[0] = 2;
  bnd[0] = 1;
  
  // 2nd point
  i1 = nic>>1;
  j1 = njc;
  k1 = nkc>>1;

  ind1 = mesh.getPos(i1, j1, k1);
  x1 = mesh.getX(ind1);
  y1 = mesh.getY(ind1);
  z1 = mesh.getZ(ind1);

  c = 1;
  distance[c] = (x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1);
  dir[c] = 2;
  bnd[c] = njc;

  for (E_Int i = c-1; i >= 0; i--)
  {
    if (distance[i] < distance[i+1])
    {
      // Swap i and i+1
      bnd1 = bnd[i+1];
      dir1 = dir[i+1];
      distance1 = distance[i+1];
      
      bnd[i+1] = bnd[i];
      bnd[i] = bnd1;
      dir[i+1] = dir[i];
      dir[i] = dir1;
      distance[i+1] = distance[i];
      distance[i] = distance1;
    }
    else
      break;
  }

  // 3rd point
  i1 = 1;
  j1 = njc>>1;
  k1 = nkc>>1;

  ind1 = mesh.getPos(i1,j1,k1);
  x1 = mesh.getX(ind1);
  y1 = mesh.getY(ind1);
  z1 = mesh.getZ(ind1);

  c = 2;
  distance[c] = (x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1);
  dir[c] = 1;
  bnd[c] = 1;

  for (E_Int i = c-1; i >= 0; i--)
  {
    if (distance[i] < distance[i+1])
    {
      // swap i and i+1
      bnd1 = bnd[i+1];
      dir1 = dir[i+1];
      distance1 = distance[i+1];
      
      bnd[i+1] = bnd[i];
      bnd[i] = bnd1;
      dir[i+1] = dir[i];
      dir[i] = dir1;
      distance[i+1] = distance[i];
      distance[i] = distance1;
    }
    else
      break;
  }

  // 4th point
  i1 = nic;
  j1 = njc>>1;
  k1 = nkc>>1;

  ind1 = mesh.getPos(i1,j1,k1);
  x1 = mesh.getX(ind1);
  y1 = mesh.getY(ind1);
  z1 = mesh.getZ(ind1);

  c = 3;
  distance[c] = (x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1);
  dir[c] = 1;
  bnd[c] = nic;

  for (E_Int i = c-1; i >= 0; i--)
  {
    if (distance[i] < distance[i+1])
    {
      // swap i and i+1
      bnd1 = bnd[i+1];
      dir1 = dir[i+1];
      distance1 = distance[i+1];
      
      bnd[i+1] = bnd[i];
      bnd[i] = bnd1;
      dir[i+1] = dir[i];
      dir[i] = dir1;
      distance[i+1] = distance[i];
      distance[i] = distance1;
    }
    else
      break;
  }

  // 5th point
  i1 = nic>>1;
  j1 = njc>>1;
  k1 = 1;

  ind1 = mesh.getPos(i1,j1,k1);
  x1 = mesh.getX(ind1);
  y1 = mesh.getY(ind1);
  z1 = mesh.getZ(ind1);

  c = 4;
  distance[c] = (x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1);
  dir[c] = 3;
  bnd[c] = 1;

  for (E_Int i = c-1; i >= 0; i--)
  {
    if (distance[i] < distance[i+1])
    {
      // swap i and i+1
      bnd1 = bnd[i+1];
      dir1 = dir[i+1];
      distance1 = distance[i+1];
      
      bnd[i+1] = bnd[i];
      bnd[i] = bnd1;
      dir[i+1] = dir[i];
      dir[i] = dir1;
      distance[i+1] = distance[i];
      distance[i] = distance1;
    }
    else
      break;
  }

  // 6th point
  i1 = nic>>1;
  j1 = njc>>1;
  k1 = nkc;

  ind1 = mesh.getPos(i1,j1,k1);
  x1 = mesh.getX(ind1);
  y1 = mesh.getY(ind1);
  z1 = mesh.getZ(ind1);

  c = 5;
  distance[c] = (x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1);
  dir[c] = 3;
  bnd[c] = nkc;

  for (E_Int i = c-1; i >= 0; i--)
  {
    if (distance[i] < distance[i+1])
    {
      // swap i and i+1
      bnd1 = bnd[i+1];
      dir1 = dir[i+1];
      distance1 = distance[i+1];
      
      bnd[i+1] = bnd[i];
      bnd[i] = bnd1;
      dir[i+1] = dir[i];
      dir[i] = dir1;
      distance[i+1] = distance[i];
      distance[i] = distance1;
    }
    else
      break;
  }
  
  // parcours le plan
  distance2 = K_CONST::E_MAX_FLOAT;
  
  for (c = 5; c >= 0; c--)
  {
    switch (dir[c])
    {
      case 1:
        i1 = bnd[c];
        i2 = bnd[c];
        j1 = 1;
        j2 = njc;
        k1 = 1;
        k2 = nkc;
        break;
      case 2:
        i1 = 1;
        i2 = nic;
        j1 = bnd[c];
        j2 = bnd[c];
        k1 = 1;
        k2 = nkc;
        break;
      case 3:
        i1 = 1;
        i2 = nic;
        j1 = 1;
        j2 = njc;
        k1 = bnd[c];
        k2 = bnd[c];
        break;
      default:
        i1 = 1;
        i2 = 1;
        j1 = 1;
        j2 = 1;
        k1 = 1;
        k2 = 1;
    }
    
    for (E_Int i = i1; i <= i2; i++)
      for (E_Int j = j1; j <= j2; j++)
        for (E_Int k = k1; k <= k2; k++)
        {
          ind1 = mesh.getPos(i,j,k);
          x1 = mesh.getX(ind1);
          y1 = mesh.getY(ind1);
          z1 = mesh.getZ(ind1);
          distance1 = (x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1);

          if (distance1 < distance2)
          {
            distance2 = distance1;
            isav = i;
            jsav = j;
            ksav = k;
          }
        }
    
    // test cell
    if (isav == nic) isav = isav-1;
    if (jsav == njc) jsav = jsav-1;
    if (ksav == nkc) ksav = ksav-1;

    if (cellNatureField.getSize() != 0)
      found = getExtrapolationCoeffForCell(x, y, z,
                                           isav, jsav, ksav, cf,
                                           cellNatureField,
                                           testNature, 0, is, js, ks, cfMax,
                                           interpType, interpMeshType);
    else
      found = getExtrapolationCoeffForCell(x, y, z,
                                           isav, jsav, ksav, cf, cfMax, interpType, interpMeshType);
    
    if (interpType == K_KINTERP::BlkInterpData::O2CF)
    {
      test = E_abs(cf[0])+E_abs(cf[1])+E_abs(cf[2])+E_abs(cf[3])+E_abs(cf[4])+
        E_abs(cf[5])+E_abs(cf[6])+E_abs(cf[7]);
    }
    else if (interpType == K_KINTERP::BlkInterpData::O3CF)
    {
      printf("Warning: extrapolation coeff not implemented for O3CF\n"); return 0;   
    }
    else
    {
      E_Int orderInterp = 2; E_Int nindi = 7;
      if (interpType == K_KINTERP::BlkInterpData::O3ABC) {orderInterp = 3; nindi = 10;}
      else {orderInterp = 5; nindi = 16;}
      test = 0;
      for (E_Int i0 = 1; i0 <= orderInterp; i0++)
        for (E_Int j0 = orderInterp+1; j0 <= 2*orderInterp;j0++)
          for (E_Int k0 = 2*orderInterp+1; k0 < nindi; k0++)
          {
            test += E_abs(cf[i0-1]*cf[j0-1]*cf[k0-1]);
          }
    }

    // if cell is OK (cfi<cfMax) return OK and cellNatureField is OK
    ic = isav; jc = jsav; kc = ksav;
    if (found == 1 && test < cfMax) return 1;

    // try neighbouring cell (is,js,ks)?
    
  }

  return 0;  
}

//=============================================================================
/* Find the best extrapolation cell for point (x,y,z) on structured meshes
   Suppose that the mesh in centers exists 
IN: (x,y,z) : coordonnees du point interpole
IN: indi : tableau permettant de retrouver les indices de la cellule d interpolation
OUT: cf : coefficient d extrapolation
IN: order : ordre de l extrapolation
IN: cellNatureField : champ de la nature Chimere des cellules sur la grille d interpolation
IN: interpType : type d interpolation (ordre 2 par tetraedre : O2CF, ordre 2,3,5, lagrangienne :  O2ABC, O3ABC, O5ABC) */
//=============================================================================
short K_KINTERP::BlkInterpWithKMesh::
getExtrapolationCellStruct(E_Float x, E_Float y, E_Float z,
                           FldArrayI& indi,
                           FldArrayF& cf,
                           E_Int order, E_Float cfMax,
                           const FldArrayF& cellNatureField,
                           K_KINTERP::BlkInterpData::InterpolationType interpType)
{
  E_Int ic, jc, kc;
  short found = searchExtrapolationCell(x, y, z, ic, jc, kc, cf, order, cfMax, cellNatureField);

  // Indices par direction des points de la molecule d'interpolation
  compStandardIndices(ic, jc, kc, indi, interpType); 

  return found;
}
//=============================================================================
/* Find the best extrapolation cell for point (x,y,z) on unstructured meshes
   Suppose that the mesh in centers exists
IN: (x,y,z) : coordonnees du point interpole
IN: indi : tableau permettant de retrouver les indices de la cellule d interpolation
OUT: cf : coefficient d extrapolation
IN: order : ordre de l extrapolation
IN: cellNatureField : champ de la nature Chimere des cellules sur la grille d interpolation
IN: interpType : type d interpolation (ordre 2 par tetraedre : O2CF, ordre 2,3,5, lagrangienne :  O2ABC, O3ABC, O5ABC) */
//=============================================================================
short K_KINTERP::BlkInterpWithKMesh::
getExtrapolationCellUnstr(E_Float x, E_Float y, E_Float z,
                          E_Int& noelt,
                          FldArrayF& cf,
                          FldArrayI& indi,
                          E_Int order,
                          const FldArrayF& cellNatureField)
{
  short found =  searchExtrapolationCell(x, y, z, noelt, cf, indi, order, cellNatureField);
  if ( found == 0 ) return found;

  indi[0] = 4; // marqueur du type de stockage
  indi[1] = noelt;
  return found;  
}

//=============================================================================
/* Compute the extrapolation coefficients (storage OiABC) )for point (x,y,z) 
   inside cell (ic,jc,kc) for structured meshes.
   It is based on the coefficients in the best tetrahedra. 
   In this version (default), the mesh in centers is supposed to be known.
IN: (x,y,z) : coordonnees du point interpole
IN: (ic,jc,kc) : indices de la cellule d'interpolation
OUT: cf : coefficient d extrapolation
IN: cellNatureField : champ de la nature Chimere des cellules sur la grille d interpolation
IN: testNature : si testNature = 0, interpolation aux centres des cellules
                 si testNature = 1, interpolation aux points EX
IN: order : ordre de l extrapolation
IN: interpType : type d interpolation (ordre 2 par tetraedre : O2CF, ordre 2,3,5, lagrangienne :  O2ABC, O3ABC, O5ABC)
IN: interpMeshType : type du maillage d interpolation (aux noeuds ou en centres etendus) */
//=============================================================================
short K_KINTERP::BlkInterpWithKMesh::
getExtrapolationCoeffForCell(E_Float x, E_Float y, E_Float z,
                             E_Int ic, E_Int jc, E_Int kc,
                             FldArrayF& cf,
                             const FldArrayI& cellNatureField,
                             E_Int testNature, E_Int order,
                             E_Int& is, E_Int& js, E_Int& ks, E_Float cfMax,
                             K_KINTERP::BlkInterpData::InterpolationType interpType,
                             K_KINTERP::BlkInterpData::InterpMeshType interpMeshType)
{
  // maillage en centres etendus ou standard selon interpMeshType
  K_KINTERP::KMesh& mesh = _mesh;
  E_Float* xl = mesh.getXVector();
  E_Float* yl = mesh.getYVector();
  E_Float* zl = mesh.getZVector();
  // indice de la cellule d interpolation
  E_Int ind;
  // vecteur pour stocker les coordonnees du centre, des sommets et des centres des faces du hexa d interpolation
  E_Float xt[15];
  E_Float yt[15];
  E_Float zt[15];
  // variables pour le calcul des coefficients tetraedriques
  E_Float xi, yi, zi;
  E_Float xp, yp, zp;
  E_Float xr, yr, zr;
  E_Float xs, ys, zs;
  E_Float xq, yq, zq;
  // dimensions du maillage
  E_Int nic, njc, nkc; // centres etendus
  E_Int im, jm, imjm;  // standard
  // nombre de coefficients, taille du tableau indi defini plus loin
  E_Int ncf, nindi;
  // pointeur sur le cellNatureField
  const E_Int* cellNp = cellNatureField.begin();// 0 : masque 1 autre 
  // ordre des interpolations
  E_Int orderInterp;
  // autres variables locales
  E_Int i;

  switch ( interpType ) 
  {
    case K_KINTERP::BlkInterpData::O2CF:
      ncf = 8; nindi = 7; orderInterp = 2;
      break;
    case K_KINTERP::BlkInterpData::O3ABC:
      ncf = 9; nindi = 10; orderInterp = 3;
      break;
    case K_KINTERP::BlkInterpData::O5ABC:
      ncf = 15; nindi = 16; orderInterp = 5;
      break;
    default:
      printf("Warning: getExtrapolationCoeffForCell: only 2nd, 3rd and 5th order interpolations are implemented.\n");
      ncf = 8; nindi = 7;  orderInterp = 2; interpType = K_KINTERP::BlkInterpData::O2CF;
      break;
  }

  // tableau de stockage des indices des cellules formant la cellule d interpolation
  FldArrayI indi(nindi);

  nic = _mesh.getIm();
  njc = _mesh.getJm();
  nkc = _mesh.getKm();

  // Passage centres etendus -> centres
  if (interpMeshType == K_KINTERP::BlkInterpData::EXT_CENTERS)
  {
    im = nic-2;
    jm = njc-2;
    fromExtendedToStandardCenters(ic, jc, kc, indi, interpType); 
  }
  else
  {
    im = nic-1;
    jm = njc-1;
    compStandardIndices(ic, jc, kc, indi, interpType); 
  }
  imjm = im*jm;

  // Index of interpolation cell (in center)
  vector<E_Int> indTab(ncf);
  E_Int isOrder2 = indi[0];
  if (interpType ==  K_KINTERP::BlkInterpData::O2CF || isOrder2 == 2 ) 
  {
    E_Int indi1 = indi[1];
    E_Int indi2 = indi[2];
    E_Int indi3im = indi[3]*im;
    E_Int indi4im = indi[4]*im;
    E_Int indi5imjm = indi[5]*imjm;
    E_Int indi6imjm = indi[6]*imjm;
    indTab[0] = indi1 + indi3im + indi5imjm;
    indTab[1] = indi2 + indi3im + indi5imjm;
    indTab[2] = indi1 + indi4im + indi5imjm;
    indTab[3] = indi2 + indi4im + indi5imjm;
    indTab[4] = indi1 + indi3im + indi6imjm;
    indTab[5] = indi2 + indi3im + indi6imjm;
    indTab[6] = indi1 + indi4im + indi6imjm;
    indTab[7] = indi2 + indi4im + indi6imjm;
  }
  else if ( interpType ==  K_KINTERP::BlkInterpData::O3CF)
  {
    E_Int c = 0;
    for (E_Int i0 = 1; i0 <= orderInterp; i0++)
      for (E_Int j0 = orderInterp+1; j0 <= 2*orderInterp;j0++)
        for (E_Int k0 = 2*orderInterp+1; k0 < nindi; k0++)
        {
          indTab[c] = indi[i0] + indi[j0]*im + indi[k0]*imjm;
          c++;
        }
  }
  else // OiABC non degrade
  {
    E_Int c = 0;
    for (E_Int i0 = 1; i0 <= orderInterp; i0++)
      for (E_Int j0 = orderInterp+1; j0 <= 2*orderInterp;j0++)
        for (E_Int k0 = 2*orderInterp+1; k0 < nindi; k0++)
        {
          indTab[c] = indi[i0] + indi[j0]*im + indi[k0]*imjm;
          c++;
        }
  }

  // Coord of interpolation cell
  //printf(SF_3D_ ": " SF_D3_ "\n",ic,jc,kc,nic,njc,nkc);
  ind = mesh.getPos(ic, jc, kc);
  E_Int icdummy,jcdummy,kcdummy;
  coordHexa(ind,nic,njc,xl,yl,zl,icdummy,jcdummy,kcdummy,xt,yt,zt);

  if (interpType ==  K_KINTERP::BlkInterpData::O2CF || isOrder2 == 2 )
  { 
    // Try to find a correct tetrahedra
    // that is a tetrahedra where cellNatureField is equal to 1.
    E_Int isom, its, inds, indq, indr, itq, ibsom;
    E_Float cf0;
    E_Int isomm;
  
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

    E_Int ibsomMin = 0;
    E_Int itsMin = 0;
    E_Int itqMin = 0;
    E_Float cfSumMin = K_CONST::E_MAX_FLOAT;
    E_Float cfSum;
    E_Float cfSumMin2 = K_CONST::E_MAX_FLOAT;
    FldArrayF cf_sav(ncf);
    cf_sav.setAllValuesAtNull();
    E_Float sum, sum2;
    E_Float thresold;
    E_Float report;
  
    xp = xt[14];
    yp = yt[14];
    zp = zt[14];

    // teste les 24 tetraedres a la recherche d'un tetraedre
    // valide pour l'interpolation (ne possedant pas de points interpoles)
    // et du tetraedre presentant le min de la somme des val. abs des coeff.
    isomm = 0;
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
        
          /* transforming to hexahedra coefficients */ 
          cf0 = K_CONST::ONE_EIGHTH*(1.-xi-yi-zi);
          for (i = 0; i < 8; i++) cf[i] = cf0;
        
          for (i = 0; i < 4; i++)
          {
            ind = indfc[i+(indq-8)*4];
            cf[ind] = cf[ind]+xi*K_CONST::ONE_FOURTH;
          }
        
          cf[indr] = cf[indr]+yi;
          cf[inds] = cf[inds]+zi;

          cfSum = E_abs(cf[0])+E_abs(cf[1])+E_abs(cf[2])+E_abs(cf[3])+
            E_abs(cf[4])+E_abs(cf[5])+E_abs(cf[6])+E_abs(cf[7]);

          // keep the tetrahedra with smallest sum of abs. coeff.
          if (cfSum < cfSumMin2)
          {
            cfSumMin2 = cfSum;
            cf_sav = cf; // copy the best extrap. coeff (best tetra)
          }

          if ((xi+yi+zi <= K_CONST::ONE)&&(xi >= 0.)&&(yi >= 0.)&&(zi >= 0.))
          {
            // indeed the best (it contains P)
            cfSumMin2 = 0.;
            cf_sav = cf;
          }

          // criterion : all cells are computed (cellN = 1)
          if (testNature == 0 && E_abs(cellNp[indTab[0]] -1)+E_abs(cellNp[indTab[1]] -1)+E_abs(cellNp[indTab[2]] -1)+E_abs(cellNp[indTab[3]] -1)+
              E_abs(cellNp[indTab[4]] -1)+E_abs(cellNp[indTab[5]] -1)+E_abs(cellNp[indTab[6]]-1)+E_abs(cellNp[indTab[7]] -1) < K_CONST::E_CUTOFF)
          {
            // keep the best valid tetra
            if (cfSum < cfSumMin)
            {
              cfSumMin = cfSum;
              ibsomMin = ibsom;
              itsMin = its;
              itqMin = itq;
            }
          }
          // criterion : all cells are computed or interpolated (cellN = 1 or 2)
          if (testNature == 1 && cellNp[indTab[0]]*cellNp[indTab[1]]*cellNp[indTab[2]]*cellNp[indTab[3]]*cellNp[indTab[4]]*cellNp[indTab[5]]*cellNp[indTab[6]]*cellNp[indTab[7]] > 0)
          {
            // keep the best valid tetra
            if (cfSum < cfSumMin)
            {
              cfSumMin = cfSum;
              ibsomMin = ibsom;
              itsMin = its;
              itqMin = itq;
            }
          }        
        }           
      }
    }
    if (order == 1) thresold = 1.e5;
    else if (order == 2) thresold = 100.; // ceci peut-etre modifie
    else thresold = 0.;

    if (cfSumMin < thresold) // OK, We found the best valid tetra (order 1)
    {
      ibsom = ibsomMin;
      its = itsMin;
      itq = itqMin;
    
      isom = indss[ibsom+4*isomm];
      indr = isom;
      xr = xt[indr];
      yr = yt[indr];
      zr = zt[indr];
    
      inds = neighbour[its+isom*3];
      xs = xt[inds];
      ys = yt[inds];
      zs = zt[inds];

      indq = center[itq+its+isom*4];
      xq = xt[indq];
      yq = yt[indq];
      zq = zt[indq];
    
      coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                       xr, yr, zr, xs, ys, zs, xi, yi, zi);
    
      /* transforming to hexahedra coefficients */ 
      cf0 = K_CONST::ONE_EIGHTH*(1.-xi-yi-zi);
      for (i = 0; i < 8; i++) cf[i] = cf0;
    
      for (i = 0; i < 4; i++)
      {
        ind = indfc[i+(indq-8)*4];
        cf[ind] = cf[ind]+xi*K_CONST::ONE_FOURTH;
      }
      cf[indr] = cf[indr]+yi;
      cf[inds] = cf[inds]+zi;

      if (cfSumMin > cfMax) return 0;
      else return 1;
    }
    else
    {
      if (order != 1)
      {
        // Array to define locally the cellNatureField
        FldArrayI locCellN(8);
        
        cf = cf_sav;
        // degenerates to order 0
        // ceci peut introduire des discontinuites de
        // traitement entre deux points consecutifs
        // -----------------
        // centre de cellule : il suffit qu'il existe au moins un point calcule pour que la molecule d extrapolation soit valide
        // point EX : il suffit qu'il existe au moins un point non masque pour que la molecule d extrapolation soit valide

        // testNature = 0: local cellNatureField : 0 (blanked, interpolated), 1 (computed)
        // testNature = 1: local cellNatureField : 0 (blanked),               1 (computed, interpolated)
        if (testNature == 0)
        {
          for (E_Int icell = 0; icell < 8; icell++)
            locCellN[icell] = 1 - E_abs(cellNp[indTab[icell]]-1);
        }
        else if (testNature == 1)
        {
          for (E_Int icell = 0; icell < 8; icell++)
            locCellN[icell] = E_min(cellNp[indTab[icell]],1);
        }

        sum = locCellN[0]+locCellN[1]+locCellN[2]+locCellN[3]+locCellN[4]+locCellN[5]+locCellN[6]+locCellN[7];

        if (sum < K_CONST::E_CUTOFF) goto next; // sum = 0 (critere d echec pour les centres de cellule et les points EX)

        sum2 = 
          (1.-locCellN[0])*cf[0]+(1.-locCellN[1])*cf[1]+
          (1.-locCellN[2])*cf[2]+(1.-locCellN[3])*cf[3]+
          (1.-locCellN[4])*cf[4]+(1.-locCellN[5])*cf[5]+
          (1.-locCellN[6])*cf[6]+(1.-locCellN[7])*cf[7];

        report = sum2 / sum;
      
        cf[0]=(cf[0]+report)*locCellN[0];
        cf[1]=(cf[1]+report)*locCellN[1];
        cf[2]=(cf[2]+report)*locCellN[2];
        cf[3]=(cf[3]+report)*locCellN[3];
        cf[4]=(cf[4]+report)*locCellN[4];
        cf[5]=(cf[5]+report)*locCellN[5];
        cf[6]=(cf[6]+report)*locCellN[6];
        cf[7]=(cf[7]+report)*locCellN[7];

        if (E_abs(cf[0])+E_abs(cf[1])+E_abs(cf[2])+E_abs(cf[3])+
            E_abs(cf[4])+E_abs(cf[5])+E_abs(cf[6])+E_abs(cf[7]) > cfMax)
          return 0;
        else return 1;
      }
    }
  }
  else if (interpType ==  K_KINTERP::BlkInterpData::O3CF)
  {
    printf("Warning : extrapolation coeff not implemented for O3CF\n"); return 0;   
  }
  else // OiABC non degrade
  {
    E_Float thresold = 1.e5;
    E_Int indint;
    // decalage de l indice de la cellule d interpolation
    // pour les frontieres max et pour le 2d
    if (interpType ==  K_KINTERP::BlkInterpData::O3ABC)
    {    
      // frontiere max
      ic = ic-1;
      jc = jc-1;
      kc = kc-1;
      // cas 2d
      if ( ic < 1) ic = 1;
      if ( jc < 1) jc = 1;
      if ( kc < 1) kc = 1;
    }
    else if (interpType ==  K_KINTERP::BlkInterpData::O5ABC)
    {
      // frontiere max
      if (ic >= nic-1)
        ic = nic-4;
      else ic = ic-2;
      
      if (jc >= njc-1)
        jc = njc-4;
      else jc = jc-2;
      
      if (kc >= nkc-1)
        kc = nkc-4;
      else kc = kc-2;
      
      // cas 2d
      if ( ic < 1) ic = 1;
      if ( jc < 1) jc = 1;
      if ( kc < 1) kc = 1;
    }
    short err = compLagrangeCoefs(x, y, z, ic, jc, kc, cf, interpType, interpMeshType);
    E_Float sum_cf=0.;
    E_Int prod_celln = 1;
    E_Int sum_celln = 0;
    E_Int c = 0;
    if (err ==0) // point non interpole
    {
      for (E_Int i0 = 1; i0 <= orderInterp; i0++)
        for (E_Int j0 = orderInterp+1; j0 <= 2*orderInterp;j0++)
          for (E_Int k0 = 2*orderInterp+1; k0 < nindi; k0++)
          {
            sum_cf +=  E_abs(cf[i0-1]*cf[j0-1]*cf[k0-1]);
            indint = indTab[c]; c++;
            sum_celln += cellNp[indint];
            prod_celln = prod_celln*cellNp[indint];
            
          }
      
      // cellule non valide (que des pts masques) ou somme des cfs trop grande
      if ((sum_celln < K_CONST::E_CUTOFF)||(sum_cf > thresold)) goto next;
      
      if (prod_celln > 0) return 1; // pas de pt masque dans la cellule d interpolation
      else goto next; // pas encore implemente : faire une redistribution des coeffs sur les pts non masques
    }
  }
  next:
  // Essai de trouver le voisin le plus probable
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
    
  is = ic;
  js = jc;
  ks = kc;
    
  if (xi >= K_CONST::ONE_HALF) is = E_min(ic+1,nic-1);
  if (xi <= -K_CONST::ONE_HALF) is = E_max(ic-1,1);
  if (yi >= K_CONST::ONE_HALF) js = E_min(jc+1,njc-1);
  if (yi <= -K_CONST::ONE_HALF) js = E_max(jc-1,1);
  if (zi >= K_CONST::ONE_HALF) ks = E_min(kc+1,nkc-1);
  if (zi <= -K_CONST::ONE_HALF) ks = E_max(kc-1,1);
    
  if ((is == ic)&&(js == jc)&&(ks == kc))
  {
    if (xi >= 0.) is = E_min(ic+1,nic-1);
    if (xi < 0.) is = E_max(ic-1,1);
    if (yi >= 0.) js = E_min(jc+1,njc-1);
    if (yi < 0.) js = E_max(jc-1,1);
    if (zi >= 0.) ks = E_min(kc+1,nkc-1);
    if (zi < 0.) ks = E_max(kc-1,1);
  }
    
  return 0; /* Cannot extrapolate from this cell */
}

//=============================================================================
/* Calcul les coefficients d extrapolation pour le pt (x,y,z) a partir des 
   indices (ic,jc,kc) du premier point de la molecule d'interpolation.
   C'est base sur les coefficients du meilleur tetraedre.
   (ic,jc,kc) correspond au KMesh _mesh 
   IN: (x,y,z) : coordonnees du point interpole
   IN: (ic,jc,kc) : indices du premier point de la cellule d interpolation
   OUT: cf : coefficients d extrapolation
   IN: interpType : type d interpolation (ordre 2 par tetraedre : O2CF, ordre 2,3,5, lagrangienne :  O2ABC, O3ABC, O5ABC)
   IN: interpMeshType : type du maillage d interpolation (aux noeuds ou en centres etendus) */
//=============================================================================
short K_KINTERP::BlkInterpWithKMesh::
getExtrapolationCoeffForCell(E_Float x, E_Float y, E_Float z,
                             E_Int ic, E_Int jc, E_Int kc,
                             FldArrayF& cf, E_Float cfMax,
                             K_KINTERP::BlkInterpData::InterpolationType interpType,
                             K_KINTERP::BlkInterpData::InterpMeshType interpMeshType)
{
  K_KINTERP::KMesh& mesh = _mesh;//const_cast<KMesh&>(_mesh);
  E_Float* xl = mesh.getXVector();
  E_Float* yl = mesh.getYVector();
  E_Float* zl = mesh.getZVector();
  E_Int ni = mesh.getIm();
  E_Int nj = mesh.getJm();
  E_Int nk = mesh.getKm();

  E_Int ind;
  E_Int i;
  E_Float xt[15];
  E_Float yt[15];
  E_Float zt[15];
  
  // nombre de coefficients
  E_Int ncf;

  switch ( interpType ) 
  {
    case K_KINTERP::BlkInterpData::O2CF:
      ncf = 8;
      break;
    case K_KINTERP::BlkInterpData::O3ABC:
      ncf = 9;
      break;
    case K_KINTERP::BlkInterpData::O5ABC:
      ncf = 15;
      break;
    default:
      printf("Warning: compAndStoreEXInterpCoefs: only 2nd, 3rd and 5th order interpolations are implemented.\n");
      ncf = 8; interpType = K_KINTERP::BlkInterpData::O2CF;
      break;
  }

  // coord of interpolation cell
  ind = mesh.getPos(ic,jc,kc);
  E_Int icdummy,jcdummy,kcdummy;
  coordHexa(ind,ni,nj,xl,yl,zl,icdummy,jcdummy,kcdummy,xt,yt,zt);
  
  if (interpType ==  K_KINTERP::BlkInterpData::O2CF)
  { 
    // Try to find the best tetrahedra for extrapolation
    E_Float xi,yi,zi;
    E_Float xp,yp,zp;
    E_Float xr,yr,zr;
    E_Float xs,ys,zs;
    E_Float xq,yq,zq;
    E_Int isom,its,inds,indq,indr,itq,ibsom;
    E_Float cf0;
    E_Int isomm;
  
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
  
    E_Float cfSumMin = K_CONST::E_MAX_FLOAT;
    E_Float cfSum;
    E_Float cfSumMin2 = K_CONST::E_MAX_FLOAT;
    FldArrayF cf_sav(ncf);
    cf_sav.setAllValuesAtNull();
  
    xp = xt[14];
    yp = yt[14];
    zp = zt[14];

    // test les 24 tetrahedres a la recherche d'un tetrahedre
    // valide pour l'interpolation (ne possedant pas de points
    // interpoles) et du tetrahedre presentant le min de la somme des val. abs des coeff.
    isomm = 0;
  
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
        
          /* transforming to hexahedra coefficients */ 
          cf0 = K_CONST::ONE_EIGHTH*(1.-xi-yi-zi);
          for (i = 0; i < 8; i++) cf[i] = cf0;
        
          for (i = 0; i < 4; i++)
          {
            ind = indfc[i+(indq-8)*4];
            cf[ind] = cf[ind]+xi*K_CONST::ONE_FOURTH;
          }
        
          cf[indr] = cf[indr]+yi;
          cf[inds] = cf[inds]+zi;

          cfSum = E_abs(cf[0])+E_abs(cf[1])+E_abs(cf[2])+E_abs(cf[3])+
            E_abs(cf[4])+E_abs(cf[5])+E_abs(cf[6])+E_abs(cf[7]);

          // keep the tetrahedra with smallest sum of abs. coeff.
          if (cfSum < cfSumMin2)
          {
            cfSumMin2 = cfSum;
            cf_sav = cf; // copy the best extrap. coeff (best tetra)
          }

          if ((xi+yi+zi <= K_CONST::ONE)&&(xi>=0.)&&(yi>=0.)&&(zi>=0.))
          {
            // indeed the best (it contains P)
            cfSumMin2 = 0.;
            cf_sav = cf;
          }
        
          // keep the best valid tetra
          if (cfSum < cfSumMin) cfSumMin = cfSum;
        }
      }
    }

    cf = cf_sav;

    // CB - check
    cfSum = E_abs(cf[0])+E_abs(cf[1])+E_abs(cf[2])+E_abs(cf[3])+
            E_abs(cf[4])+E_abs(cf[5])+E_abs(cf[6])+E_abs(cf[7]);
    if (cfSum < cfMax) return 1;
    else return 0;
  }
  else if (interpType ==  K_KINTERP::BlkInterpData::O3CF)
  {
    printf("Warning: extrapolation coeff not implemented for O3CF.\n"); return 0;   
  }
  else // OiABC non degrade
  {
    // decalage de l indice de la cellule d interpolation
    // pour les frontieres max et pour le 2d
    if (interpType ==  K_KINTERP::BlkInterpData::O3ABC)
    {    
      // frontiere max
      ic = ic-1;
      jc = jc-1;
      kc = kc-1;
      // cas 2d
      if ( ic < 1) ic = 1;
      if ( jc < 1) jc = 1;
      if ( kc < 1) kc = 1;
    }
    else if (interpType ==  K_KINTERP::BlkInterpData::O5ABC)
    {
      // frontiere max
      if (ic >= ni-1)
        ic = ni-4;
      else ic = ic-2;
      
      if (jc >= nj-1)
        jc = nj-4;
      else jc = jc-2;
      
      if (kc >= nk-1)
        kc = nk-4;
      else kc = kc-2;
      
      // cas 2d
      if ( ic < 1) ic = 1;
      if ( jc < 1) jc = 1;
      if ( kc < 1) kc = 1;
    }
    short err = compLagrangeCoefs(x, y, z, ic, jc, kc, cf, interpType, interpMeshType);
    if (err ==0) return 0; // point non interpole
    return 1;
  }
}

//=============================================================================
/* A partir du pt (ic,jc,kc) reconstruit les tableaux indi, indj, indk des 
   indices en i,j,k des points de la molecule d interpolation
   retourne 1 si ic,jc,kc est dans la grille en centres
   retourne 0 si ic,jc,kc est sur les bords (la ou la solution est extrapolee)
*/
//=============================================================================
E_Int K_KINTERP::BlkInterpWithKMesh::
fromExtendedToStandardCenters(E_Int ic, E_Int jc, E_Int kc,
                              FldArrayI& indTab, 
                              K_KINTERP::BlkInterpData::InterpolationType interpType)
{
  E_Int imr = _mesh.getIm()-2;
  E_Int jmr = _mesh.getJm()-2;
  E_Int kmr = _mesh.getKm()-2;
  E_Int inci1, inci2, incj1, incj2, inck1, inck2;
  E_Int inci3, inci4, incj3, incj4, inck3, inck4;
  E_Int ret = 1;

  switch (interpType)
  {
    case O2CF:
    case O2ABC:
      
      ic = ic-2; jc = jc-2; kc = kc-2;

      indTab[0] = 2;
      indTab[1] = ic;
      indTab[2] = ic+1;
      indTab[3] = jc;
      indTab[4] = jc+1;
      indTab[5] = kc;
      indTab[6] = kc+1;

      //direction i : bords
      if (ic == -1)
      { indTab[1] = 0; ret = 0; }
      
      else if (ic == imr-1)
      { indTab[2] = ic; ret = 0; }
 
      // direction j : bords
      if (jc == -1)
      { indTab[3] = 0; ret = 0; }

      else if (jc == jmr-1)
      { indTab[4] = jc; ret = 0; }

      //direction k : bords
      if (kc == -1)
      { indTab[5] = 0; ret = 0; }
      else if (kc == kmr-1)
      { indTab[6] = kc; ret = 0; }      
      break;

    case O3CF:
    case O3ABC:
      
      ic = ic-2; jc = jc-2; kc = kc-2;
      inci1 = 1; inci2 = 2;
      incj1 = 1; incj2 = 2;
      inck1 = 1; inck2 = 2;

      // Dim  = 2
      if (imr == 1 ) 
      {
        ic = 0;
        inci1 = 0;
        inci2 = 0;
      }
     
      else if (jmr == 1)
      {
        jc = 0; incj1 = 0; incj2 = 0;
      }
      else if (kmr == 1)
      {
        kc = 0; inck1 = 0; inck2 = 0;  
      }

      indTab[0] = 3;
      indTab[1] = ic;
      indTab[2] = ic+inci1;
      indTab[3] = ic+inci2;
      indTab[4] = jc;
      indTab[5] = jc+incj1;
      indTab[6] = jc+incj2;
      indTab[7] = kc;
      indTab[8] = kc+inck1;
      indTab[9] = kc+inck2; 

      //direction i : bords 
      if (ic == -1)
      { indTab[1] = 0; ret = 0; }
      
      else if (ic == imr-2)
      { indTab[3] = ic+1; ret = 0; }
 
      else if (ic == imr-1)
      {
        indTab[2] = ic;
        indTab[3] = ic;
        ret = 0;
      }

      //direction j : bords 
      if (jc == -1)
      { indTab[4] = 0; ret = 0; }
      
      else if (jc == jmr-2)
      { indTab[6] = jc+1; ret = 0; }
 
      else if (jc == jmr-1)
      {
        indTab[5] = jc;
        indTab[6] = jc;
        ret = 0;
      }

      //direction k : bords 
      if (kc == -1)
      { indTab[7] = 0; ret = 0; }
      
      else if (kc == kmr-2)
      { indTab[9] = kc+1; ret = 0; }
 
      else if (kc == kmr-1)
      {
        indTab[8] = kc;
        indTab[9] = kc;
        ret = 0;
      }
      break;

    case O5ABC: 
      
      ic = ic-2; jc = jc-2; kc = kc-2;
      inci1 = 1; incj1 = 1; inck1 = 1;
      inci2 = 2; incj2 = 2; inck2 = 2;
      inci3 = 3; incj3 = 3; inck3 = 3;
      inci4 = 4; incj4 = 4; inck4 = 4;

      // Dimension 2
      if (imr == 1) 
      {
        ic = 0;
        inci1 = 0;
        inci2 = 0;
        inci3 = 0;
        inci4 = 0;
      }
     
      else if (jmr == 1)
      {
        jc = 0;
        incj1 = 0;
        incj2 = 0;
        incj3 = 0;
        incj4 = 0;
      }
      else if (kmr == 1)
      {
        kc = 0;   
        inck1 = 0;
        inck2 = 0;
        inck3 = 0;
        inck4 = 0;
      }

      indTab[0] = 5;
      indTab[1] = ic;
      indTab[2] = ic+inci1;
      indTab[3] = ic+inci2;
      indTab[4] = ic+inci3;
      indTab[5] = ic+inci4;
      indTab[6] = jc;
      indTab[7] = jc+incj1;
      indTab[8] = jc+incj2;
      indTab[9] = jc+incj3;
      indTab[10] = jc+incj4;
      indTab[11] = kc;
      indTab[12] = kc+inck1;
      indTab[13] = kc+inck2; 
      indTab[14] = kc+inck3;
      indTab[15] = kc+inck4;
 
      //direction i : bords 
      if (ic == -1)
      { indTab[1] = 0; ret = 0; }
      
      else if (ic == imr-4)
      { indTab[5] = imr-1; ret = 0; }

      else if (ic == imr-3)
      {     
        indTab[4] = imr-1;
        indTab[5] = imr-1;
        ret = 0;
      }
      else if (ic == imr-2)
      {
        indTab[3] = imr-1;
        indTab[4] = imr-1;
        indTab[5] = imr-1;
        ret = 0;
      }
      else if (ic == imr-1)
      {
        indTab[2] = ic;
        indTab[3] = ic;
        indTab[4] = ic;
        indTab[5] = ic;
        ret = 0;
      }
      
      //direction j : bords 
      if (jc == -1)
      { indTab[6] = 0; ret = 0; }
      
      else if (jc == jmr-4)
      { indTab[10] = jmr-1;ret = 0; }

      else if (jc == jmr-3)
      {     
        indTab[9] = jmr-1;
        indTab[10] = jmr-1;
        ret = 0;
      }
      else if (jc == jmr-2)
      {
        indTab[8] = jmr-1;
        indTab[9] = jmr-1;
        indTab[10] = jmr-1;
        ret = 0;
      }
      else if (jc == jmr-1)
      {
        indTab[7] = jc;
        indTab[8] = jc;
        indTab[9] = jc;
        indTab[10] = jc;
        ret = 0;
      }
      
      // direction k : bords 
      if (kc == -1)
      { indTab[11] = 0; ret = 0; }
      
      else if (kc == kmr-4)
      { indTab[15] = kmr-1; ret = 0; }

      else if (kc == kmr-3)
      {     
        indTab[14] = kmr-1;
        indTab[15] = kmr-1;
        ret = 0;
      }
      else if (kc == kmr-2)
      {
        indTab[13] = kmr-1;
        indTab[14] = kmr-1;
        indTab[15] = kmr-1;
        ret = 0;
      }
      else if (kc == kmr-1)
      {
        indTab[12] = kc;
        indTab[13] = kc;
        indTab[14] = kc;
        indTab[15] = kc;
        ret = 0;
      }
      break;

    default:
      printf("The interpolation type is unknown.\n");
      exit(0);
  }
  return ret;
}
//=============================================================================
// A partir du pt (ic,jc,kc) reconstruit les tableaux indi, indj, indk des 
// indices en i,j,k des points de la molecule d interpolation
// Version tableau : points traites entre istart et iend-1 inclus
//=============================================================================
void K_KINTERP::BlkInterpWithKMesh::
fromExtendedToStandardCentersv(E_Int istart, E_Int iend,
                               FldArrayI& ic, FldArrayI& jc, FldArrayI& kc,
                               FldArrayI& indTab, FldArrayI& extrap,
                               K_KINTERP::BlkInterpData::InterpolationType interpType)
{
  E_Int imr = _mesh.getIm()-2;
  E_Int jmr = _mesh.getJm()-2;
  E_Int kmr = _mesh.getKm()-2;

  E_Int inci1, inci2, incj1, incj2, inck1, inck2;
  E_Int inci3, inci4, incj3, incj4, inck3, inck4;

  E_Int ic0, jc0, kc0;
  E_Int* indTab1 = indTab.begin(1);
  E_Int* indTab2 = indTab.begin(2);
  E_Int* indTab3 = indTab.begin(3);
  E_Int* indTab4 = indTab.begin(4);
  E_Int* indTab5 = indTab.begin(5);
  E_Int* indTab6 = indTab.begin(6);
  E_Int* indTab7 = indTab.begin(7);
  E_Int* indTab8 = NULL;
  E_Int* indTab9 = NULL;
  E_Int* indTab10 = NULL;
  E_Int* indTab11 = NULL;
  E_Int* indTab12 = NULL;
  E_Int* indTab13 = NULL;
  E_Int* indTab14 = NULL;
  E_Int* indTab15 = NULL;
  E_Int* indTab16 = NULL;

  switch (interpType)
  {
    case O2CF:
    case O2ABC:
    
      for (E_Int p = istart; p < iend; p++)
      {
        ic0 = ic[p]-2; jc0 = jc[p]-2; kc0 = kc[p]-2;

        indTab1[p] = 2;
        indTab2[p] = ic0;
        indTab3[p] = ic0+1;
        indTab4[p] = jc0;
        indTab5[p] = jc0+1;
        indTab6[p] = kc0;
        indTab7[p] = kc0+1;

        //direction i : bords
        if (ic0 == -1)
        {indTab2[p] = 0;  extrap[p] = 0;}
        else if (ic0 == imr-1)
        {indTab3[p] = ic0; extrap[p] = 0;}
        
        // direction j : bords
        if (jc0 == -1)
        {indTab4[p] = 0; extrap[p] = 0;}
        else if (jc0 == jmr-1)
        {indTab5[p] = jc0; extrap[p] = 0;}

        //direction k : bords
        if (kc0 == -1)
        {indTab6[p] = 0; extrap[p] = 0;}
        else if (kc0 == kmr-1)
        {indTab7[p] = kc0; extrap[p] = 0;}
      }
      break;
      
    case O3CF:
    case O3ABC:
      indTab8 = indTab.begin(8);
      indTab9 = indTab.begin(9);
      indTab10 = indTab.begin(10);
      for (E_Int p = istart; p < iend; p++)
      {
        ic0 = ic[p]-2; jc0 = jc[p]-2; kc0 = kc[p]-2;
        inci1 = 1; incj1 = 1; inck1 = 1;
        inci2 = 2; incj2 = 2; inck2 = 2;

        // Dim  = 2
        if ( imr == 1 ) 
        {ic0 = 0;  inci1 = 0; inci2 = 0;}
        
        else if ( jmr == 1 )
        {jc0 = 0; incj1 = 0; incj2 = 0;}
        else if ( kmr == 1 )
        {kc0 = 0; inck1 = 0; inck2 = 0;}

        indTab1[p] = 3;
        indTab2[p] = ic0;
        indTab3[p] = ic0+inci1;
        indTab4[p] = ic0+inci2;
        indTab5[p] = jc0;
        indTab6[p] = jc0+incj1;
        indTab7[p] = jc0+incj2;
        indTab8[p] = kc0;
        indTab9[p] = kc0+inck1;
        indTab10[p] = kc0+inck2; 

        //direction i : bords 
        if (ic0 == -1)
        {indTab2[p] = 0; extrap[p] = 0;}
        
        else if (ic0 == imr-2)
        {indTab4[p] = ic0+1;extrap[p] = 0;}
        
        else if ( ic0 == imr-1)
        {indTab3[p] = ic0; indTab4[p] = ic0; extrap[p] = 0;}
        
        //direction j : bords 
        if (jc0 == -1)
        {indTab5[p] = 0; extrap[p] = 0;}
        
        else if (jc0 == jmr-2)
        {indTab7[p] = jc0+1; extrap[p] = 0;}
        
        else if ( jc0 == jmr-1)
        {
          indTab6[p] = jc0;
          indTab7[p] = jc0;
          extrap[p] = 0;
        }
        
        //direction k : bords 
        if (kc0 == -1)
        {indTab8[p] = 0; extrap[p] = 0;}
        
        else if (kc0 == kmr-2)
        {indTab10[p] = kc0+1; extrap[p] = 0;}
        
        else if ( kc0 == kmr-1)
        {
          indTab9[p] = kc0;
          indTab10[p] = kc0;
          extrap[p] = 0;
        }
      }
      break;

    case O5ABC:      
      indTab8 = indTab.begin(8);
      indTab9 = indTab.begin(9);
      indTab10 = indTab.begin(10);
      indTab11 = indTab.begin(11);
      indTab12 = indTab.begin(12);
      indTab13 = indTab.begin(13);
      indTab14 = indTab.begin(14);
      indTab15 = indTab.begin(15);
      indTab16 = indTab.begin(16);

      for (E_Int p = istart; p < iend; p++)
      {
        ic0 = ic[p]-2; jc0 = jc[p]-2; kc0 = kc[p]-2;      
        inci1 = 1; incj1 = 1; inck1 = 1;
        inci2 = 2; incj2 = 2; inck2 = 2;
        inci3 = 3; incj3 = 3; inck3 = 3;
        inci4 = 4; incj4 = 4; inck4 = 4;

        // Dimension 2
        if ( imr == 1 ) 
        {
          ic0 = 0;
          inci1 = 0;
          inci2 = 0;
          inci3 = 0;
          inci4 = 0;
        }
     
        else if ( jmr == 1 )
        {
          jc0 = 0;
          incj1 = 0;
          incj2 = 0;
          incj3 = 0;
          incj4 = 0;
        }
        else if ( kmr == 1 )
        {
          kc0 = 0;   
          inck1 = 0;
          inck2 = 0;  
          inck3 = 0;
          inck4 = 0;
        }

        indTab1[p] = 5;
        indTab2[p] = ic0;
        indTab3[p] = ic0+inci1;
        indTab4[p] = ic0+inci2;
        indTab5[p] = ic0+inci3;
        indTab6[p] = ic0+inci4;
        indTab7[p] = jc0;
        indTab8[p] = jc0+incj1;
        indTab9[p] = jc0+incj2;
        indTab10[p] = jc0+incj3;
        indTab11[p] = jc0+incj4;
        indTab12[p] = kc0;
        indTab13[p] = kc0+inck1;
        indTab14[p] = kc0+inck2;
        indTab15[p] = kc0+inck3;
        indTab16[p] = kc0+inck4;

        //direction i : bords 
        if (ic0 == -1)
        { indTab2[p] = 0; extrap[p] = 0;}
       
        else if ( ic0 == imr-4)
        {indTab6[p] = imr-1; extrap[p] = 0;}

        else if ( ic0 == imr-3)
        {     
          indTab5[p] = imr-1;
          indTab6[p] = imr-1;
          extrap[p] = 0;
        }
        else if ( ic0 == imr-2)
        {
          indTab4[p] = imr-1;
          indTab5[p] = imr-1;
          indTab6[p] = imr-1;
          extrap[p] = 0;
        }
        else if ( ic0 == imr-1)
        {
          indTab3[p] = ic0;
          indTab4[p] = ic0;
          indTab5[p] = ic0;
          indTab6[p] = ic0;
          extrap[p] = 0;
        }
      
        //direction j : bords 
        if (jc0 == -1)
        {indTab7[p] = 0; extrap[p] = 0;}
        
        else if (jc0 == jmr-4)
        {indTab11[p] = jmr-1; extrap[p] = 0;}

        else if (jc0 == jmr-3)
        {     
          indTab10[p] = jmr-1;
          indTab11[p] = jmr-1;
          extrap[p] = 0;
        }
        else if (jc0 == jmr-2)
        {
          indTab9[p] = jmr-1;
          indTab10[p] = jmr-1;
          indTab11[p] = jmr-1;
          extrap[p] = 0;
        }
        else if ( jc0 == jmr-1)
        {
          indTab8[p] = jc0;
          indTab9[p] = jc0;
          indTab10[p] = jc0;
          indTab11[p] = jc0;
          extrap[p] = 0;
        }
        
        //direction k : bords 
        if (kc0 == -1)
        {indTab12[p] = 0; extrap[p] = 0;}
        
        else if ( kc0 == kmr-4)
        {indTab16[p] = kmr-1; extrap[p] = 0;}

        else if ( kc0 == kmr-3)
        {     
          indTab15[p] = kmr-1;
          indTab16[p] = kmr-1;
          extrap[p] = 0;
        }
        else if ( kc0 == kmr-2)
        {
          indTab14[p] = kmr-1;
          indTab15[p] = kmr-1;
          indTab16[p] = kmr-1;
          extrap[p] = 0;
        }
        else if ( kc0 == kmr-1)
        {
          indTab13[p] = kc0;
          indTab14[p] = kc0;
          indTab15[p] = kc0;
          indTab16[p] = kc0;
          extrap[p] = 0;
        }
      }
      break;
    default:
      printf("The interpolation type is unknown.\n");
      exit(0);
  }
}
//===========================================================================
// compute Lagrange polynomials : si x,y,z mal approxime alors compLagrangeCoefs
// retourne 1, sinon si bonne approximation retourne 0
//===========================================================================
short K_KINTERP::BlkInterpWithKMesh::
compLagrangeCoefs(E_Float x, E_Float y, E_Float z,
                  E_Int ic, E_Int jc, E_Int kc,
                  FldArrayF& cf,
                  K_KINTERP::BlkInterpData::InterpolationType interpType,
                  K_KINTERP::BlkInterpData::InterpMeshType interpMeshType)
{
  const E_Float* xl = _mesh.getXVector();
  const E_Float* yl = _mesh.getYVector();
  const E_Float* zl = _mesh.getZVector();
  E_Int ni = _mesh.getIm();
  E_Int nj = _mesh.getJm();
  E_Int npts = ni*nj* _mesh.getKm();

  E_Float ksi, eta, zeta; // coordonnees de (x,y,z) dans element de ref
  // calcul des coordonnees ksi, eta, zeta
  E_Int npts_interp_1D, npts_interp_3D;

  switch ( interpType )
  {
    case O2ABC :
      npts_interp_1D = 2;
      npts_interp_3D = 8;
      break;
    case O3ABC : 
      npts_interp_1D = 3;
      npts_interp_3D = 27;
      break;
    case O5ABC :
      npts_interp_1D = 5;
      npts_interp_3D = 125;
      break;
    default :
      printf("Error: compLagrangeCoefs: not a valid interpolation type.\n"); 
      return 1;
  }

  // declaration des tableaux utilises
  FldArrayF x_interp(npts_interp_3D);
  FldArrayF y_interp(npts_interp_3D);
  FldArrayF z_interp(npts_interp_3D);
  FldArrayF base(npts_interp_1D);
  FldArrayF A1(npts_interp_3D,npts_interp_3D);
  FldArrayF B1(npts_interp_3D,3);
  FldArrayI indxc(npts_interp_3D);
  FldArrayI indxr(npts_interp_3D);
  FldArrayI ipiv(npts_interp_3D);
  E_Int err = 0;//vaut 1 si approximation de x,y,z mauvaise
  
  compinterpolatedptinrefelt_(
    xl, yl, zl, npts, ic, jc, kc, ni, nj, 
    x, y, z, npts_interp_1D, npts_interp_3D, 
    x_interp.begin(), y_interp.begin(), z_interp.begin(), 
    base.begin(), A1.begin(), B1.begin(),
    indxc.begin(), indxr.begin(), ipiv.begin(),  
    ksi, eta, zeta, err);
  if ( err == 1)  return 1; // point interpole mal calcule 

  E_Float* cfp = cf.begin();

  //------------------------------
  // point interpole bien calcule 
  //------------------------------ 
  // Calcul des polynomes correspondants
  E_Float ksip = 1.+ksi;
  E_Float ksim = 1.-ksi;
  E_Float etap = 1.+eta;
  E_Float etam = 1.-eta;
  E_Float zetap = 1.+zeta;
  E_Float zetam = 1.-zeta;
  
  E_Float inv24 = 0.041666666667;
  E_Float inv6 =  0.166666666667;
  E_Float inv4 = K_CONST::ONE_FOURTH;

  E_Float ksip2 = 2.+ksi;
  E_Float ksim2 = 2.-ksi;
  E_Float etap2 = 2.+eta;
  E_Float etam2 = 2.-eta;
  E_Float zetap2 = 2.+zeta;
  E_Float zetam2 = 2.-zeta;

  switch (interpType)
  {
    case O2ABC:
      cfp[0] = ksim;
      cfp[1] = ksi;
      cfp[2] = etam;
      cfp[3] = eta;
      cfp[4] = zetam;
      cfp[5] = zeta;
      break;

    case O3ABC:
      cfp[0] = -K_CONST::ONE_HALF * ksi * ksim; //alpha0
      cfp[1] =  ksip * ksim; //alpha1
      cfp[2] =  K_CONST::ONE_HALF * ksi * ksip; //alpha2
      cfp[3] = -K_CONST::ONE_HALF * eta * etam; //beta0
      cfp[4] =  etap * etam; //beta1
      cfp[5] =  K_CONST::ONE_HALF * eta * etap; //beta2
      cfp[6] = -K_CONST::ONE_HALF * zeta * zetam; //gamma0
      cfp[7] =  zetap * zetam; //gamma1
      cfp[8] =  K_CONST::ONE_HALF * zeta * zetap; //gamma2
      break;

    case O5ABC:
      cfp[0] = inv24 * ksi*ksip*ksim*ksim2; //alpha1
      cfp[1] = -inv6 * ksi*ksip2*ksim*ksim2; //alpha2
      cfp[2] =  inv4 * ksim2*ksim*ksip*ksip2;//alpha3
      cfp[3] = inv6 * ksi*ksim2*ksip*ksip2; //alpha4
      cfp[4] = -inv24 * ksip2*ksip*ksi*ksim; //alpha5

      cfp[5] = inv24 * eta*etap*etam*etam2; //beta1
      cfp[6] = -inv6 * eta*etap2*etam*etam2 ; //beta2
      cfp[7] =  inv4 * etam2*etam*etap*etap2; //beta3
      cfp[8] = inv6 * eta*etam2*etap*etap2; //beta4
      cfp[9] = -inv24 * etap2*etap*eta*etam; //beta5

      cfp[10] = inv24 * zeta*zetap*zetam*zetam2; //gamma1
      cfp[11] = -inv6 * zeta*zetap2*zetam*zetam2; //gamma2
      cfp[12] =  inv4 * zetam2*zetam*zetap*zetap2; //gamma3
      cfp[13] = inv6 * zeta*zetam2*zetap*zetap2; //gamma4
      cfp[14] = -inv24 * zetap2*zetap*zeta*zetam; //gamma5
      break;
      
    default:
      printf("Error: compLagrangeCoefs: not a valid interpolation type.\n"); 
      return 1;
  }
  return 0;
}

//=============================================================================
/* Calcul des indices i,i+1... de la molecule d' interpolation 
   Le premier element de indi : indi[0] donne l info sur l ordre 
   d interpolation */
//=============================================================================
void K_KINTERP::BlkInterpWithKMesh::
compStandardIndices(E_Int ic, E_Int jc, E_Int kc, 
                    FldArrayI& indi, 
                    K_KINTERP::BlkInterpData::InterpolationType interpType)
{
  switch (interpType)
  {
    case O2CF :
    case O2ABC :
      indi[0] = 2;
      indi[1] = ic-1;
      indi[2] = ic;
      indi[3] = jc-1;
      indi[4] = jc;   
      indi[5] = kc-1;
      indi[6] = kc;
      break;
    case O3ABC:
    case O3CF:
      indi[0] = 3;
      indi[1] = ic-1;
      indi[2] = ic;
      indi[3] = ic+1;
      indi[4] = jc-1;
      indi[5] = jc;
      indi[6] = jc+1;
      indi[7] = kc-1;
      indi[8] = kc;
      indi[9] = kc+1;
      break;
    case O5ABC:
      indi[0] = 5;
      indi[1] = ic-1;
      indi[2] = ic;
      indi[3] = ic+1;
      indi[4] = ic+2;
      indi[5] = ic+3;
      indi[6] = jc-1;
      indi[7] = jc;
      indi[8] = jc+1;
      indi[9] = jc+2;
      indi[10] = jc+3;
      indi[11] = kc-1;
      indi[12] = kc;
      indi[13] = kc+1;
      indi[14] = kc+2;
      indi[15] = kc+3;
      break;
    default:
      printf("compStandardIndices: unknown interpolation type.\n");
      exit(0);
  }
}
//=============================================================================
/* Calcul des indices i,i+1... des molecules d' interpolation 
   comprises entre istart et iend-1 inclus
   Les premiers elemts de indiv :indiv(i,1) sont deja calcules avant
   et valent l ordre d interpolation local au pt 
*/
//=============================================================================
void K_KINTERP::BlkInterpWithKMesh::
compStandardIndicesv(E_Int istart, E_Int iend,
                     FldArrayI& ic, FldArrayI& jc, FldArrayI& kc, 
                     FldArrayI& indi, 
                     K_KINTERP::BlkInterpData::InterpolationType interpType)
{
  switch (interpType)
  {
    case O2CF:
    case O2ABC:
      for (E_Int i = istart; i < iend; i++)
      {
        indi(i,1) = 2;
        indi(i,2) = ic[i]-1;
        indi(i,3) = ic[i];
        indi(i,4) = jc[i]-1;
        indi(i,5) = jc[i];   
        indi(i,6) = kc[i]-1;
        indi(i,7) = kc[i];
      }
      break;
    case O3CF:
    case O3ABC:
      for (E_Int i = istart; i < iend; i++)
      {
        indi(i,1) = 3;
        indi(i,2) = ic[i]-1;
        indi(i,3) = ic[i];
        indi(i,4) = ic[i]+1;
        indi(i,5) = jc[i]-1;
        indi(i,6) = jc[i];
        indi(i,7) = jc[i]+1;   
        indi(i,8) = kc[i]-1;
        indi(i,9) = kc[i];
        indi(i,10) = kc[i]+1;
      }    
      break;
    case O5ABC:
      for (E_Int i = istart; i < iend; i++)
      {
        indi(i,1) = 5;
        indi(i,2) = ic[i]-1;
        indi(i,3) = ic[i];
        indi(i,4) = ic[i]+1;
        indi(i,5) = ic[i]+2;
        indi(i,6) = ic[i]+3;
        indi(i,7) = jc[i]-1;
        indi(i,8) = jc[i];
        indi(i,9) = jc[i]+1;
        indi(i,10) = jc[i]+2;
        indi(i,11) = jc[i]+3;
        indi(i,12) = kc[i]-1;
        indi(i,13) = kc[i];
        indi(i,14) = kc[i]+1;
        indi(i,15) = kc[i]+2;
        indi(i,16) = kc[i]+3;
      }
      break;
    default:
      printf("compStandardIndices: unknown interpolation type.\n");
      exit(0);
  }
}

// ============= Interp/BlkInterpWithKMesh.cpp ===============
