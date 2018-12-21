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
# include "Interp/Interp.h"
# include <list>
using namespace std;
using namespace K_FLD;
E_Float _EPS_DET   = 1.e-16;
E_Float _EPS_TETRA = 1.e-4;
E_Float _EPS_GEOM  = K_CONST::E_GEOM_CUTOFF;
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

//=============================================================================
/* a est la connectivite cEEN dans le cas non structure, NULL si pas besoin
   IN: x,y,z: point a extrapoler
   IN: interpData: l'interpData du maillage donneur
   IN: meshType: type du maillage donneur (2: Non-Structure, 1: structure)
   IN: ni, nj, nk: taille du maillage donneur
   IN: xl, yl, zl: maillage donneur 
   IN: cellN: cellN du maillage donneur
   IN: interpType: type d'extrapolation
   IN: nature=O aucun cellN=0 ds la molecule donneuse 
       nature=1 que des cellN=1 ds la molecule donneuse
   IN: constraint: contrainte sur la valeur absolu des coeffs
   IN: extrapOrder: 0 (ordre 0), 1 (lineaire)
*/
//=============================================================================
E_Int K_INTERP::getExtrapolationData(
  E_Float x, E_Float y, E_Float z,
  InterpAdt* interpData, void* c1,  void* c2, 
  E_Int meshtype, E_Int ni, E_Int nj, E_Int nk,
  E_Float* xl, E_Float* yl, E_Float* zl, E_Float* cellN,
  E_Int& isBorder, E_Int& type, FldArrayI& indi, FldArrayF& cf, 
  K_INTERP::InterpAdt::InterpolationType interpType,
  E_Int nature, E_Float constraint, E_Int extrapOrder)
{ 
  E_Int dim = 3;
  if (nk == 1) dim = 2;
  short found = 0;
  E_Int ic, jc, kc, noet;
  isBorder = 0;
  // O2CF et nature = 0 (pas de pts masques)
  switch (interpType)
  {
    case K_INTERP::InterpAdt::O2CF:
      if (meshtype == 2) 
      {
        FldArrayI& cEV = *(FldArrayI*)c1;
        found = searchExtrapolationCellUnstruct(
          xl, yl, zl, cellN, cEV, interpData, 
          x, y, z, noet, cf, nature, extrapOrder, constraint);
        if (found < 1) return found;
        //if (c2 != NULL) 
        //{
        //  vector< vector<E_Int> > cEEN;
        //  cEEN = *(vector< vector<E_Int> >*)c2;
        //  if (cEEN[noet].size() < 4) isBorder = 1; 
        //}
        type = 4; indi[0] = noet;
      }
      else 
      {
        found = searchExtrapolationCellStruct(
          ni, nj, nk, xl, yl, zl, cellN, interpData, 
          x, y, z, ic, jc, kc, cf, nature, extrapOrder, constraint);

        if (found < 1) return found;
        ic = ic-1; jc = jc-1; kc = kc-1; // pour demarrer a 0        
        if ( dim == 3 ) 
        {
          if (ic == 0 || ic == ni-2 || jc == 0 || jc == nj-2 || kc == 0 || kc == nk-2) isBorder = 1;
          type = 2; indi[0] = ic + jc*ni + kc*ni*nj;
        }
        else 
        {
          if (ic == 0 || ic == ni-2 || jc == 0 || jc == nj-2 ) isBorder = 1;
          type = 22; indi[0] = ic + jc*ni;
          //report des coefs
          for (E_Int nocf = 0; nocf < 4; nocf++)
          {
            cf[nocf] += cf[nocf+4];
            cf[nocf+4] = 0.;
          }
        }
      }
      break;
      
    case K_INTERP::InterpAdt::O3ABC:
    case K_INTERP::InterpAdt::O5ABC:
      //if (meshtype == 2) {printf("Error: getExtrapolationData: interpType ABC not valid for unstructured zones.\n"); return -1;}
      found = searchExtrapolationCellStruct(
        ni, nj, nk, xl, yl, zl, cellN, interpData, 
        x, y, z, ic, jc, kc, cf, nature, extrapOrder, constraint);
      if (found < 1) return found;
      ic = ic-1; jc = jc-1; kc = kc-1;// pour demarrer à 0        
      if (ic == 0 || ic == ni-2 || jc == 0 || jc == nj-2 || kc == 0 || kc == nk-2) isBorder = 1;
      type = 2; indi[0] = ic + jc*ni + kc*ni*nj;
      break;

    default:
      return -1;
  }
  return found;
}
//=============================================================================
/* a est la connectivite cEEN dans le cas non structure, NULL si pas besoin  
   IN: xl, yl, zl: maillage sur lequel les coefs d interpolation sont calcules 
   IN: ni, nj, nk: taille du maillage ci dessus
   IN: nature=O aucun cellN=0 ds la molecule donneuse 
       nature=1 que des cellN=1 ds la molecule donneuse */
//=============================================================================
E_Int K_INTERP::getInterpolationData(
  E_Float x, E_Float y, E_Float z,
  InterpAdt* interpData, void* c1, void* c2, 
  E_Int meshtype, E_Int ni, E_Int nj, E_Int nk,
  E_Float* xl, E_Float* yl, E_Float* zl, E_Float* cellN,
  E_Int& isBorder, E_Int& type, FldArrayI& indi, FldArrayF& cf, 
  K_INTERP::InterpAdt::InterpolationType interpType,
  E_Int nature)
{
  E_Int dim = 3;
  if (nk == 1) dim = 2;
  E_Float geomCutOff = K_CONST::E_GEOM_CUTOFF;
  short corr = 0; // correction d'ordre 2 -> corr=1 
  short found = 0;
  E_Int ics, jcs, kcs, ic, jc, kc, noet, ind;
  isBorder = 0;
  E_Float cellN0;
  E_Float val; E_Int d;
  //2D case : z must be equal to _zmin 
  if (dim == 2) 
  {
    if (K_FUNC::E_abs(z-interpData->_zmin)>geomCutOff) return 0;
  }
  switch (interpType)
  {
    case K_INTERP::InterpAdt::O2CF:
      if (meshtype == 2)
      {
        FldArrayI& cEV = *(FldArrayI*)c1;
        found = searchInterpolationCellUnstruct(
          xl, yl, zl, cEV, interpData, x, y, z, noet, cf);
        if (found < 1) return found;
        //if (c2 != NULL)
        //{
        //  vector< vector<E_Int> > cEEN;
        //  cEEN = *(vector< vector<E_Int> >*)c2;
        //  if (cEEN[noet].size() < 4) isBorder = 1; 
        //}
        type = 4; indi[0] = noet;
        if (cellN != NULL) 
        {
          if (nature == 0) // pas de pt masque ds la molecule donneuse
          {
            val = 1.;
            for (E_Int nov = 1; nov <= 4; nov++)
            {
              ind = cEV(noet, nov); 
              cellN0 = cellN[ind];
              d = cf[nov-1] > 1.e-12;
              val *= cellN0-d+1;
            }
            if (K_FUNC::fEqualZero(val, geomCutOff) == true) return 0;// pas interpolable 
          }
          else // pas de pt masque ou interpole dans la cellule donneuse 
          {
            val = 0.;
            for (E_Int nov = 1; nov <= 4; nov++)
            {
              ind = cEV(noet, nov); 
              cellN0 = cellN[ind];
              val += K_FUNC::E_abs(cf[nov-1])*K_FUNC::E_abs(1.-cellN0);
            }
            if (K_FUNC::fEqualZero(val, geomCutOff) == false) return 0;// pas interpolable 
          }
        }
      }
      else
      { 
        found = searchInterpolationCellStruct(
          ni, nj, nk, xl, yl, zl, interpData, x, y, z, ic, jc, kc, cf); 
        if (found < 1) return found;

        ic = ic-1; jc = jc-1; kc = kc-1; // pour demarrer a 0   
        if ( dim == 3 ) 
        {
          if (ic == 0 || ic == ni-2 || jc == 0 || jc == nj-2 || kc == 0 || kc == nk-2) isBorder = 1;

          type = 2; indi[0] = ic + jc*ni + kc*ni*nj;
          if (cellN != NULL)
          {
            if (nature == 0) // pas de pt masque ds la molecule donneuse
            {
              val = 1.;
              for (E_Int kk = 0; kk < 2; kk++)
                for (E_Int jj = 0; jj < 2; jj++)
                  for (E_Int ii = 0; ii < 2; ii++)
                  {
                    ind = (ic+ii) + (jc+jj)*ni + (kc+kk)*ni*nj;
                    cellN0 = cellN[ind];
                    d = cf[ii+jj*2+kk*4] > 1.e-12;
                    val *= cellN0-d+1;
                  }
              if (K_FUNC::fEqualZero(val, geomCutOff) == true) return 0; // pas interpolable 
            }
            else // pas de pt masque ou interpole dans la cellule donneuse
            {
              val = 0.;
              for (E_Int kk = 0; kk < 2; kk++)
                for (E_Int jj = 0; jj < 2; jj++)
                  for (E_Int ii = 0; ii < 2; ii++)
                  {
                    ind = (ic+ii) + (jc+jj)*ni + (kc+kk)*ni*nj;
                    cellN0 = cellN[ind];
                    val += K_FUNC::E_abs(cf[ii+jj*2+kk*4])*K_FUNC::E_abs(1.-cellN0);
                  }
              if (K_FUNC::fEqualZero(val, geomCutOff) == false) return 0;// pas interpolable 
            }
          }
        }//fin 3D
        else 
        {
          if (ic == 0 || ic == ni-2 || jc == 0 || jc == nj-2) isBorder = 1;
          type = 22; indi[0] = ic + jc*ni; 
          //report des coefs
          E_Float sumCf = 0.;
          for (E_Int nocf = 0; nocf < 4; nocf++)
          {
            cf[nocf] += cf[nocf+4];
            cf[nocf+4] = 0.;
            sumCf += K_FUNC::E_abs(cf[nocf]);
          }
          if (cellN != NULL)
          {
            if (nature == 0) // pas de pt masque ds la molecule donneuse
            {
              val = 1.;
              for (E_Int jj = 0; jj < 2; jj++)
                for (E_Int ii = 0; ii < 2; ii++)
                {
                  ind = (ic+ii) + (jc+jj)*ni;
                  cellN0 = cellN[ind];
                  d = cf[ii+jj*2] > 1.e-12;
                  val *= cellN0-d+1;
                }
              if (K_FUNC::fEqualZero(val, geomCutOff) == true) return 0;// pas interpolable 
            }
            else // pas de pt masque ou interpole dans la cellule donneuse
            {
              val = 0.;           
              for (E_Int jj = 0; jj < 2; jj++)
                for (E_Int ii = 0; ii < 2; ii++)
                {
                  ind = (ic+ii) + (jc+jj)*ni;
                  cellN0 = cellN[ind];
                  val += K_FUNC::E_abs(cf[ii+jj*2])*K_FUNC::E_abs(1.-cellN0);                  
                }              
              if (K_FUNC::fEqualZero(val, geomCutOff) == false) return 0;// pas interpolable 
            }
          }
        }//fin 2D
      }
      break;
      
    case K_INTERP::InterpAdt::O3ABC:
      //if (meshtype == 2) {printf("Error: getInterpolationData: interpType O3ABC not valid for unstructured zones.\n"); return -1;}
      if (ni < 3 || nj < 3 || nk < 3)
      { 
        //printf("Error: getInterpolationData: 3rd order interpolation requires at least 3 points per direction.\n");
        return -1;
      }
      found = searchInterpolationCellStruct(ni, nj, nk, xl, yl, zl, interpData,x, y, z, ic, jc, kc, cf); 
      if (found < 1) return found;
      ic = ic-1; jc = jc-1; kc = kc-1;
      //pour les cas 2d
      if (ic < 1) ic = 1;
      if (jc < 1) jc = 1;
      if (kc < 1) kc = 1;
      
      ics = ic; jcs = jc; kcs = kc;// pour le Lagrange: decalage de 1
      ic = ic-1; jc = jc-1; kc = kc-1;

      type = 3; indi[0] = ic + jc * ni + kc * ni*nj;
      if (ic == 0 || ic == ni-3 || jc == 0 || jc == nj-3 || kc == 0 || kc == nk-3) isBorder = 1;
      if (cellN != NULL)
      {
        if (nature == 0) // pas de pt masque ds la molecule donneuse
        {
          val = 1.;
          for (E_Int kk = 0; kk < 3; kk++)
            for (E_Int jj = 0; jj < 3; jj++)
              for (E_Int ii = 0; ii < 3; ii++)
              {
                ind = (ic + ii) + (jc + jj ) *ni + (kc+kk)*ni*nj;
                cellN0 = cellN[ind];
                d = cf[ii]*cf[jj+3]*cf[kk+6] > 1.e-12;
                val *= cellN0-d+1;
              }
          if (K_FUNC::fEqualZero(val, geomCutOff) == true) return 0;// pas interpolable 
        }
        else // pas de pt masque ou interpole dans la cellule donneuse       
        {
          val = 0.;
          for (E_Int kk = 0; kk < 3; kk++)
            for (E_Int jj = 0; jj < 3; jj++)
              for (E_Int ii = 0; ii < 3; ii++)
              {
                ind = (ic+ii) + (jc+jj)*ni + (kc+kk)*ni*nj;
                cellN0 = cellN[ind];
                val += K_FUNC::E_abs(cf[ii]*cf[jj+3]*cf[kk+6])*K_FUNC::E_abs(1.-cellN0);
              }
          if (K_FUNC::fEqualZero(val,geomCutOff) == false) return 0;// pas interpolable 
        }
      }
      corr = compLagrangeCoefs(x, y, z, ics, jcs, kcs, ni, nj, nk, xl, yl, zl, 
                               cf, interpType);
      if (corr == 1) // mauvaise approx de (x,y,z) -> ordre 2 type O2CF
      {
        found = searchInterpolationCellStruct(
          ni, nj, nk, xl, yl, zl, interpData, x, y, z, ic, jc, kc, cf); 
        if (found < 1) return found;
	ic = ic-1; jc = jc-1; kc = kc-1;
        if (ic == 0 || ic == ni-2 || jc == 0 || jc == nj-2 || kc == 0 || kc == nk-2) isBorder = 1;
        else isBorder = 0;
        type = 2; indi[0] = ic + jc*ni + kc*ni*nj;
      }
      break;
      
    case K_INTERP::InterpAdt::O5ABC:
      //if (meshtype == 2) {printf("Error: getInterpolationData: interpType O5ABC not valid for unstructured zones.\n"); return -1;}
      if (ni < 5 || nj < 5 || nk < 5)
      { 
        //printf("Error: getInterpolationData: 5th order interpolation requires at least 5 points per direction.\n");
        return -1;
      }
      found = searchInterpolationCellStruct(
        ni, nj, nk, xl, yl, zl, interpData,x, y, z, ic, jc, kc, cf); 
      if (found < 1) return found;
      //decalage pour frontiere max
      if (ic >= ni-1) ic = ni-4;
      else ic = ic-2;
    
      if (jc >= nj-1) jc = nj-4;
      else jc = jc-2;
      
      if (kc >= nk-1) kc = nk-4;
      else kc = kc-2;

      // pour les cas 2d
      if (ic < 1) ic = 1;
      if (jc < 1) jc = 1;
      if (kc < 1) kc = 1;

      ics = ic; jcs = jc; kcs = kc;//sauvegarde pour compLagrange pour qui les indices demarrent a 1
      ic = ic-1; jc = jc-1; kc = kc-1;
      type = 5; indi[0] = ic + jc * ni + kc * ni*nj;

      if (ic == 0 || ic == ni-5 || jc == 0 || jc == nj-5 || kc == 0 || kc == nk-5) isBorder = 1;

      if (cellN != NULL)
      {
        if (nature == 0) // traitement simple: pas de pt masque ds la molecule donneuse
        {
          val = 1.;
          for (E_Int kk = 0; kk < 5; kk++)
            for (E_Int jj = 0; jj < 5; jj++)
              for (E_Int ii = 0; ii < 5; ii++)
              {
                ind = (ic + ii) + (jc + jj)*ni + (kc+kk)*ni*nj;
                cellN0 = cellN[ind];
                d = cf[ii]*cf[jj+5]*cf[kk+10] > 1.e-12;
                val *= cellN0-d+1;
              }
          if (K_FUNC::fEqualZero(val, geomCutOff) == true) return 0;// pas interpolable 
        }
        else // nature=1 pas de pt masque ou interpole dans la cellule
        {
          val = 0.;
          for (E_Int kk = 0; kk < 5; kk++)
            for (E_Int jj = 0; jj < 5; jj++)
              for (E_Int ii = 0; ii < 5; ii++)
              {
                ind = (ic+ii) + (jc+jj)*ni + (kc+kk)*ni*nj;
                cellN0 = cellN[ind];
                val += K_FUNC::E_abs(cf[ii]*cf[jj+5]*cf[kk+10])*K_FUNC::E_abs(1.-cellN0);
              }
          if (K_FUNC::fEqualZero(val, geomCutOff) == false) return 0;// pas interpolable 
        } 
      }
      corr = compLagrangeCoefs(x, y, z, ics, jcs, kcs, ni, nj, nk, xl, yl, zl, 
                               cf, interpType);
      if (corr == 1) // mauvaise approx de (x,y,z) -> ordre 2
      {
        found = searchInterpolationCellStruct(
          ni, nj, nk, xl, yl, zl, interpData, x, y, z, ic, jc, kc, cf); 
        if (found < 1) return found;
        ic = ic-1; jc = jc-1; kc = kc-1;
        if (ic == 0 || ic == ni-2 || jc == 0 || jc == nj-2 || kc == 0 || kc == nk-2) isBorder = 1;
        else isBorder = 0;
        type = 2; indi[0] = ic + jc*ni + kc*ni*nj;
      }
      break;

    default:
      return -1;
  }
  return found;
}
// ============================================================================
/* Find the interpolation cell for point (x,y,z) for a tetra kmesh
   le numero de l'elt trouve est retourne dans noelt */
// ============================================================================
short K_INTERP::searchExtrapolationCellUnstruct(
  E_Float* xt, E_Float* yt, E_Float* zt,
  E_Float* cellNp,
  FldArrayI& connectEV,
  K_INTERP::InterpAdt* interpData,
  E_Float x, E_Float y, E_Float z,
  E_Int& noelt, FldArrayF& cf,
  E_Int nature, E_Int extrapOrder, E_Float constraint)
{
  // datas for interpolation cell (tetrahedra)
  E_Float xp, yp, zp;
  E_Float xq, yq, zq;
  E_Float xr, yr, zr;
  E_Float xs, ys, zs;
  E_Int indp, indq, indr, inds;
  E_Float xi, yi, zi;

  E_Int foundInDomain = 0;
  E_Int noeltsave = 0;
  E_Int ncf = 8;
  FldArrayF cf_sav(ncf);

  // Init cf
  cf.setAllValuesAtNull();

  // search list of candidate cells
  list<E_Int> listOfCandidateCells;
  E_Float alphaTol = K_FUNC::E_max( (2*constraint-5.)*0.1, 0.);
  E_Int found = interpData->getListOfCandidateCells(
    x, y, z, 
    listOfCandidateCells, alphaTol);
  if (found == 0) return 0; // listOfCandidateCells empty

  /* Find the right cell among candidate cells */
  list<E_Int>::iterator itr;
  E_Float sum;
  E_Float sum_coef;
  E_Float saved_sum_coef = K_CONST::E_MAX_FLOAT;
  E_Float max_coef;
  E_Int* cn1 = connectEV.begin(1);
  E_Int* cn2 = connectEV.begin(2);
  E_Int* cn3 = connectEV.begin(3);
  E_Int* cn4 = connectEV.begin(4);

  // Loop on all candidate cells : computation of coefficients cf
  for (itr = listOfCandidateCells.begin();
       itr != listOfCandidateCells.end();
       itr++)
  {
    // index of candidate cell
    noelt = *itr;
    
    // compute interpolation coefficients for tetrahedra
    indp = cn1[noelt]-1;
    indq = cn2[noelt]-1;
    indr = cn3[noelt]-1;
    inds = cn4[noelt]-1;
    
    xp = xt[indp]; yp = yt[indp]; zp = zt[indp];
    xq = xt[indq]; yq = yt[indq]; zq = zt[indq];
    xr = xt[indr]; yr = yt[indr]; zr = zt[indr];
    xs = xt[inds]; ys = yt[inds]; zs = zt[inds];
    coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                     xr, yr, zr, xs, ys, zs, xi, yi, zi);
    
    // cellN for points of interpolation cell
    E_Float cellNmp = 1.;
    E_Float cellNmq = 1.;
    E_Float cellNmr = 1.;
    E_Float cellNms = 1.;
    
    if (cellNp != NULL)
    {
      cellNmp = cellNp[indp];
      cellNmq = cellNp[indq];
      cellNmr = cellNp[indr];
      cellNms = cellNp[inds];
    }
    // extrapolation - order 0
    cf.setAllValuesAtNull();
    if (extrapOrder == 0)
    {
      for (E_Int i = 0; i < cf.getSize(); i++) cf[i]=0.;
      // test which point of the tetrahedron will be used for extrapolation
      max_coef = K_FUNC::E_max(xi,yi); max_coef = K_FUNC::E_max(max_coef,zi); 
      if (max_coef <= 0.5) cf[0] = 1.; // point P
      else if (K_FUNC::fEqualZero(max_coef,xi)) cf[1] = 1.; // point Q
      else if (K_FUNC::fEqualZero(max_coef,yi)) cf[2] = 1.; // point R
      else cf[3] = 1.; // point S         
    }
    else if (extrapOrder == 1) // extrapolation - order 1
    {
      cf[0] = 1-xi-yi-zi; 
      cf[1] = xi;
      cf[2] = yi;
      cf[3] = zi;
      sum_coef = K_FUNC::E_abs(cf[1])+K_FUNC::E_abs(cf[2])+K_FUNC::E_abs(cf[3]);
      // if sum_coef exceeds constraint, degenerate to order 0
      if (sum_coef >= constraint)
      {
        for (E_Int i = 0; i < cf.getSize(); i++) cf[i] = 0.;
        // test will point of the tetrahedra will be used for extrapolation
        max_coef = K_FUNC::E_max(xi,yi); max_coef = K_FUNC::E_max(max_coef,zi); 
        if (max_coef <= 0.5) cf[0] = 1.; // point P
        else if (K_FUNC::fEqualZero(max_coef,xi)) cf[1] = 1.; // point Q
        else if (K_FUNC::fEqualZero(max_coef,yi)) cf[2] = 1.; // point R
        else cf[3] = 1.; // point S         
      }
    }
    // Keep only valid extrapolation cell (no blanked cell in tetrahedron)
    sum = cellNmp*cellNmq*cellNmr*cellNms;
    if (sum > K_CONST::E_CUTOFF)
    {
      // Keep extrapolation cell with minimum distance (~ minimum sum of absolute coefficients) 
      // from the interpolated point
      sum_coef = K_FUNC::E_abs(cf[1])+K_FUNC::E_abs(cf[2])+K_FUNC::E_abs(cf[3]);
      if (sum_coef < saved_sum_coef)
      {
        if (sum_coef < constraint) foundInDomain = 1;
        else foundInDomain = 0;
        saved_sum_coef = sum_coef;
        noeltsave = noelt;
        cf_sav = cf;
      }
    }
  }

  noelt = noeltsave;
  cf = cf_sav;
  return foundInDomain;
}
// ============================================================================
/* Find the interpolation cell for point (x,y,z) in case of a structured donor
   zone */
// ============================================================================
short K_INTERP::searchExtrapolationCellStruct(
  E_Int ni, E_Int nj, E_Int nk, 
  E_Float* xl, E_Float* yl, E_Float* zl,
  E_Float* cellNp,
  K_INTERP::InterpAdt* interpData,
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

  // Init cf
  cf.setAllValuesAtNull();

  // search list of candidate cells
  list<E_Int> listOfCandidateCells;
  E_Float alphaTol = K_FUNC::E_max( (2*constraint-5.)*0.1, 0.);
  E_Int found = interpData->getListOfCandidateCells(
    x, y, z, 
    listOfCandidateCells, alphaTol);
  
  if (found == 0) return 0; // listOfCandidateCells empty

  /* Find the right cell among candidate cells */
  list<E_Int>::iterator itr;
  E_Float sum_coef;
  E_Float saved_sum_coef = K_CONST::E_MAX_FLOAT;
  E_Float saved_max_diff_coef = K_CONST::E_MAX_FLOAT;
  E_Float diff_coeff = K_CONST::E_MAX_FLOAT;
  E_Int ind;
  // parameters for extrapolation routine
  E_Int is, js, ks, ret;
  E_Int nij;

  // Loop on all candidate cells: computation of coefficients cf
  for (itr = listOfCandidateCells.begin();
       itr != listOfCandidateCells.end();
       itr++)
  {
    // 1D-index of interpolation cell
    ind = *itr;
 
    // (i,j,k)-indices of interpolation cell
    nij = ni*nj;
    k = ind/nij; 
    j = (ind-k*nij)/ni;
    i = ind-j*ni-k*nij;
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
    //printf("xyz: %f %f %f\n", x, y, z);
    //printf("%f %f %f %f %f %f %f %f\n", cfp[0],cfp[1],cfp[2],cfp[3],cfp[4],cfp[5],cfp[6],cfp[7]);
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
// ============================================================================
/* Find the interpolation cell for point (x,y,z) for a tetra kmesh
   le numero de l'elt trouve est retourne dans noelt */
// ============================================================================
short K_INTERP::searchInterpolationCellUnstruct(
  E_Float* xt, E_Float* yt, E_Float* zt,
  FldArrayI& connectEV,
  K_INTERP::InterpAdt* interpData,
  E_Float x, E_Float y, E_Float z,
  E_Int& noelt, FldArrayF& cf)
{ 
  cf.setAllValuesAtNull();

  // recherche de la liste des cellules candidates
  list<E_Int> listOfCandidateCells; 
  E_Int found = interpData->getListOfCandidateCells(
    x, y, z, listOfCandidateCells);
  if (found == 0) return 0; // listOfCandidateCells vide

  const E_Float EPS = _EPS_TETRA;
  E_Float xp, yp, zp;
  E_Float xq, yq, zq;
  E_Float xr, yr, zr;
  E_Float xs, ys, zs;
  E_Int indp, indq, indr, inds;
  E_Float xi, yi, zi, sum;
  E_Int et;
  list<E_Int>::iterator itr;
  E_Int* cn1 = connectEV.begin(1);
  E_Int* cn2 = connectEV.begin(2);
  E_Int* cn3 = connectEV.begin(3);
  E_Int* cn4 = connectEV.begin(4);

  for (itr = listOfCandidateCells.begin(); 
       itr != listOfCandidateCells.end();
       itr++)
  {
    et = *itr;
    indp = cn1[et]-1;
    indq = cn2[et]-1;
    indr = cn3[et]-1;
    inds = cn4[et]-1;

    xp = xt[indp]; yp = yt[indp]; zp = zt[indp];
    xq = xt[indq]; yq = yt[indq]; zq = zt[indq];
    xr = xt[indr]; yr = yt[indr]; zr = zt[indr];
    xs = xt[inds]; ys = yt[inds]; zs = zt[inds];

    coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                     xr, yr, zr, xs, ys, zs, xi, yi, zi);

    sum = xi+yi+zi;
    if (xi > -EPS && yi > -EPS && zi > -EPS && sum < K_CONST::ONE+3*EPS)
    {
      cf[0] = 1-sum; 
      cf[1] = xi; cf[2] = yi; cf[3] = zi;
      noelt = et;
      return 1;
    }
  }
  noelt = -1;
  return 0;
}

// ============================================================================
/* Find the interpolation cell for point (x,y,z) in case of a structured donor
   zone */
// ============================================================================
short K_INTERP::searchInterpolationCellStruct(
  E_Int ni, E_Int nj, E_Int nk, 
  E_Float* xl, E_Float* yl, E_Float* zl,
  K_INTERP::InterpAdt* interpData,
  E_Float x, E_Float y, E_Float z,
  E_Int& ic, E_Int& jc, E_Int& kc,
  FldArrayF& cf)
{
  // recherche de la liste des cellules candidates
  list<E_Int> listOfCandidateCells; 
  E_Int found = interpData->getListOfCandidateCells(x, y, z, 
                                                    listOfCandidateCells);
  /* Find the right cell among candidate cells */
  
  if (found == 0) return 0; // listOfCandidateCells vide
  
  list<E_Int>::iterator itr;
  E_Boolean JUMP;
  E_Int ind, i, isomm, is, js, ks;
  E_Float xi, yi, zi;
  E_Float xt[15], yt[15], zt[15];
  E_Int nij = ni*nj;
  
  /* For 8 cells, apply the technique of jump to find interpolation cell */
  JUMP = false;
  itr = listOfCandidateCells.begin();
  ind = *itr;
  for (i = 0; i < 8; i++)
  { 
    coordHexa(ind, ni, nj, nk,
              xl, yl, zl,
              ic, jc, kc,
              xt, yt, zt);
    
    if (getCellJump(x, y, z,
                    xt, yt, zt,              
                    isomm,
                    xi, yi, zi) == false)
    {
      if (i == 7)
      {
        JUMP =  false;
        break;
      }
      /* Apply a technique of jump */
      is = K_FUNC::E_min(8, E_Int((xi+K_FUNC::E_sign(xi))*K_CONST::ONE_HALF));
      is = K_FUNC::E_max(-8, is);
      js = K_FUNC::E_min(8, E_Int((yi+K_FUNC::E_sign(yi))*K_CONST::ONE_HALF));
      js = K_FUNC::E_max(-8, js);
      ks = K_FUNC::E_min(8, E_Int((zi+K_FUNC::E_sign(zi))*K_CONST::ONE_HALF));
      ks = K_FUNC::E_max(-8, ks);
      ind = ind+is+js*ni+ks*nij;

      kc = ind/nij;
      E_Int kcnij = kc*nij;
      jc = (ind-kcnij)/ni;
      E_Int jcni = jc*ni;
      ic = ind-jcni-kcnij;

      if (ic<0 || ic>ni-2 || jc<0 || jc>nj-2 || kc<0 || kc>nk-2)
      {
        JUMP =  false;
        break;
      }
    }
    else
    {        
      for (itr = listOfCandidateCells.begin();
           itr != listOfCandidateCells.end();
           itr++)
      {
        if (ind == *itr)
        {
          JUMP = true;
          goto saut;
        }
      }
      JUMP =  false;
      break;
    }
  }

  /* If the technique of jump succeeds, find the interpolation coefficients
     in the cell by cut it in 24 tetrahedras */
  saut:
  if (JUMP)  
  {
    if (getCoeffInterpHexa(x, y, z,
                          isomm,
                          xi, yi, zi, 
                          xt, yt, zt,
                          cf) == true)
    {      
      return 1;
    }
    else
    {
      JUMP =  false;
      listOfCandidateCells.erase(itr);
    }
  }
 
  /* If the technique of jump fails, we test all candidate cells */
  if (JUMP ==  false) 
  {    
    for (itr = listOfCandidateCells.begin();
         itr != listOfCandidateCells.end();
         itr++)
    {
      ind = (*itr);
          
      coordHexa(ind, ni, nj, nk,
                xl, yl, zl,
                ic, jc, kc,
                xt, yt, zt);
      if (coeffInterpHexa(x, y, z,
                          xt, yt, zt,
                          cf) == true)
      {
        return 1;
      }
    }
  }
  return 0;
}

// ============================================================================
/* Get the hexahedra coordinates+coord of centers of faces+center
   This is the reference for the cell edges numerotation. */
// ============================================================================
void K_INTERP::coordHexa(E_Int ind, E_Int ni, E_Int nj, E_Int nk,
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
E_Boolean K_INTERP::coeffInterpHexa(E_Float x, E_Float y, E_Float z,
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
  is1 = E_Int((K_CONST::ONE + K_FUNC::E_sign(xi-yi)) * K_CONST::ONE_HALF);
  is2 = E_Int((K_CONST::ONE + K_FUNC::E_sign(yi-zi)) * K_CONST::ONE_HALF);
  is3 = E_Int((K_CONST::ONE + K_FUNC::E_sign(zi-xi)) * K_CONST::ONE_HALF);
  
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
    cf0 = K_CONST::ONE_EIGHT*(K_CONST::ONE-xi-yi-zi);
    
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
          cf0 = K_CONST::ONE_EIGHT*(1.-xi-yi-zi);
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
E_Boolean K_INTERP::getCellJump(E_Float x, E_Float y, E_Float z,
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
void K_INTERP::coeffInterpTetra(E_Float x, E_Float y, E_Float z,
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
E_Boolean K_INTERP::getCoeffInterpHexa(E_Float x, E_Float y, E_Float z,
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
  E_Int is1 = static_cast<E_Int>((K_CONST::ONE + K_FUNC::E_sign(xi-yi)) * K_CONST::ONE_HALF);
  E_Int is2 = static_cast<E_Int>((K_CONST::ONE + K_FUNC::E_sign(yi-zi)) * K_CONST::ONE_HALF);
  E_Int is3 = static_cast<E_Int>((K_CONST::ONE + K_FUNC::E_sign(zi-xi)) * K_CONST::ONE_HALF);
  
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
      
    cf0 = K_CONST::ONE_EIGHT*(K_CONST::ONE-xi-yi-zi);
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
          cf0 = K_CONST::ONE_EIGHT*(1.-xi-yi-zi);
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
// compute Lagrange polynomials: si x,y,z mal approxime alors compLagrangeCoefs
// retourne 1, sinon si bonne approximation retourne 0
/* IN: xl, yl, zl: maillage sur lequel les coefs d interpolation sont calcules 
   IN: ni, nj, nk: taille du maillage ci dessus
 */
//=============================================================================
short K_INTERP::compLagrangeCoefs(E_Float x, E_Float y, E_Float z,
                                  E_Int ic, E_Int jc, E_Int kc,
                                  E_Int ni, E_Int nj, E_Int nk,
                                  E_Float* xl, E_Float* yl, E_Float* zl,
                                  FldArrayF& cf, 
                                  K_INTERP::InterpAdt::InterpolationType interpType)
{
  E_Int npts = ni*nj*nk;
  E_Float ksi, eta, zeta; // coordonnees de (x,y,z) dans element de ref
  // calcul des coordonnees ksi, eta, zeta
  E_Int npts_interp_1D, npts_interp_3D;

  switch (interpType)
  {
    case K_INTERP::InterpAdt::O2ABC:
      npts_interp_1D = 2;
      npts_interp_3D = 8;
      break;
    case K_INTERP::InterpAdt::O3ABC: 
      npts_interp_1D = 3;
      npts_interp_3D = 27;
      break;
    case K_INTERP::InterpAdt::O5ABC:
      npts_interp_1D = 5;
      npts_interp_3D = 125;
      break;
    default :
      printf("Error: compLagrangeCoefs: not a valid interpolation type.\n"); 
      return -1;
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
  // ICI ic,jc,kc demarrent a 1
  compinterpolatedptinrefelt_(
    xl, yl, zl, npts, ic, jc, kc, ni, nj, 
    x, y, z, npts_interp_1D, npts_interp_3D, 
    x_interp.begin(), y_interp.begin(), z_interp.begin(), 
    base.begin(), A1.begin(), B1.begin(),
    indxc.begin(), indxr.begin(), ipiv.begin(),  
    ksi, eta, zeta, err);
  if (err == 1)  return 1; // point interpole mal calcule 

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
    case K_INTERP::InterpAdt::O2ABC:
      cfp[0] = ksim;
      cfp[1] = ksi;
      cfp[2] = etam;
      cfp[3] = eta;
      cfp[4] = zetam;
      cfp[5] = zeta;
      break;

    case K_INTERP::InterpAdt::O3ABC:
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

    case K_INTERP::InterpAdt::O5ABC:
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
      return -1;
  }
  return 0;
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
short K_INTERP::getExtrapolationCoeffForCell(
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
  E_Boolean valid = false;
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
        cf0 = K_CONST::ONE_EIGHT*(K_CONST::ONE-xi-yi-zi);
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
      if ( dim == 3)
      {
        for (E_Int kk = 0; kk < 2; kk++)
          for (E_Int jj = 0; jj < 2; jj++)
            for (E_Int ii = 0; ii < 2; ii++)
            {
              E_Int ind = (ic + ii) + (jc + jj)*ni + (kc+kk)*ni*nj;
              if (nature == 1) locCellN[icell] = 1.-(cellNp[ind]-1.);
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
            // if (nature == 1) locCellN[icell] = 1.-(cellNp[ind]-1.);
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
