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
# include "Interp/Interp.h"
# include <list>
using namespace std;
using namespace K_FLD;
// E_Float _EPS_DET   = 1.e-16;
// E_Float _EPS_TETRA = 1.e-4;
// E_Float _EPS_GEOM  = K_CONST::E_GEOM_CUTOFF;
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

# define CELLN(order)                        \
  if (nature == 0)                           \
  {                                          \
      val = 1.;                              \
      for (E_Int kk = 0; kk < order; kk++)       \
        for (E_Int jj = 0; jj < order; jj++)     \
          for (E_Int ii = 0; ii < order; ii++)   \
            {                                \
              ind = (ic + ii) + (jc + jj) *ni + (kc+kk)*ni*nj;\
              cellN0 = cellN[ind];           \
              d = cf[ii]*cf[jj+order]*cf[kk+order*2] > 1.e-12;\
              val *= cellN0-d+1;             \
            }                                \
      if (K_FUNC::fEqualZero(val, geomCutOff) == true) return 0;\
  }                                          \
  else                                       \
  {                                          \
    val = 0.;                                \
    for (E_Int kk = 0; kk < order; kk++)         \
      for (E_Int jj = 0; jj < order; jj++)       \
        for (E_Int ii = 0; ii < order; ii++)     \
          {                                  \
            ind = (ic+ii) + (jc+jj)*ni + (kc+kk)*ni*nj;\
            cellN0 = cellN[ind];             \
            val += K_FUNC::E_abs(cf[ii]*cf[jj+order]*cf[kk+order*2])*K_FUNC::E_abs(1.-cellN0);\
          }                                  \
    if (K_FUNC::fEqualZero(val,geomCutOff) == false) return 0;\
  }

//=============================================================================
/* a est la connectivite cEEN dans le cas non structure, NULL si pas besoin
   IN: x,y,z: point a extrapoler
   IN: InterpData: l'InterpData du maillage donneur
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
  InterpData* InterpData, void* c1,  void* c2, 
  E_Int meshtype, E_Int ni, E_Int nj, E_Int nk,
  E_Float* xl, E_Float* yl, E_Float* zl, E_Float* cellN,
  E_Int& isBorder, E_Int& type, FldArrayI& indi, FldArrayF& cf, 
  K_INTERP::InterpData::InterpolationType interpType,
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
    case K_INTERP::InterpData::O2CF:
    if (InterpData->_topology == 2) 
    {
      FldArrayI& cEV = *(FldArrayI*)c1;
      found = InterpData->searchExtrapolationCellUnstruct(xl, yl, zl, cellN, cEV, 
                                                         x, y, z, noet, cf, nature, extrapOrder, constraint);
      if (found < 1) return found;
      type = 4; indi[0] = noet;
    }
    else 
    {
      if (InterpData->_topology == 1) 
        found = InterpData->searchExtrapolationCellStruct(ni, nj, nk, xl, yl, zl, cellN, 
                                                          x, y, z, ic, jc, kc, cf, nature, extrapOrder, constraint);
      else if (InterpData->_topology == 0) 
        found = InterpData->searchExtrapolationCellCart(ni, nj, nk, xl, yl, zl, cellN, 
                                                        x, y, z, ic, jc, kc, cf, nature, 
                                                        extrapOrder, constraint);
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
      
    case K_INTERP::InterpData::O3ABC:
    case K_INTERP::InterpData::O5ABC:
    if (InterpData->_topology == 2) 
    {
      return -1;
    }
    else  
    {
      found = InterpData->searchExtrapolationCellStruct(ni, nj, nk, xl, yl, zl, cellN, 
                                                       x, y, z, ic, jc, kc, cf, nature, extrapOrder, constraint);
      if (found < 1) return found;
      ic = ic-1; jc = jc-1; kc = kc-1;// pour demarrer � 0        

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
  InterpData* InterpData, void* c1, void* c2, 
  E_Int meshtype, E_Int ni, E_Int nj, E_Int nk,
  E_Float* xl, E_Float* yl, E_Float* zl, E_Float* cellN,
  E_Int& isBorder, E_Int& type, FldArrayI& indi, FldArrayF& cf, 
  K_INTERP::InterpData::InterpolationType interpType,
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
  E_Int nij = ni*nj;
  //2D case : z must be equal to _zmin 
  if (dim == 2) 
  {
    if (K_FUNC::E_abs(z-InterpData->_zmin)>geomCutOff) return 0;
  }
  
  switch (interpType)
  {
    case K_INTERP::InterpData::O2CF:
      if (InterpData->_topology==2)
      {
        FldArrayI& cEV = *(FldArrayI*)c1;
        found = InterpData->searchInterpolationCellUnstruct(xl, yl, zl, cEV, x, y, z, noet, cf);
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
              ind = cEV(noet, nov)-1;
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
              ind = cEV(noet, nov)-1; 
              cellN0 = cellN[ind];
              val += K_FUNC::E_abs(cf[nov-1])*K_FUNC::E_abs(1.-cellN0);
            }
            if (K_FUNC::fEqualZero(val, geomCutOff) == false) return 0;// pas interpolable 
          }
        }
      }
      else
      { 
        if (InterpData->_topology == 1)
          found = InterpData->searchInterpolationCellStruct(ni, nj, nk, xl, yl, zl, 
                                                            x, y, z, ic, jc, kc, cf); 
        else if (InterpData->_topology == 0) //cart
          found = InterpData->searchInterpolationCellCartO2(ni, nj, nk, x, y, z, ic, jc, kc, cf);

        if (found < 1) return found;
        
        ic = ic-1; jc = jc-1; kc = kc-1; // pour demarrer a 0   
        if (dim == 3)
        {
          if (ic == 0 || ic == ni-2 || jc == 0 || jc == nj-2 || kc == 0 || kc == nk-2) isBorder = 1;

          type = 2; indi[0] = ic + jc*ni + kc*nij;
          
          if (cellN != NULL)
          {
            if (nature == 0) // pas de pt masque ds la molecule donneuse
            {
              val = 1.;
              for (E_Int kk = 0; kk < 2; kk++)
                for (E_Int jj = 0; jj < 2; jj++)
                  for (E_Int ii = 0; ii < 2; ii++)
                  {
                    ind = (ic+ii) + (jc+jj)*ni + (kc+kk)*nij;
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
                    ind = (ic+ii) + (jc+jj)*ni + (kc+kk)*nij;
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
      
    case K_INTERP::InterpData::O3ABC:
      if (InterpData->_topology == 2) return -1;

      if (ni < 3 || nj < 3 || nk < 3)
      { 
        //printf("Error: getInterpolationData: 3rd order interpolation requires at least 3 points per direction.\n");
        return -1;
      }

      if (InterpData->_topology == 1) //ADT sur maillage curviligne
      { 
        found = InterpData->searchInterpolationCellStruct(ni, nj, nk, xl, yl, zl, x, y, z, ic, jc, kc, cf); 
    
        if (found < 1) return found;
        
        ic = ic-1; jc = jc-1; kc = kc-1;
        //pour les cas 2d
        if (ic < 1) ic = 1;
        if (jc < 1) jc = 1;
        if (kc < 1) kc = 1;

        ics = ic; jcs = jc; kcs = kc;// pour le Lagrange: decalage de 1
        ic = ic-1; jc = jc-1; kc = kc-1;

        if (ic == 0 || ic == ni-3 || jc == 0 || jc == nj-3 || kc == 0 || kc == nk-3) isBorder = 1;
        corr = compLagrangeCoefs(x, y, z, ics, jcs, kcs, ni, nj, nk, xl, yl, zl, 
                                 cf, interpType);
        type = 3; indi[0] = ic + jc*ni + kc*nij;
        if (cellN != NULL) { CELLN(3); }

        if (corr == 1) // mauvaise approx de (x,y,z) -> ordre 2 type O2CF
        {
          found = InterpData->searchInterpolationCellStruct(ni, nj, nk, xl, yl, zl, x, y, z, ic, jc, kc, cf); 
          if (found < 1) return found;
          ic = ic-1; jc = jc-1; kc = kc-1;
          if (ic == 0 || ic == ni-2 || jc == 0 || jc == nj-2 || kc == 0 || kc == nk-2) isBorder = 1;
          else isBorder = 0;
          type = 2; indi[0] = ic + jc*ni + kc*nij;
        }
      }
      else if (InterpData->_topology == 0) // CART sur maillage cartesien
      {
        found = InterpData->searchInterpolationCellCartO3(ni, nj, nk, x, y, z, ic, jc, kc, cf);
        if (found < 1) return found;

        ic = ic-1; jc = jc-1; kc = kc-1;// indices demarrent a 0
        type = 3; indi[0] = ic + jc*ni + kc*nij;
        if (ic == 0 || ic == ni-3 || jc == 0 || jc == nj-3 || kc == 0 || kc == nk-3) isBorder = 1;
        if (cellN != NULL) { CELLN(3); }
      }
      break;
      
    case K_INTERP::InterpData::O4ABC:
      if (InterpData->_topology == 2) return -1;

      if (ni < 4 || nj < 4 || nk < 4)
      {
        //printf("Error: getInterpolationData: 4rd order interpolation requires at least 4 points per direction.\n");
        return -1;
      }


      if (InterpData->_topology == 1) //ADT sur maillage curviligne
      {
        found = InterpData->searchInterpolationCellStruct(ni, nj, nk, xl, yl, zl, x, y, z, ic, jc, kc, cf);

        if (found < 1) return found;

        if (ic   > 1 ) ic -=1;          //on decale le donneur de 1
        if (ic+3 > ni) ic -= ic+3- ni;  // on decentre la cellule donneuse sitrop proche bord
        if (jc   > 1 ) jc -=1;
        if (jc+3 > nj) jc -= jc+3-nj;
        if (kc   > 1 ) kc -=1;
        if (kc+3 > nk) kc -= kc+3-nk;

        corr = compLagrangeCoefs(x, y, z, ic, jc, kc, ni, nj, nk, xl, yl, zl, cf, interpType);
        //ics = ic; jcs = jc; kcs = kc;// pour le Lagrange: decalage de 1

        ic = ic-1; jc = jc-1; kc = kc-1;  // indices demarrent a 0 en c

        if (ic == 0 || ic == ni-4 || jc == 0 || jc == nj-4 || kc == 0 || kc == nk-4) isBorder = 1;

        type = 44; indi[0] = ic + jc*ni + kc*nij;
        if (cellN != NULL) { CELLN(4); }

        if (corr == 1) // mauvaise approx de (x,y,z) -> ordre 2 type O2CF
        {
          found = InterpData->searchInterpolationCellStruct(ni, nj, nk, xl, yl, zl, x, y, z, ic, jc, kc, cf);
          if (found < 1) return found;
          ic = ic-1; jc = jc-1; kc = kc-1;
          if (ic == 0 || ic == ni-2 || jc == 0 || jc == nj-2 || kc == 0 || kc == nk-2) isBorder = 1;
          else isBorder = 0;
          type = 2; indi[0] = ic + jc*ni + kc*nij;
        }
      }
      else if (InterpData->_topology == 0) // CART sur maillage cartesien
      {
        found = InterpData->searchInterpolationCellCartO4(ni, nj, nk, x, y, z, ic, jc, kc, cf);
        if (found < 1) return found;

        ic = ic-1; jc = jc-1; kc = kc-1;// indices demarrent a 0
        type = 44; indi[0] = ic + jc*ni + kc*nij;
        if (ic == 0 || ic == ni-4 || jc == 0 || jc == nj-4 || kc == 0 || kc == nk-4) isBorder = 1;
        if (cellN != NULL) { CELLN(4); }
      }
      break;


    case K_INTERP::InterpData::O5ABC:

      if (InterpData->_topology == 2) return -1;
      if (ni < 5 || nj < 5 || nk < 5)
      { 
        //printf("Error: getInterpolationData: 5th order interpolation requires at least 5 points per direction.\n");
        return -1;
      }
      if (InterpData->_topology == 1)
      {
        found = InterpData->searchInterpolationCellStruct(ni, nj, nk, xl, yl, zl, x, y, z, ic, jc, kc, cf); 
      }
      else if (InterpData->_topology==0)
      {
        printf(" INTERP O5ABC NOT YET IMPLEMENTED FOR CARTESIAN GRIDS \n"); return -1;
      }
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
      type = 5; indi[0] = ic + jc * ni + kc * nij;

      if (ic == 0 || ic == ni-5 || jc == 0 || jc == nj-5 || kc == 0 || kc == nk-5) isBorder = 1;

      if (cellN != NULL) {  CELLN(5); }
      corr = compLagrangeCoefs(x, y, z, ics, jcs, kcs, ni, nj, nk, xl, yl, zl, 
                               cf, interpType);
      if (corr == 1) // mauvaise approx de (x,y,z) -> ordre 2
      {
        found = InterpData->searchInterpolationCellStruct(
          ni, nj, nk, xl, yl, zl, x, y, z, ic, jc, kc, cf); 
        if (found < 1) return found;
        ic = ic-1; jc = jc-1; kc = kc-1;
        if (ic == 0 || ic == ni-2 || jc == 0 || jc == nj-2 || kc == 0 || kc == nk-2) isBorder = 1;
        else isBorder = 0;
        type = 2; indi[0] = ic + jc*ni + kc*nij;
      }
      break;

    default:
      return -1;
  }
  return found;
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
                                  K_INTERP::InterpData::InterpolationType interpType)
{
  E_Int npts = ni*nj*nk;
  E_Float ksi, eta, zeta; // coordonnees de (x,y,z) dans element de ref
  // calcul des coordonnees ksi, eta, zeta
  E_Int npts_interp_1D, npts_interp_3D;

  switch (interpType)
  {
    case K_INTERP::InterpData::O2ABC:
      npts_interp_1D = 2;
      npts_interp_3D = 8;
      break;
    case K_INTERP::InterpData::O3ABC:
      npts_interp_1D = 3;
      npts_interp_3D = 27;
      break;
    case K_INTERP::InterpData::O4ABC:
      npts_interp_1D = 4;
      npts_interp_3D = 64;
      break;
    case K_INTERP::InterpData::O5ABC:
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

  //printf("interpTyp %d  %f %f %f %f  %f \n",interpType, ksim2 , ksim,  ksi, ksip, ksip2);

  switch (interpType)
  {
    case K_INTERP::InterpData::O2ABC:
      cfp[0] = ksim;
      cfp[1] = ksi;
      cfp[2] = etam;
      cfp[3] = eta;
      cfp[4] = zetam;
      cfp[5] = zeta;
      break;

    case K_INTERP::InterpData::O3ABC:
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

    case K_INTERP::InterpData::O4ABC:
      cfp[0] = -inv6              * ksim2  * ksim  * ksi ;  //alpha4  p
      cfp[1] =  K_CONST::ONE_HALF * ksim2  * ksim  * ksip;  //alpha3  i
      cfp[2] =  K_CONST::ONE_HALF * ksim2  * ksi   * ksip;  //alpha2  m
      cfp[3] = -inv6              * ksim   * ksi   * ksip;  //alpha1  m2

      cfp[4] = -inv6              * etam2  * etam  * eta ;  //beta4   p
      cfp[5] =  K_CONST::ONE_HALF * etam2  * etam  * etap;  //beta3   i
      cfp[6] =  K_CONST::ONE_HALF * etam2  * eta   * etap;  //beta2   m
      cfp[7] = -inv6              * etam   * eta   * etap;  //beta1   m2

      cfp[8] = -inv6              * zetam2 * zetam * zeta ; //zeta4   p
      cfp[9] =  K_CONST::ONE_HALF * zetam2 * zetam * zetap; //zeta3   i
      cfp[10]=  K_CONST::ONE_HALF * zetam2 * zeta  * zetap; //zeta2   m
      cfp[11]= -inv6              * zetam  * zeta  * zetap; //zeta1   m2
      break;

    case K_INTERP::InterpData::O5ABC:
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
