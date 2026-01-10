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

# include <stdlib.h>
# include "BlkInterp.h"
# include "Def/DefCplusPlusConst.h"

using namespace K_FUNC;
using namespace std;
using namespace K_FLD;

// ============================================================================
/* Constructor */
// ============================================================================
K_KINTERP::BlkInterpData::BlkInterpData() :
  _EPS_DET(1.e-16), _EPS_TETRA(1.e-4), _EPS_GEOM(K_CONST::E_GEOM_CUTOFF)
{
}

// ============================================================================
/* Destructor */
// ============================================================================
K_KINTERP::BlkInterpData::~BlkInterpData()
{
}

// ============================================================================
/* Search and return interpolation cells with coefficients. 
   Cas  non structure : indi retourne 4 et le no de l elt donneur */
// ============================================================================
short K_KINTERP::BlkInterpData::
getInterpolationCellUnstruct( E_Float x, E_Float y, E_Float z,
                              E_Int& noet,
                              FldArrayI& indi, FldArrayF& cf )
{
  E_Int res =  searchInterpolationCellUnstruct(x, y, z, noet, indi, cf);
  indi[0] = 4; indi[1] = noet;
  return res;
}
// ============================================================================
/* Search and return interpolation cells with coefficients. 
   Cas  structure*/
// ============================================================================
short K_KINTERP::BlkInterpData::
getInterpolationCellStruct(E_Float x, E_Float y, E_Float z,
                           E_Int& ic, E_Int& jc, E_Int& kc,
                           FldArrayF& cf,   
                           K_KINTERP::BlkInterpData::InterpMeshType interpMeshType,
                           K_KINTERP::BlkInterpData::InterpolationType interpType)
{
  short found;
  switch ( interpType )
  {
    case O2CF :
      found = searchInterpolationCellO2CF( x, y, z, ic, jc, kc, cf );
      if ( found < 1) return found;      
      break;
      
    case O3CF :
      found = searchInterpolationCellO3CF( x, y, z, ic, jc, kc, cf );
      if ( found < 1) return found;
      break;   

    case O3ABC :
      found = searchInterpolationCellO3ABC( x, y, z, ic, jc, kc);
      if ( found < 1) return found;   
      break;

    case O5ABC :
      found = searchInterpolationCellO5ABC( x, y, z, ic, jc, kc);
      if ( found < 1) return found;    
      break;
  
    default:
      printf("getInterpolationCellStruct : unknown type of interpolation.\n");
      return -1;
      break;
  }
  return found;
}
// ============================================================================
/* Search and return interpolation cells with coefficients. 
   Cas  structure*/
// ============================================================================
short K_KINTERP::BlkInterpData::
getInterpolationCellStruct(E_Float x, E_Float y, E_Float z,
                           FldArrayI& indi,
                           FldArrayF& cf,   
                           K_KINTERP::BlkInterpData::InterpMeshType interpMeshType,
                           K_KINTERP::BlkInterpData::InterpolationType interpType)
{
  short found;
  E_Int ic, jc, kc;
  short corr = 0; // correction d ordre 2 -> corr=1 
  switch ( interpType )
  {
    case O2CF :
      if ( indi.getSize() != 7 )
      {
        printf("Error: getInterpolationCellStruct: wrong size of indi.\n");
        exit(0);  
      }
      found = searchInterpolationCellO2CF( x, y, z, ic, jc, kc, cf );
      if ( found < 1) return found;      
      break;
      
    case O3CF :
      if ( indi.getSize() != 10 )
      {
        printf("Error: getInterpolationCellStruct: wrong size of indi.\n");
        exit(0);
      }
      found = searchInterpolationCellO3CF( x, y, z, ic, jc, kc, cf );
      if ( found < 1) return found;
      break;   

    case O3ABC :
      if ( indi.getSize() != 10 )
      {
        printf("Error: getInterpolationCellStruct: wrong size of indi.\n");
        exit(0);
      }
      found = searchInterpolationCellO3ABC( x, y, z, ic, jc, kc);
      if (found < 1) return found;
      corr =  compLagrangeCoefs(x, y, z, ic, jc, kc, cf, 
                                interpType, interpMeshType);
      if (corr == 1) // mauvaise approx de (x,y,z) -> ordre 2
        return correctInterpToO2CF(x, y, z, indi, cf, interpMeshType);    
  
      break;

    case O5ABC :
      if ( indi.getSize() != 16 )
      {
        printf("Error: getInterpolationCellStruct: wrong size of indi.");
        exit(0);  
      } 
      found = searchInterpolationCellO5ABC( x, y, z, ic, jc, kc);
      if (found < 1) return found;
      
      corr = compLagrangeCoefs(x, y, z, ic, jc, kc, cf, 
                               interpType, interpMeshType);
      if (corr == 1) // mauvaise approx de (x,y,z) -> ordre 2
        return correctInterpToO2CF(x, y, z, indi, cf, interpMeshType);    
     
      break;
  
    default:
      printf("Error: getInterpolationCellStruct: unknown type of interpolation.\n");
      return -1;
      break;
  }

  // Indices par direction des points de la molecule d'interpolation
  if (interpMeshType == EXT_CENTERS)
    fromExtendedToStandardCenters(ic, jc, kc, indi, interpType);
  else 
    compStandardIndices(ic, jc, kc, indi, interpType); 

  return found;
}
// ============================================================================
void K_KINTERP::BlkInterpData::
getInterpolationCellStructv(FldArrayF& coord, 
                            E_Int istart, E_Int iend,
                            FldArrayI& indi,
                            FldArrayF& cf, FldArrayIS& found, 
                            FldArrayI& extrap,
                            K_KINTERP::BlkInterpData::InterpMeshType interpMeshType,
                            K_KINTERP::BlkInterpData::InterpolationType interpType)
{
  E_Int size = coord.getSize();
  FldArrayI icv(size); 
  FldArrayI jcv(size);
  FldArrayI kcv(size);
  FldArrayIS corr(size);
  switch ( interpType )
  {
    case O2CF :
      if (indi.getNfld() != 7 )
      {
        printf("Error: getInterpolationCellStructv: wrong size of indi.\n");
        exit(0);
      }
      searchInterpolationCellO2CFv(coord, istart, iend, icv, jcv, kcv, 
                                   cf, found);
      // Indices par direction des points de la molecule d'interpolation
      if ( interpMeshType == EXT_CENTERS)
        fromExtendedToStandardCentersv(
          istart, iend, icv, jcv, kcv, indi, extrap, interpType);
      else 
        compStandardIndicesv(istart, iend, icv, jcv, kcv, indi, interpType);      
      break;

    case O3CF:
      if (indi.getNfld() != 10 )
      {
        printf("Error: getInterpolationCellStructv: wrong size of indi.\n");
        exit(0);
      }
      searchInterpolationCellO3CFv(coord, istart, iend, icv, jcv, kcv, 
                                   cf, found);

      // Indices par direction des points de la molecule d'interpolation
      if ( interpMeshType == EXT_CENTERS)
        fromExtendedToStandardCentersv(
          istart, iend, icv, jcv, kcv, indi, extrap, interpType);
      else 
        compStandardIndicesv(istart, iend, icv, jcv, kcv, indi, interpType);
      
      break;
   
    case O3ABC:
      if (indi.getNfld() != 10 )
      {
        printf("getInterpolationCellStructv : wrong size of indi.\n");
        exit(0);
      }
      searchInterpolationCellO3ABCv(coord, istart, iend, 
                                    icv, jcv, kcv, found);
      compLagrangeCoefsv(coord, istart, iend, icv, jcv, kcv, found, cf, 
                         corr, interpType, interpMeshType);

      // Indices par direction des points de la molecule d'interpolation
      if ( interpMeshType == EXT_CENTERS)
        fromExtendedToStandardCentersv(
          istart, iend, icv, jcv, kcv, indi, extrap, interpType);
      else 
        compStandardIndicesv(istart, iend, icv, jcv, kcv, indi, interpType);
    
      
      // modification des coefficients d interpolation : ordre 2 O2CF 
      // pour les points mal approximes 
      // modif locale de indi
      correctInterpToO2CFv(coord, istart, iend, indi,
                           found, cf, corr, interpMeshType);
      break;

    case O5ABC :
      if (indi.getNfld() != 16)
      {
        printf("Error: getInterpolationCellStructv: wrong size of indi.\n");
        exit(0);
      }
      searchInterpolationCellO5ABCv(coord, istart, iend, 
                                    icv, jcv, kcv, found);
      compLagrangeCoefsv(coord, istart, iend, icv, jcv, kcv, found, cf, 
                         corr, interpType, interpMeshType);
      
      // Indices par direction des points de la molecule d'interpolation
      if ( interpMeshType == EXT_CENTERS)
        fromExtendedToStandardCentersv(
          istart, iend, icv, jcv, kcv, indi, extrap, interpType);
      else 
        compStandardIndicesv(istart, iend, icv, jcv, kcv, indi, interpType);
    
      
      // modification des coefficients d'interpolation : ordre 2 O2CF 
      // pour les points mal approximes 
      // modif locale de indi
      correctInterpToO2CFv(coord, istart, iend, indi,
                           found, cf, corr, interpMeshType);
      break;

    default:
      printf("Error: getInterpolationCellStructv: unknown type of interpolation.\n");
      break;
  }  
}

// ============================================================================
void K_KINTERP::BlkInterpData::
getInterpolationCellStructv(FldArrayF& coord, 
                            E_Int istart, E_Int iend,
                            FldArrayI& indi, FldArrayI& icv, FldArrayI& jcv, FldArrayI& kcv,
                            FldArrayF& cf, FldArrayIS& found, 
                            FldArrayI& extrap,
                            K_KINTERP::BlkInterpData::InterpMeshType interpMeshType,
                            K_KINTERP::BlkInterpData::InterpolationType interpType)
{
  E_Int size = coord.getSize();
  FldArrayIS corr(size);

  switch ( interpType )
  {
    case O2CF :
      if (indi.getNfld() != 7 )
      {
        printf("Error: getInterpolationCellStructv: wrong size of indi.\n");
        exit(0);
      }
      searchInterpolationCellO2CFv(coord, istart, iend, icv, jcv, kcv, 
                                   cf, found);

      // Indices par direction des points de la molecule d'interpolation
      if ( interpMeshType == EXT_CENTERS)
        fromExtendedToStandardCentersv(
          istart, iend, icv, jcv, kcv, indi, extrap, interpType);
      else 
        compStandardIndicesv(istart, iend, icv, jcv, kcv, indi, interpType);      
      break;

    case O3CF:
      if (indi.getNfld() != 10 )
      {
        printf("Error: getInterpolationCellStructv: wrong size of indi.\n");
        exit(0);
      }
      searchInterpolationCellO3CFv(coord, istart, iend, icv, jcv, kcv, 
                                   cf, found);

      // Indices par direction des points de la molecule d'interpolation
      if ( interpMeshType == EXT_CENTERS)
        fromExtendedToStandardCentersv(
          istart, iend, icv, jcv, kcv, indi, extrap, interpType);
      else 
        compStandardIndicesv(istart, iend, icv, jcv, kcv, indi, interpType);
      
      break;
   
    case O3ABC:
      if (indi.getNfld() != 10 )
      {
        printf("getInterpolationCellStructv: wrong size of indi.\n");
        exit(0);
      }
      searchInterpolationCellO3ABCv(coord, istart, iend, 
                                    icv, jcv, kcv, found);
      compLagrangeCoefsv(coord, istart, iend, icv, jcv, kcv, found, cf, 
                         corr, interpType, interpMeshType);

      // Indices par direction des points de la molecule d'interpolation
      if (interpMeshType == EXT_CENTERS)
        fromExtendedToStandardCentersv(
          istart, iend, icv, jcv, kcv, indi, extrap, interpType);
      else 
        compStandardIndicesv(istart, iend, icv, jcv, kcv, indi, interpType);
    
      
      // modification des coefficients d interpolation : ordre 2 O2CF 
      // pour les points mal approximes 
      // modif locale de indi
      correctInterpToO2CFv(coord, istart, iend, indi,
                           found, cf, corr, interpMeshType);
      break;

    case O5ABC :
      if (indi.getNfld() != 16)
      {
        printf("Error: getInterpolationCellStructv: wrong size of indi.\n");
        exit(0);
      }
      searchInterpolationCellO5ABCv(coord, istart, iend, 
                                    icv, jcv, kcv, found);
      compLagrangeCoefsv(coord, istart, iend, icv, jcv, kcv, found, cf, 
                         corr, interpType, interpMeshType);
      
      // Indices par direction des points de la molecule d'interpolation
      if ( interpMeshType == EXT_CENTERS)
        fromExtendedToStandardCentersv(
          istart, iend, icv, jcv, kcv, indi, extrap, interpType);
      else 
        compStandardIndicesv(istart, iend, icv, jcv, kcv, indi, interpType);
    
      
      // modification des coefficients d interpolation : ordre 2 O2CF 
      // pour les points mal approximes 
      // modif locale de indi
      correctInterpToO2CFv(coord, istart, iend, indi,
                           found, cf, corr, interpMeshType);
      break;

    default:
      printf("Error: getInterpolationCellStructv: unknown type of interpolation.\n");
      break;
  }  
}

//=============================================================================
/* Find the interp. coeff of point (x,y,z) in the given tetrahedra.
   Taken from FLU3M. */
//=============================================================================
void K_KINTERP::BlkInterpData::
coeffInterpTetra(E_Float x, E_Float y, E_Float z,
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
  if (E_abs(det) < EPS)
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

// ============= Interp/BlkInterpData.cpp ====================================
