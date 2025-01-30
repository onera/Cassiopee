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

# include "CompGeom/compGeom.h"
# include <vector>
# include <math.h>

using namespace K_FLD;
using namespace std;

extern "C"
{
  void k6boundboxofstructcell_(
    E_Int* ind, const E_Int& ni, 
    const E_Float* x, const E_Float* y, const E_Float* z,
    E_Float& xmin, E_Float& xmax, E_Float& ymin,
    E_Float& ymax, E_Float& zmin, E_Float& zmax);

  void k6boundbox_(const E_Int& im, const E_Int& jm, const E_Int& km,
                   const E_Float* x, const E_Float* y, const E_Float* z,
                   E_Float& xmax, E_Float& ymax, E_Float& zmax, 
                   E_Float& xmin, E_Float& ymin, E_Float& zmin );

  void k6boundboxunstr_(const E_Int& npts, 
                        const E_Float* x, const E_Float* y, const E_Float* z, 
                        E_Float& xmax, E_Float& ymax, E_Float& zmax, 
                        E_Float& xmin, E_Float& ymin, E_Float& zmin);
}

//===========================================================================
/* Calcul de l'aire d un triangle ABC a partir des longueurs de ses 3 cotes
   par la formule de Heron */
//===========================================================================
E_Float K_COMPGEOM::compTriangleArea(E_Float a, E_Float b, E_Float c)
{
  E_Float ps2 = (a+b+c)*0.5;
  return sqrt(ps2*(ps2-a)*(ps2-b)*(ps2-c));
}
//===========================================================================
/* Calcul de la bounding box d'une grille structuree ou non structuree. 
   k6boundboxunstr calcule la bounding box a partir de tous les pts et 
   pas les frontieres uniquement*/ 
//===========================================================================
void K_COMPGEOM::boundingBox(E_Int npts, E_Float* xt, E_Float* yt, E_Float* zt,
                             E_Float& xmin, E_Float& ymin, E_Float& zmin,
                             E_Float& xmax, E_Float& ymax, E_Float& zmax)
{
  k6boundboxunstr_(npts, xt, yt, zt, xmax, ymax, zmax, xmin, ymin, zmin);
}
//===========================================================================
/* Calcul de la bounding box d'une grille structuree ou non structuree. 
   k6boundboxunstr calcule la bounding box a partir de tous les pts et 
   pas les frontieres uniquement
*/ 
//===========================================================================
void K_COMPGEOM::boundingBox(E_Int posx, E_Int posy, E_Int posz,
                             FldArrayF& field,
                             E_Float& xmin, E_Float& ymin, E_Float& zmin,
                             E_Float& xmax, E_Float& ymax, E_Float& zmax)
{
  E_Int size = field.getSize();
  k6boundboxunstr_( 
    size, field.begin(posx), field.begin(posy), field.begin(posz),
    xmax, ymax, zmax, xmin, ymin, zmin);
}
//=============================================================================
// Calcul de la bounding box d'un array structure
//=============================================================================
void K_COMPGEOM::boundingBox(E_Int im, E_Int jm, E_Int km, 
                             E_Int posx, E_Int posy, E_Int posz,
                             FldArrayF& field,
                             E_Float& xmin, E_Float& ymin, E_Float& zmin,
                             E_Float& xmax, E_Float& ymax, E_Float& zmax)
{
  k6boundbox_(im, jm, km, 
              field.begin(posx), field.begin(posy), field.begin(posz),
              xmax, ymax, zmax, xmin, ymin, zmin );
}
//===========================================================================
/* Calcul de la bounding box d'un ensemble de grilles structurees et/ou 
   non structurees.
*/ 
//===========================================================================
void K_COMPGEOM::globalBoundingBox(
  vector<E_Int>& posxt, vector<E_Int>& posyt, vector<E_Int>& poszt,
  vector<FldArrayF*>& listOfFields,
  E_Float& xmin, E_Float& ymin, E_Float& zmin,
  E_Float& xmax, E_Float& ymax, E_Float& zmax)
{
  xmin = +K_CONST::E_MAX_FLOAT;
  ymin = +K_CONST::E_MAX_FLOAT;
  zmin = +K_CONST::E_MAX_FLOAT;
  xmax = -K_CONST::E_MAX_FLOAT;
  ymax = -K_CONST::E_MAX_FLOAT;
  zmax = -K_CONST::E_MAX_FLOAT;

  E_Int nd = listOfFields.size();
  E_Float xminl, yminl, zminl, xmaxl, ymaxl, zmaxl;

  for (E_Int i = 0; i < nd; i++)
  {
    FldArrayF* field = listOfFields[i];
    E_Int size = field->getSize();
    E_Int posx = posxt[i]; 
    E_Int posy = posyt[i];
    E_Int posz = poszt[i];

    k6boundboxunstr_( 
      size, field->begin(posx), field->begin(posy), field->begin(posz),
      xmaxl, ymaxl, zmaxl, xminl, yminl, zminl);

    xmin = K_FUNC::E_min(xmin, xminl);
    ymin = K_FUNC::E_min(ymin, yminl);
    zmin = K_FUNC::E_min(zmin, zminl);
    xmax = K_FUNC::E_max(xmax, xmaxl);
    ymax = K_FUNC::E_max(ymax, ymaxl);
    zmax = K_FUNC::E_max(zmax, zmaxl);
  }
}
//==========================================================================
/* Bounding box de toutes les cellules d'une grille structuree
   IN: im, jm, km: dimensions de l'array definissant la grille
   IN: coord: coordonnees de la grille
   OUT: bbox(ncells, 6): xmin, ymin, zmin, xmax, ymax, zmax
*/
//==========================================================================
void K_COMPGEOM::boundingBoxOfStructCells(E_Int im, E_Int jm, E_Int km, 
                                          K_FLD::FldArrayF& coord,
                                          K_FLD::FldArrayF& bbox)
{
  
  E_Float xmin, ymin, zmin, xmax, ymax, zmax;
  E_Int im1 = K_FUNC::E_max(1, im-1);
  E_Int jm1 = K_FUNC::E_max(1, jm-1);
  E_Int km1 = K_FUNC::E_max(1, km-1);
  E_Int imjm = im*jm;
  E_Int im1jm1 = im1*jm1;

  FldArrayI indtab(8);
  if (bbox.getSize() == 0) bbox.malloc(im1jm1*km1, 6);
  
  E_Int stepi = 1;
  E_Int stepj = im;
  E_Int stepk = im*jm;
  if (im == 1) stepi = 0;
  if (jm == 1) stepj = 0;
  if (km == 1) stepk = 0;
  
  E_Int indcell;
  for (E_Int k = 0; k < km1; k++)
    for (E_Int j = 0; j < jm1; j++)
      for (E_Int i = 0; i < im1; i++)
      {
        indcell = i + j * im1 + k * im1jm1;
        indtab[0] = i+j*im+k*imjm;
        indtab[1] = indtab[0] + stepi;
        indtab[2] = indtab[0] + stepj;
        indtab[3] = indtab[2] + stepi;
        indtab[4] = indtab[0] + stepk;
        indtab[5] = indtab[1] + stepk;
        indtab[6] = indtab[2] + stepk;
        indtab[7] = indtab[3] + stepk;

        k6boundboxofstructcell_( 
          indtab.begin(), coord.getSize(),
          coord.begin(1), coord.begin(2), coord.begin(3), 
          xmin, xmax, ymin, ymax, zmin, zmax);
        bbox(indcell,1) = xmin;
        bbox(indcell,2) = ymin;
        bbox(indcell,3) = zmin;
        bbox(indcell,4) = xmax;
        bbox(indcell,5) = ymax;
        bbox(indcell,6) = zmax;
      }
}
//======================================================================
/* Bounding box de toutes les cellules d'une grille non structuree
   IN: connect: connectivite de la grille
   IN: coord: coordonnees de la grille
   OUT: bbox(nelts, 6): xmin, ymin, zmin, xmax, ymax, zmax
   bbox est alloue ici. */
//======================================================================
void K_COMPGEOM::boundingBoxOfUnstrCells(K_FLD::FldArrayI& connect,
                                         E_Float* xt, E_Float* yt, E_Float* zt,
                                         K_FLD::FldArrayF& bbox)
{
  E_Int nelts = connect.getSize();
  if (bbox.getSize() == 0) bbox.malloc(nelts, 6);
  E_Float* bbox1 = bbox.begin(1);
  E_Float* bbox2 = bbox.begin(2);
  E_Float* bbox3 = bbox.begin(3);
  E_Float* bbox4 = bbox.begin(4);
  E_Float* bbox5 = bbox.begin(5);
  E_Float* bbox6 = bbox.begin(6);

#pragma omp parallel default(shared)
  {
    E_Float xmin, ymin, zmin, xmax, ymax,zmax;
    #pragma omp for
    for (E_Int et = 0; et < nelts; et++)
    {
      boundingBoxOfUnstrCell(et, connect, xt, yt, zt,
                             xmax, ymax, zmax, xmin, ymin, zmin);
      bbox1[et] = xmin; bbox2[et] = ymin; bbox3[et] = zmin;
      bbox4[et] = xmax; bbox5[et] = ymax;  bbox6[et] = zmax;
    }
  }
}
//=============================================================================
/* Find the bounding box of a cell of a structured array - issu de KMesh */
//=============================================================================
void K_COMPGEOM::boundingBoxOfCell(E_Int im, E_Int jm, E_Int km, E_Int ind,
                                   FldArrayF& coord,
                                   E_Float& xmax, E_Float& ymax, E_Float& zmax, 
                                   E_Float& xmin, E_Float& ymin, E_Float& zmin) 
{
  E_Int imjm = im*jm;
  E_Int k = ind/imjm;
  E_Int j = (ind - k * imjm) / im;
  E_Int i = ind - j * im + k * imjm;

  E_Int alpha = 1;
  E_Int beta  = 1;
  E_Int gamma = 1;          
  if (i == im-1) alpha = -1;
  if (j == jm-1) beta = -1;
  if (k == km-1) gamma = -1;
  if (im == 1) alpha = 0;
  if (jm == 1) beta = 0;
  if (km == 1) gamma = 0;
  
  E_Int indtab[8];
  indtab[0] = ind;
  indtab[1] = (i+alpha) + j*im + k*imjm;
  indtab[2] = (i+alpha) + (j+beta)*im + k*imjm;
  indtab[3] = i + (j+beta)*im + k*imjm;
  indtab[4] = i + j*im + (k+gamma)*imjm;
  indtab[5] = (i+alpha) + j*im + (k+gamma)*imjm;
  indtab[6] = (i+alpha) + (j+beta)*im + (k+gamma)*imjm;  
  indtab[7] = i + (j+beta)*im + (k+gamma)*imjm;
  
  k6boundboxofstructcell_(indtab, coord.getSize(), 
                          coord.begin(1), coord.begin(2), coord.begin(3), 
                          xmin, xmax, ymin, ymax, zmin, zmax);
}

//=============================================================================
/* Calcul de la bounding box d'une cellule non structuree */
//=============================================================================
void K_COMPGEOM::boundingBoxOfUnstrCell(
  E_Int noet, FldArrayI& connect, 
  E_Float* xt, E_Float* yt, E_Float* zt,
  E_Float& xmax, E_Float& ymax, E_Float& zmax, 
  E_Float& xmin, E_Float& ymin, E_Float& zmin) 
{

  xmin = K_CONST::E_MAX_FLOAT;
  ymin = K_CONST::E_MAX_FLOAT;
  zmin = K_CONST::E_MAX_FLOAT;
  xmax = -K_CONST::E_MAX_FLOAT;
  ymax = -K_CONST::E_MAX_FLOAT;
  zmax = -K_CONST::E_MAX_FLOAT;

  E_Int nvert = connect.getNfld();
  
  for (E_Int vert = 1; vert <= nvert; vert++)
  {
    E_Int ind = connect(noet,vert)-1;
    xmin = K_FUNC::E_min(xmin,xt[ind]);
    ymin = K_FUNC::E_min(ymin,yt[ind]);
    zmin = K_FUNC::E_min(zmin,zt[ind]);
    xmax = K_FUNC::E_max(xmax,xt[ind]);
    ymax = K_FUNC::E_max(ymax,yt[ind]);
    zmax = K_FUNC::E_max(zmax,zt[ind]);
  }
}
//=============================================================================
// Intersection de bbox de 2 grilles structurees
// retourne 1 si les bounding boxes s intersectent, 0 sinon
// similaire a BlkOverlapData::testBBIntersection
//=============================================================================
E_Int 
K_COMPGEOM::compBoundingBoxIntersection(E_Int ni1, E_Int nj1, E_Int nk1, 
                                        E_Int posx1, E_Int posy1, 
                                        E_Int posz1, FldArrayF& f1, 
                                        E_Int ni2, E_Int nj2, E_Int nk2, 
                                        E_Int posx2, E_Int posy2, 
                                        E_Int posz2, FldArrayF& f2, 
                                        E_Float& xmin1, E_Float& xmax1, 
                                        E_Float& ymin1, E_Float& ymax1, 
                                        E_Float& zmin1, E_Float& zmax1, 
                                        E_Float& xmin2, E_Float& xmax2, 
                                        E_Float& ymin2, E_Float& ymax2, 
                                        E_Float& zmin2, E_Float& zmax2,
                                        E_Float tol)
{
  // bbox1 ds repere absolu
  k6boundbox_(ni1, nj1, nk1, 
              f1.begin(posx1), f1.begin(posy1), f1.begin(posz1), 
              xmax1, ymax1, zmax1, xmin1, ymin1, zmin1 );
  
  //bbox2 ds repere absolu
  k6boundbox_(ni2, nj2, nk2, 
              f2.begin(posx2), f2.begin(posy2), f2.begin(posz2), 
              xmax2, ymax2, zmax2, xmin2, ymin2, zmin2 );

  if ( xmin1  <=  xmax2+tol && xmax1  >=  xmin2-tol &&
       ymin1  <=  ymax2+tol && ymax1  >=  ymin2-tol &&
       zmin1  <=  zmax2+tol && zmax1  >=  zmin2-tol )
    return 1;
  else return 0;
}
//=============================================================================
// Intersection de 2 bbox donnees
// retourne 1 si les bounding boxes s'intersectent, 0 sinon
//=============================================================================
E_Int 
K_COMPGEOM::compBoundingBoxIntersection(E_Float xmin1, E_Float xmax1, 
                                        E_Float ymin1, E_Float ymax1, 
                                        E_Float zmin1, E_Float zmax1, 
                                        E_Float xmin2, E_Float xmax2, 
                                        E_Float ymin2, E_Float ymax2, 
                                        E_Float zmin2, E_Float zmax2,
                                        E_Float tol)
{
  if (xmin1 > xmax2+tol) return 0;
  if (xmax1 < xmin2-tol) return 0;
  if (ymin1 > ymax2+tol) return 0;
  if (ymax1 < ymin2-tol) return 0;
  if (zmin1 > zmax2+tol) return 0;
  if (zmax1 < zmin2-tol) return 0;
  return 1;
}
//=============================================================================
/* Recherche si les CEBB  de 2 grilles s intersectent */
//=============================================================================
E_Int 
K_COMPGEOM::compCEBBIntersection(E_Int ni1, E_Int nj1, E_Int nk1, 
                                 E_Int posx1, E_Int posy1, E_Int posz1,
                                 FldArrayF& f1, 
                                 E_Int ni2, E_Int nj2, E_Int nk2,
                                 E_Int posx2, E_Int posy2, E_Int posz2,
                                 FldArrayF& f2, E_Float tol)
{
  // 1- test des bbox globales
  E_Float xmin1, xmax1, ymin1, ymax1, zmin1, zmax1;
  E_Float xmin2, xmax2, ymin2, ymax2, zmin2, zmax2;

  E_Int isIntersect = compBoundingBoxIntersection(
    ni1, nj1, nk1, posx1, posy1, posz1, f1,
    ni2, nj2, nk2, posx2, posy2, posz2, f2, 
    xmin1, xmax1, ymin1, ymax1, zmin1, zmax1,
    xmin2, xmax2, ymin2, ymax2, zmin2, zmax2, tol);

  if ( isIntersect == 0 ) return 0;

  // 2 - test des CEBB 2 a 2
  FldArrayF cartMin1;
  FldArrayF cartMax1;
  FldArrayF cartMin2;
  FldArrayF cartMax2;

  // creation des tableaux de coordonnees
  E_Int npts1 = ni1*nj1*nk1;
  FldArrayF coord1(npts1, 3);
  E_Float* c1x = coord1.begin(1);
  E_Float* c1y = coord1.begin(2);
  E_Float* c1z = coord1.begin(3);
  E_Float* f1x = f1.begin(posx1);
  E_Float* f1y = f1.begin(posy1);
  E_Float* f1z = f1.begin(posz1);
  for (E_Int i1 = 0; i1 < npts1; i1++)
  {
    c1x[i1] = f1x[i1];
    c1y[i1] = f1y[i1];
    c1z[i1] = f1z[i1];
  }
  
  E_Int npts2 = ni2*nj2*nk2;
  FldArrayF coord2(npts2, 3);
  E_Float* c2x = coord2.begin(1);
  E_Float* c2y = coord2.begin(2);
  E_Float* c2z = coord2.begin(3);
  E_Float* f2x = f2.begin(posx2);
  E_Float* f2y = f2.begin(posy2);
  E_Float* f2z = f2.begin(posz2);
  for (E_Int i2 = 0; i2 < npts2; i2++)
  {
    c2x[i2] = f2x[i2];
    c2y[i2] = f2y[i2];
    c2z[i2] = f2z[i2];
  }
  
  FldArrayF cartEltArray1;
  FldArrayF cartEltArray2;
  E_Boolean isDegenerated = true;

  for (E_Int dir1 = 1; dir1 <= 3; dir1++)
  {
    short isok1 = compCartEltsArray(
      dir1, ni1, nj1, nk1, xmin1, ymin1, zmin1, 
      xmax1, ymax1, zmax1, coord1, cartEltArray1);
    
    if ( isok1 == 1 ) // direction possible pour les elts cart 
    { 
      for (E_Int dir2 = 1; dir2 <= 3; dir2++)
      {
        short isok2 = compCartEltsArray(
          dir2, ni2, nj2, nk2, xmin2, ymin2, zmin2, 
          xmax2, ymax2, zmax2, coord2, cartEltArray2);

        if (isok2 == 1)
        {
          isDegenerated = false;
          // test des intersection des cebboxes
          for (E_Int elt1 = 0; elt1 < cartEltArray1.getSize(); elt1++)
          {
            E_Float xmin11 = cartEltArray1(elt1,1);
            E_Float ymin11 = cartEltArray1(elt1,2);
            E_Float zmin11 = cartEltArray1(elt1,3);
            E_Float xmax11 = cartEltArray1(elt1,4);
            E_Float ymax11 = cartEltArray1(elt1,5);
            E_Float zmax11 = cartEltArray1(elt1,6);
            
            for (E_Int elt2 = 0; elt2 < cartEltArray2.getSize(); elt2++)
            {
              E_Float xmin21 = cartEltArray2(elt2,1);
              E_Float ymin21 = cartEltArray2(elt2,2);
              E_Float zmin21 = cartEltArray2(elt2,3);
              E_Float xmax21 = cartEltArray2(elt2,4);
              E_Float ymax21 = cartEltArray2(elt2,5);
              E_Float zmax21 = cartEltArray2(elt2,6);
            
              // test intersection
              if ( xmin11 <=  xmax21+tol && xmax11 >=  xmin21-tol &&
                   ymin11 <=  ymax21+tol && ymax11 >=  ymin21-tol &&
                   zmin11 <=  zmax21+tol && zmax11 >=  zmin21-tol )
                return 1;
            }
          }
        }
      }
    } 
  }
  // Dans le cas ou on ne peut pas calculer la CEBB, on renvoit
  // le resultat de l'intersection des BB
  if (isDegenerated == true) return 1;
  return 0;
}
//============================ CompGeom/compGeom.cpp =======================
