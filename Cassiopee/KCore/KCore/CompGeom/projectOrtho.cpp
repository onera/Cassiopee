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

# include "CompGeom/compGeom.h"
# include "Nuga/include/KdTree.h"
# include "Nuga/include/BbTree.h"
# include "Nuga/include/ArrayAccessor.h"
# include <stdio.h>

using namespace K_FLD;
using namespace std;

//=============================================================================
/* Projete un point (x,y,z) sur un surface array (TRI ou BAR) orthogonalement
   ou projete sur le sommet du plus proche triangle
   IN: x,y,z: coord du point a projeter
   IN: fx2, fy2, fz2, cn2: coord et connectivite de la surface
   OUT: xo, yo, zo: les coord. du point projete.
   Retourne le no du triangle sur lequel s'effectue la projection
   Retourne -1 si impossible de projeter */
//=============================================================================
E_Int K_COMPGEOM::projectOrtho(E_Float x, E_Float y, E_Float z,
                               E_Float* fx2, E_Float* fy2, E_Float* fz2,
                               K_FLD::FldArrayI& cn2, 
                               E_Float& xo, E_Float& yo, E_Float& zo,
                               E_Float* p0, E_Float* p1, E_Float* p2, E_Float* p)
{
  E_Int noet = -1;
  E_Int ret; 
  p[0] = x; p[1] = y; p[2] = z;
  
  E_Float distc = 1.e6;
  E_Float dist2; E_Float sigma0, sigma1;
  E_Float xp, yp, zp;
  E_Int ind1, ind2, ind3;
  E_Bool in;
  E_Int* cn2p1 = cn2.begin(1);
  E_Int* cn2p2 = cn2.begin(2);
  E_Int nvert = cn2.getNfld();
  xo = x; yo = y; zo = z;

  if (nvert == 2) // bar
  {
    for (E_Int e = 0; e < cn2.getSize(); e++)
    {
      ind1 = cn2p1[e]-1; ind2 = cn2p2[e]-1;     
      p0[0] = fx2[ind1]; p0[1] = fy2[ind1]; p0[2] = fz2[ind1];
      p1[0] = fx2[ind2]; p1[1] = fy2[ind2]; p1[2] = fz2[ind2];
      ret = K_COMPGEOM::distanceToBar(p0, p1, p, 1, xp, yp, zp, in, dist2);
      
      if (ret == 0)
      {
        if (dist2 < distc) 
        {xo = xp; yo = yp; zo = zp; distc = dist2; noet = e;}
      }  
    }
  }
  else if (nvert == 3) // tri
  {
    E_Int* cn2p3 = cn2.begin(3);
    for (E_Int e = 0; e < cn2.getSize(); e++)
    {
      ind1 = cn2p1[e]-1; ind2 = cn2p2[e]-1;  ind3 = cn2p3[e]-1;
      p0[0] = fx2[ind1]; p0[1] = fy2[ind1]; p0[2] = fz2[ind1];
      p1[0] = fx2[ind2]; p1[1] = fy2[ind2]; p1[2] = fz2[ind2];
      p2[0] = fx2[ind3]; p2[1] = fy2[ind3]; p2[2] = fz2[ind3];        
      ret = K_COMPGEOM::distanceToTriangle(p0, p1, p2, p, 2, 
                                           dist2, in, 
                                           xp, yp, zp,
                                           sigma0, sigma1);
      if (ret == 0)
      {
        //if (in == false) dist2 = dist2 * 1.1;
        if (dist2 < distc)
        { xo = xp; yo = yp; zo = zp; distc = dist2; noet = e; }
      }
    }
  }
  else 
  {
    printf("Warning: projectOrtho: only valid for BAR, TRI elements.\n");
    return -1;
  }
  return noet;
}

//=============================================================================
/* Projete un point (x,y,z) sur un surface array (TRI ou BAR) orthogonalement
   ou projete sur le sommet du plus proche triangle
   IN: x,y,z: coord du point a projeter
   IN: fx2, fy2, fz2, cn2, coord et connect de la surface
   IN: indices est la liste des triangles a tester
   OUT: xo, yo, zo: les coord. du point projete.
   Retourne le no du triangle sur lequel s'effectue la projection
   Retourne -1 si impossible de projeter */
//=============================================================================
E_Int K_COMPGEOM::projectOrthoPrecond(
  E_Float x, E_Float y, E_Float z,
  E_Float* fx2, E_Float* fy2, E_Float* fz2,
  std::vector<E_Int>& indices, K_FLD::FldArrayI& cn2, 
  E_Float& xo, E_Float& yo, E_Float& zo,
  E_Float* p0, E_Float* p1, E_Float* p2, E_Float* p)
{
  E_Int noet = -1;
  E_Int ret;
  p[0] = x; p[1] = y; p[2] = z;
  
  E_Float distc = 1.e6;
  E_Float dist2; E_Float sigma0, sigma1;
  E_Float xp, yp, zp;
  E_Int ind1, ind2, ind3;
  E_Bool in;
  E_Int* cn2p1 = cn2.begin(1);
  E_Int* cn2p2 = cn2.begin(2);
  E_Int nvert = cn2.getNfld();
  xo = x; yo = y; zo = z;
  E_Int nbb = indices.size();
  E_Int et = 0;
  if (nvert == 2) // bar
  {
    for (E_Int noe = 0; noe < nbb; noe++)
    {
      et = indices[noe];
      ind1 = cn2p1[et]-1; ind2 = cn2p2[et]-1;     
      p0[0] = fx2[ind1]; p0[1] = fy2[ind1]; p0[2] = fz2[ind1];
      p1[0] = fx2[ind2]; p1[1] = fy2[ind2]; p1[2] = fz2[ind2];
      ret = K_COMPGEOM::distanceToBar(p0, p1, p, 1, xp, yp, zp, in, dist2);
      
      if (ret == 0)
      {
        if (dist2 < distc) 
        {xo = xp; yo = yp; zo = zp; distc = dist2; noet = et;}
      }  
    }
  }
  else if (nvert == 3) // tri
  {
    E_Int* cn2p3 = cn2.begin(3);
    for (E_Int noe = 0; noe < nbb; noe++)
    {
      et = indices[noe];
      ind1 = cn2p1[et]-1; ind2 = cn2p2[et]-1; ind3 = cn2p3[et]-1;
      p0[0] = fx2[ind1]; p0[1] = fy2[ind1]; p0[2] = fz2[ind1];
      p1[0] = fx2[ind2]; p1[1] = fy2[ind2]; p1[2] = fz2[ind2];
      p2[0] = fx2[ind3]; p2[1] = fy2[ind3]; p2[2] = fz2[ind3];        
      ret = K_COMPGEOM::distanceToTriangle(p0, p1, p2, p, 2, 
                                           dist2, in, 
                                           xp, yp, zp,
                                           sigma0, sigma1);
      if (ret == 0)
      {
        if (dist2 < distc) 
        {xo = xp; yo = yp; zo = zp; distc = dist2; noet = et;}
      }
    }
  }
  else 
  {
    printf("Warning: projectOrthoPrecond: only valid for BAR, TRI elements.\n");
    return -1;
  }
  return noet;
}

//=============================================================================
// Algorithme de projection sans preconditionnement
// IN: npts: nombre de pts a projeter
// IN: fx, fy, fz: coord des points a projeter
// IN: fx2, fy2, fz2, cn2: coord et connect de la surface
// OUT: fx, fy, fz modifie (contenant les points projetes)
//=============================================================================
void K_COMPGEOM::projectOrthoWithoutPrecond(
  E_Int npts, FldArrayI& cn2, 
  E_Float* fx2, E_Float* fy2, E_Float* fz2,
  E_Float* fx, E_Float* fy, E_Float* fz)
{
  #pragma omp parallel
  {
    E_Int ret = 0;
    E_Float xo, yo, zo;
    E_Float p0[3]; E_Float p1[3]; E_Float p2[3]; E_Float p[3];

    #pragma omp for
    for (E_Int ind = 0; ind < npts; ind++)
    {
      ret = projectOrtho(fx[ind], fy[ind], fz[ind], 
                         fx2, fy2, fz2, cn2, xo, yo, zo, p0, p1, p2, p);
      if (ret != -1) {fx[ind] = xo; fy[ind] = yo; fz[ind] = zo;}
    }
  }
}
//=============================================================================
/* Algorithme de projection avec preconditionnement.
   IN: fx2, fy2, fz2, cn2: coord et connect de la surface de projection.
   IN/OUT: fields: vecteur des surfaces que l'on projete sur fx2.
*/
//=============================================================================
void K_COMPGEOM::projectOrthoWithPrecond(
  E_Int posx2, E_Int posy2, E_Int posz2, 
  FldArrayI& cn2, FldArrayF& f2, vector<E_Int>& sizet,
  vector<E_Float*>& fxt, vector<E_Float*>& fyt, vector<E_Float*>& fzt) 
{
  E_Int nelts2 = cn2.getSize();
  E_Float* fx2 = f2.begin(posx2);
  E_Float* fy2 = f2.begin(posy2);
  E_Float* fz2 = f2.begin(posz2);
  // Creation du kdtree
  K_FLD::ArrayAccessor<FldArrayF> coordAcc(f2, posx2, posy2, posz2);
  K_SEARCH::KdTree<FldArrayF> kdt(coordAcc);
  // Creation du bboxtree
  typedef K_SEARCH::BoundingBox<3>  BBox3DType; 
  vector<BBox3DType*> boxes(nelts2);// liste des bbox de ts les elements de a2
  K_FLD::FldArrayF bbox(nelts2, 6);// xmin, ymin, zmin, xmax, ymax, zmax
  K_COMPGEOM::boundingBoxOfUnstrCells(cn2, fx2, fy2, fz2, bbox);

  E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(4);
  E_Float* yminp = bbox.begin(2); E_Float* ymaxp = bbox.begin(5);
  E_Float* zminp = bbox.begin(3); E_Float* zmaxp = bbox.begin(6);
  
  #pragma omp parallel
  {
    E_Float minB[3]; E_Float maxB[3];
    #pragma omp for
    for (E_Int et = 0; et < nelts2; et++)
    {
      minB[0] = xminp[et]; minB[1] = yminp[et]; minB[2] = zminp[et];
      maxB[0] = xmaxp[et]; maxB[1] = ymaxp[et]; maxB[2] = zmaxp[et]; 
      boxes[et] = new BBox3DType(minB, maxB);
    }
  }
    
  // Build the box tree.
  K_SEARCH::BbTree3D bbtree(boxes);
  
  // projection des pts de f sur f2
  E_Int nzones = fxt.size();

  #pragma omp parallel
  {
    E_Float xo, yo, zo; E_Float pt[3];
    E_Float p0[3]; E_Float p1[3]; E_Float p2[3]; E_Float p[3];
    E_Int ret=0; E_Int indp=0; 
    E_Float rx, ry, rz, rad;
    E_Float minB[3]; E_Float maxB[3];
    vector<E_Int> indicesBB; // liste des indices des facettes intersectant la bbox
  
    for (E_Int v = 0; v < nzones; v++)
    {
      E_Int npts = sizet[v];
      E_Float* fx = fxt[v];
      E_Float* fy = fyt[v];
      E_Float* fz = fzt[v];
    
      #pragma omp for
      for (E_Int ind = 0; ind < npts; ind++)
      {
        // recherche du pt le plus proche P' de P
        pt[0] = fx[ind]; pt[1] = fy[ind]; pt[2] = fz[ind];
        indp = kdt.getClosest(pt);
      
        // calcul de la bounding box de la sphere de rayon PP'
        rx = pt[0]-fx2[indp]; ry = pt[1]-fy2[indp]; rz = pt[2]-fz2[indp];
        rad = sqrt(rx*rx+ry*ry+rz*rz);
        minB[0] = pt[0]-rad; minB[1] = pt[1]-rad; minB[2] = pt[2]-rad;
        maxB[0] = pt[0]+rad; maxB[1] = pt[1]+rad; maxB[2] = pt[2]+rad;
        bbtree.getOverlappingBoxes(minB, maxB, indicesBB);
      
        ret = projectOrthoPrecond(fx[ind], fy[ind], fz[ind], 
                                  fx2, fy2, fz2, indicesBB, cn2, xo, yo, zo,
                                  p0, p1, p2, p);
        if (ret != -1) {fx[ind] = xo; fy[ind] = yo; fz[ind] = zo;}
        indicesBB.clear();
      }
    }
  }

  // delete boxes
  E_Int size = boxes.size();
  for (E_Int v = 0; v < size; v++) delete boxes[v];    
}

//=============================================================================
// Algorithme de projection avec preconditionnement
// IN: npts: nombre de pts a projeter
// IN: fx, fy, fz: coord des points a projeter
// IN: fx2, fy2, fz2, cn2: coord et connect de la surface de projection.
// OUT: fx, fy, fz modifie (contenant les points projetes)
//=============================================================================
void K_COMPGEOM::projectOrthoWithPrecond(
  E_Int npts, FldArrayI& cn2, 
  E_Int posx2, E_Int posy2, E_Int posz2,
  FldArrayF& f2, E_Float* fx, E_Float* fy, E_Float* fz)
{
  E_Int nelts2 = cn2.getSize();
  E_Float* fx2 = f2.begin(posx2);
  E_Float* fy2 = f2.begin(posy2);
  E_Float* fz2 = f2.begin(posz2);
 
  // Creation du kdtree
  K_FLD::ArrayAccessor<FldArrayF> coordAcc(f2, posx2, posy2, posz2);
  K_SEARCH::KdTree<FldArrayF> kdt(coordAcc);

  // Creation du bboxtree
  typedef K_SEARCH::BoundingBox<3>  BBox3DType; 
  vector<BBox3DType*> boxes(nelts2);// liste des bbox de ts les elements de a2
  K_FLD::FldArrayF bbox(nelts2,6);// xmin, ymin, zmin, xmax, ymax, zmax
  K_COMPGEOM::boundingBoxOfUnstrCells(cn2, fx2, fy2, fz2, bbox);

  E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(4);
  E_Float* yminp = bbox.begin(2); E_Float* ymaxp = bbox.begin(5);
  E_Float* zminp = bbox.begin(3); E_Float* zmaxp = bbox.begin(6);
  
  #pragma omp parallel
  {
    E_Float minB[3];  E_Float maxB[3];
    #pragma omp for
    for (E_Int et = 0; et < nelts2; et++)
    {
      minB[0] = xminp[et]; minB[1] = yminp[et]; minB[2] = zminp[et];
      maxB[0] = xmaxp[et]; maxB[1] = ymaxp[et]; maxB[2] = zmaxp[et]; 
      boxes[et] = new BBox3DType(minB, maxB);
    }
  }
  // Build the box tree.
  K_SEARCH::BbTree3D bbtree(boxes);
  
  // projection des pts de f sur f2
  #pragma omp parallel
  {
    E_Float xo, yo, zo; E_Float pt[3];
    E_Float p0[3]; E_Float p1[3]; E_Float p2[3]; E_Float p[3];
    E_Float minB[3];  E_Float maxB[3];
    E_Int ret=0; E_Int indp=0; 
    E_Float rx, ry, rz, rad;
    vector<E_Int> indicesBB; // liste des indices des facettes intersectant la bbox

    #pragma omp for
    for (E_Int ind = 0; ind < npts; ind++)
    {
      // recherche du pt le plus proche P' de P
      pt[0] = fx[ind]; pt[1] = fy[ind]; pt[2] = fz[ind];
      indp = kdt.getClosest(pt);

      // calcul de la bounding box de la sphere de rayon PP'
      rx = pt[0]-fx2[indp]; ry = pt[1]-fy2[indp]; rz = pt[2]-fz2[indp];
      rad = sqrt(rx*rx+ry*ry+rz*rz);
      minB[0] = pt[0]-rad; minB[1] = pt[1]-rad; minB[2] = pt[2]-rad;
      maxB[0] = pt[0]+rad; maxB[1] = pt[1]+rad; maxB[2] = pt[2]+rad;
      bbtree.getOverlappingBoxes(minB, maxB, indicesBB);

      ret = projectOrthoPrecond(fx[ind], fy[ind], fz[ind], 
                                fx2, fy2, fz2, indicesBB, cn2, xo, yo, zo,
                                p0, p1, p2, p);
      if (ret != -1) {fx[ind] = xo; fy[ind] = yo; fz[ind] = zo;}
      indicesBB.clear();
    }
  }
  // delete boxes
  E_Int size = boxes.size();
  for (E_Int v = 0; v < size; v++) delete boxes[v];
}
