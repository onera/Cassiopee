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
# include <stdio.h>
# include "Nuga/include/BbTree.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* Projette un point sur un surface array2 (TRI ou BAR) suivant une 
   direction donnee.
   Retourne le no du triangle sur lequel s'effectue la projection 
   -1 si impossible de projeter */
//=============================================================================
E_Int K_COMPGEOM::projectDir(E_Float x, E_Float y, E_Float z,
                             E_Float nx, E_Float ny, E_Float nz,
                             E_Float* fx2, E_Float* fy2, E_Float* fz2,
                             FldArrayI& cn2, 
                             E_Float& xo, E_Float& yo, E_Float& zo,
                             E_Int oriented)
{
  E_Int precond = 1;
  if (precond == 0) 
    return projectOneDirWithoutPrecond(x, y, z, nx, ny, nz, 
                                       fx2, fy2, fz2, cn2, xo ,yo, zo, oriented);
  else
    return projectOneDirWithPrecond(x, y, z, nx, ny, nz, 
                                    fx2, fy2, fz2, cn2, xo ,yo, zo, oriented);
}
//=============================================================================
/* Projette un point sur un surface array2 (TRI ou BAR) suivant une 
   direction donnee a partir d'une liste de triangles indices.
   Retourne le no du triangle sur lequel s'effectue la projection 
   -1 si impossible de projeter */
//=============================================================================
E_Int K_COMPGEOM::projectDir(E_Float x, E_Float y, E_Float z,
                             E_Float nx, E_Float ny, E_Float nz,
                             E_Float* fx2, E_Float* fy2, E_Float* fz2,
                             std::vector<E_Int>& indices, FldArrayI& cn2, 
                             E_Float& xo, E_Float& yo, E_Float& zo, E_Int oriented)
{ 
  E_Int noet = -1;
  E_Int ret;
  E_Float p0[3]; E_Float p1[3]; E_Float p2[3]; E_Float p[3];
  E_Float pr1[3]; E_Float pr2[3]; E_Float pi[3];
  p[0] = x; p[1] = y; p[2] = z;
  pr1[0] = p[0]; pr1[1] = p[1]; pr1[2] = p[2];
  pr2[0] = p[0] + nx; pr2[1] = p[1] + ny; pr2[2] = p[2] + nz;

  E_Float distc = 1.e6;
  E_Float dist;
  E_Int ind1, ind2, ind3;
  E_Int nvert = cn2.getNfld();
  xo = x; yo = y; zo = z; 

  E_Int nbboxes = indices.size();
  E_Int et;
  E_Float dx, dy, dz, normp, ps;
  
  if (nvert == 3) // tri
  {
    E_Int* cn2p1 = cn2.begin(1);
    E_Int* cn2p2 = cn2.begin(2);
    E_Int* cn2p3 = cn2.begin(3);
    for (E_Int noe = 0; noe < nbboxes; noe++)   
    {
      et = indices[noe];
      ind1 = cn2p1[et]-1; ind2 = cn2p2[et]-1; ind3 = cn2p3[et]-1;
      p0[0] = fx2[ind1]; p0[1] = fy2[ind1]; p0[2] = fz2[ind1];
      p1[0] = fx2[ind2]; p1[1] = fy2[ind2]; p1[2] = fz2[ind2];
      p2[0] = fx2[ind3]; p2[1] = fy2[ind3]; p2[2] = fz2[ind3];        
      ret = K_COMPGEOM::intersectRayTriangle(p0, p1, p2, 
                                             pr1, pr2, pi);
      if (ret == 1)
      {
        dx = pi[0]-p[0]; dy = pi[1]-p[1]; dz = pi[2]-p[2];
        dist = dx*dx + dy*dy + dz*dz;
        if (oriented != 0)
        {
          normp = sqrt(nx*nx+ny*ny+nz*nz);//normale
          ps = nx*dx+ny*dy+nz*dz/(dist*normp);
          if ( ps > 0. ) 
          {
            if (dist < distc) 
            {xo = pi[0]; yo = pi[1]; zo = pi[2]; distc = dist; noet = et;}
          }
        }
        else 
        {
          if (dist < distc) 
          {xo = pi[0]; yo = pi[1]; zo = pi[2]; distc = dist; noet = et;}
        }
      }
    }
  }
  else 
  {
    printf("Warning: projectDirPrecond: only valid for TRI elements.\n");
    return -1;
  }
  return noet;
}
//=============================================================================
E_Int K_COMPGEOM::projectOneDirWithPrecond(
  E_Float x, E_Float y, E_Float z,
  E_Float nx, E_Float ny, E_Float nz,
  E_Float* fx2, E_Float* fy2, E_Float* fz2,
  FldArrayI& cn2, 
  E_Float& xo, E_Float& yo, E_Float& zo, E_Int oriented)
{
  E_Float tol = 1.e-6;
  typedef K_SEARCH::BoundingBox<3>  BBox3DType;

  E_Int nelts2 = cn2.getSize();
  // Creation de la bboxtree
  vector<BBox3DType*> boxes(nelts2);// liste des bbox de ts les elements de a2
  K_FLD::FldArrayF bbox(nelts2,6);// xmin, ymin, zmin, xmax, ymax, zmax
  K_COMPGEOM::boundingBoxOfUnstrCells(cn2, fx2, fy2, fz2, bbox);
  E_Float minB[3];  E_Float maxB[3];
  E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(4);
  E_Float* yminp = bbox.begin(2); E_Float* ymaxp = bbox.begin(5);
  E_Float* zminp = bbox.begin(3); E_Float* zmaxp = bbox.begin(6);
  for (E_Int et = 0; et < nelts2; et++)
  {
    minB[0] = xminp[et]; minB[1] = yminp[et]; minB[2] = zminp[et];
    maxB[0] = xmaxp[et]; maxB[1] = ymaxp[et]; maxB[2] = zmaxp[et]; 
    boxes[et] = new BBox3DType(minB, maxB);
  }
  // Build the box tree.
  K_SEARCH::BbTree3D bbtree(boxes);
  
  E_Int noet = -1;
  E_Float p0[3]; E_Float p1[3]; E_Float p2[3];
  E_Int ret; 
  E_Float p[3]; E_Float pr1[3]; E_Float pr2[3]; E_Float pi[3];
  p[0] = x; p[1] = y; p[2] = z;
  pr1[0] = p[0]; pr1[1] = p[1]; pr1[2] = p[2];
  pr2[0] = p[0] + nx; pr2[1] = p[1] + ny; pr2[2] = p[2] + nz;

  E_Float distc = 1.e6;
  E_Float dist;
  E_Int ind1, ind2, ind3;
  E_Int nvert = cn2.getNfld();
  xo = x; yo = y; zo = z; 

  vector<E_Int> indicesBB;
  bbtree.getIntersectingBoxes(pr1, pr2, indicesBB, tol);
  E_Int nbboxes = indicesBB.size();
  E_Int et;
  E_Float dx, dy, dz, normp, ps;

  if (nvert == 3) // tri
  {
    E_Int* cn2p1 = cn2.begin(1);
    E_Int* cn2p2 = cn2.begin(2);
    E_Int* cn2p3 = cn2.begin(3);
    for (E_Int noe = 0; noe < nbboxes; noe++)   
    {
      et = indicesBB[noe];
      ind1 = cn2p1[et]-1; ind2 = cn2p2[et]-1;  ind3 = cn2p3[et]-1;
      p0[0] = fx2[ind1]; p0[1] = fy2[ind1]; p0[2] = fz2[ind1];
      p1[0] = fx2[ind2]; p1[1] = fy2[ind2]; p1[2] = fz2[ind2];
      p2[0] = fx2[ind3]; p2[1] = fy2[ind3]; p2[2] = fz2[ind3];        
      ret = K_COMPGEOM::intersectRayTriangle(p0, p1, p2, 
                                             pr1, pr2, pi);
      if (ret == 1)
      {
        dx = pi[0]-p[0]; dy = pi[1]-p[1]; dz = pi[2]-p[2];
        dist = dx*dx + dy*dy + dz*dz;
        if (oriented != 0)
        {
          normp = sqrt(nx*nx+ny*ny+nz*nz);//normale
          ps = nx*dx+ny*dy+nz*dz/(dist*normp);
          if (ps > 0.) 
          {
            if (dist < distc) 
            {xo = pi[0]; yo = pi[1]; zo = pi[2]; distc = dist; noet = et;}
          }
        }
        else 
        {
          if (dist < distc) 
          {xo = pi[0]; yo = pi[1]; zo = pi[2]; distc = dist; noet = et;}
        }
      }
    }
  }
  else 
  {
    E_Int boxesSize = boxes.size();
    for (E_Int v = 0; v < boxesSize; v++) delete boxes[v];
    printf("Warning: projectDir: only valid for TRI elements.\n");
    return -1;
  }
  E_Int boxesSize = boxes.size();
  for (E_Int v = 0; v < boxesSize; v++) delete boxes[v];
  return noet;
}

//=============================================================================
E_Int K_COMPGEOM::projectOneDirWithoutPrecond(
  E_Float x, E_Float y, E_Float z,
  E_Float nx, E_Float ny, E_Float nz,
  E_Float* fx2, E_Float* fy2, E_Float* fz2,
  FldArrayI& cn2, 
  E_Float& xo, E_Float& yo, E_Float& zo, E_Int oriented)    
{
  E_Int noet = -1;
  E_Float p0[3]; E_Float p1[3]; E_Float p2[3];
  E_Int ret; 
  E_Float p[3]; E_Float pr1[3]; E_Float pr2[3]; E_Float pi[3];
  p[0] = x; p[1] = y; p[2] = z;
  pr1[0] = p[0]; pr1[1] = p[1]; pr1[2] = p[2];
  pr2[0] = p[0] + nx; pr2[1] = p[1] + ny; pr2[2] = p[2] + nz;
  
  E_Float distc = 1e6;
  E_Float dist, dx, dy, dz, normp, ps;
  E_Int ind1, ind2, ind3;
  E_Int nvert = cn2.getNfld();
  xo = x; yo = y; zo = z; 

  if (nvert == 3) // tri
  {
    E_Int* cn2p1 = cn2.begin(1);
    E_Int* cn2p2 = cn2.begin(2);
    E_Int* cn2p3 = cn2.begin(3);
    for (E_Int e = 0; e < cn2.getSize(); e++)
    {
      ind1 = cn2p1[e]-1; ind2 = cn2p2[e]-1;  ind3 = cn2p3[e]-1;
      p0[0] = fx2[ind1]; p0[1] = fy2[ind1]; p0[2] = fz2[ind1];
      p1[0] = fx2[ind2]; p1[1] = fy2[ind2]; p1[2] = fz2[ind2];
      p2[0] = fx2[ind3]; p2[1] = fy2[ind3]; p2[2] = fz2[ind3];        
      ret = K_COMPGEOM::intersectRayTriangle(p0, p1, p2, 
                                             pr1, pr2, pi);
      if (ret == 1)
      {
        dx = pi[0]-p[0]; dy = pi[1]-p[1]; dz = pi[2]-p[2];
        dist = dx*dx + dy*dy + dz*dz;
        if ( oriented != 0 )
        {
          normp = sqrt(nx*nx+ny*ny+nz*nz);//normale
          ps = nx*dx+ny*dy+nz*dz/(dist*normp);
          if ( ps > 0. ) 
          {
            if (dist < distc) 
            {xo = pi[0]; yo = pi[1]; zo = pi[2]; distc = dist; noet = e;}
          }
        }
        else
        { 
          if (dist < distc) 
          {xo = pi[0]; yo = pi[1]; zo = pi[2]; distc = dist; noet = e;}
        }
      }
    }
  }
  else 
  {
    printf("Warning: projectDir: only valid for TRI elements.\n");
    return -1;
  }
  return noet;
}

//=============================================================================
/* projection avec precond par bboxtrees */
//=============================================================================
void K_COMPGEOM::projectDirWithPrecond(
  E_Float nx, E_Float ny, E_Float nz,
  E_Int npts, E_Int nelts2, FldArrayI& cn2,
  E_Float* fx2, E_Float* fy2, E_Float* fz2,
  E_Float* fx, E_Float* fy, E_Float* fz, E_Int oriented)
{
  E_Float tol = 1.e-6;
  typedef K_SEARCH::BoundingBox<3> BBox3DType;

  E_Int* cn2p1 = cn2.begin(1); E_Int* cn2p2 = cn2.begin(2); E_Int* cn2p3 = cn2.begin(3);

  // Creation de la bboxtree
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
    
  // Algorithme de projection
  #pragma omp parallel
  {
    E_Float p[3]; E_Float pr1[3]; E_Float pr2[3]; E_Float pi[3];
    E_Float p0[3]; E_Float p1[3]; E_Float p2[3];
    E_Float dist; E_Float distc; 
    E_Int ret; E_Int ind1, ind2, ind3;
    E_Float dx, dy, dz, normp, ps;
    vector<E_Int> indicesBB;
    E_Int nbboxes, et;
    
    #pragma omp for
    for (E_Int ind = 0; ind < npts; ind++)
    {
      p[0] = fx[ind]; p[1] = fy[ind]; p[2] = fz[ind];
      pr1[0] = p[0]; pr1[1] = p[1]; pr1[2] = p[2];
      pr2[0] = p[0] + nx; pr2[1] = p[1] + ny; pr2[2] = p[2] + nz;
      bbtree.getIntersectingBoxes(pr1, pr2, indicesBB, tol);
    
      distc = 1.e6;
      nbboxes = indicesBB.size();
      for (E_Int noe = 0; noe < nbboxes; noe++)
      {
        et = indicesBB[noe];
        ind1 = cn2p1[et]-1; ind2 = cn2p2[et]-1; ind3 = cn2p3[et]-1;
        p0[0] = fx2[ind1]; p0[1] = fy2[ind1]; p0[2] = fz2[ind1];
        p1[0] = fx2[ind2]; p1[1] = fy2[ind2]; p1[2] = fz2[ind2];
        p2[0] = fx2[ind3]; p2[1] = fy2[ind3]; p2[2] = fz2[ind3];        
        ret = K_COMPGEOM::intersectRayTriangle(p0, p1, p2,
                                               pr1, pr2, pi);
        if (ret == 1)
        {
          dx = pi[0]-p[0]; dy = pi[1]-p[1]; dz = pi[2]-p[2];
          dist = dx*dx + dy*dy + dz*dz;
          if (oriented != 0)
          {
            normp = sqrt(nx*nx+ny*ny+nz*nz);//normale
            ps = nx*dx+ny*dy+nz*dz/(dist*normp);
            if (ps > 0.) 
            {
              if (dist < distc) 
              {fx[ind] = pi[0]; fy[ind] = pi[1]; fz[ind] = pi[2]; distc = dist;}
            }
          }
          else 
          {
            if (dist < distc) 
            {fx[ind] = pi[0]; fy[ind] = pi[1]; fz[ind] = pi[2]; distc = dist;}
          }
        }
      }
      indicesBB.clear();
    }
  }

  // delete
  for (E_Int et = 0; et < nelts2; et++) delete boxes[et];
  boxes.clear();
}

//=============================================================================
/* projection avec precond par bboxtrees 
   effectue sur une liste de surfaces a projeter sur la meme surface */
//=============================================================================
void K_COMPGEOM::projectDirWithPrecond(
  E_Float nx, E_Float ny, E_Float nz,
  FldArrayI& cn2, E_Float* fx2, E_Float* fy2, E_Float* fz2,
  vector<E_Int>& sizet,
  vector<E_Float*>& fxt, vector<E_Float*>& fyt, vector<E_Float*>& fzt, E_Int oriented) 
{
  E_Float tol = 1.e-6;
  typedef K_SEARCH::BoundingBox<3> BBox3DType;
  E_Int nelts2 = cn2.getSize();
  E_Int* cn2p1 = cn2.begin(1); E_Int* cn2p2 = cn2.begin(2); E_Int* cn2p3 = cn2.begin(3);

  // Creation de la bboxtree
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
    
  // Algorithme de projection
  E_Int nzones = sizet.size();

  #pragma omp parallel
  {
    E_Float p[3]; E_Float pr1[3]; E_Float pr2[3]; E_Float pi[3];
    E_Float p0[3]; E_Float p1[3]; E_Float p2[3];
    E_Float dist; E_Float distc; 
    E_Int ret; E_Int ind1, ind2, ind3;
    E_Float dx, dy, dz, normp, ps;
    vector<E_Int> indicesBB;
    E_Int nbboxes, et;
  
    for (E_Int v = 0; v < nzones; v++)
    {
      E_Int npts = sizet[v];
      E_Float* fx = fxt[v];
      E_Float* fy = fyt[v];
      E_Float* fz = fzt[v];

      #pragma omp for
      for (E_Int ind = 0; ind < npts; ind++)
      {
        p[0] = fx[ind]; p[1] = fy[ind]; p[2] = fz[ind];
        pr1[0] = p[0]; pr1[1] = p[1]; pr1[2] = p[2];
        pr2[0] = p[0] + nx; pr2[1] = p[1] + ny; pr2[2] = p[2] + nz;
        bbtree.getIntersectingBoxes(pr1, pr2, indicesBB, tol);
      
        distc = 1e6;
        nbboxes = indicesBB.size();
        for (E_Int noe = 0; noe < nbboxes; noe++)
        {
          et = indicesBB[noe];
          ind1 = cn2p1[et]-1; ind2 = cn2p2[et]-1; ind3 = cn2p3[et]-1;
          p0[0] = fx2[ind1]; p0[1] = fy2[ind1]; p0[2] = fz2[ind1];
          p1[0] = fx2[ind2]; p1[1] = fy2[ind2]; p1[2] = fz2[ind2];
          p2[0] = fx2[ind3]; p2[1] = fy2[ind3]; p2[2] = fz2[ind3];        
          ret = K_COMPGEOM::intersectRayTriangle(p0, p1, p2,
                                                pr1, pr2, pi);
          if (ret == 1)
          {
            dx = pi[0]-p[0]; dy = pi[1]-p[1]; dz = pi[2]-p[2];
            dist = dx*dx + dy*dy + dz*dz;
            if (oriented != 0)
            {
              normp = sqrt(nx*nx+ny*ny+nz*nz);//normale
              ps = nx*dx+ny*dy+nz*dz/(dist*normp);
              if (ps > 0.) 
              {
                if (dist < distc) 
                {fx[ind] = pi[0]; fy[ind] = pi[1]; fz[ind] = pi[2]; distc = dist;}
              }
            }
            else 
            {
              if (dist < distc) 
              {fx[ind] = pi[0]; fy[ind] = pi[1]; fz[ind] = pi[2]; distc = dist;}
            }           
          }
        }
        indicesBB.clear();
      }
    }
  }
  // delete
  for (E_Int et = 0; et < nelts2; et++) delete boxes[et];
  boxes.clear();
}

//=============================================================================
/* projection avec precond par bboxtrees 
   effectue sur une liste de surfaces a projeter sur la meme surface
   La direction est variable pour chaque point a projeter.
*/
//=============================================================================
void K_COMPGEOM::projectDirWithPrecond(
  vector<E_Float*>& nxt, vector<E_Float*>& nyt, vector<E_Float*>& nzt,
  FldArrayI& cn2, E_Float* fx2, E_Float* fy2, E_Float* fz2,
  vector<E_Int>& sizet,
  vector<E_Float*>& fxt, vector<E_Float*>& fyt, vector<E_Float*>& fzt, E_Int oriented) 
{
  E_Float tol = 1.e-6;
  typedef K_SEARCH::BoundingBox<3> BBox3DType;
  E_Int nelts2 = cn2.getSize();
  E_Int* cn2p1 = cn2.begin(1); E_Int* cn2p2 = cn2.begin(2); E_Int* cn2p3 = cn2.begin(3);

  // Creation de la bboxtree
  vector<BBox3DType*> boxes(nelts2);// liste des bbox de ts les elements de a2
  K_FLD::FldArrayF bbox(nelts2,6);// xmin, ymin, zmin, xmax, ymax, zmax
  K_COMPGEOM::boundingBoxOfUnstrCells(cn2, fx2, fy2, fz2, bbox);
  E_Float minB[3];  E_Float maxB[3];
  E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(4);
  E_Float* yminp = bbox.begin(2); E_Float* ymaxp = bbox.begin(5);
  E_Float* zminp = bbox.begin(3); E_Float* zmaxp = bbox.begin(6);
  for (E_Int et = 0; et < nelts2; et++)
  {
    minB[0] = xminp[et]; minB[1] = yminp[et]; minB[2] = zminp[et];
    maxB[0] = xmaxp[et]; maxB[1] = ymaxp[et]; maxB[2] = zmaxp[et]; 
    boxes[et] = new BBox3DType(minB, maxB);
  }
  // Build the box tree.
  K_SEARCH::BbTree3D bbtree(boxes);
    
  // Algorithme de projection
  E_Float p[3]; E_Float pr1[3]; E_Float pr2[3]; E_Float pi[3];
  E_Float p0[3]; E_Float p1[3]; E_Float p2[3];
  E_Float dist; E_Float distc; 
  E_Int ret; E_Int ind1, ind2, ind3;
  E_Float dx, dy, dz, normp, ps;
  vector<E_Int> indicesBB;
  E_Int nbboxes, et;
  E_Int nzones = sizet.size();
  E_Float nxl, nyl, nzl;
  for (E_Int v = 0; v < nzones; v++)
  {
    E_Int npts = sizet[v];
    E_Float* fx = fxt[v];
    E_Float* fy = fyt[v];
    E_Float* fz = fzt[v];
    E_Float* nx = nxt[v];
    E_Float* ny = nyt[v];
    E_Float* nz = nzt[v];
    for (E_Int ind = 0; ind < npts; ind++)
    {
      nxl = nx[ind]; nyl = ny[ind]; nzl = nz[ind];
      p[0] = fx[ind]; p[1] = fy[ind]; p[2] = fz[ind];
      pr1[0] = p[0]; pr1[1] = p[1]; pr1[2] = p[2];
      pr2[0] = p[0] + nxl; pr2[1] = p[1] + nyl; pr2[2] = p[2] + nzl;
      bbtree.getIntersectingBoxes(pr1, pr2, indicesBB, tol);
      
      distc = 1e6;
      nbboxes = indicesBB.size();
      for (E_Int noe = 0; noe < nbboxes; noe++)
      {
        et = indicesBB[noe];
        ind1 = cn2p1[et]-1; ind2 = cn2p2[et]-1; ind3 = cn2p3[et]-1;
        p0[0] = fx2[ind1]; p0[1] = fy2[ind1]; p0[2] = fz2[ind1];
        p1[0] = fx2[ind2]; p1[1] = fy2[ind2]; p1[2] = fz2[ind2];
        p2[0] = fx2[ind3]; p2[1] = fy2[ind3]; p2[2] = fz2[ind3];        
        ret = K_COMPGEOM::intersectRayTriangle(p0, p1, p2,
                                               pr1, pr2,
                                               pi);
        
        if (ret == 1)
        {
          dx = pi[0]-p[0]; dy = pi[1]-p[1]; dz = pi[2]-p[2];
          dist = dx*dx + dy*dy + dz*dz;
          if ( oriented != 0 )
          {
            normp = sqrt(nxl*nxl+nyl*nyl+nzl*nzl);//normale
            ps = nxl*dx+nyl*dy+nzl*dz/(dist*normp);
            if ( ps > 0. ) 
            {
              if (dist < distc) 
              {fx[ind] = pi[0]; fy[ind] = pi[1]; fz[ind] = pi[2]; distc = dist;}
            }
          }
          else 
          {
            if (dist < distc) 
            {fx[ind] = pi[0]; fy[ind] = pi[1]; fz[ind] = pi[2]; distc = dist;}
          }          
        }
      }
      indicesBB.clear();
    }
  }
  for (E_Int et = 0; et < nelts2; et++) delete boxes[et];
  boxes.clear();
}

//=============================================================================
/* projection sans preconditionnement: retourne fx, fy, fz les coordonnees de
   la zone projetee */
//=============================================================================
void K_COMPGEOM::projectDirWithoutPrecond(
  E_Float nx, E_Float ny, E_Float nz,
  E_Int npts, E_Int nelts2, FldArrayI& cn2,
  E_Float* fx2, E_Float* fy2, E_Float* fz2,
  E_Float* fx, E_Float* fy, E_Float* fz, E_Int oriented)
{
  E_Int* cn2p1 = cn2.begin(1); 
  E_Int* cn2p2 = cn2.begin(2); 
  E_Int* cn2p3 = cn2.begin(3);
  // Algorithme de projection
  E_Float p[3]; E_Float pr1[3]; E_Float pr2[3]; E_Float pi[3];
  E_Float p0[3]; E_Float p1[3]; E_Float p2[3];
  E_Float dist; E_Float distc; 
  E_Int ret; E_Int ind1, ind2, ind3;
  E_Float dx, dy, dz, normp, ps;
  for (E_Int ind = 0; ind < npts; ind++)
  {
    p[0] = fx[ind]; p[1] = fy[ind]; p[2] = fz[ind];
    pr1[0] = p[0]; pr1[1] = p[1]; pr1[2] = p[2];
    pr2[0] = p[0] + nx; pr2[1] = p[1] + ny; pr2[2] = p[2] + nz;
    distc = 1e6;
    for (E_Int e = 0; e < nelts2; e++)
    {
      ind1 = cn2p1[e]-1;
      ind2 = cn2p2[e]-1;
      ind3 = cn2p3[e]-1;
      p0[0] = fx2[ind1]; p0[1] = fy2[ind1]; p0[2] = fz2[ind1];
      p1[0] = fx2[ind2]; p1[1] = fy2[ind2]; p1[2] = fz2[ind2];
      p2[0] = fx2[ind3]; p2[1] = fy2[ind3]; p2[2] = fz2[ind3];        
     
      ret = K_COMPGEOM::intersectRayTriangle(p0, p1, p2,
                                             pr1, pr2,
                                             pi);
      if (ret == 1)
      {
        dx = pi[0]-p[0]; dy = pi[1]-p[1]; dz = pi[2]-p[2];
        dist = dx*dx + dy*dy + dz*dz;
        if ( oriented != 0 )
        {
          normp = sqrt(nx*nx+ny*ny+nz*nz);//normale
          ps = nx*dx+ny*dy+nz*dz/(dist*normp);
          if ( ps > 0. ) 
            {
              if (dist < distc) 
              {fx[ind] = pi[0]; fy[ind] = pi[1]; fz[ind] = pi[2]; distc = dist;}
            }
        }
        else 
        {
          if (dist < distc) 
          {fx[ind] = pi[0]; fy[ind] = pi[1]; fz[ind] = pi[2]; distc = dist;}
        }   
      }
    }
  }
}
