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
# include "CompGeom/compGeom.h"

using namespace std;
using namespace K_FLD;

//===========================================================================
/* Triangulation de Delaunay */
//===========================================================================
void K_COMPGEOM::delaunay(E_Float coefa, E_Float coefb, E_Float coefc, 
                          E_Float coefd,
                          FldArrayF& coord, FldArrayI& connect, 
                          E_Int keepBB)
{
  E_Float eps = 1.e-12;
  list<Triangle*> triangles;
  list<Triangle*> compTriangles;
  list<Triangle*>::iterator itrT;
  list<Triangle*>::iterator itrT2;
  list<Edge*> edges;
  E_Int indA, indB, indC;

  E_Float dmax2;
  E_Bool isInCC;
  FldArrayF norm2; // SQUARED norm of points in field
  // Sort coordinates with respect to their norm
  K_SORT::sortCoordinates(coord, norm2);
  // Remove the points that have identical coordinates in field
  K_SORT::removeDoublePoints(coord, norm2);

  E_Float* norm2p = norm2.begin();
  E_Int sizeIni = coord.getSize();
  E_Int nfld = coord.getNfld();

  /* Compute the bbox and make 2 triangles with it */
  compAndTileBoundingBox(coefa, coefb, coefc, coefd, 
                         coord, triangles);
  
  /*---------------*/
  /* Insert points */
  /*---------------*/
  for (E_Int ind = 0; ind < sizeIni; ind++)
  {
    // 1-iter on triangles
    itrT  = triangles.begin(); //init

    while (itrT != triangles.end())
    {
      (*itrT)->getVertexIndices(indA, indB,  indC);
      dmax2 = (*itrT)->getDistMaxOfCC();
      
      if (norm2p[ind] > dmax2 + eps)
      { 
        compTriangles.push_back(*itrT);
        itrT2 = itrT;
        itrT++;
        triangles.erase(itrT2);
      }
      else // test if pt is in circum circle
      {
        isInCC = (*itrT)->isPointInCircumCircle(ind, coord);
        
        if (isInCC == true)
        {
          insertTriangleEdgesInList(*itrT, edges);
          itrT2 = itrT;
          itrT++;
          delete *itrT2;
          triangles.erase(itrT2);
        }
        else itrT++;
      }
    }

    /* Create new triangles with list of triangles and current point
       and reinit edges list */
    if (edges.size() != 0)
    {
      insertNewTriangles(coefa, coefb, coefc, coefd, 
                         ind, coord, edges, triangles);
    }
  }

  // Fusion des deux listes
  triangles.splice(triangles.begin(), compTriangles);

  /* Remove triangles using bbox vertices */
  if (keepBB == 0) 
  {
    removeTrianglesWithBBoxVertices(sizeIni, triangles);
    // Remove bbox vertices from field array
    coord.reAllocMat(sizeIni, nfld);
  }
  /*--------------------*/
  /* Build connectivity */
  /*--------------------*/
  buildConnectivity(triangles, connect);

  /* Nettoyage triangles */
  for (itrT = triangles.begin(); itrT != triangles.end(); itrT++)
    delete *itrT;
}

//=============================================================================
/* Tile the bounding box of the cloud of points in triangles */
//=============================================================================
void K_COMPGEOM::compAndTileBoundingBox(E_Float coefa, E_Float coefb, 
                                        E_Float coefc, E_Float coefd,
                                        FldArrayF& field, 
                                        list<Triangle*>& triangles)
{ 
  E_Float xmin, xmax, ymin, ymax, zmin, zmax;
  E_Float inv, pt1, pt2, pt3, pt4;
  E_Int npts = field.getSize();

  // Compute the bbox
  boundingBoxUnstruct(npts, field.begin(1), field.begin(2), field.begin(3), xmin, ymin, zmin, xmax, ymax, zmax);

  // Realloc field array to put 4 new points
  E_Int sizeIni = field.getSize();
  E_Int nfld = field.getNfld();
  field.reAllocMat(sizeIni+4, nfld);
  
  E_Float* xt = field.begin(1);
  E_Float* yt = field.begin(2);
  E_Float* zt = field.begin(3);

  // Extend bounding box from delta
  E_Float delta = 0.1; //10
  E_Float deltax = delta * (xmax - xmin);
  E_Float deltay = delta * (ymax - ymin);
  E_Float deltaz = delta * (zmax - zmin);

  xmin = xmin - deltax;
  xmax = xmax + deltax;
  ymin = ymin - deltay;
  ymax = ymax + deltay;
  zmin = zmin - deltaz;
  zmax = zmax + deltaz;

  // Compute bounding box in plane
  if (K_FUNC::fEqualZero(coefc, 1.e-12) == false)
  {
    inv = 1.0 / coefc;
    pt1 = -(coefa*xmin+coefb*ymin+coefd) * inv;
    pt2 = -(coefa*xmax+coefb*ymin+coefd) * inv;
    pt3 = -(coefa*xmax+coefb*ymax+coefd) * inv;
    pt4 = -(coefa*xmin+coefb*ymax+coefd) * inv;

    xt[sizeIni] = xmin;
    yt[sizeIni] = ymin;
    zt[sizeIni] = pt1;
    
    xt[sizeIni+1] = xmax;
    yt[sizeIni+1] = ymin;
    zt[sizeIni+1] = pt2;

    xt[sizeIni+2] = xmax;
    yt[sizeIni+2] = ymax;
    zt[sizeIni+2] = pt3;
  
    xt[sizeIni+3] = xmin;
    yt[sizeIni+3] = ymax;
    zt[sizeIni+3] = pt4;
  }
  else if (K_FUNC::fEqualZero(coefb, 1.e-12) == false)
  {
    inv = 1. / coefb;
    pt1 = -(coefa*xmin+coefc*zmin+coefd) * inv;
    pt2 = -(coefa*xmax+coefc*zmin+coefd) * inv;
    pt3 = -(coefa*xmax+coefc*zmax+coefd) * inv;
    pt4 = -(coefa*xmin+coefc*zmax+coefd) * inv;
 
    xt[sizeIni] = xmin;
    yt[sizeIni] = pt1;
    zt[sizeIni] = zmin;

    xt[sizeIni+1] = xmax;
    yt[sizeIni+1] = pt2;
    zt[sizeIni+1] = zmin;
    
    xt[sizeIni+2] = xmax;
    yt[sizeIni+2] = pt3;
    zt[sizeIni+2] = zmax;
  
    xt[sizeIni+3] = xmin;
    yt[sizeIni+3] = pt4;
    zt[sizeIni+3] = zmax;
  }
  else
  {
    inv = 1. / coefa;
    pt1 = -(coefb*ymin+coefc*zmin+coefd) * inv;
    pt2 = -(coefb*ymax+coefc*zmin+coefd) * inv;
    pt3 = -(coefb*ymax+coefc*zmax+coefd) * inv;
    pt4 = -(coefb*ymin+coefc*zmax+coefd) * inv;

    xt[sizeIni] = pt1;
    yt[sizeIni] = ymin;
    zt[sizeIni] = zmin;
  
    xt[sizeIni+1] = pt2;
    yt[sizeIni+1] = ymax;
    zt[sizeIni+1] = zmin;

    xt[sizeIni+2] = pt3;
    yt[sizeIni+2] = ymax;
    zt[sizeIni+2] = zmax;

    xt[sizeIni+3] = pt4;
    yt[sizeIni+3] = ymin;
    zt[sizeIni+3] = zmax;
  }

  // Cut with respect to the diagonal AC
  Triangle* t1 = 
    new Triangle(coefa, coefb, coefc, coefd, 
                 sizeIni, sizeIni+1, sizeIni+2, field);
  Triangle* t2 = 
    new Triangle(coefa, coefb, coefc, coefd,
                 sizeIni, sizeIni+2, sizeIni+3, field);
  triangles.push_back(t1);
  triangles.push_back(t2);

}

//=============================================================================
/* Insert the triangle edges AB, BC and CA in the list of edges 
   If one of the edges is already in the list, erase it */
//=============================================================================
void K_COMPGEOM::insertTriangleEdgesInList(Triangle* oneTriangle, 
                                           list<Edge*>& edges)
{
  E_Int vertex1, vertex2;
  E_Int indA, indB, indC;
  oneTriangle->getVertexIndices(indA, indB, indC);
  
  E_Bool found1 = false;
  E_Bool found2 = false;
  E_Bool found3 = false;

  list<Edge*>::iterator itr;
  list<Edge*>::iterator itr2;

  itr = edges.begin();

  while (itr != edges.end())
  {
    vertex1 = (*itr)->s1;
    vertex2 = (*itr)->s2;
    // check for AB
    if ( (vertex1 == indA && vertex2 == indB) ||
         (vertex1 == indB && vertex2 == indA))
    {
      itr2 = itr;
      itr++;
      free(*itr2);
      edges.erase(itr2);
      found1 = true;
    }
    // check for BC
    else if ((vertex1 == indB && vertex2 == indC) ||
             (vertex1 == indC && vertex2 == indB))
    {
      itr2 = itr;
      itr++;
      free(*itr2);
      edges.erase(itr2);
      found2 = true;
    }
    // check for CA
    else if ((vertex1 == indC && vertex2 == indA) ||
             (vertex1 == indA && vertex2 == indC))
    {
      itr2 = itr;
      itr++;
      free(*itr2);
      edges.erase(itr2);
      found3 = true;
    }    
    else itr++;
  }

  if ( found1 == false )
  {
    Edge* edge = (Edge*)malloc( sizeof(Edge) );
    edge->s1 = indA;
    edge->s2 = indB;
    edges.push_back(edge);
  }  
  if ( found2 == false )
  {
    Edge* edge = (Edge*)malloc( sizeof(Edge) );
    edge->s1 = indB;
    edge->s2 = indC;
    edges.push_back(edge);
  } 
  if ( found3 == false )
  {
    Edge* edge = (Edge*)malloc( sizeof(Edge) );
    edge->s1 = indC;
    edge->s2 = indA;
    edges.push_back(edge);
  }
}

//=============================================================================
/* Insert new triangles made of pt H(field(ind,.)) and points from edges list
 and remove from the edge list at the same time each edge */
//=============================================================================
void K_COMPGEOM::insertNewTriangles(E_Float coefa, E_Float coefb, 
                                    E_Float coefc, E_Float coefd,
                                    E_Int ind, FldArrayF& field,
                                    list<Edge*>& edges, 
                                    list<Triangle*>& triangles)
{
  E_Int ind1, ind2;
  list<Edge*>::iterator itr;
  list<Edge*>::iterator itr2;
  itr = edges.begin();

  while (itr != edges.end())
  {
    ind1 = (*itr)->s1;
    ind2 = (*itr)->s2;
    Triangle* t = new Triangle(coefa, coefb, coefc, coefd, 
                               ind1, ind2, ind, field);
    triangles.push_back(t); 
    itr2 = itr;
    itr++;
    free(*itr2);
  }
  edges.clear();
}
//=============================================================================
/* Remove triangles using the bounding box vertices  */
//=============================================================================
void K_COMPGEOM::removeTrianglesWithBBoxVertices(
  E_Int sizeIni, list<Triangle*>& triangles)
{
  E_Int ind1, ind2, ind3;
  list<Triangle*>::iterator itr = triangles.begin();
  list<Triangle*>::iterator itr2;
  while (itr != triangles.end())
  {
    (*itr)->getVertexIndices(ind1, ind2, ind3);

    // The triangle contains one bbox vertex
    if (ind1 >= sizeIni || ind2 >= sizeIni || ind3 >= sizeIni)
    {
      itr2 = itr;
      itr++;
      delete *itr2;
      triangles.erase(itr2);
    }
    else itr++;
  }
}

//=============================================================================
/* Computes connectivity for triangles  */
//=============================================================================
void K_COMPGEOM::buildConnectivity(list<Triangle*>& triangles,
                                   FldArrayI& connect)
{
  E_Int nelts = triangles.size(); 
  E_Int ind1, ind2, ind3;
  connect.malloc(nelts, 3);
  E_Int n = 0;

  E_Int* connect1 = connect.begin(1);
  E_Int* connect2 = connect.begin(2);
  E_Int* connect3 = connect.begin(3);
  
  for ( list<Triangle*>::iterator itr = triangles.begin(); 
        itr != triangles.end(); itr++)
  {
    (*itr)->getVertexIndices(ind1, ind2, ind3);
    connect1[n] = ind1+1;
    connect2[n] = ind2+1;
    connect3[n] = ind3+1;
    n++;
  }
}
