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

# include "generator.h"
# include "Search/BbTree.h"

using namespace K_FLD;
using namespace std;

namespace K_GENERATOR
{
//=============================================================================
/* Calcul le point interpole dans un quad.
   Utilise une interpolation bilineaire si le quad est convexe.
   Utilise un decoupage en triangles sinon.
   IN: les 4 points du quad
   IN: alpha, beta: alpha et beta le long des lignes i et j
   IN: edge: le style de decoupage (0,1 ou 2)
   OUT: Coordonnees du point(alpha, beta)
*/
#define _point_(x1,x2,x3,xo,alpha,beta) \
  A = (1.-alpha)*x1 + alpha*x2;         \
  B = (1.-beta)*x1 + beta*x3;           \
  C = (1.-beta)*x2 + beta*x3;           \
  D = (1.-alpha-beta)*x1 + (alpha+beta)*x2;\
  E = (1.-alpha-beta)*x1 + (alpha+beta)*x3;\
  F = (1.-alpha)*x3 + alpha*x2;         \
  xo = (1.-alpha-beta)*( A + B - x1) +\
    alpha*( C + D - x2 ) +\
    beta*( E + F - x3);
//=============================================================================
void buildPoint(E_Float x1, E_Float y1, E_Float z1,
                E_Float x2, E_Float y2, E_Float z2,
                E_Float x3, E_Float y3, E_Float z3,
                E_Float x4, E_Float y4, E_Float z4,
                E_Int edge,
                E_Float alpha, E_Float beta,
                E_Float& xo, E_Float& yo, E_Float& zo)
{
  E_Float A, B, C, D, E, F;
  E_Float alpha1 = 1.-alpha;
  E_Float beta1 = 1.-beta;
  if (edge == 0) // le quad est convexe, l'interpolation bilineaire est OK
  {
    A = alpha1*x1 + alpha*x2;
    B = alpha1*x4 + alpha*x3;
    xo = beta1*A + beta*B;
    A = alpha1*y1 + alpha*y2;
    B = alpha1*y4 + alpha*y3;
    yo = beta1*A + beta*B;
    A = alpha1*z1 + alpha*z2;
    B = alpha1*z4 + alpha*z3;
    zo = beta1*A + beta*B;
  }
  else if (edge == 1) // decoupage 1-3
  {
    if (beta <= alpha) // triangle 1-2-3
    {
      _point_(x2,x1,x3,xo,alpha1,beta);
      _point_(y2,y1,y3,yo,alpha1,beta);
      _point_(z2,z1,z3,zo,alpha1,beta);
    }
    else // triangle 1-3-4
    {
      _point_(x4,x3,x1,xo,alpha,beta1);
      _point_(y4,y3,y1,yo,alpha,beta1);
      _point_(z4,z3,z1,zo,alpha,beta1);
    }
  }
  else // edge = 2 - decoupage 2-4
  {
    if (alpha+beta <= 1) // triangle 1-2-4
    {
      _point_(x1,x2,x4,xo,alpha,beta);
      _point_(y1,y2,y4,yo,alpha,beta);
      _point_(z1,z2,z4,zo,alpha,beta);
    }
    else // triangle 2-3-4
    {
      _point_(x3,x4,x2,xo,alpha1,beta1);
      _point_(y3,y4,y2,yo,alpha1,beta1);
      _point_(z3,z4,z2,zo,alpha1,beta1);
    }
  }
  
}
//=============================================================================
/* Engendre un ensemble de grilles structurees */
//=============================================================================
PyObject* front2Struct(PyObject* self, PyObject* args)
{
  PyObject *front, *surface, *distrib;
  E_Int Vmin;
  E_Float dist;
#if defined E_DOUBLEREAL && defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "OOOld", &front, &surface, &distrib, &Vmin, 
                        &dist)) return NULL;
#elif defined E_DOUBLEREAL && !defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "OOOid", &front, &surface, &distrib, &Vmin, 
                        &dist)) return NULL;
#elif !defined E_DOUBLEREAL && defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "OOOlf", &front, &surface, &distrib, &Vmin, 
                        &dist)) return NULL;
#else
  if (!PyArg_ParseTuple(args, "OOOif", &front, &surface, &distrib, &Vmin, 
                        &dist)) return NULL;
#endif

  // Check front: must be quad
  E_Int ni1, nj1, nk1;
  FldArrayF* f1; FldArrayI* cn1;
  char* varString1; char* eltType1;
  E_Int res1 = K_ARRAY::getFromArray(front, varString1, 
                                     f1, ni1, nj1, nk1, cn1, eltType1);
  if (res1 != 2 || strcmp(eltType1, "QUAD") != 0 ) 
  {
    if ( res1 == 1 ) delete f1;
    if ( res1 == 2 ) {delete f1; delete cn1;}
    PyErr_SetString(PyExc_TypeError, 
                    "front2Struct: surface must be a QUAD array.");
    return NULL;
  }
  E_Int posx1 = K_ARRAY::isCoordinateXPresent(varString1);
  E_Int posy1 = K_ARRAY::isCoordinateYPresent(varString1);
  E_Int posz1 = K_ARRAY::isCoordinateZPresent(varString1);
  if ( posx1 == -1 || posy1 == -1 || posz1 == -1 )
  {
    delete f1; delete cn1;
    PyErr_SetString(PyExc_TypeError, 
                    "front2Struct: coords must be present in front.");
    return NULL;
  }

  // Check projection surface : must be TRI
  E_Int ni2, nj2, nk2;
  FldArrayF* f2;
  char* varString2;
  char* eltType2;
  FldArrayI* cn2;
  E_Int res2 = K_ARRAY::getFromArray(surface, varString2, 
                                     f2, ni2, nj2, nk2, cn2, eltType2);
  if ( res2 != 2 || strcmp(eltType2, "TRI") != 0 ) 
  {
    delete f1; delete cn1;
    if ( res2 == 1 ) delete f2;
    if ( res2 == 2 ) {delete f2; delete cn2;}
    PyErr_SetString(PyExc_TypeError, 
                    "front2Struct: surface must be a TRI array.");
    return NULL;
  }
  E_Int posx2 = K_ARRAY::isCoordinateXPresent(varString2);
  E_Int posy2 = K_ARRAY::isCoordinateYPresent(varString2);
  E_Int posz2 = K_ARRAY::isCoordinateZPresent(varString2);
  if ( posx2 == -1 || posy2 == -1 || posz2 == -1 )
  {
    delete f1; delete cn1; delete f2; delete cn2;
    PyErr_SetString(PyExc_TypeError, 
                    "front2Struct: coords must be present in surface.");
    return NULL;
  }
  
  // Check distrib
  E_Int nid, njd, nkd;
  FldArrayF* fd;
  char* varStringd;
  char* eltTyped;
  FldArrayI* cnd;
  E_Int resd = K_ARRAY::getFromArray(distrib, varStringd, 
                                     fd, nid, njd, nkd, cnd, eltTyped);
  if ( resd != 1 )
  {
    delete f1; delete cn1; delete f2; delete cn2;
    delete fd; if ( res2 == 2 )  delete cnd;
    PyErr_SetString(PyExc_TypeError, 
                    "front2Struct: distrib must be a structured array.");
    return NULL;
  }
  if ( nid < 2 || njd != 1 || nkd != 1) 
  {
    delete f1; delete cn1; delete f2; delete cn2; delete fd;
    PyErr_SetString(PyExc_TypeError, 
                    "front2Struct: distrib must be an i-array.");
    return NULL;
  }  
  E_Int posxd = K_ARRAY::isCoordinateXPresent(varStringd);
  E_Int posyd = K_ARRAY::isCoordinateYPresent(varStringd);
  E_Int poszd = K_ARRAY::isCoordinateZPresent(varStringd);
  if ( posxd == -1 || posyd == -1 || poszd == -1 )
  {
     delete f1; delete cn1; delete f2; delete cn2; delete fd;
    PyErr_SetString(PyExc_TypeError, 
                    "front2Struct: coords must be present in distrib.");
    return NULL;
  }
  posx1++; posy1++; posz1++; 
  posx2++; posy2++; posz2++; 
  posxd++; posyd++; poszd++;

  // Check vmin
  if (Vmin < 2) Vmin = 2;

  E_Float* xt1 = f1->begin(posx1);
  E_Float* yt1 = f1->begin(posy1);
  E_Float* zt1 = f1->begin(posz1);
  E_Float* xt2 = f2->begin(posx2);
  E_Float* yt2 = f2->begin(posy2);
  E_Float* zt2 = f2->begin(posz2);

  // Projection du front
  E_Int npts = f1->getSize();
  FldArrayF proj(npts, 3); proj.setAllValuesAtNull();
  proj.setOneField(*f1, posx1, 1);
  proj.setOneField(*f1, posy1, 2);
  proj.setOneField(*f1, posz1, 3);
  E_Float* px = proj.begin(1);
  E_Float* py = proj.begin(2);
  E_Float* pz = proj.begin(3);
  
  // precondt par bbox/kdtree
  K_COMPGEOM::projectOrthoWithPrecond(npts, *cn2, posx2, posy2, posz2, *f2, 
                                      px, py, pz);

  // Info sur la projection orthogonale
  E_Float dMin = 1.e6;
  E_Float dMax = -1.e6;
  E_Float lx, ly, lz, ll;
  for (E_Int i = 0; i < npts; i++)
  {
    lx = px[i] - xt1[i];
    ly = py[i] - yt1[i];
    lz = pz[i] - zt1[i];
    ll = sqrt(lx*lx + ly*ly + lz*lz);
    dMin = K_FUNC::E_min(dMin, ll);
    dMax = K_FUNC::E_max(dMax, ll);
  }
  printf("Info: proj ortho dist min %f.\n", dMin);
  printf("Info: proj ortho dist max %f.\n", dMax);

  // Precond pour projectDir
  typedef K_SEARCH::BoundingBox<3> BBox3DType;
  E_Int nelts2 = cn2->getSize();
  vector<BBox3DType*> boxes(nelts2);
  K_FLD::FldArrayF bbox(nelts2, 6);// xmin, ymin, zmin, xmax, ymax, zmax
  K_COMPGEOM::boundingBoxOfUnstrCells(*cn2, xt2, yt2, zt2, bbox);
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

  // Compteurs
  E_Int nQuads = cn1->getSize();
  E_Int concaveQuads = 0;
  E_Int projectedPoints = 0;

  // Passage pour chaque quad
  E_Float x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;
  E_Float xp1, yp1, zp1, xp2, yp2, zp2, xp3, yp3, zp3, xp4, yp4, zp4;
  FldArrayI& cnf = *cn1;
  E_Int Nk = fd->getSize();
  E_Float* fdx = fd->begin(1);
  PyObject* meshes = PyList_New(0);
  PyObject* tpl;
  E_Int ind, ret1, ret2;
  E_Int ret = 0; 
  E_Int edge = 0;
  E_Float xo, yo, zo;
  E_Float alpha, beta, gamma, lp;
  E_Float l=1.e6;
  E_Float Cx, Fx;
  E_Float Cy, Fy;
  E_Float Cz, Fz;
  E_Float Vmin1 = 1./(Vmin-1.);
  E_Float pr1[3]; E_Float pr2[3];
  vector<E_Int> indicesBB;
  E_Float p1[3]; E_Float p2[3]; E_Float p3[3]; E_Float p4[3];

  for (E_Int n = 0; n < nQuads; n++)
  {
    ind = cnf(n,1)-1;
    x1 = xt1[ind]; y1 = yt1[ind]; z1 = zt1[ind];
    xp1 = px[ind]; yp1 = py[ind]; zp1 = pz[ind];

    ind = cnf(n,2)-1;
    x2 = xt1[ind]; y2 = yt1[ind]; z2 = zt1[ind];
    xp2 = px[ind]; yp2 = py[ind]; zp2 = pz[ind];

    ind = cnf(n,3)-1;
    x3 = xt1[ind]; y3 = yt1[ind]; z3 = zt1[ind];
    xp3 = px[ind]; yp3 = py[ind]; zp3 = pz[ind];

    ind = cnf(n,4)-1;
    x4 = xt1[ind]; y4 = yt1[ind]; z4 = zt1[ind];
    xp4 = px[ind]; yp4 = py[ind]; zp4 = pz[ind];
 
    // Verifie la convexite des quad
    p1[0] = x1; p1[1] = y1; p1[2] = z1;
    p2[0] = x2; p2[1] = y2; p2[2] = z2;
    p3[0] = x3; p3[1] = y3; p3[2] = z3;
    p4[0] = x4; p4[1] = y4; p4[2] = z4;
    ret1 = K_COMPGEOM::checkQuadConvexity(p1, p2, p3, p4);
    p1[0] = xp1; p1[1] = yp1; p1[2] = zp1;
    p2[0] = xp2; p2[1] = yp2; p2[2] = zp2;
    p3[0] = xp3; p3[1] = yp3; p3[2] = zp3;
    p4[0] = xp4; p4[1] = yp4; p4[2] = zp4;
    ret2 = K_COMPGEOM::checkQuadConvexity(p1, p2, p3, p4);
    if (ret2 != 0) concaveQuads++;
    if (ret1 == 0 && ret2 == 0) edge = 0;
    else if (ret1 == 0 && ret2 != 0) edge = ret2;
    else if (ret1 != 0 && ret2 == 0) edge = ret1;
    else if (ret1 != 0 && ret2 != 0)
    {
      if (ret1 == ret2) edge = ret1;
      else
      {printf("Warning: front2Struct: mesh may contain folded cells.\n");
        edge = ret1;}
    }
    FldArrayF* f = new FldArrayF(Vmin*Vmin*Nk, 3);
    E_Float* fx = f->begin(1);
    E_Float* fy = f->begin(2);
    E_Float* fz = f->begin(3);

    l = 1.e6;

    for (E_Int j = 0; j < Vmin; j++)
      for (E_Int i = 0; i < Vmin; i++)
      {
        alpha = i*Vmin1;
        beta = j*Vmin1;
        //alpha1 = 1.-alpha;
        //beta1 = 1.-beta;

        buildPoint(x1, y1, z1, x2, y2, z2,
                   x3, y3, z3, x4, y4, z4, edge, alpha, beta,
                   Cx, Cy, Cz);
        buildPoint(xp1, yp1, zp1, xp2, yp2, zp2,
                   xp3, yp3, zp3, xp4, yp4, zp4, edge, alpha, beta, 
                   Fx, Fy, Fz);

        pr1[0] = Fx; pr1[1] = Fy; pr1[2] = Fz;
        pr2[0] = Cx; pr2[1] = Cy; pr2[2] = Cz;
        bbtree.getIntersectingBoxes(pr1, pr2, indicesBB, 1.e-6);

        ret = 
          K_COMPGEOM::projectDir(
            Fx, Fy, Fz,
            Fx-Cx, Fy-Cy, Fz-Cz,
            xt2, yt2, zt2, indicesBB,
            *cn2, xo, yo, zo);
        indicesBB.clear();

        // On force la non projection pour l'instant
        ret = -1;

        if (ret != -1)
        {
          // Le point se projete a une distance lp
          lp = sqrt( (xo-Fx)*(xo-Fx) + (yo-Fy)*(yo-Fy) +(zo-Fz)*(zo-Fz) );
          if (lp > dist)
          {
            // Projection limitee a dist
            Fx = Fx + dist*(xo-Fx)/lp;
            Fy = Fy + dist*(yo-Fy)/lp;
            Fz = Fz + dist*(zo-Fz)/lp;
            //Fx = xo; Fy = yo; Fz = zo;
          }
          else // projection sur la surface
          {
            Fx = xo; Fy = yo; Fz = zo;
            projectedPoints++;
          }
        }
        else
        {
          // Pas de projection : on prend l'interpolation (deja dans Fx)
        }
        l = sqrt( (Fx-Cx)*(Fx-Cx) + (Fy-Cy)*(Fy-Cy) +(Fz-Cz)*(Fz-Cz) );
        for (E_Int k = 0; k < Nk; k++)
        {
          ind = i + j*Vmin + k*Vmin*Vmin;
          //gamma = k*1. / (Nk-1.);
          gamma = fdx[k] / l;
          fx[ind] = Fx + gamma*(Cx-Fx);
          fy[ind] = Fy + gamma*(Cy-Fy);
          fz[ind] = Fz + gamma*(Cz-Fz);
        }
    }
    tpl = K_ARRAY::buildArray(*f, "x,y,z", Vmin, Vmin, Nk);
    delete f;
    PyList_Append(meshes, tpl);
    Py_DECREF(tpl);
  }

  // Resume
  printf("Info: concave quads: %d over %d quads.\n", concaveQuads, nQuads);
  printf("Info: projected points: %d over %d points.\n", projectedPoints, 
         int(nQuads*Vmin*Vmin));

  // Delete les boxes
  E_Int boxesSize = boxes.size();
  for (E_Int v = 0; v < boxesSize; v++) delete boxes[v];

  // Sortie
  delete f1; delete cn1; delete f2; delete cn2; delete fd;  
  return meshes;
}

//=============================================================================
/* Engendre un ensemble de grilles structurees */
//=============================================================================
PyObject* fillWithStruct( PyObject* self, PyObject* args )
{
  PyObject *mesh;
  E_Int Vmin;
#ifdef E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "Ol", &mesh, &Vmin)) return NULL;
#else
  if (!PyArg_ParseTuple(args, "Oi", &mesh, &Vmin)) return NULL;
#endif
  // Check mesh : must be quad pour l'instant
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res1 = K_ARRAY::getFromArray(mesh, varString, 
                                     f, ni, nj, nk, cn, eltType);
  if (res1 != 2 || strcmp(eltType, "QUAD") != 0) 
  {
    if (res1 == 1) delete f;
    if (res1 == 2) {delete f; delete cn;}
    PyErr_SetString(PyExc_TypeError, 
                    "fillWithStruct: mesh must be a QUAD array.");
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if ( posx == -1 || posy == -1 || posz == -1 )
  {
    delete f; delete cn;
    PyErr_SetString(PyExc_TypeError, 
                    "fillWithStruct: coords must be present in mesh.");
    return NULL;
  }
  posx++; posy++; posz++; 

  // Check vmin
  if (Vmin < 2) Vmin = 2;

  E_Float* xt = f->begin(posx);
  E_Float* yt = f->begin(posy);
  E_Float* zt = f->begin(posz);

  // Compteurs
  E_Int nQuads = cn->getSize();
  E_Int concaveQuads = 0;
  //E_Int projectedPoints = 0;

  // Passage pour chaque quad
  E_Float x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;
  //E_Float xp1, yp1, zp1, xp2, yp2, zp2, xp3, yp3, zp3, xp4, yp4, zp4;
  FldArrayI& cnf = *cn;
  PyObject* meshes = PyList_New(0);
  PyObject* tpl;
  E_Int ind, ret1;
  //E_Int ret = 0; 
  E_Int edge = 0;
  //E_Float xo, yo, zo;
  E_Float alpha, beta;
  E_Float Cx;
  E_Float Cy;
  E_Float Cz;
  E_Float Vmin1 = 1./(Vmin-1.);
  E_Float p1[3], p2[3], p3[3], p4[3];

  for (E_Int n = 0; n < nQuads; n++)
  {
    ind = cnf(n,1)-1;
    x1 = xt[ind]; y1 = yt[ind]; z1 = zt[ind];

    ind = cnf(n,2)-1;
    x2 = xt[ind]; y2 = yt[ind]; z2 = zt[ind];

    ind = cnf(n,3)-1;
    x3 = xt[ind]; y3 = yt[ind]; z3 = zt[ind];

    ind = cnf(n,4)-1;
    x4 = xt[ind]; y4 = yt[ind]; z4 = zt[ind];
 
    // Verifie la convexite des quad
    p1[0] = x1; p1[1] = y1; p1[2] = z1;
    p2[0] = x2; p2[1] = y2; p2[2] = z2;
    p3[0] = x3; p3[1] = y3; p3[2] = z3;
    p4[0] = x4; p4[1] = y4; p4[2] = z4;
    ret1 = K_COMPGEOM::checkQuadConvexity(p1, p2, p3, p4);

    if (ret1 != 0) concaveQuads++;
    edge = ret1;
    
    FldArrayF* fl = new FldArrayF(Vmin*Vmin, 3);
    E_Float* fx = fl->begin(1);
    E_Float* fy = fl->begin(2);
    E_Float* fz = fl->begin(3);

    ind = 0;
    for (E_Int j = 0; j < Vmin; j++)
      for (E_Int i = 0; i < Vmin; i++)
      {
        alpha = i*Vmin1;
        beta = j*Vmin1;
        //alpha1 = 1.-alpha;
        //beta1 = 1.-beta;

        buildPoint(x1, y1, z1, x2, y2, z2,
                   x3, y3, z3, x4, y4, z4, edge, alpha, beta,
                   Cx, Cy, Cz);

        fx[ind] = Cx;
        fy[ind] = Cy;
        fz[ind] = Cz; ind++;
      }
    
    tpl = K_ARRAY::buildArray(*fl, "x,y,z", Vmin, Vmin, 1);
    delete fl;
    PyList_Append(meshes, tpl);
    Py_DECREF(tpl);
  }

  // Resume
  printf("Info: concave quads : %d over %d quads.\n", concaveQuads, nQuads);

  // Sortie
  delete f; delete cn; 
  return meshes;
}


}
