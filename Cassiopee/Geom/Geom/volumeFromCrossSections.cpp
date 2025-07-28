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

// Volume from cross sections
# include "geom.h"

using namespace K_FLD;
using namespace std;
using namespace K_CONST;
using namespace K_SEARCH;

//=============================================================================
PyObject* K_GEOM::volumeFromCrossSections(PyObject* self,
                                          PyObject* args )
{
  PyObject* array1; PyObject* array2;
  PyObject* contour1; PyObject* contour2;

  if (!PyArg_ParseTuple(args, "OOOO", &array1, &array2, 
                        &contour1, &contour2)) return NULL;

  // Extraction des 2 contours
  E_Int imc1, jmc1, kmc1;
  FldArrayF* fc1;
  FldArrayI* cnpoly1;
  char* varStringc1; char* eltTypec1;
  E_Int imc2, jmc2, kmc2;
  FldArrayF* fc2;
  FldArrayI* cnpoly2;
  char* varStringc2; char* eltTypec2;
  E_Int resc1 = 
    K_ARRAY::getFromArray3(contour1, varStringc1, fc1, imc1, jmc1, kmc1, 
                           cnpoly1, eltTypec1);  
  E_Int resc2 = 
    K_ARRAY::getFromArray3(contour2, varStringc2, fc2, imc2, jmc2, kmc2, 
                           cnpoly2, eltTypec2);
  // check contours
  E_Int posxc1, posyc1, poszc1;
  E_Int posxc2, posyc2, poszc2;

  if (resc1 == 2 && resc2 == 2)
  {
    if ( strcmp(eltTypec1, "BAR") != 0 || strcmp(eltTypec2, "BAR") != 0)
    {
      RELEASESHAREDU(contour1, fc1, cnpoly1);
      RELEASESHAREDU(contour2, fc2, cnpoly2);
      PyErr_SetString(PyExc_ValueError,
                      "volumeFromCrossSections: contours must be BAR-arrays.");
      return NULL;
    }
   
    posxc1 = K_ARRAY::isCoordinateXPresent(varStringc1);
    posyc1 = K_ARRAY::isCoordinateYPresent(varStringc1);
    poszc1 = K_ARRAY::isCoordinateZPresent(varStringc1);
    posxc2 = K_ARRAY::isCoordinateXPresent(varStringc2);
    posyc2 = K_ARRAY::isCoordinateYPresent(varStringc2);
    poszc2 = K_ARRAY::isCoordinateZPresent(varStringc2);
    
    if (posxc1 == -1 || posyc1 == -1 || poszc1 == -1 ||
        posxc2 == -1 || posyc2 == -1 || poszc2 == -1 )
    {
      RELEASESHAREDU(contour1, fc1, cnpoly1);
      RELEASESHAREDU(contour2, fc2, cnpoly2);
      PyErr_SetString(PyExc_ValueError,
                      "volumeFromCrossSections: coordinates not found in contours.");
      return NULL;
    }
    posxc1++; posyc1++; poszc1++; posxc2++; posyc2++; poszc2++;
  }
  else 
  {
    if ( resc1 == 1 ) RELEASESHAREDS(contour1, fc1);
    if ( resc2 == 1 ) RELEASESHAREDS(contour2, fc2);
    PyErr_SetString(PyExc_ValueError,
                    "volumeFromCrossSections: contours must be BAR-arrays.");
    return NULL;
  }
  // FldArrayF coordpoly1(fc1->getSize(),3);
  // FldArrayF coordpoly2(fc2->getSize(),3);
  // coordpoly1.setOneField(*fc1, posxc1, 1);
  // coordpoly1.setOneField(*fc1, posyc1, 2);
  // coordpoly1.setOneField(*fc1, poszc1, 3);
  // coordpoly2.setOneField(*fc2, posxc2, 1);
  // coordpoly2.setOneField(*fc2, posyc2, 2);
  // coordpoly2.setOneField(*fc2, poszc2, 3);

  // Extraction des triangles de Delaunay des 2 cross-sections
  E_Int im1, jm1, km1;
  FldArrayF* f1; FldArrayI* cn1;
  char* varString1; char* eltType1;
  E_Int im2, jm2, km2;
  FldArrayF* f2; FldArrayI* cn2;
  char* varString2; char* eltType2;

  E_Int res1 = 
    K_ARRAY::getFromArray3(array1, varString1, f1, im1, jm1, km1, cn1, 
                           eltType1); 
  E_Int res2 = 
    K_ARRAY::getFromArray3(array2, varString2, f2, im2, jm2, km2, cn2, 
                           eltType2);
  // Check args
  if (res1 != 2 || strcmp(eltType1, "TRI") != 0)
  {
    RELEASESHAREDB(res1, array1, f1, cn1);
    RELEASESHAREDB(res2, array2, f2, cn2);    
    delete fc1; delete fc2; delete cnpoly1; delete cnpoly2;

          
    PyErr_SetString(PyExc_ValueError,
                    "volumeFromCrossSections: array1 must be a TRI-array.");
    return NULL;
  }
  if (res2 != 2 || strcmp(eltType2, "TRI") != 0)
  {
    RELEASESHAREDB(res1, array1, f1, cn1);
    RELEASESHAREDB(res2, array2, f2, cn2);    
    RELEASESHAREDU(contour1, fc1, cnpoly1);
    RELEASESHAREDU(contour2, fc2, cnpoly2);
    PyErr_SetString(PyExc_ValueError,
                    "volumeFromCrossSections: array2 must be a TRI-array.");
    return NULL;
  }
  
  E_Int posx1 = K_ARRAY::isCoordinateXPresent( varString1);
  E_Int posy1 = K_ARRAY::isCoordinateYPresent( varString1);
  E_Int posz1 = K_ARRAY::isCoordinateZPresent( varString1);
  if (posx1 == -1 || posy1 == -1 || posz1 == -1)
  {
    RELEASESHAREDU(array1, f1, cn1);
    RELEASESHAREDU(array2, f2, cn2);
    RELEASESHAREDU(contour1, fc1, cnpoly1);
    RELEASESHAREDU(contour2, fc2, cnpoly2);
    PyErr_SetString(PyExc_TypeError,
                    "volumeFromCrossSections: can't find coordinates in array.");
    return NULL;
  }
  posx1++; posy1++; posz1++;

  E_Int posx2 = K_ARRAY::isCoordinateXPresent( varString2 );
  E_Int posy2 = K_ARRAY::isCoordinateYPresent( varString2 );
  E_Int posz2 = K_ARRAY::isCoordinateZPresent( varString2 );
  if (posx2 == -1 || posy2 == -1 || posz2 == -1)
  {
    RELEASESHAREDU(array1, f1, cn1);
    RELEASESHAREDU(array2, f2, cn2);
    RELEASESHAREDU(contour1, fc1, cnpoly1);
    RELEASESHAREDU(contour2, fc2, cnpoly2);
    PyErr_SetString(PyExc_TypeError,
                    "volumeFromCrossSections: can't find coordinates in array.");
    return NULL;
  }
  posx2++; posy2++; posz2++;

  // perturbation des triangulations
  FldArrayF coord1(f1->getSize(), 3);
  FldArrayF coord2(f2->getSize(), 3);
  perturbate2D(posx1, posy1, posz1, *f1, posx2, posy2, posz2, *f2, 
               coord1, coord2);

  // Determine les hauteurs de chaque coupe
  // On suppose que les coupes sont en z=cte
  //E_Float z1 = 0.;
  //E_Float z2 = 1.;
  //if (coord1.getSize() > 0) z1 = coord1(0, posz1);
  //if (coord2.getSize() > 0) z2 = coord2(0, posz2);

  // calcul des centres des cercles circonscrits a t1 et t2
  FldArrayF coordCC1; 
  FldArrayF coordCC2;
  compCCCenters(coord1.begin(1), coord1.begin(2),coord1.begin(3),*cn1, coordCC1);
  compCCCenters(coord2.begin(1), coord2.begin(2),coord2.begin(3),*cn2, coordCC2);

  // Creation des tableaux resultat
  E_Int nv = 0;
  E_Int nvmax = 1000;
  E_Int ne = 0;
  E_Int nemax = 1000;
  FldArrayF* an = new FldArrayF(nvmax, 3);
  FldArrayF& coord = *an;
  FldArrayI* cni = new FldArrayI(nemax, 4);
  FldArrayI& cn = *cni;
  
  // calcul des tetraedres de type T1 et T2
  FldArrayI type1(cn1->getSize()); //type des triangles de T1
  FldArrayI type2(cn2->getSize()); //type des triangles de T2
  type1.setAllValuesAtNull(); type2.setAllValuesAtNull();

  FldArrayI type(nemax); // type des tetraedres crees
  compTetraType1(1, posxc1, posyc1, poszc1, *fc1, *cnpoly1, 
                 *cn1, coordCC1, coord1, coord2, 
                 nv, ne, type1, coord, cn, type);

  compTetraType1(2, posxc2, posyc2, poszc2, *fc2, *cnpoly2, 
                 *cn2, coordCC2, coord2, coord1, 
                 nv, ne, type2, coord, cn, type);

  // Creation des tetradres de type T12 : issus de 2 aretes 
  // des diagrammes de Voronoi
  compTetraType12(posxc1, posyc1, poszc1, *fc1, *cnpoly1, 
                  *cn1, coordCC1, coord1, 
                  posxc2, posyc2, poszc2, *fc2,  *cnpoly2,
                  *cn2, coordCC2, coord2, type1, type2,
                  nv, ne, coord, cn, type);
  cn.reAllocMat(ne, 4);
  coord.reAllocMat(nv, 3);
  
  // Nettoyage de la connectivite et des points (doublons...)
  K_CONNECT::cleanConnectivity(1, 2, 3, 1.e-10, eltType1, coord, cn);

  // check tetraedres de type t1/t2 : si pas de voisin de type t12, 
  // t1 ou t2 elimine
  checkTetrahedra(type, cn, coord);

  PyObject* tpl = K_ARRAY::buildArray(coord, "x,y,z", cn, 4, "TETRA");
  delete an; delete cni;
  RELEASESHAREDU(array1, f1, cn1);
  RELEASESHAREDU(array2, f2, cn2); 
  RELEASESHAREDU(contour1, fc1, cnpoly1);
  RELEASESHAREDU(contour2, fc2, cnpoly2);
  return tpl;
}

//===========================================================================
/* Verification des tetraedres : elimination de ceux de type t1 ou t2 qui 
   n'ont aucun voisin de type t12 */
//===========================================================================
void K_GEOM::checkTetrahedra(FldArrayI& type, FldArrayI& cn, FldArrayF& coord)
{
  // Determination de la connectivite tetra-tetra voisins
  E_Int nelts = cn.getSize();
  vector< vector<E_Int> > cEEN(nelts);
  K_CONNECT::connectEV2EENbrs("TETRA", coord.getSize(), cn, cEEN);

  E_Int cnt = 0;
  FldArrayI temp(nelts, 4);
  
  // parcours des elements
  for (E_Int et1 = 0; et1 < nelts; et1++)
  {
    if ( type[et1] < 10 )
    { 
      vector<E_Int>& eltVoisins = cEEN[et1];//liste des tetras voisins de et1
      E_Int nvoisins = eltVoisins.size();

      E_Int found = 0; // un voisin est il t12 ?
      for (E_Int n = 0; n < nvoisins; n++ )
      {
        E_Int et2 = eltVoisins[n]; // indice du tetra voisins 
        
        if ( type[et2] == 12 )
        {
          found = 1;
          break;
        }
      }
      
      if ( found == 1 )
      {
        for (E_Int v = 1; v <= 4; v++)
          temp(cnt, v) = cn(et1,v);
        cnt++;
      }
    }
    else 
    {
      for (E_Int v = 1; v <= 4; v++)
        temp(cnt, v) = cn(et1,v);
      cnt++; 
    }
  }
  temp.reAllocMat(cnt, 4);
 
  cn = temp;
}
//===========================================================================
/* Calcul des centres des cercles circonscrits de tous les triangles 
   coordCC est alloue par la routine
*/
//===========================================================================
void K_GEOM::compCCCenters(E_Float* xt, E_Float* yt, E_Float* zt, FldArrayI& cn, FldArrayF& coordCC)
{
  E_Int nelts = cn.getSize();
  coordCC.malloc(nelts,3);
  E_Float p1[3];
  E_Float p2[3];
  E_Float p3[3];
  E_Float pc[3];
  E_Float R;
  E_Float* xCC = coordCC.begin(1);
  E_Float* yCC = coordCC.begin(2);
  E_Float* zCC = coordCC.begin(3);
  
  E_Int* cn1 = cn.begin(1);
  E_Int* cn2 = cn.begin(2);
  E_Int* cn3 = cn.begin(3);
  for (E_Int et = 0; et < nelts; et++)
  {
    E_Int ind1 = cn1[et]-1;
    E_Int ind2 = cn2[et]-1;
    E_Int ind3 = cn3[et]-1;
    p1[0] = xt[ind1]; p1[1] = yt[ind1]; p1[2] = zt[ind1];
    p2[0] = xt[ind2]; p2[1] = yt[ind2]; p2[2] = zt[ind2];
    p3[0] = xt[ind3]; p3[1] = yt[ind3]; p3[2] = zt[ind3];
    K_COMPGEOM::circumCircle(p1, p2, p3, pc, R);
 
    xCC[et] = pc[0];
    yCC[et] = pc[1];
    zCC[et] = pc[2];
  }
}

//===========================================================================
/* Construction des tetraedres formes par une base tri et un vertex de 
   type Ti */
//===========================================================================
void K_GEOM::compTetraType1(E_Int type0, 
                            E_Int posxc1, E_Int posyc1, E_Int poszc1, 
                            FldArrayF& coordpoly1, FldArrayI& cnpoly1,
                            FldArrayI& cn1, FldArrayF& coordCC1, 
                            FldArrayF& coord1, FldArrayF& coord2,
                            E_Int& nv, E_Int& ne, FldArrayI& type1,
                            FldArrayF& coord, FldArrayI& cn,
                            FldArrayI& type)
{
  E_Int ne1 = coordCC1.getSize();
  E_Float P[3];
  E_Int nvmax = coord.getSize();
  E_Int nemax = cn.getSize();
  
  E_Float* xCC1 = coordCC1.begin(1);
  E_Float* yCC1 = coordCC1.begin(2);
  E_Float* zCC1 = coordCC1.begin(3);
  E_Int* cn11 = cn1.begin(1);
  E_Int* cn12 = cn1.begin(2);
  E_Int* cn13 = cn1.begin(3);

  E_Float* xt1 = coord1.begin(1);
  E_Float* yt1 = coord1.begin(2);
  E_Float* zt1 = coord1.begin(3);
  E_Float* xt2 = coord2.begin(1);
  E_Float* yt2 = coord2.begin(2);
  E_Float* zt2 = coord2.begin(3);

  E_Float* xt = coord.begin(1);
  E_Float* yt = coord.begin(2);
  E_Float* zt = coord.begin(3);
  E_Int* cn01 = cn.begin(1);
  E_Int* cn02 = cn.begin(2);
  E_Int* cn03 = cn.begin(3);
  E_Int* cn04 = cn.begin(4);

  // calcul des bounding box du contour 1 
  E_Float xmin1, xmax1, ymin1, ymax1, zmin1, zmax1;
  K_COMPGEOM::boundingBox(posxc1, posyc1, poszc1, coordpoly1, 
                          xmin1, ymin1,zmin1, xmax1, ymax1, zmax1);

  // Extend bounding box from delta
  E_Float delta = 0.05;//10% ds la triangulation de Delaunay
  E_Float deltax1 = delta * (xmax1 - xmin1);
  E_Float deltay1 = delta * (ymax1 - ymin1);
  E_Float deltaz1 = delta * (zmax1 - zmin1);
  xmin1 = xmin1 - deltax1; xmax1 = xmax1 + deltax1;
  ymin1 = ymin1 - deltay1; ymax1 = ymax1 + deltay1;
  zmin1 = zmin1 - deltaz1; zmax1 = zmax1 + deltaz1;

  ArrayAccessor<FldArrayF> coordAcc(coord2, 1,2,3);
  KdTree<FldArrayF> kdt(coordAcc, E_EPSILON);

  E_Float Q[3]; 
  for (E_Int et = 0; et < ne1; et++) //parcours des triangles
  {
    P[0] = xCC1[et];
    P[1] = yCC1[et];
    P[2] = zCC1[et];
    E_Int ind = kdt.getClosest(P);
    
    if ( ind != -1 )
    {
      E_Int ind1 = cn11[et]-1;
      E_Int ind2 = cn12[et]-1;
      E_Int ind3 = cn13[et]-1;
      
      // test si un des sommets est exterieur au contour : 
      // creation du tetraedre
      Q[0] = xt1[ind1]; Q[1] = yt1[ind1]; Q[2] = zt1[ind1];
      E_Int in1 = 
        K_COMPGEOM::pointInBB(xmin1, ymin1, zmin1, xmax1, ymax1, zmax1, Q);
      Q[0] = xt1[ind2]; Q[1] = yt1[ind2]; Q[2] = zt1[ind2];
      E_Int in2 = 
        K_COMPGEOM::pointInBB(xmin1, ymin1, zmin1, xmax1, ymax1, zmax1, Q);
      Q[0] = xt1[ind3]; Q[1] = yt1[ind3]; Q[2] = zt1[ind3];
      E_Int in3 = 
        K_COMPGEOM::pointInBB(xmin1, ymin1, zmin1, xmax1, ymax1, zmax1, Q);

      if ( in1 == 1 && in2 == 1 && in3 == 1 )
      {
        E_Int in = testIfEdgesCentersInContour(ind1, ind2, ind3, coord1, 
                                               posxc1, posyc1, poszc1,
                                               coordpoly1, cnpoly1);

        if ( in == 1 )
        {
          type1[et] = type0; //le triangle est ds un tetraedre de type type0

          xt[nv] = xt1[ind1];
          yt[nv] = yt1[ind1];
          zt[nv] = zt1[ind1];
          xt[nv+1] = xt1[ind2];
          yt[nv+1] = yt1[ind2];
          zt[nv+1] = zt1[ind2];
          xt[nv+2] = xt1[ind3];
          yt[nv+2] = yt1[ind3];
          zt[nv+2] = zt1[ind3];
          xt[nv+3] = xt2[ind];
          yt[nv+3] = yt2[ind];
          zt[nv+3] = zt2[ind];

          cn01[ne] = nv+1; cn02[ne] = nv+2; cn03[ne] = nv+3; cn04[ne] = nv+4;
          type[ne] = type0;
        
          ne++;

          if (ne >= nemax-1) 
          {
            nemax = nemax + 1000;
            cn.reAllocMat(nemax, 4);
            cn01 = cn.begin(1);
            cn02 = cn.begin(2);
            cn03 = cn.begin(3);
            cn04 = cn.begin(4);
            type.reAlloc(nemax);
          }
          nv = nv + 4;
          if ( nv >= nvmax-4 ) 
          {
            nvmax = nvmax + 1000;
            coord.reAllocMat(nvmax, 3);
            xt = coord.begin(1);
            yt = coord.begin(2);
            zt = coord.begin(3);
          }
        }
      }
    }
  }
}

//===========================================================================
/* Creation du tableau des aretes du diagramme de Voronoi d'une 
   triangulation */
//===========================================================================
void K_GEOM::compVoronoiEdges(E_Int nvertex, FldArrayI& cn, FldArrayI& vedges)
{
  char* eltType = new char[10];
  strcpy(eltType,"TRI");
  E_Int nelts = cn.getSize();
  vedges.malloc(3*nelts,2); // allocation max
  vector< vector<E_Int> > cEEN(nelts);
  K_CONNECT::connectEV2EENbrs(eltType, nvertex, cn, cEEN);
  
  E_Int nv = 0;
  FldArrayIS dejaVu(nelts, nelts);
  dejaVu.setAllValuesAtNull();
  
  E_Int* vedges1 = vedges.begin(1);
  E_Int* vedges2 = vedges.begin(2);

  for (E_Int et1 = 0; et1 < nelts; et1++)
  {
    vector<E_Int>& eltVoisins = cEEN[et1];
    E_Int nvoisins = eltVoisins.size();
    for (E_Int et2 = 0; et2 < nvoisins; et2++)
    {
      E_Int ev = eltVoisins[et2];
      if ( dejaVu(et1, ev+1) == 0 || dejaVu(ev,et1+1) == 0 )
      {
        vedges1[nv] = et1;
        vedges2[nv] = ev;
        dejaVu(et1, ev+1) = 1;
        dejaVu(ev, et1+1) = 1;
        nv++;
      }
    }
  }
  vedges.reAllocMat(nv, 2);
  delete [] eltType;
}

//===========================================================================
/* Calcul les tetraedres de type T12 : issus de 2 aretes intersectantes
   des diagrammes de Voronoi de T1 et T2 */
//===========================================================================
void K_GEOM::compTetraType12(E_Int posxc1, E_Int posyc1, E_Int poszc1, 
                             FldArrayF& fc1, FldArrayI& cnpoly1,
                             FldArrayI& cn1,
                             FldArrayF& coordCC1, FldArrayF& coord1,
                             E_Int posxc2, E_Int posyc2, E_Int poszc2, 
                             FldArrayF& fc2,FldArrayI& cnpoly2,
                             FldArrayI& cn2, FldArrayF& coordCC2,
                             FldArrayF& coord2,
                             FldArrayI& type1, FldArrayI& type2, 
                             E_Int& nv, E_Int& ne,
                             FldArrayF& coord, FldArrayI& cn, FldArrayI& type)
{
  E_Int nvmax = coord.getSize();
  E_Int nemax = cn.getSize();

  //calcul des bounding box des contours1 et 2
  E_Float xmin1, xmax1, ymin1, ymax1, zmin1, zmax1;
  E_Float xmin2, xmax2, ymin2, ymax2, zmin2, zmax2;
  K_COMPGEOM::boundingBox(posxc1, posyc1, poszc1, fc1, 
                          xmin1, ymin1,zmin1, xmax1, ymax1, zmax1);
  K_COMPGEOM::boundingBox(posxc2, posyc2, poszc2, fc2, 
                          xmin2, ymin2, zmin2, xmax2, ymax2, zmax2);

  // Extend bounding box from delta
  E_Float delta = 0.05;//10% ds la triangulation de Delaunay
  E_Float deltax1 = delta * (xmax1 - xmin1);
  E_Float deltay1 = delta * (ymax1 - ymin1);
  E_Float deltaz1 = delta * (zmax1 - zmin1);
  E_Float deltax2 = delta * (xmax2 - xmin2);
  E_Float deltay2 = delta * (ymax2 - ymin2);
  E_Float deltaz2 = delta * (zmax2 - zmin2);
  xmin1 = xmin1 - deltax1; xmin2 = xmin2 - deltax2;
  xmax1 = xmax1 + deltax1; xmax2 = xmax2 + deltax2;
  ymin1 = ymin1 - deltay1; ymin2 = ymin2 - deltay2;
  ymax1 = ymax1 + deltay1; ymax2 = ymax2 + deltay2;
  zmin1 = zmin1 - deltaz1; zmin2 = zmin2 - deltaz2;
  zmax1 = zmax1 + deltaz1; zmax2 = zmax2 + deltaz2;

  FldArrayI vedges1; //vedges(.,1) : indice du 1er elt, vedges(.,2) : 2eme elt 
  FldArrayI vedges2;

  // Calcul des aretes du diagramme de Voronoi de T1
  compVoronoiEdges(coord1.getSize(), cn1, vedges1);

  // Calcul des aretes du diagramme de Voronoi de T2
  compVoronoiEdges(coord2.getSize(), cn2, vedges2);

  E_Int v1Size = vedges1.getSize();
  E_Int v2Size = vedges2.getSize();
  E_Float ps1a[3]; E_Float ps1b[3];
  E_Float ps2a[3]; E_Float ps2b[3];
  E_Float pi0[3]; E_Float pi1[3];
  E_Float P[3];
  FldArrayI ind1(3); FldArrayI ind2(3);

  E_Float* xCC1 = coordCC1.begin(1);
  E_Float* yCC1 = coordCC1.begin(2);
  E_Float* xCC2 = coordCC2.begin(1);
  E_Float* yCC2 = coordCC2.begin(2);
  E_Float* xt1 = coord1.begin(1);
  E_Float* yt1 = coord1.begin(2);
  E_Float* zt1 = coord1.begin(3);
  E_Float* xt2 = coord2.begin(1);
  E_Float* yt2 = coord2.begin(2);
  E_Float* zt2 = coord2.begin(3);

  E_Float* xt = coord.begin(1);
  E_Float* yt = coord.begin(2);
  E_Float* zt = coord.begin(3);
  E_Int* cn01 = cn.begin(1);
  E_Int* cn02 = cn.begin(2);
  E_Int* cn03 = cn.begin(3);
  E_Int* cn04 = cn.begin(4);

  // test intersection des vedges1 et 2
  for (E_Int v1 = 0; v1 < v1Size; v1++)
  {
    E_Int et1a = vedges1(v1, 1);
    E_Int et1b = vedges1(v1, 2);
    
    if (type1[et1a] != 0 || type1[et1b] != 0 ) //de type type[et]-solide
    {      
      ps1a[0] = xCC1[et1a];
      ps1a[1] = yCC1[et1a];
      ps1a[2] = 0.; //coordCC1(et1a,3);   
      
      ps1b[0] = xCC1[et1b];
      ps1b[1] = yCC1[et1b];
      ps1b[2] = 0.; //coordCC1(et1b,3);
      
      for (E_Int v2 = 0; v2 < v2Size; v2++)
      {
        E_Int et2a = vedges2(v2, 1);
        E_Int et2b = vedges2(v2, 2);
        if ( type2[et2a] != 0 || type2[et2b] != 0 ) //is i-Solid
        {
          ps2a[0] = xCC2[et2a];
          ps2a[1] = yCC2[et2a];
          ps2a[2] = 0.; //coordCC2(et2a,3);   
      
          ps2b[0] = xCC2[et2b];
          ps2b[1] = yCC2[et2b];
          ps2b[2] = 0.; //coordCC2(et2b,3);

          E_Int intersect = K_COMPGEOM::intersect2Segments(ps1a, ps1b, 
                                                           ps2a, ps2b,
                                                           pi0, pi1);
          if ( intersect == 1 )
          {
            getDelaunayVertices(et1a, et1b, cn1, ind1);
            getDelaunayVertices(et2a, et2b, cn2, ind2);

            E_Int ind1a = ind1[0]-1; E_Int ind1b = ind1[1]-1;
            E_Int ind2a = ind2[0]-1; E_Int ind2b = ind2[1]-1;

            P[0] = xt1[ind1a]; P[1] = yt1[ind1a]; P[2] = zt1[ind1a];
            E_Int in1 = 
              K_COMPGEOM::pointInBB(xmin1, ymin1, zmin1, xmax1, ymax1, zmax1,P);

            P[0] = xt1[ind1b]; P[1] = yt1[ind1b]; P[2] = zt1[ind1b];
            E_Int in2 = 
              K_COMPGEOM::pointInBB(xmin1, ymin1, zmin1, xmax1, ymax1, zmax1,P);

            P[0] = xt2[ind2a]; P[1] = yt2[ind2a]; P[2] = zt2[ind2a];
            E_Int in3 = 
              K_COMPGEOM::pointInBB(xmin2, ymin2, zmin2, xmax2, ymax2, zmax2,P);

            P[0] = xt2[ind2b]; P[1] = yt2[ind2b]; P[2] = zt2[ind2b];
            E_Int in4 = 
              K_COMPGEOM::pointInBB(xmin2, ymin2, zmin2, xmax2, ymax2, zmax2,P);

            if (in1 == 1 && in2 == 1 && in3 == 1 && in4 == 1)
            {
              xt[nv] = xt1[ind1a];
              yt[nv] = yt1[ind1a];
              zt[nv] = zt1[ind1a]; 
              xt[nv+1] = xt1[ind1b];
              yt[nv+1] = yt1[ind1b]; 
              zt[nv+1] = zt1[ind1b]; 
              xt[nv+2] = xt2[ind2a];
              yt[nv+2] = yt2[ind2a];
              zt[nv+2] = zt2[ind2a]; 
              xt[nv+3] = xt2[ind2b];   
              yt[nv+3] = yt2[ind2b];
              zt[nv+3] = zt2[ind2b];
              cn01[ne] = nv+1; cn02[ne] = nv+2; 
              cn03[ne] = nv+3; cn04[ne] = nv+4;
              type[ne] = 12;
              ne++;

              if (ne >= nemax-1) 
              {
                nemax = nemax + 1000;
                cn.reAllocMat(nemax, 4);
                cn01 = cn.begin(1);
                cn02 = cn.begin(2);
                cn03 = cn.begin(3);
                cn04 = cn.begin(4);
                type.reAlloc(nemax);
              }
              nv = nv+4;
              if ( nv >= nvmax-4 ) 
              {
                nvmax = nvmax + 1000;
                coord.reAllocMat(nvmax, 3);
                xt = coord.begin(1);
                yt = coord.begin(2);
                zt = coord.begin(3);
              }
            }
          }
        }
      }
    }
  }
}

//===========================================================================
/* Determine les indices des vertices dans la triangulation de Delaunay 
   communs aux 2 elements et1 et et2. Retourne les indices des sommets 
   ind demarre a 1 */
//===========================================================================
void K_GEOM::getDelaunayVertices(E_Int et1, E_Int et2, FldArrayI& cEV, 
                                 FldArrayI& ind)
{
  E_Int indi1[3]; E_Int indi2[3];

  //indices des 3 sommets des triangles
  indi1[0] = cEV(et1,1); indi1[1] = cEV(et1,2); indi1[2] = cEV(et1,3);
  indi2[0] = cEV(et2,1); indi2[1] = cEV(et2,2); indi2[2] = cEV(et2,3);
  
  ind.setAllValuesAtNull();

  E_Int n = 0;
  for (E_Int i1 = 0; i1 < 3; i1++)
  {
    for (E_Int i2 = 0; i2 < 3; i2++)
    {
      if ( indi1[i1] == indi2[i2] ) 
      {
        ind[n] = indi1[i1];
        n++;
        break;
      }
    }
    if ( n == 2 ) break;
  }
  if ( n != 2 ) 
    printf("Warning: getDelaunayVertices: a common edge not found.\n");
}
//==========================================================================
/* Pertubation d'un des points si identiques dans le plan (x,y) 
   coord1, coord2 : 
*/
//==========================================================================
void 
K_GEOM::perturbate2D(E_Int posx1, E_Int posy1, E_Int posz1, FldArrayF& f1, 
                     E_Int posx2, E_Int posy2, E_Int posz2, FldArrayF& f2,
                     FldArrayF& coord1, FldArrayF& coord2)
{
  E_Float perturb = 1.e-12;
  E_Float dmin = 1.e-12;

  E_Int size1 = f1.getSize();
  E_Int size2 = f2.getSize();

  coord1.setOneField(f1, posx1, 1);
  coord1.setOneField(f1, posy1, 2);
  coord1.setOneField(f1, posz1, 3);
  coord2.setOneField(f2, posx2, 1);
  coord2.setOneField(f2, posy2, 2);
  coord2.setOneField(f2, posz2, 3);
  E_Float* xt1 = coord1.begin(1);
  E_Float* yt1 = coord1.begin(2);
  E_Float* xt2 = coord2.begin(1);
  E_Float* yt2 = coord2.begin(2);
  E_LONG idum = -1;

  for (E_Int ind1 = 0; ind1 < size1; ind1++)
  {
    for (E_Int ind2 = 0; ind2 < size2; ind2++)
    {
      E_Float dx = xt2[ind2]-xt1[ind1];
      E_Float dy = yt2[ind2]-yt1[ind1];
      
      if ( K_FUNC::E_abs(dx) < dmin && K_FUNC::E_abs(dy) < dmin)
      {
        coord1(ind1,1) = coord1(ind1,1) + K_NOISE::stdRand(&idum)*perturb;
        coord1(ind1,2) = coord1(ind1,2) + K_NOISE::stdRand(&idum)*perturb;
        break;
      }
    }
  }
}
//===========================================================================
/* 
   Test si le centre des 3 aretes d'un triangle defini par les 3 indices
   ind1, ind2, ind3 de coord sont toutes internes au 
   contour definit par (coordc, cnc).
   Une tolerance est ajoutee de facon a ce que la routine retourne 1
   meme si les centres sont sur le contour lui-meme.
   Retourne 1 si le centre des 3 aretes est a l'interieur du contour.
*/
//===========================================================================
E_Int 
K_GEOM::testIfEdgesCentersInContour(E_Int ind1, E_Int ind2, E_Int ind3,
                                    FldArrayF& coord, 
                                    E_Int posxc1, E_Int posyc1, E_Int poszc1, 
                                    FldArrayF& coordc, FldArrayI& cnc)
{
  E_Float x1 = coord(ind1,1); E_Float y1 = coord(ind1,2);
  E_Float x2 = coord(ind2,1); E_Float y2 = coord(ind2,2); 
  E_Float x3 = coord(ind3,1); E_Float y3 = coord(ind3,2);
  E_Float z = coord(ind1,3);

  E_Float dx12 = K_FUNC::E_abs(x2-x1); E_Float dy12 = K_FUNC::E_abs(y2-y1);
  E_Float dx23 = K_FUNC::E_abs(x2-x3); E_Float dy23 = K_FUNC::E_abs(y2-y3);
  E_Float dx31 = K_FUNC::E_abs(x3-x1); E_Float dy31 = K_FUNC::E_abs(y3-y1);
  E_Float dx = K_FUNC::E_max(dx12, dx31); 
  E_Float dy = K_FUNC::E_max(dy12, dy31);
  dx = K_FUNC::E_max(dx, dx23); dy = K_FUNC::E_max(dy, dy23);

  E_Float epsx = 0.01 * dx;
  E_Float epsy = 0.01 * dy;
  E_Float eps = K_FUNC::E_max(epsx, epsy);

  E_Float P[3];
  P[0] = 0.5*(x1+x2); P[1] = 0.5*(y1+y2); P[2] = z;

  E_Int isIn = K_COMPGEOM::pointInPolygon2D(coordc.begin(posxc1), coordc.begin(posyc1), coordc.begin(poszc1), 
                                            cnc, P, eps, 
                                            NULL, NULL);
  if ( isIn == 0 ) return 0;
  
  P[0] = 0.5*(x3+x2); P[1] = 0.5*(y3+y2); P[2] = z;
  isIn = K_COMPGEOM::pointInPolygon2D(coordc.begin(posxc1), coordc.begin(posyc1), coordc.begin(poszc1), 
                                      cnc, P, eps, 
                                      NULL, NULL);
  if ( isIn == 0 ) return 0;

  P[0] = 0.5*(x1+x3); P[1] = 0.5*(y1+y3); P[2] = z;
  isIn = K_COMPGEOM::pointInPolygon2D(coordc.begin(posxc1), coordc.begin(posyc1), coordc.begin(poszc1), 
                                      cnc, P, eps, 
                                      NULL, NULL);
  if ( isIn == 0 ) return 0;
  return 1;
}
