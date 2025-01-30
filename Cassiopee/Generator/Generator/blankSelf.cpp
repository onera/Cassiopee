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

# include "generator.h"

using namespace std;

struct facet
{
  E_Float x0, y0, z0; // points de la face
  E_Float x1, y1, z1;
  E_Float x2, y2, z2;
  E_Float x3, y3, z3;
  E_Float xmin, ymin, zmin;
  E_Float xmax, ymax, zmax;
  E_Int cell; // indice de la cellule voisine
  E_Int ind0, ind1, ind2, ind3; // indices des sommets de la face
};

// Cree la face noFace de la cellule indCell
// IN: cn: connectivite elts-vertex, ff: table des faces
facet* createFace(E_Int indCell, E_Int noFace, FldArrayI& cn, E_Int* ff, 
                  E_Float* fx, E_Float* fy, E_Float* fz)
{
  facet* face = new facet;
      
  face->ind0 = cn(indCell, ff[noFace+6*0])-1;
  face->x0 = fx[face->ind0];
  face->y0 = fy[face->ind0];
  face->z0 = fz[face->ind0];
        
  face->ind1 = cn(indCell, ff[noFace+6*1])-1;
  face->x1 = fx[face->ind1];
  face->y1 = fy[face->ind1];
  face->z1 = fz[face->ind1];
        
  face->ind2 = cn(indCell, ff[noFace+6*2])-1;
  face->x2 = fx[face->ind2];
  face->y2 = fy[face->ind2];
  face->z2 = fz[face->ind2];
        
  face->ind3 = cn(indCell, ff[noFace+6*3])-1;
  face->x3 = fx[face->ind3];
  face->y3 = fy[face->ind3];
  face->z3 = fz[face->ind3];
  
  face->xmin = K_FUNC::E_min(face->x0, face->x1);
  face->xmin = K_FUNC::E_min(face->xmin, face->x2);
  face->xmin = K_FUNC::E_min(face->xmin, face->x3);
  face->ymin = K_FUNC::E_min(face->y0, face->y1);
  face->ymin = K_FUNC::E_min(face->ymin, face->y2);
  face->ymin = K_FUNC::E_min(face->ymin, face->y3);
  face->zmin = K_FUNC::E_min(face->z0, face->z1);
  face->zmin = K_FUNC::E_min(face->zmin, face->z2);
  face->zmin = K_FUNC::E_min(face->zmin, face->z3);
        
  face->cell = indCell;
  return face;
}

// Retourne 1 si les deux faces s'intersectent
E_Int intersect(facet* face1, facet* face2)
{
  E_Float P0[3]; E_Float P1[3]; E_Float P2[3];
  E_Float Q0[3]; E_Float Q1[3]; E_Float Q2[3];
  
  E_Int mul = 0;
  E_Int ind0, ret;
  ind0 = face1->ind0;
  if (ind0 == face2->ind0 || ind0 == face2->ind1 || ind0 == face2->ind2 || ind0 == face2->ind3) mul++;
  ind0 = face1->ind1;
  if (ind0 == face2->ind0 || ind0 == face2->ind1 || ind0 == face2->ind2 || ind0 == face2->ind3) mul++;
  ind0 = face1->ind2;
  if (ind0 == face2->ind0 || ind0 == face2->ind1 || ind0 == face2->ind2 || ind0 == face2->ind3) mul++;
  ind0 = face1->ind3;
  if (ind0 == face2->ind0 || ind0 == face2->ind1 || ind0 == face2->ind2 || ind0 == face2->ind3) mul++;
  if (mul < 2)
  {
    // calcul les intersections TRI/TRI du decoupage de faces
    P0[0] = face1->x0; P0[1] = face1->y0; P0[2] = face1->z0;
    P1[0] = face1->x1; P1[1] = face1->y1; P1[2] = face1->z1;
    P2[0] = face1->x2; P2[1] = face1->y2; P2[2] = face1->z2;
          
    Q0[0] = face2->x0; Q0[1] = face2->y0; Q0[2] = face2->z0;
    Q1[0] = face2->x1; Q1[1] = face2->y1; Q1[2] = face2->z1;
    Q2[0] = face2->x2; Q2[1] = face2->y2; Q2[2] = face2->z2;
          
    ret = K_COMPGEOM::trianglesIntersection(P0,P1,P2,Q0,Q1,Q2,1.e-10);
    if (ret != 0) return 1;
           
    Q0[0] = face2->x0; Q0[1] = face2->y0; Q0[2] = face2->z0;
    Q1[0] = face2->x2; Q1[1] = face2->y2; Q1[2] = face2->z2;
    Q2[0] = face2->x3; Q2[1] = face2->y3; Q2[2] = face2->z3;

    ret = K_COMPGEOM::trianglesIntersection(P0,P1,P2,Q0,Q1,Q2,1.e-10);
    if (ret != 0) return 1;
          
    P0[0] = face1->x0; P0[1] = face1->y0; P0[2] = face1->z0;
    P1[0] = face1->x2; P1[1] = face1->y2; P1[2] = face1->z2;
    P2[0] = face1->x3; P2[1] = face1->y3; P2[2] = face1->z3;
          
    ret = K_COMPGEOM::trianglesIntersection(P0,P1,P2,Q0,Q1,Q2,1.e-10);
    if (ret != 0) return 1;
            
    Q0[0] = face2->x0; Q0[1] = face2->y0; Q0[2] = face2->z0;
    Q1[0] = face2->x1; Q1[1] = face2->y1; Q1[2] = face2->z1;
    Q2[0] = face2->x2; Q2[1] = face2->y2; Q2[2] = face2->z2;
          
    ret = K_COMPGEOM::trianglesIntersection(P0,P1,P2,Q0,Q1,Q2,1.e-10);
    if (ret != 0) return 1;  
  }
  return 0;
}

//=============================================================================
/* 
   blank self. 
   IN: hexa mesh + cellN in center
   OUT: cellN modifie in-place
 */
//=============================================================================
PyObject* K_GENERATOR::blankSelf(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* cellNObject;
  if (!PyArg_ParseTuple(args, "OO", &array, &cellNObject))
  {
    return NULL;
  }

  // Check array: must be HEXA mesh
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray2(array, varString, f, ni, nj, nk, 
                                     cn, eltType);
  if (res != 2 || strcmp(eltType, "HEXA") != 0) 
  {
    if (res != 0) RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "blankSelf: a must be HEXA.");
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError, 
                    "front2Hexa: coords must be present in a1.");
    return NULL;       
  }
  E_Float* fx = f->begin(posx+1);
  E_Float* fy = f->begin(posy+1);
  E_Float* fz = f->begin(posz+1);
  
  // Check cellN (like array but in center)
  E_Int ni1, nj1, nk1;
  FldArrayF* f1; FldArrayI* cn1;
  char* varString1; char* eltType1;
  E_Int res1 = K_ARRAY::getFromArray2(cellNObject, varString1, f1, ni1, nj1, nk1, 
                                      cn1, eltType1);
  
  E_Int posCellN = K_ARRAY::isCellNatureField1Present(varString1);
  if (posCellN == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    RELEASESHAREDB(res1, cellNObject, f1, cn1);
    PyErr_SetString(PyExc_TypeError,
                    "blankSelf: cellN must be present in second array.");
    return NULL;
  }
  E_Int nv = f->getSize();
  E_Int ne = cn->getSize();

  // tableau des faces pour un HEXA
  E_Int ff[24];
  ff[0 + 0*6] = 1; ff[0 + 1*6] = 4; ff[0 + 2*6] = 3; ff[0 + 3*6] = 2;
  ff[1 + 0*6] = 1; ff[1 + 1*6] = 2; ff[1 + 2*6] = 6; ff[1 + 3*6] = 5;
  ff[2 + 0*6] = 2; ff[2 + 1*6] = 3; ff[2 + 2*6] = 7; ff[2 + 3*6] = 6;
  ff[3 + 0*6] = 3; ff[3 + 1*6] = 4; ff[3 + 2*6] = 8; ff[3 + 3*6] = 7;
  ff[4 + 0*6] = 1; ff[4 + 1*6] = 5; ff[4 + 2*6] = 8; ff[4 + 3*6] = 4;
  ff[5 + 0*6] = 5; ff[5 + 1*6] = 6; ff[5 + 2*6] = 7; ff[5 + 3*6] = 8;
  
  // Construction de la connectivite elt/elt voisins
  vector< vector<E_Int> > cEEN(ne);
  vector< vector<E_Int> > commonFace(ne);
  K_CONNECT::connectEV2EENbrs(eltType, nv, *cn, cEEN, commonFace);
  
  // Extraction des faces
  E_Int cellNi; E_Int ind; E_Int cellNv;
  struct facet* face; struct facet* face1; struct facet* face2; 
  E_Int ret; 
  E_Float* cellN = f1->begin(posCellN+1);
  vector<struct facet*> faces;
  
  for (E_Int i = 0; i < ne; i++)
  {
    cellNi = cellN[i];
    for (size_t n = 0; n < cEEN[i].size(); n++)
    {
      ind = cEEN[i][n];
      cellNv = cellN[ind];
      if ((cellNi == 1 && cellNv == 0) || (cellNi == 0 && cellNv == 1))
      {
        //printf("%d %d\n", cellNi, cellNv);
        // face entre deux cellules cellN=0 et cellN=1
        E_Int noFace = commonFace[i][n];
        face = createFace(i, noFace, *cn, ff, fx, fy, fz);
        faces.push_back(face);
        // faces horizontales
        if (cellN[i] == 1) 
        { face = createFace(i, 0, *cn, ff, fx, fy, fz); faces.push_back(face); }
        if (cellN[i] == 1) 
        { face = createFace(i, 5, *cn, ff, fx, fy, fz); faces.push_back(face); }
        if (cellN[ind] == 1)
        { face = createFace(ind, 0, *cn, ff, fx, fy, fz); faces.push_back(face); }
      }
    }
  }
  
  // calcul des intersections des faces deux a deux
  for (size_t i = 0; i < faces.size(); i++)
  {
    face1 = faces[i];
    for (size_t j = 0; j < faces.size(); j++)
    {
      face2 = faces[j];
      if (i != j)
      {
        ret = intersect(face1, face2);
        if (ret == 1) { cellN[face1->cell] = 0; cellN[face2->cell] = 0; } // disymetrique volontaire
      }
    }
  }
  
  // nettoyage des facettes
  for (size_t i = 0; i < faces.size(); i++) delete faces[i];
  
  RELEASESHAREDB(res, array, f, cn);
  RELEASESHAREDB(res1, cellNObject, f1, cn1);
  
  Py_INCREF(Py_None);
  return Py_None;
}


//=============================================================================
/* 
   blank first. 
   Prend les faces opposee de chaque cellule, verifie si elles
   s'intersectent, si oui, blank
   IN: hexa mesh + cellN in center
   OUT: cellN modifie in-place
 */
//=============================================================================
PyObject* K_GENERATOR::blankFirst(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* cellNObject;
  if (!PyArg_ParseTuple(args, "OO", &array, &cellNObject))
  {
    return NULL;
  }

  // Check array: must be HEXA mesh
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray2(array, varString, f, ni, nj, nk, 
                                     cn, eltType);
  if (res != 2 || strcmp(eltType, "HEXA") != 0) 
  {
    if (res != 0) RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "blankFist: a must be HEXA.");
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError, 
                    "blankFirst: coords must be present in a1.");
    return NULL;       
  }
  E_Float* fx = f->begin(posx+1);
  E_Float* fy = f->begin(posy+1);
  E_Float* fz = f->begin(posz+1);

  // Check cellN (like array but in center)
  E_Int ni1, nj1, nk1;
  FldArrayF* f1; FldArrayI* cn1;
  char* varString1; char* eltType1;
  E_Int res1 = K_ARRAY::getFromArray2(cellNObject, varString1, f1, ni1, nj1, nk1, 
                                      cn1, eltType1);
  
  E_Int posCellN = K_ARRAY::isCellNatureField1Present(varString1);
  if (posCellN == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    RELEASESHAREDB(res1, cellNObject, f1, cn1);
    PyErr_SetString(PyExc_TypeError,
                    "blankFirst: cellN must be present in second array.");
    return NULL;
  }
  E_Int ne = cn->getSize();
  E_Float* cellN = f1->begin(posCellN+1);
  
  E_Int ff[24];
  ff[0 + 0*6] = 1; ff[0 + 1*6] = 4; ff[0 + 2*6] = 3; ff[0 + 3*6] = 2;
  ff[1 + 0*6] = 1; ff[1 + 1*6] = 2; ff[1 + 2*6] = 6; ff[1 + 3*6] = 5;
  ff[2 + 0*6] = 2; ff[2 + 1*6] = 3; ff[2 + 2*6] = 7; ff[2 + 3*6] = 6;
  ff[3 + 0*6] = 3; ff[3 + 1*6] = 4; ff[3 + 2*6] = 8; ff[3 + 3*6] = 7;
  ff[4 + 0*6] = 1; ff[4 + 1*6] = 5; ff[4 + 2*6] = 8; ff[4 + 3*6] = 4;
  ff[5 + 0*6] = 5; ff[5 + 1*6] = 6; ff[5 + 2*6] = 7; ff[5 + 3*6] = 8;
  
  facet* face1; facet* face2; E_Int ret;
  for (E_Int e = 0; e < ne; e++)
  {
    // faces opposees 1/3 et 2/4
    face1 = createFace(e, 0, *cn, ff, fx, fy, fz);  
    face2 = createFace(e, 3, *cn, ff, fx, fy, fz);
    ret = intersect(face1, face2);
    if (ret == 1) cellN[e] = 0;
    face1 = createFace(e, 2, *cn, ff, fx, fy, fz);  
    face2 = createFace(e, 4, *cn, ff, fx, fy, fz);
    ret = intersect(face1, face2);
    if (ret == 1) cellN[e] = 0;
    delete face1; delete face2;
  }

  RELEASESHAREDB(res, array, f, cn);
  RELEASESHAREDB(res1, cellNObject, f1, cn1);
  
  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
/* 
   blank Ext. 
   Prend les faces externes, verifie si elles
   s'intersectent, si oui, blank
   IN: hexa mesh + cellN in center
   OUT: cellN modifie in-place
 */
//=============================================================================
PyObject* K_GENERATOR::blankExt(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* cellNObject;
  if (!PyArg_ParseTuple(args, "OO", &array, &cellNObject))
  {
    return NULL;
  }

  // Check array: must be HEXA mesh
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray2(array, varString, f, ni, nj, nk, 
                                     cn, eltType);
  if (res != 2 || strcmp(eltType, "HEXA") != 0) 
  {
    if (res != 0) RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "blankExt: a must be HEXA.");
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError, 
                    "blankExt: coords must be present in a1.");
    return NULL;       
  }
  E_Float* fx = f->begin(posx+1);
  E_Float* fy = f->begin(posy+1);
  E_Float* fz = f->begin(posz+1);

  // Check cellN (like array but in center)
  E_Int ni1, nj1, nk1;
  FldArrayF* f1; FldArrayI* cn1;
  char* varString1; char* eltType1;
  E_Int res1 = K_ARRAY::getFromArray2(cellNObject, varString1, f1, ni1, nj1, nk1, 
                                      cn1, eltType1);
  
  E_Int posCellN = K_ARRAY::isCellNatureField1Present(varString1);
  if (posCellN == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    RELEASESHAREDB(res1, cellNObject, f1, cn1);
    PyErr_SetString(PyExc_TypeError,
                    "blankExt: cellN must be present in second array.");
    return NULL;
  }
  E_Int ne = cn->getSize();
  E_Float* cellN = f1->begin(posCellN+1);
  
  E_Int ff[24];
  ff[0 + 0*6] = 1; ff[0 + 1*6] = 4; ff[0 + 2*6] = 3; ff[0 + 3*6] = 2;
  ff[1 + 0*6] = 1; ff[1 + 1*6] = 2; ff[1 + 2*6] = 6; ff[1 + 3*6] = 5;
  ff[2 + 0*6] = 2; ff[2 + 1*6] = 3; ff[2 + 2*6] = 7; ff[2 + 3*6] = 6;
  ff[3 + 0*6] = 3; ff[3 + 1*6] = 4; ff[3 + 2*6] = 8; ff[3 + 3*6] = 7;
  ff[4 + 0*6] = 1; ff[4 + 1*6] = 5; ff[4 + 2*6] = 8; ff[4 + 3*6] = 4;
  ff[5 + 0*6] = 5; ff[5 + 1*6] = 6; ff[5 + 2*6] = 7; ff[5 + 3*6] = 8;
  
  facet* face; facet* face1; facet* face2; E_Int ret;
  vector<struct facet*> faces;

  for (E_Int e = 0; e < ne; e++)
  {
    // faces externe
    if (cellN[e] == 1) 
    {
      face = createFace(e, 5, *cn, ff, fx, fy, fz); 
      faces.push_back(face); 
    }
  }
  
   // calcul des intersections des faces deux a deux
  for (size_t i = 0; i < faces.size(); i++)
  {
    face1 = faces[i];
    for (size_t j = 0; j < faces.size(); j++)
    {
      face2 = faces[j];
      if (i != j)
      {
        ret = intersect(face1, face2);
        if (ret == 1) { cellN[face1->cell] = 0; cellN[face2->cell] = 0; } // disymetrique volontaire
      }
    }
  }
  
  // nettoyage des facettes
  for (size_t i = 0; i < faces.size(); i++) delete faces[i];

  RELEASESHAREDB(res, array, f, cn);
  RELEASESHAREDB(res1, cellNObject, f1, cn1);
  
  Py_INCREF(Py_None);
  return Py_None;
}



//=============================================================================
/* 
   blank prev.
   Prend le bloc des couches precedentes, trouve les faces limites cellN=0, cellN=1
   Prend la couche actuelle, trouve les faces limites cellN=0, cellN=1
   Blank les cellules du bloc courant qui intersectent celles des couches precedentes
   IN: hexa mesh + cellN in center couche courante et couches precedentes
   OUT: cellN de couche courante modifie in-place
 */
//=============================================================================
PyObject* K_GENERATOR::blankPrev(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* cellNObject;
  PyObject* arrayPrev;
  PyObject* cellNObjectPrev;
  if (!PyArg_ParseTuple(args, "OOOO", &array, &cellNObject, &arrayPrev, &cellNObjectPrev))
  {
    return NULL;
  }

  // Check array: must be HEXA mesh (couche courante)
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray2(array, varString, f, ni, nj, nk, 
                                     cn, eltType);
  if (res != 2 || strcmp(eltType, "HEXA") != 0) 
  {
    if (res != 0) RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "blankPrev: a must be HEXA.");
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError, 
                    "blankPrev: coords must be present in a1.");
    return NULL;       
  }
  E_Float* fx = f->begin(posx+1);
  E_Float* fy = f->begin(posy+1);
  E_Float* fz = f->begin(posz+1);

  // Check cellN (like array but in center)
  E_Int ni1, nj1, nk1;
  FldArrayF* f1; FldArrayI* cn1;
  char* varString1; char* eltType1;
  E_Int res1 = K_ARRAY::getFromArray2(cellNObject, varString1, f1, ni1, nj1, nk1, 
                                      cn1, eltType1);
  
  E_Int posCellN = K_ARRAY::isCellNatureField1Present(varString1);
  if (posCellN == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    RELEASESHAREDB(res1, cellNObject, f1, cn1);
    PyErr_SetString(PyExc_TypeError,
                    "blankPrev: cellN must be present in second array.");
    return NULL;
  }
  E_Float* cellN = f1->begin(posCellN+1);

  // Check arrayPrev: must be HEXA mesh (couches precedentes)
  E_Int nip, njp, nkp;
  FldArrayF* fp; FldArrayI* cnp;
  char* varStringp; char* eltTypep;
  E_Int resp = K_ARRAY::getFromArray2(arrayPrev, varStringp, fp, nip, njp, nkp, 
                                      cnp, eltTypep);
  if (resp != 2 || strcmp(eltTypep, "HEXA") != 0) 
  {
    if (resp != 0)
    { 
      RELEASESHAREDB(res, array, f, cn);
      RELEASESHAREDB(res1, cellNObject, f1, cn1);
      RELEASESHAREDB(resp, arrayPrev, fp, cnp);
    }
    PyErr_SetString(PyExc_TypeError,
                    "blankPrev: a must be HEXA.");
    return NULL;
  }
  posx = K_ARRAY::isCoordinateXPresent(varStringp);
  posy = K_ARRAY::isCoordinateYPresent(varStringp);
  posz = K_ARRAY::isCoordinateZPresent(varStringp);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    RELEASESHAREDB(res1, cellNObject, f1, cn1);
    RELEASESHAREDB(resp, arrayPrev, fp, cnp);
    PyErr_SetString(PyExc_TypeError, 
                    "blankPrev: coords must be present in a1.");
    return NULL;       
  }
  E_Float* fxp = fp->begin(posx+1);
  E_Float* fyp = fp->begin(posy+1);
  E_Float* fzp = fp->begin(posz+1);

  // Check cellN (like array but in center)
  E_Int nip1, njp1, nkp1;
  FldArrayF* fp1; FldArrayI* cnp1;
  char* varStringp1; char* eltTypep1;
  E_Int resp1 = K_ARRAY::getFromArray2(cellNObjectPrev, varStringp1, fp1, nip1, njp1, nkp1, 
                                       cnp1, eltTypep1);
  
  posCellN = K_ARRAY::isCellNatureField1Present(varStringp1);
  if (posCellN == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    RELEASESHAREDB(res1, cellNObject, f1, cn1);
    RELEASESHAREDB(resp, arrayPrev, fp, cnp);
    RELEASESHAREDB(resp1, cellNObjectPrev, fp1, cnp1);

    PyErr_SetString(PyExc_TypeError,
                    "blankPrev: cellN must be present in second array.");
    return NULL;
  }
  E_Float* cellNp = fp1->begin(posCellN+1);
  
  E_Int ff[24];
  ff[0 + 0*6] = 1; ff[0 + 1*6] = 4; ff[0 + 2*6] = 3; ff[0 + 3*6] = 2;
  ff[1 + 0*6] = 1; ff[1 + 1*6] = 2; ff[1 + 2*6] = 6; ff[1 + 3*6] = 5;
  ff[2 + 0*6] = 2; ff[2 + 1*6] = 3; ff[2 + 2*6] = 7; ff[2 + 3*6] = 6;
  ff[3 + 0*6] = 3; ff[3 + 1*6] = 4; ff[3 + 2*6] = 8; ff[3 + 3*6] = 7;
  ff[4 + 0*6] = 1; ff[4 + 1*6] = 5; ff[4 + 2*6] = 8; ff[4 + 3*6] = 4;
  ff[5 + 0*6] = 5; ff[5 + 1*6] = 6; ff[5 + 2*6] = 7; ff[5 + 3*6] = 8;

  // Construction de la connectivite elt/elt voisins
  E_Int ne = cn->getSize();
  E_Int nv = f->getSize();
  vector< vector<E_Int> > cEEN(ne);
  vector< vector<E_Int> > commonFace(ne);
  K_CONNECT::connectEV2EENbrs(eltType, nv, *cn, cEEN, commonFace);
      
  // Pool des faces limites de array
  E_Int cellNi; E_Int ind; E_Int cellNv;
  struct facet* face; struct facet* face1; struct facet* face2; 
  E_Int ret; 
  vector<struct facet*> faces;
  
  for (E_Int i = 0; i < ne; i++)
  {
    cellNi = cellN[i];
    for (size_t n = 0; n < cEEN[i].size(); n++)
    {
      ind = cEEN[i][n];
      cellNv = cellN[ind];
      if ((cellNi == 1 && cellNv == 0) || (cellNi == 0 && cellNv == 1))
      {
        // face entre deux cellules cellN=0 et cellN=1
        E_Int noFace = commonFace[i][n];
        face = createFace(i, noFace, *cn, ff, fx, fy, fz);
        faces.push_back(face);
        // faces externes aussi a ajouter
      }
    }
  }
     
  // Construction de la connectivite elt/elt voisins
  ne = cnp->getSize();
  nv = fp->getSize();
  cEEN.clear(); cEEN.resize(ne);
  commonFace.clear(); commonFace.resize(ne);
  K_CONNECT::connectEV2EENbrs(eltTypep, nv, *cnp, cEEN, commonFace); 
  
  // Pool des faces limites de arrayPrev
  vector<struct facet*> facesPrev;
  
  for (E_Int i = 0; i < ne; i++)
  {
    cellNi = cellNp[i];
    for (size_t n = 0; n < cEEN[i].size(); n++)
    {
      ind = cEEN[i][n];
      cellNv = cellNp[ind];
      if ((cellNi == 1 && cellNv == 0) || (cellNi == 0 && cellNv == 1))
      {
        // face entre deux cellules cellN=0 et cellN=1
        E_Int noFace = commonFace[i][n];
        face = createFace(i, noFace, *cnp, ff, fxp, fyp, fzp);
        facesPrev.push_back(face);
      }
    }
  }

  // calcul des intersections des faces deux a deux
  for (size_t i = 0; i < faces.size(); i++)
  {
    face1 = faces[i];
    for (size_t j = 0; j < facesPrev.size(); j++)
    {
      face2 = facesPrev[j];
      ret = intersect(face1, face2);
      if (ret == 1) { cellN[face1->cell] = 0; cellN[face2->cell] = 0; break; } 
    }
  }
  
  // nettoyage des facettes
  for (size_t i = 0; i < faces.size(); i++) delete faces[i];
  for (size_t i = 0; i < facesPrev.size(); i++) delete facesPrev[i];


  RELEASESHAREDB(res, array, f, cn);
  RELEASESHAREDB(res1, cellNObject, f1, cn1);
  RELEASESHAREDB(resp, arrayPrev, fp, cnp);
  RELEASESHAREDB(resp1, cellNObjectPrev, fp1, cnp1);
  
  Py_INCREF(Py_None);
  return Py_None;
}
