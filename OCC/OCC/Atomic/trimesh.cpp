/*    
    Copyright 2013-2023 Onera.

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

#include "occ.h"
#include "OCCSurface.h"
#include "TopoDS.hxx"
#include "TopoDS_Face.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "Nuga/include/SurfaceMesher.h"

// trimesh : mesh with TRI on a CAD patch
// IN: contour UV, no de la face dans la CAD
PyObject* K_OCC::trimesh(PyObject* self, PyObject* args)
{
    PyObject* hook;
    E_Int faceNo; // no de la face
    PyObject* arrayUV;
    E_Float hmax, hausd;
    if (!PYPARSETUPLE(args, "OOldd", "OOidd", "OOlff", "OOiff", 
                      &hook, &arrayUV, &faceNo, &hmax, &hausd)) return NULL;  
    void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
    packet = (void**) PyCObject_AsVoidPtr(hook);
#else
    packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif
    TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];
    TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];
  
    // Cree la OCCSurface
    const TopoDS_Face& F = TopoDS::Face(surfaces(faceNo));
    const OCCSurface& occ_surf = K_OCC::OCCSurface(F, edges, 0);
    
    
    // Get from array
    FldArrayF* fi; E_Int ni, nj, nk;
    char* varString; FldArrayI* ci; char* eltType;
    E_Int ret = K_ARRAY::getFromArray2(arrayUV, varString, fi, ni, nj, nk, ci, eltType);
    if (ret != 2)
    {
        PyErr_SetString(PyExc_TypeError, "trimesh: invalid array.");
        return NULL;
    }

    // Cree le mailleur
    DELAUNAY::SurfaceMesher<OCCSurface> mesher;

    // Recuperation des donnees
    E_Int n = fi->getSize();
    K_FLD::FloatArray pos3D(3, n); // pos3D: les coords reelles
    for (E_Int i = 0; i < n; i++) pos3D(0,i) = (*fi)(i,1);
    for (E_Int i = 0; i < n; i++) pos3D(1,i) = (*fi)(i,2);
    for (E_Int i = 0; i < n; i++) pos3D(2,i) = (*fi)(i,3);
    
    K_FLD::FloatArray UVcontour(2, n);
    for (E_Int i = 0; i < n; i++) UVcontour(0,i) = (*fi)(i,4);
    for (E_Int i = 0; i < n; i++) UVcontour(1,i) = (*fi)(i,5);
    
    K_FLD::IntArray connectB(*ci); // connectivite
    E_Int ne = ci->getSize();
    for (E_Int i = 0; i < ne; i++) connectB(0,i) -= 1;
    for (E_Int i = 0; i < ne; i++) connectB(1,i) -= 1;

    /*
    printf("UVContour\n");
    printf("rows=%d cols=%d\n", UVcontour.rows(), UVcontour.cols());
    for (E_Int j = 0; j < UVcontour.cols(); j++)
    printf("%d = %f %f\n", j, UVcontour(0,j), UVcontour(1,j));

    printf("pos3D\n");
    printf("rows=%d cols=%d\n", pos3D.rows(), pos3D.cols());
    for (E_Int j = 0; j < pos3D.cols(); j++)
    printf("%d = %f %f %f\n", j, pos3D(0,j), pos3D(1,j), pos3D(2,j));

    printf("connectB\n");
    printf("rows=%d cols=%d\n", connectB.rows(), connectB.cols());
    for (E_Int j = 0; j < connectB.cols(); j++)
    printf("%d = %d %d\n", j, connectB(0,j), connectB(1,j));
    */
    DELAUNAY::SurfaceMeshData<OCCSurface> data(UVcontour, pos3D, connectB, occ_surf);

    DELAUNAY::SurfaceMesherMode mode;
    mode.chordal_error = hausd; // chordal error set
    //if (aniso) mode.metric_mode = mode.ANISO;
    //else mode.metric_mode = mode.ISO_CST;
    // ANISO mode ne marche pas!!!!!!!!!!
    mode.metric_mode = mode.ISO_CST;
    mode.symmetrize = false;
    mode.hmax = hmax; // h moyen
    mode.growth_ratio = -1;
    //mode.ignore_coincident_nodes = true; // pour bypasser les pbs d'insertion 
    mesher.mode = mode;

    E_Int err = 0;
    err = mesher.run(data);
    if (err || (data.connectM.cols() == 0))
    {
        // connectM doit etre la sortie
        RELEASESHAREDB(ret, arrayUV, fi, ci);
        PyErr_SetString(PyExc_TypeError, "trimesh: mesher has failed.");
        return NULL;
    }
    
    // recupere la sortie    
    FldArrayF* coords = new FldArrayF;
    data.pos3D.convert(*coords);
    FldArrayI* cn = new FldArrayI;
    data.connectM.convert(*cn, 1);
    PyObject* tpl = K_ARRAY::buildArray(*coords, "x,y,z", *cn, -1, "TRI");
    delete coords; delete cn;

    RELEASESHAREDB(ret, arrayUV, fi, ci);
    return tpl;
}
