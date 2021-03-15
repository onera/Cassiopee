/*    
    Copyright 2013-2021 Onera.

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
    const OCCSurface& occ_surf = K_OCC::OCCSurface(F, edges, faceNo-1);
    
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
    K_FLD::FloatArray pos3D(3, n); // pos3D: les coords reelles?
    for (E_Int i = 0; i < n; i++) pos3D(0,i) = (*fi)(i,1);
    for (E_Int i = 0; i < n; i++) pos3D(1,i) = (*fi)(i,2);
    for (E_Int i = 0; i < n; i++) pos3D(2,i) = (*fi)(i,3);
    
    K_FLD::FloatArray UVcontour(2, n);
    for (E_Int i = 0; i < n; i++) UVcontour(0,i) = (*fi)(i,4);
    for (E_Int i = 0; i < n; i++) UVcontour(1,i) = (*fi)(i,5);
    
    K_FLD::IntArray connectB(*ci); // connectivite
    printf("connect %d %d\n", connectB.rows(), connectB.cols());

    DELAUNAY::SurfaceMeshData<OCCSurface> data(UVcontour, pos3D, connectB, occ_surf);

    DELAUNAY::SurfaceMesherMode mode;
    mode.chordal_error = hausd; // chordal error set
    //if (aniso) mode.metric_mode = mode.ANISO;
    //else mode.metric_mode = mode.ISO

    // old mode (mais c'est le meilleur)
    mode.symmetrize = false;
    mode.hmax = hmax; // h moyen
    mode.growth_ratio = 0.;
    //mode.ignore_coincident_nodes = true; // pour bypasser les pb d'insertion 
    mesher.mode = mode;

    E_Int err = 0;
    err = mesher.run(data);
    if (err || (data.connectM.cols() == 0))
    {
        // connectM doit etre la sortie
        RELEASESHAREDB(ret, arrayUV, fi, ci);
        PyErr_SetString(PyExc_TypeError, "trimesh: meshes has failed.");
        return NULL;
    }
    
    // recupere la sortie
    //FldArrayI cn(data.connectM);
    //FldArrayF coords(data.pos3D);
    
    RELEASESHAREDB(ret, arrayUV, fi, ci);
    Py_INCREF(Py_None);
    return Py_None;
}
