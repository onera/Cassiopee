/*    
    Copyright 2013-2024 Onera.

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
  E_Float hmin, hmax, hausd, grading;
  if (!PYPARSETUPLE_(args, OO_ I_ RRRR_, 
                    &hook, &arrayUV, &faceNo, &hmin, &hmax, &hausd, &grading)) return NULL;  
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

  // Return null if input edges doesnt contain any bar
  if (ne == 0)
  {
    printf("Error: trimesh: input edge has " SF_D_ " nodes and " SF_D_ " bars.\n", n, ne);
    return NULL;
  }
  /*
  printf("UVContour\n");
  printf("rows=%d cols=%d\n", UVcontour.rows(), UVcontour.cols());
  for (E_Int j = 0; j < UVcontour.cols(); j++)
  printf("%d = " SF_F2_ "\n", j, UVcontour(0,j), UVcontour(1,j));

  printf("pos3D\n");
  printf("rows=%d cols=%d\n", pos3D.rows(), pos3D.cols());
  for (E_Int j = 0; j < pos3D.cols(); j++)
  printf("%d = " SF_F2_ " " SF_F_ "\n", j, pos3D(0,j), pos3D(1,j), pos3D(2,j));

  printf("connectB\n");
  printf("rows=%d cols=%d\n", connectB.rows(), connectB.cols());
  for (E_Int j = 0; j < connectB.cols(); j++)
  printf("%d = %d %d\n", j, connectB(0,j), connectB(1,j));
  */
  DELAUNAY::SurfaceMeshData<OCCSurface> data(UVcontour, pos3D, connectB, occ_surf);

  DELAUNAY::SurfaceMesherMode mode;

  E_Int aniso = false;

  if ( (hausd < 0 && hmax > 0) ||
       (hausd > 0 && std::abs(hmax-hmin) < 1.e-12) ) // iso hmax
  {
    // mode pure hmax
    E_Float dx = (hmax-hmin)/hmax;
    // uniform h
    mode.metric_mode = mode.ISO_CST;
    mode.hmax = hmax; // hmax
    mode.hmin = hmax; // hmax
    mode.chordal_error = 20000.; // not used
    mode.growth_ratio = 1.; // grading forced
    mode.nb_smooth_iter = 0; // iter de lissage de la metrique
    mode.symmetrize = false;
    if (dx > 0.2) mode.growth_ratio = 1.1;
    //printf("trimesh uniform hmin=" SF_F_ " hmax=" SF_F_ " grading=" SF_F_ "\n", mode.hmin, mode.hmax, mode.growth_ratio);      
  }
  else if (hausd > 0 && hmax > 0 && hmin >= 0 && aniso == true) // aniso mix
  {
    // mode pure hausd
    mode.metric_mode = mode.ANISO; //ISO_RHO impose la courbure minimum dans les deux directions
    mode.hmax = hmax; // h moyen
    mode.hmin = hmin; // h moyen
    mode.chordal_error = hausd; // chordal error set
    mode.growth_ratio = grading; // grading forced
    mode.nb_smooth_iter = 2; // iter de lissage pour assurer le grading
    mode.symmetrize = true;
  }
  else if (hausd > 0 && hmax > 0 && hmin >= 0 && aniso == false) // iso rho mix
  {
    // mode mix hmin/hmax/hausd
    mode.metric_mode = mode.ISO_RHO; //ISO_RHO impose la courbure minimum dans les deux directions;
    mode.hmax = hmax; // h max
    mode.hmin = hmin; // h min
    //mode.hmin = K_CONST::E_MAX_FLOAT; // hmin as in landier
    mode.chordal_error = hausd; // chordal error set
    //mode.growth_ratio = grading; // grading forced (pas coherent avec hausd?)
    mode.growth_ratio = 1.1; // grading ne sert pas si pas de lissage
    mode.nb_smooth_iter = 0; // iter de lissage de la metrique
    mode.symmetrize = false;
  }
  else 
  {
    PyErr_SetString(PyExc_TypeError, "trimesh: invalid hmin,hmax,hausd.");
    return NULL;
  }

  if (mode.metric_mode == mode.ISO_CST) printf("trimesh mode=ISO_CST ");
  else if (mode.metric_mode == mode.ISO_RHO) printf("trimesh mode=ISO_RHO ");
  else printf("trimesh mode=ANISO ");
  printf("hmin=%g hmax=%g hausd=%g sym=%d grading=%g  smooth=" SF_D_ "\n", 
          mode.hmin, mode.hmax, mode.chordal_error, 
          mode.symmetrize, mode.growth_ratio, mode.nb_smooth_iter);

  mode.metric_interpol_type = mode.LINEAR;
  //mode.metric_interpol_type = mode.GEOMETRIC;
    
  //mode.ignore_coincident_nodes = true; // pour bypasser les pbs d'insertion 
  mesher.mode = mode;

  E_Int err = 0;
  mesher.clear(); // landier
  err = mesher.run(data);
  if (err || (data.connectM.cols() == 0))
  {
      // connectM doit etre la sortie
      printf("error = " SF_D_ "\n", err);
      printf("cols = " SF_D_ "\n", data.connectM.cols());
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
