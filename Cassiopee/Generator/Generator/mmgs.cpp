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

// Remaillage surfacique avec mmgs

#include "generator.h"
#include "MMGS/libmmgs.h"

using namespace std;
using namespace K_FLD;
using namespace K_FUNC; 

// ============================================================================
/* MMGS
   IN: maillage TRI
   IN: eventuellement metric ou solution
   IN: 
   IN: 
   OUT: maillage TRI remaille. */
// ============================================================================
PyObject* K_GENERATOR::mmgs(PyObject* self, PyObject* args)
{
  E_Float ridgeAngle=45.; // angle detection
  E_Float hmin = 0.;
  E_Float hmax = 0.;
  E_Float hausd = 0.;
  E_Float hgrad = -1.;
  E_Int anisotropy = 0;
  E_Int optim = 0;
  PyObject* array;
  PyObject* fixedNodes; // imposed nodes
  PyObject* fixedEdges; // imposed edge
  if (!PYPARSETUPLE_(args, O_ RRRR_ R_ II_ OO_,
                    &array, &ridgeAngle, &hmin, &hmax, &hausd, &hgrad, 
                    &anisotropy, &optim, &fixedNodes, &fixedEdges))
  {
    return NULL;
  }
  
  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, 
                                     cn, eltType);
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "mmgs: invalid array.");
    return NULL;
  }

  if (strcmp(eltType, "TRI") != 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "mmgs: array must be TRI.");
    RELEASESHAREDB(res, array, f, cn);
    return NULL; 
  }

  // Check posx,posy,posz
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "mmgs: coordinates not found in array.");
    RELEASESHAREDB(res, array, f, cn);
    return NULL;
  }
  posx++; posy++; posz++;

  // Check sizemap field (if any) -> isotropic remeshing
  E_Int posf = K_ARRAY::isNamePresent("sizemap", varString);
  posf++;

  // check sizemap tensor (if any) -> anistropic remeshing
  E_Int posf11 = K_ARRAY::isNamePresent("sizemap11", varString);
  posf11++;
  E_Int posf12 = K_ARRAY::isNamePresent("sizemap12", varString);
  posf12++;
  E_Int posf13 = K_ARRAY::isNamePresent("sizemap13", varString);
  posf13++;
  E_Int posf22 = K_ARRAY::isNamePresent("sizemap22", varString);
  posf22++;
  E_Int posf23 = K_ARRAY::isNamePresent("sizemap23", varString);
  posf23++;
  E_Int posf33 = K_ARRAY::isNamePresent("sizemap33", varString);
  posf33++;

  /* update Info */
  MMG5_Info info;
  info.imprim = -99;
  info.ddebug = 0;
  info.mem = -1;
  
  /* ridge angle */
  info.dhd = ridgeAngle;
  info.dhd = max(0.0, min(180.0,info.dhd));
  info.dhd = cos(info.dhd*3.14159265359/180.0);
  
  info.hmax = hmax;
  info.hmin = hmin;
  info.hausd = hausd;
  if (hgrad < 0.0) hgrad = -1;
  else hgrad = log(hgrad);
  info.hgrad = hgrad;

  // Structure pour stocker le maillage et la metrique
  MMG5_pMesh mesh;
  MMG5_pSol  metric;

  //====================================================================
  // Comment faire tout en memoire (d'elephant)
  //====================================================================
  
  mesh = NULL;
  metric = NULL;
  MMGS_Init_mesh(MMG5_ARG_start,
                 MMG5_ARG_ppMesh, &mesh, MMG5_ARG_ppMet, &metric,
                 MMG5_ARG_end);

  E_Int nargs; E_Int argc; char** vals; char** argv;

  if (optim == 1) // optimisation du maillage, pas de settings
  {
    nargs = 0; vals = NULL;
    argc = 1;
    argv = new char* [argc];
    argv[0] = new char [20]; strcpy(argv[0], "mmgs_O3");
    //argv[1] = new char [20]; strcpy(argv[1], "-optim");
  }
  else if (posf == 0 && 
           posf11 == 0 && posf12 == 0 && posf13 == 0 &&
           posf22 == 0 && posf23 == 0 && posf33 == 0) // pas de metric field, hmin,hmax,hausd,hgrad,ridgeAngle can be set
  {
    nargs = 5;
    vals = new char* [nargs];
    vals[0] = new char [20]; sprintf(vals[0], "%g", hausd);
    vals[1] = new char [20]; sprintf(vals[1], "%g", ridgeAngle);
    vals[2] = new char [20]; sprintf(vals[2], "%g", hmin);
    vals[3] = new char [20]; sprintf(vals[3], "%g", hmax);
    vals[4] = new char [20]; sprintf(vals[4], "%g", hgrad);
  
    argc = 2*nargs+1;
    if (anisotropy == 1) argc += 1;
    argv = new char* [argc];
    argv[0] = new char [20]; strcpy(argv[0], "mmgs_O3");
    argv[1] = new char [20]; strcpy(argv[1], "-hausd");
    argv[2] = new char [20]; strcpy(argv[2], vals[0]);
    argv[3] = new char [20]; strcpy(argv[3], "-ar");
    argv[4] = new char [20]; strcpy(argv[4], vals[1]);
    argv[5] = new char [20]; strcpy(argv[5], "-hmin");
    argv[6] = new char [20]; strcpy(argv[6], vals[2]);
    argv[7] = new char [20]; strcpy(argv[7], "-hmax");
    argv[8] = new char [20]; strcpy(argv[8], vals[3]);
    argv[9] = new char [20]; strcpy(argv[9], "-hgrad");
    argv[10] = new char [20]; strcpy(argv[10], vals[4]);
    if (anisotropy == 1) { argv[11] = new char [20]; strcpy(argv[11], "-A"); }
    printf("INFO: hausd=%s hmin=%s hmax=%s hgrad=%s\n", vals[0], vals[2], vals[3], vals[4]);
  }
  else // metric field, ridgeAngle, hgrad can be set
  {
    nargs = 2;
    vals = new char* [nargs];
    vals[0] = new char [20]; sprintf(vals[0], "%g", ridgeAngle);
    vals[1] = new char [20]; sprintf(vals[1], "%g", hgrad);
    argc = 2*nargs+1;
    if (anisotropy == 1) argc += 1;
    argv = new char* [argc];
    argv[0] = new char [20]; strcpy(argv[0], "mmgs_O3");
    argv[1] = new char [20]; strcpy(argv[1], "-ar");
    argv[2] = new char [20]; strcpy(argv[2], vals[0]);
    argv[3] = new char [20]; strcpy(argv[3], "-hgrad");
    argv[4] = new char [20]; strcpy(argv[4], vals[1]);
    if (anisotropy == 1) { argv[5] = new char [20]; strcpy(argv[5], "-A"); }
  }
  
  MMGS_parsar(argc, argv, mesh, metric, NULL);

  // Conversion de l'array vers mesh
  // dimensionne nbre de vertex, nbre de triangles, nbre d'edges
  E_Int ne = cn->getSize(); // nbre d'elements
  E_Int np = f->getSize(); // nbre de vertex
  E_Int na = 0; // nbre d'edges
  if (fixedEdges != Py_None)
  {
    // we add imposed edges in the mesh -> increase na
    E_Int ns = PyList_Size(fixedEdges);
    for (E_Int l = 0; l < ns; l++)
    {
      PyObject* o = PyList_GetItem(fixedEdges, l); // BAR connect numpy
      E_Int* ptr; E_Int nelts; E_Int nfld;
      K_NUMPY::getFromNumpyArray(o, ptr, nelts, nfld);
      na += nelts;
      Py_DECREF(o);
    }
  }
  MMGS_Set_meshSize(mesh, np, ne, na);

  E_Float* fx = f->begin(posx);
  E_Float* fy = f->begin(posy);
  E_Float* fz = f->begin(posz);
  for (E_Int i = 0; i < np; i++)
    MMGS_Set_vertex(mesh, fx[i], fy[i], fz[i], 0, i+1);

  E_Int* c1 = cn->begin(1);
  E_Int* c2 = cn->begin(2);
  E_Int* c3 = cn->begin(3);
  E_Int stride = cn->getStride(); 
  for (E_Int i = 0; i < ne; i++)
    MMGS_Set_triangle(mesh, c1[i*stride], c2[i*stride], c3[i*stride], 3, i+1);

  if (posf >= 1)
  {
    E_Float* ff = f->begin(posf);
    MMGS_Set_solSize(mesh, metric, MMG5_Vertex, np, MMG5_Scalar);
    for (E_Int k = 1; k <= np; k++)
    {
      MMGS_Set_scalarSol(metric, ff[k-1], k); // existe aussi set_TensorSol
    }
  }

  if (posf11 >=1 && posf12 && posf13 >=1 && posf22 >=1 && posf23>=1 && posf33>=1)
  {
    E_Float* ff11 = f->begin(posf11);
    E_Float* ff12 = f->begin(posf12);
    E_Float* ff13 = f->begin(posf13);
    E_Float* ff22 = f->begin(posf22);
    E_Float* ff23 = f->begin(posf23);
    E_Float* ff33 = f->begin(posf33);
    MMGS_Set_solSize(mesh, metric, MMG5_Vertex, np, MMG5_Tensor);
    for (E_Int k = 1; k <= np; k++)
    {
      MMGS_Set_tensorSol(metric, ff11[k-1], ff12[k-1], ff13[k-1],
        ff22[k-1], ff23[k-1], ff33[k-1], k);
    }
  }

  // Use fixedNodes if any
  if (fixedNodes != Py_None)
  {
    E_Int ns = PyList_Size(fixedNodes);
    for (E_Int l = 0; l < ns; l++)
    {
      PyObject* o = PyList_GetItem(fixedNodes, l);
      E_Int* ptr; E_Int npts; E_Int nfld;
      K_NUMPY::getFromNumpyArray(o, ptr, npts, nfld);
      for (E_Int i = 0; i < npts; i++) MMGS_Set_requiredVertex(mesh, ptr[i]); // +1?
      Py_DECREF(o);
    }
  }

  // Used fixedEdges if any (through a BAR connectivity)
  if (fixedEdges != Py_None)
  {
    E_Int ns = PyList_Size(fixedEdges);
    for (E_Int l = 0; l < ns; l++)
    {
      PyObject* o = PyList_GetItem(fixedEdges, l); // BAR connect numpy
      E_Int* ptr; E_Int nelts; E_Int nfld;
      K_NUMPY::getFromNumpyArray(o, ptr, nelts, nfld);
      for (E_Int i = 0; i < nelts; i++)
      {
        E_Int v0 = ptr[i];
        E_Int v1 = ptr[i+nelts];
        MMGS_Set_edge(mesh, v0, v1, 2, i+1); // define edge i+1
        MMGS_Set_requiredEdge(mesh, i+1);
      }
      Py_DECREF(o);
    }
  }

  // Remesh
  MMGS_mmgslib(mesh, metric);

  // Save the mesh
  //MMGS_saveMesh(mesh, "mineOut.mesh");

  // 1) Manually get the mesh 
  // get the size of the mesh: vertices, tetra, triangles, edges
  E_Int npo, nto, nao; // flaws
  MMGS_Get_meshSize(mesh, &npo, &nto, &nao);
  printf("INFO: output mesh has " SF_D_ " vertices, " SF_D_ " triangles.\n", npo, nto);

  // Table to know if a vertex is corner
  int* corner = new int [npo+1];
  
  // Table to know if a vertex/tetra/tria/edge is required
  E_Int size = std::max(npo, nto); size = std::max(size, nao);
  int* required = new int [ size+1 ];
  
  // Table to know if a component is corner and/or required
  E_Int* ridge = new E_Int [nao+1];

  // Allocate array2
  PyObject* o = K_ARRAY::buildArray3(3, "x,y,z", npo, nto, "TRI", false, f->getApi());
  FldArrayF* fo; FldArrayI* co;
  K_ARRAY::getFromArray3(o, fo, co);
  E_Float* fox = fo->begin(1);
  E_Float* foy = fo->begin(2);
  E_Float* foz = fo->begin(3);
  E_Int* co1 = co->begin(1);
  E_Int* co2 = co->begin(2);
  E_Int* co3 = co->begin(3);
  stride = co->getStride();

  E_Int ref; E_Float px,py,pz;
  for (E_Int k=1; k <= npo; k++) 
  {
    MMGS_Get_vertex(mesh, &px, &py, &pz,
                    &ref,&(corner[k]),&(required[k]));
    fox[k-1] = px;
    foy[k-1] = py;
    foz[k-1] = pz;
  }

  E_Int ind1,ind2,ind3;
  for (E_Int k=1; k <= nto; k++) 
  {
    MMGS_Get_triangle(mesh, &ind1,&ind2,&ind3,
                      &ref, &(required[k]));
    co1[(k-1)*stride] = ind1;
    co2[(k-1)*stride] = ind2;
    co3[(k-1)*stride] = ind3;
    //printf("%d %d %d\n", ind1,ind2,ind3);
  }
  
  // Free all
  for (E_Int i = 0; i < nargs; i++) delete [] vals[i];
  delete [] vals;
  for (E_Int i = 0; i < argc; i++) delete [] argv[i];
  delete [] argv;

  delete [] corner;
  delete [] required;
  delete [] ridge;

  MMGS_Free_all(MMG5_ARG_start,
                MMG5_ARG_ppMesh, &mesh, MMG5_ARG_ppMet, &metric,
                MMG5_ARG_end);

  /* Analysis */
  //_MMG5_analys(&mesh);
  
  RELEASESHAREDB(res, array, f, cn);
  RELEASESHAREDU(o, fo, co);

  return o;
}
