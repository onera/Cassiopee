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
// convertCAD2Arrays (algo=2)

//#define DEBUG_CAD_READER

#include <stdio.h>
#include <string.h>
#include "occ.h"
#include "CADviaOCC.h"

using namespace std;
using namespace K_FLD;

#include <iostream>

namespace K_OCC
{
E_Int CADread2
(char* file, char* fileFmt, E_Float h, E_Float chordal_err, E_Float gr, 
 E_Float merge_tol,
 vector<FldArrayF*>& unstructField,
 vector<FldArrayI*>& connect,
 vector<E_Int>& eltType,
 vector<char*>& zoneNames,
 bool do_join);
}

// ============================================================================
/* Convert CAD file to arrays using OpenCascade */
// ============================================================================
PyObject* K_OCC::convertCAD2Arrays2(PyObject* self, PyObject* args)
{
  char* fileName; char* fileFmt;
  E_Float h, chordal_err, gr(-1.), merge_tol(1.e-6);
  E_Int do_join(1);

  if (!PYPARSETUPLE_(args, SS_ RRRR_ I_, &fileName, &fileFmt, &h, &chordal_err,
                    &gr, &merge_tol, &do_join)) return NULL;

  // Check recognised formats
  if (K_STRING::cmp(fileFmt, "fmt_iges") != 0 && K_STRING::cmp(fileFmt, "fmt_step") != 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertCAD2Arrays: unknown file format.");
    return NULL;
  }
  
  vector<FldArrayI*> c;
  vector<FldArrayF*> ufield;  // field read for each zone
  vector<E_Int> et;
  vector<char*> zoneNames;

  printf("Reading %s (%s)...", fileName, fileFmt);
  fflush(stdout);
  
  E_Int ret = CADread2(fileName, fileFmt, h, chordal_err, gr, merge_tol, ufield, c, et, zoneNames, do_join);

  if (ret == 1)
  {
    char error[256];
    sprintf(error, "convertCAD2Arrays: fail to read %s (%s).", 
            fileName, fileFmt);
    PyErr_SetString(PyExc_IOError, error);
    return NULL;
  }
  printf("done.\n");

  // Building numpy arrays
  PyObject* tpl;
    
  if (c.size() == 0)
  {
    printf("Warning: convertCAD2Arrays: no block in file.\n");
    return PyList_New(0);
  }
  
  PyObject* l = PyList_New(0);
  
  E_Int n = ufield.size();    
  for (E_Int i = 0; i < n; i++)
  {
    tpl = K_ARRAY::buildArray(*ufield[i], "x,y,z",
                              *c[i], et[i]);
    delete ufield[i]; delete c[i];
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  // build zoneNames list. Les fonctions de lecture ont alloue un char* par
  // zone lue dans l'ordre (zones structurees, non structurees)
  for (size_t i = 0; i < zoneNames.size(); i++) delete [] zoneNames[i];
  
  return l;
}

// sub routine
E_Int K_OCC::CADread2
(char* file, char* fileFmt, E_Float h, E_Float chordal_err, 
 E_Float gr, E_Float merge_tol,
 vector<FldArrayF*>& unstructField,
 vector<FldArrayI*>& connect,
 vector<E_Int>& eltType,
 vector<char*>& zoneNames,
 bool do_join)
{
   std::vector<K_FLD::FloatArray> crds;
   std::vector<K_FLD::IntArray> connectMs;
   bool aniso = false; // fixme : currently aniso mode mixed with growth ratio does not work
      
#ifdef DEBUG_CAD_READER
  std::cout << "import_OCC_CAD_wrapper::import_cad..." << std::endl;
  std::cout << "import_cad..." << std::endl;
#endif
    
  // CAD --> OCC Shape with associated homemade graph to link flat storage ids between faces and edges.
  CADviaOCC reader;
  E_Int err = reader.import_cad(file, fileFmt, h, chordal_err, gr);
  if (err) return err;
  
#ifdef DEBUG_CAD_READER
  std::cout << "import_cad done." << std::endl;
  std::cout << "mesh_edges..." << std::endl;
#endif
    
  // Mesh the edges.
  std::vector<K_FLD::IntArray> connectEs;
  K_FLD::FloatArray coords;
  err = reader.mesh_edges(coords, connectEs);
  if (err) return err;
  
#ifdef DEBUG_CAD_READER
  std::cout << "mesh_edges done." << std::endl;
  std::cout << "build_loops..." << std::endl;
#endif
  
  // Prepare loops.
  std::vector<K_FLD::IntArray> connectBs;
  err = reader.build_loops(coords, connectEs, connectBs, merge_tol);
  if (err) return err;
  
#ifdef DEBUG_CAD_READER
  std::cout << "build_loops done." << std::endl;
  std::cout << "mesh_faces..." << std::endl;
#endif
    
  // Mesh the surfaces.
  connectMs.clear();
  err = reader.mesh_faces2(coords, connectBs, crds, connectMs, aniso, do_join);
  
  if (err) return err;

  E_Int nmeshes = connectMs.size();

  unstructField.resize(nmeshes, nullptr);
  connect.resize(nmeshes, nullptr);
  eltType.resize(nmeshes, 0);
  zoneNames.resize(nmeshes, nullptr);
 
  for (size_t i=0; i < connectMs.size(); i++)
  {
    if (connectMs[i].cols() == 0) continue; //failed to mesh it
     
    FldArrayF* crd = new FldArrayF;
    crds[i].convert(*crd);
    unstructField[i] = crd;
    FldArrayI* cnt = new FldArrayI;
    connectMs[i].convert(*cnt,1/*shift*/);
    connect[i] = cnt;
    
    char* zoneName = new char [128];
    sprintf(zoneName, "Zone%zu",i);

    E_Int row = connectMs[i].rows();

    zoneNames[i] = zoneName;

    eltType[i] = row -1 ; // 2->TRI or 3->QUAD
  }

  return 0;
}
