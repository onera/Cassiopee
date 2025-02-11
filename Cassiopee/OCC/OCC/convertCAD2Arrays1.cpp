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
// convertCAD2Arrays (algo=1)

#include <stdio.h>
#include <string.h>
#include "occ.h"
#include "import_OCC_CAD_wrapper.h"

using namespace std;
using namespace K_FLD;

#include <iostream>

namespace K_OCC 
{
  E_Int CADread(
  char* file, char* fileFmt, E_Float h, E_Float chordal_err, E_Float gr, char*& varString,
  std::vector<K_FLD::FldArrayF*>& unstructField,
  std::vector<K_FLD::FldArrayI*>& connect,
  std::vector<E_Int>& eltType,
  std::vector<char*>& zoneNames,
  bool do_join);
}

// ============================================================================
/* Convert CAD file to arrays using OpenCascade */
// ============================================================================
PyObject* K_OCC::convertCAD2Arrays1(PyObject* self, PyObject* args)
{
  char* fileName; char* fileFmt;
  E_Float h, chordal_err, gr(-1.);
  E_Int do_join(1);

  if (!PYPARSETUPLE_(args, SS_ RRR_ I_, &fileName, &fileFmt, &h, &chordal_err,
                    &gr, &do_join)) return NULL;

  // Check recognised formats
  if (K_STRING::cmp(fileFmt, "fmt_iges") != 0 && K_STRING::cmp(fileFmt, "fmt_step") != 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertCAD2Arrays: unknown file format.");
    return NULL;
  }
  
  char* varString = NULL;
  vector<FldArrayI*> c;
  vector<FldArrayF*> ufield;  // field read for each zone
  vector<E_Int> et;
  vector<char*> zoneNames;

  printf("Reading %s (%s)...", fileName, fileFmt);
  fflush(stdout);
  
  E_Int ret = CADread(fileName, fileFmt, h, chordal_err, gr, varString, ufield, c, et, zoneNames, do_join);

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
    tpl = K_ARRAY::buildArray(*ufield[i], varString,
                              *c[i], et[i]);
    delete ufield[i]; delete c[i];
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  // build zoneNames list. Les fonctions de lecture ont alloue un char* par
  // zone lue dans l'ordre (zones structurees, non structurees)
  E_Int znsize = zoneNames.size();
  for (E_Int i = 0; i < znsize; i++) delete [] zoneNames[i];

  delete [] varString;

  return l;
}

E_Int K_OCC::CADread
(char* file, char* fileFmt, E_Float h, E_Float chordal_err, E_Float gr, char*& varString,
 vector<FldArrayF*>& unstructField,
 vector<FldArrayI*>& connect,
 vector<E_Int>& eltType,
 vector<char*>& zoneNames,
 bool do_join)
{
   std::vector<K_FLD::FloatArray> crds;
   std::vector<K_FLD::IntArray> connectMs;
   bool aniso = false; // fixme : currently aniso mode mixed with growth ratio does not work
   E_Int err = import_OCC_CAD_wrapper::import_cad(file, fileFmt, crds, connectMs, h, chordal_err, gr, aniso, do_join);
   if (err) return err;

   int nmeshes = connectMs.size();

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
     sprintf(zoneName, "Zone%zu", i);

    E_Int row = connectMs[i].rows();

    zoneNames[i] = zoneName;

    eltType[i] = row -1 ; // 2->TRI or 3->QUAD
   }

  varString = new char [8];
  strcpy(varString, "x,y,z");
  
  return 0;
}
