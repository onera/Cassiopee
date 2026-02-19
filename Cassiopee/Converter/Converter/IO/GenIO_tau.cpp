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

// Binary NETCDF3 TAU (.grid) file support
// Netcdf3 can NOT store int64

# include "GenIO.h"
#include "netcdf.h"
#include "kcore.h"

/* ------------------------------------------------------------------------- */
// return new index when input is global orig index in connectivities
/* ------------------------------------------------------------------------- */
E_Int newIndex(E_Int index, 
               std::vector<E_Int>& oldStart, std::vector<E_Int>& oldEnd,
               std::vector<E_Int>& newStart)
{
  E_Int b = -1;
  for (size_t i = 0; i < oldStart.size(); i++)
  {
    if (index >= oldStart[i] && index <= oldEnd[i]) { b = i; break; }
  }
  if (b == -1)
  {
    printf("Error: tauwrite: PL index out of range: " SF_D_ "\n", index);
    return 0;
  }

  E_Int offset = index-oldStart[b];
  //printf("%d %d %d %d\n", index, oldStart[i], newStart[i], offset);
  //printf("%d %d\n", index, newStart[i]+offset);
  return newStart[b]+offset-1;
}

/* ------------------------------------------------------------------------- */
// Create GridElement node for basic elements connectivity
/* ------------------------------------------------------------------------- */
E_Int createGridElements4Tau(E_Int eltType, const char* name, E_Int ncells, E_Int istart, 
  PyArrayObject* rc, PyObject*& GE)
{
  // GridElements
  std::vector<npy_intp> npy_dim_vals(2);
  npy_dim_vals[0] = 2;
#ifdef E_DOUBLEINT
  PyArrayObject* r1 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT64, 1);
  int64_t* pp1 = (int64_t*)PyArray_DATA(r1);
#else
  PyArrayObject* r1 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);
  int32_t* pp1 = (int32_t*)PyArray_DATA(r1);
#endif
  pp1[0] = eltType; 
  if (eltType == 5 || eltType == 7) pp1[1] = 1; // tag as BC connect
  else pp1[1] = 0; 
  PyObject* children1 = PyList_New(0);
  GE = Py_BuildValue("[sOOs]", name, r1, children1, "Elements_t");

  // Element range
#ifdef E_DOUBLEINT
  PyArrayObject* r2 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT64, 1);
  int64_t* pp2 = (int64_t*)PyArray_DATA(r2);
#else
  PyArrayObject* r2 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);
  int32_t* pp2 = (int32_t*)PyArray_DATA(r2);
#endif
  pp2[0] = istart; pp2[1] = istart+ncells-1;
  PyObject* children2 = PyList_New(0);
  PyObject* er = Py_BuildValue("[sOOs]", "ElementRange", r2, children2, "IndexRange_t");
  PyList_Append(children1, er); Py_DECREF(er);

  // Element connectivity
  PyObject* children3 = PyList_New(0);
  PyObject* ec = Py_BuildValue("[sOOs]", "ElementConnectivity", rc, children3, "DataArray_t");
  PyList_Append(children1, ec); Py_DECREF(ec);

  return 0;
}

//=============================================================================
/* 
   tauread
*/
//=============================================================================
E_Int K_IO::GenIO::tauread(char* file, PyObject*& tree)
{
  tree = Py_None;

  /* Open file */
  int ncid, ret;
  ret = nc_open(file, NC_NOWRITE, &ncid);
  if (ret != NC_NOERR)
  {
    printf("Warning: tauread: can not open file %s.\n", file);
    return 1;
  }

  // All numpy to set
  IMPORTNUMPY;
  PyArrayObject* xc = NULL;
  PyArrayObject* yc = NULL;
  PyArrayObject* zc = NULL;
  PyArrayObject* hexa = NULL;
  PyArrayObject* tetra = NULL;
  PyArrayObject* penta = NULL;
  PyArrayObject* pyra = NULL;
  PyArrayObject* tri = NULL;
  PyArrayObject* quad = NULL;
  int32_t* bctag = NULL;
  E_Int nhexa=0; E_Int ntetra=0; E_Int npenta=0; E_Int npyra=0; E_Int ntri=0; E_Int nquad=0;
  E_Int nvertex=0; E_Int ncells=0;
  std::vector<npy_intp> npy_dim_vals(2);
        
  // Get the number of datasets in the file
  int nd;
  ret = nc_inq_nvars(ncid, &nd);
  //printf("Number of datasets: %d\n", nd);

  // Iterate over all datasets
  nc_type dtype;
  for (int id = 0; id < nd; id++) 
  {
    // Get dataset name
    char name[NC_MAX_NAME + 1];
    nc_inq_varname(ncid, id, name);
    //printf(">> DataSet Name: %s\n", name);

    // Get number of dimensions of dataset
    int ndims;
    nc_inq_varndims(ncid, id, &ndims);
    //printf("Number of dimensions: %d\n", ndims);

    // Retrieve dataset dims
    int dimsid[10]; dimsid[0] = 0; dimsid[1] = 0; dimsid[2] = 0;
    size_t dims[10];
    nc_inq_vardimid(ncid, id, dimsid);
    for (E_Int i = 0; i < ndims; i++) 
    {
      nc_inq_dim(ncid, dimsid[i], NULL, &dims[i]);
    }

    size_t start[10]; start[0] = 0; start[1] = 0; start[2] = 0;
    size_t count[10]; count[0] = 0; count[1] = 0; count[2] = 0;
    size_t size = 0;

    // Get dataset type
    nc_inq_vartype(ncid, id, &dtype);

    if (strcmp(name, "points_xc") == 0)
    {
      if (dtype == NC_DOUBLE)
      {
        npy_dim_vals[0] = dims[0];
        count[0] = dims[0];
        xc = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_FLOAT64, 1);
        E_Float* pxc = (E_Float*)PyArray_DATA(xc);
        nc_get_vara_double(ncid, id, start, count, pxc);
        nvertex = dims[0];
      }
    }
    else if (strcmp(name, "points_yc") == 0)
    {
      if (dtype == NC_DOUBLE)
      {
        npy_dim_vals[0] = dims[0];
        count[0] = dims[0];
        yc = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_FLOAT64, 1);
        E_Float* pyc = (E_Float*)PyArray_DATA(yc);
        nc_get_vara_double(ncid, id, start, count, pyc);
      }
    }
    else if (strcmp(name, "points_zc") == 0)
    {
      if (dtype == NC_DOUBLE)
      {
        npy_dim_vals[0] = dims[0];
        count[0] = dims[0];
        zc = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_FLOAT64, 1);
        E_Float* pzc = (E_Float*)PyArray_DATA(zc);
        nc_get_vara_double(ncid, id, start, count, pzc);
      }
    }
    else if (strcmp(name, "points_of_hexaeders") == 0) // HEXA
    {
      if (dtype == NC_INT)
      {
        size = dims[0]*dims[1];
        npy_dim_vals[0] = size;
        nhexa = dims[0];
        count[0] = dims[0]; count[1] = dims[1];
#ifdef E_DOUBLEINT
        hexa = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT64, 1);  
        int64_t* pp3 = (int64_t*)PyArray_DATA(hexa);
        int32_t* local = new int32_t [size];
        nc_get_vara_int(ncid, id, start, count, local);
        for (size_t i = 0; i < size; i++) pp3[i] = local[i]+1;
        delete [] local;
#else
        hexa = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);  
        int32_t* pp3 = (int32_t*)PyArray_DATA(hexa);
        nc_get_vara_int(ncid, id, start, count, pp3);
        for (size_t i = 0; i < size; i++) pp3[i] += 1;
#endif
      }
    }
    else if (strcmp(name, "points_of_tetraeders") == 0) // TETRA
    {
      if (dtype == NC_INT)
      {
        size = dims[0]*dims[1];
        npy_dim_vals[0] = size;
        ntetra = dims[0];
        count[0] = dims[0]; count[1] = dims[1];
#ifdef E_DOUBLEINT
        tetra = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT64, 1);  
        int64_t* pp3 = (int64_t*)PyArray_DATA(tetra);
        int32_t* local = new int32_t [size];
        nc_get_vara_int(ncid, id, start, count, local);
        for (size_t i = 0; i < size; i++) pp3[i] = local[i]+1;
        delete [] local;
#else
        tetra = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);  
        int32_t* pp3 = (int32_t*)PyArray_DATA(tetra);
        nc_get_vara_int(ncid, id, start, count, pp3);
        for (size_t i = 0; i < size; i++) pp3[i] += 1;
#endif
      }
    }
    else if (strcmp(name, "points_of_prisms") == 0) // PENTA
    {
      if (dtype == NC_INT)
      {
        size = dims[0]*dims[1];
        npy_dim_vals[0] = size;
        npenta = dims[0];
        count[0] = dims[0]; count[1] = dims[1];
#ifdef E_DOUBLEINT
        penta = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT64, 1);  
        int64_t* pp3 = (int64_t*)PyArray_DATA(penta);
        int32_t* local = new int32_t [size];
        nc_get_vara_int(ncid, id, start, count, local);
        for (size_t i = 0; i < size; i++) pp3[i] = local[i]+1;
        delete [] local;
#else
        penta = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);
        int32_t* pp3 = (int32_t*)PyArray_DATA(penta);
        nc_get_vara_int(ncid, id, start, count, pp3);
        for (size_t i = 0; i < size; i++) pp3[i] += 1;
#endif
      }
    }
    else if (strcmp(name, "points_of_pyramids") == 0) // PYRA
    {
      if (dtype == NC_INT)
      {
        size = dims[0]*dims[1];
        npy_dim_vals[0] = size;
        npyra = dims[0];
        count[0] = dims[0]; count[1] = dims[1];
#ifdef E_DOUBLEINT
        pyra = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT64, 1);  
        int64_t* pp3 = (int64_t*)PyArray_DATA(pyra);
        int32_t* local = new int32_t [size];
        nc_get_vara_int(ncid, id, start, count, local);
        for (size_t i = 0; i < size; i++) pp3[i] = local[i]+1;
        delete [] local;
#else
        pyra = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);  
        int32_t* pp3 = (int32_t*)PyArray_DATA(pyra);
        nc_get_vara_int(ncid, id, start, count, pp3);
        for (size_t i = 0; i < size; i++) pp3[i] += 1;
#endif
      }
    }
    else if (strcmp(name, "points_of_surfacetriangles") == 0) // TRI
    {
      if (dtype == NC_INT)
      {
        size = dims[0]*dims[1];
        npy_dim_vals[0] = size;
        ntri = dims[0];
        count[0] = dims[0]; count[1] = dims[1];
#ifdef E_DOUBLEINT
        tri = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT64, 1);  
        int64_t* pp3 = (int64_t*)PyArray_DATA(tri);
        int32_t* local = new int32_t [size];
        nc_get_vara_int(ncid, id, start, count, local);
        for (size_t i = 0; i < size; i++) pp3[i] = local[i]+1;
        delete [] local;
#else
        tri = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);  
        int32_t* pp3 = (int32_t*)PyArray_DATA(tri);
        nc_get_vara_int(ncid, id, start, count, pp3);
        for (size_t i = 0; i < size; i++) pp3[i] += 1;
#endif
      }
    }
    else if (strcmp(name, "points_of_surfacequadrilaterals") == 0) // QUADS
    {
      if (dtype == NC_INT)
      {
        size = dims[0]*dims[1];
        npy_dim_vals[0] = size;
        nquad = dims[0];
        count[0] = dims[0]; count[1] = dims[1];
#ifdef E_DOUBLEINT
        quad = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT64, 1);  
        int64_t* pp3 = (int64_t*)PyArray_DATA(quad);
        int32_t* local = new int32_t [size];
        nc_get_vara_int(ncid, id, start, count, local);
        for (size_t i = 0; i < size; i++) pp3[i] = local[i]+1;
        delete [] local;
#else
        quad = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);  
        int32_t* pp3 = (int32_t*)PyArray_DATA(quad);
        nc_get_vara_int(ncid, id, start, count, pp3);
        for (size_t i = 0; i < size; i++) pp3[i] += 1;
#endif
      }
    }
    else if (strcmp(name, "boundarymarker_of_surfaces") == 0) // MARKERS
    {
      if (dtype == NC_INT)
      {
        size = dims[0];
        count[0] = dims[0];
        bctag = new int32_t [size];
        nc_get_vara_int(ncid, id, start, count, bctag);
      }
    }

    /* all netcdf types
    switch (dtype)
    {
      case NC_BYTE:
        printf("Data Type: NC_BYTE\n");
        break;
      case NC_CHAR:
        printf("Data Type: NC_CHAR\n");
        break;
      case NC_SHORT:
        printf("Data Type: NC_SHORT\n");
        break;
      case NC_INT:
        printf("Data Type: NC_INT\n");
        break;
      case NC_FLOAT:
        printf("Data Type: NC_FLOAT\n");
        break;
      case NC_DOUBLE:
        printf("Data Type: NC_DOUBLE\n");
        break;
      default:
        printf("Data Type: Unknown\n");
        break;
    }
    */
  }

  // Create tree
  PyObject* children1 = PyList_New(0);
  tree = Py_BuildValue("[sOOs]", "tree", Py_None, children1, "CGNSTree");
  
  // Create version
  npy_dim_vals[0] = 1;
  PyObject* children2 = PyList_New(0);
  PyArrayObject* r2 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_FLOAT32, 1);
  float* pf2 = (float*)PyArray_DATA(r2);
  pf2[0] = 4.0;
  PyObject* version = Py_BuildValue("[sOOs]", "CGNSLibraryVersion", (PyObject*)r2, children2, "CGNSLibraryVersion_t");
  PyList_Append(children1, version); Py_DECREF(version);
  
  // Create Base
  npy_dim_vals[0] = 2;
#ifdef E_DOUBLEINT
  PyArrayObject* r3 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT64, 1);
  int64_t* pp3 = (int64_t*)PyArray_DATA(r3);
#else
  PyArrayObject* r3 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);
  int32_t* pp3 = (int32_t*)PyArray_DATA(r3);
#endif
  pp3[0] = 3; pp3[1] = 3;
  PyObject* children3 = PyList_New(0);
  PyObject* base = Py_BuildValue("[sOOs]", "Base", r3, children3, "CGNSBase_t");
  PyList_Append(children1, base); Py_DECREF(base);
  
  // Create Zone
  npy_dim_vals[0] = 1; npy_dim_vals[1] = 3;
#ifdef E_DOUBLEINT
  PyArrayObject* r4 = (PyArrayObject*)PyArray_EMPTY(2, &npy_dim_vals[0], NPY_INT64, 1);
  int64_t* pp4 = (int64_t*)PyArray_DATA(r4);
#else
  PyArrayObject* r4 = (PyArrayObject*)PyArray_EMPTY(2, &npy_dim_vals[0], NPY_INT32, 1);
  int32_t* pp4 = (int32_t*)PyArray_DATA(r4);
#endif
  pp4[0] = 0; pp4[1] = 0; pp4[2] = 0;
  PyObject* children4 = PyList_New(0);
  PyObject* zone = Py_BuildValue("[sOOs]", "Zone", r4, children4, "Zone_t");
  PyList_Append(children3, zone); Py_DECREF(zone);
  
  // Create ZoneType
  PyObject* children9 = PyList_New(0);
  npy_dim_vals[0] = 12;
  PyArrayObject* r9 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_STRING, 1);
  char* pp9 = (char*)PyArray_DATA(r9);
  K_STRING::cpy(pp9, "Unstructured", 12, false);
  PyObject* zoneType = Py_BuildValue("[sOOs]", "ZoneType", r9, children9, "ZoneType_t");
  PyList_Append(children4, zoneType); Py_DECREF(zoneType);
  
  // Create GridCoordinates
  PyObject* children5 = PyList_New(0);
  PyObject* GC = Py_BuildValue("[sOOs]", "GridCoordinates", Py_None, children5, "GridCoordinates_t");
  PyList_Append(children4, GC); Py_DECREF(GC);
  
  // CoordinateX
  PyObject* children6 = PyList_New(0);
  PyObject* nxc = Py_BuildValue("[sOOs]", "CoordinateX", xc, children6, "DataArray_t");
  PyList_Append(children5, nxc); Py_DECREF(nxc);

  // CoordinateY
  PyObject* children7 = PyList_New(0);
  PyObject* nyc = Py_BuildValue("[sOOs]", "CoordinateY", yc, children7, "DataArray_t");
  PyList_Append(children5, nyc); Py_DECREF(nyc);

  // CoordinateZ
  PyObject* children8 = PyList_New(0);
  PyObject* nzc = Py_BuildValue("[sOOs]", "CoordinateZ", zc, children8, "DataArray_t");
  PyList_Append(children5, nzc); Py_DECREF(nzc);

  // Create GridElements
  PyObject* GE; E_Int istart=1;
  if (hexa != NULL) 
  {
    createGridElements4Tau(17, "HEXA", nhexa, istart, hexa, GE);
    PyList_Append(children4, GE); Py_DECREF(GE);
    istart += nhexa;
    ncells += nhexa;
  }
  if (tetra != NULL) 
  {
    createGridElements4Tau(10, "TETRA", ntetra, istart, tetra, GE);
    PyList_Append(children4, GE); Py_DECREF(GE);
    istart += ntetra;
    ncells += ntetra;
  }
  if (penta != NULL) 
  {
    createGridElements4Tau(14, "PENTA", npenta, istart, penta, GE);
    PyList_Append(children4, GE); Py_DECREF(GE);
    istart += npenta;
    ncells += npenta;
  }
  if (pyra != NULL) 
  {
    createGridElements4Tau(12, "PYRA", npyra, istart, pyra, GE);
    PyList_Append(children4, GE); Py_DECREF(GE);
    istart += npyra;
    ncells += npyra;
  }
  if (tri != NULL) 
  {
    createGridElements4Tau(5, "TRI", ntri, istart, tri, GE);
    PyList_Append(children4, GE); Py_DECREF(GE);
    istart += ntri;
  }
  if (quad != NULL) 
  {
    createGridElements4Tau(7, "QUAD", nquad, istart, quad, GE);
    PyList_Append(children4, GE); Py_DECREF(GE);
    istart += nquad;
  }
  if (bctag != NULL)
  {
    // Create ZoneBC
    PyObject* children9 = PyList_New(0);
    PyObject* nzbc = Py_BuildValue("[sOOs]", "ZoneBC", Py_None, children9, "ZoneBC_t");
    PyList_Append(children4, nzbc); Py_DECREF(nzbc);

    // Build PL corresponding to tags
    size_t size = ntri+nquad;
    size_t nvol = ntetra + nhexa + npenta + npyra;
    std::map<E_Int, E_Int> tagmap; // tag -> nbre elements tagged
    for (size_t i = 0; i < size; i++)
    {
      if (tagmap.find(bctag[i]) == tagmap.end()) tagmap[bctag[i]] = 1;
      else tagmap[bctag[i]] += 1;
    }
    char bcname[256]; char bctype[256];
    for (const auto& pair : tagmap) // for each tag
    {
      E_Int tag = pair.first; // tag
      E_Int nfaces = pair.second; // nbre de faces pour ce tag
      //printf("tag " SF_D_ " is set " SF_D_ " times.\n", tag, nfaces);
      if (nfaces == 0) continue;
      // Create BC_t for tag
      PyObject* children10 = PyList_New(0);
      strcpy(bctype, "FamilySpecified");
      sprintf(bcname, "BC" SF_D_, tag);
      npy_dim_vals[0] = strlen(bctype);
      PyArrayObject* r10 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_STRING, 1);
      char* pp10 = (char*)PyArray_DATA(r10);
      K_STRING::cpy(pp10, bctype, strlen(bctype), false);
      PyObject* bc = Py_BuildValue("[sOOs]", bcname, r10, children10, "BC_t");
      PyList_Append(children9, bc); Py_DECREF(bc);
      // Create point list
      npy_dim_vals[0] = 1;
      npy_dim_vals[1] = nfaces;
#ifdef E_DOUBLEINT
      PyArrayObject* r11 = (PyArrayObject*)PyArray_EMPTY(2, &npy_dim_vals[0], NPY_INT64, 1);
      int64_t* pp11 = (int64_t*)PyArray_DATA(r11);
#else
      PyArrayObject* r11 = (PyArrayObject*)PyArray_EMPTY(2, &npy_dim_vals[0], NPY_INT32, 1);
      int32_t* pp11 = (int32_t*)PyArray_DATA(r11);
#endif
      E_Int c = 0;
      for (size_t i = 0; i < size; i++)
      {
        if (bctag[i] == tag) { pp11[c] = i+1+nvol; c++; }
      }
      PyObject* children11 = PyList_New(0);
      PyObject* pl = Py_BuildValue("[sOOs]", "PointList", r11, children11, "IndexArray_t");
      PyList_Append(children10, pl); Py_DECREF(pl);
      // Create GridLocation
      npy_dim_vals[0] = 10;
      PyArrayObject* r12 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_STRING, 1);
      char* pp12 = (char*)PyArray_DATA(r12);
      K_STRING::cpy(pp12, "FaceCenter", 10, false);
      PyObject* children12 = PyList_New(0);
      PyObject* gl = Py_BuildValue("[sOOs]", "GridLocation", r12, children12, "GridLocation_t");
      PyList_Append(children10, gl); Py_DECREF(gl);
      // Create BC FamilyName
      PyObject* children13 = PyList_New(0);
      sprintf(bctype, "BCType" SF_D_, tag);
      npy_dim_vals[0] = strlen(bctype);
      PyArrayObject* r13 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_STRING, 1);
      char* pp13 = (char*)PyArray_DATA(r13);
      K_STRING::cpy(pp13, bctype, strlen(bctype), false);
      PyObject* famName = Py_BuildValue("[sOOs]", "FamilyName", r13, children13, "FamilyName_t");
      PyList_Append(children10, famName); Py_DECREF(famName);
    }
  }

  // Modify the zone attribute
  pp4[0] = nvertex; pp4[1] = ncells;

  // close file
  nc_close(ncid);
  
  if (bctag != NULL) delete [] bctag;
  
  return 0;
}

//=============================================================================
/* 
   tauwrite
*/
//=============================================================================
E_Int K_IO::GenIO::tauwrite(char* file, PyObject* tree)
{
  /* Open file for writing */
  int ncid, ret;
  ret = nc_create(file, NC_CLOBBER, &ncid);
  if (ret != NC_NOERR)
  {
    printf("Warning: tauwrite: can not open file %s for writing.\n", file);
    return 1;
  }

  // Check only for one zone
  std::vector<PyObject*> bases;
  K_PYTREE::getNodesFromType1(tree, "CGNSBase_t", bases);

  std::vector<PyObject*> zones;
  for (size_t i = 0; i < bases.size(); i++)
  {
    K_PYTREE::getNodesFromType1(bases[i], "Zone_t", zones);
  }
  if (zones.size() > 1)
  {
    printf("Warning: tauwrite: only first zone is written.\n");
  }
  if (zones.size() == 0) return 0;
  
  PyObject* zone = zones[0];

  // Write file attributes
  nc_put_att_text(ncid, NC_GLOBAL, "type", 24, "Primary Grid: Tau Format");

  // Write connectivities
  std::vector<PyObject*> connects;
  K_PYTREE::getNodesFromType1(zone, "Elements_t", connects);

  std::vector<PyArrayObject*> hook;
  int ni_tetra, ni_hexa, ni_penta, ni_pyra, ni_tri, ni_quad;
  int nj_tetra, nj_hexa, nj_penta, nj_pyra, nj_tri, nj_quad;
  int no_of_points, no_of_markers;
  int no_of_elements, no_of_surfaceelements;
  int var_tetra, var_hexa, var_penta, var_pyra, var_tri, var_quad;
  int var_xc, var_yc, var_zc, var_markers;
  E_Int nhexa=0; E_Int ntetra=0; E_Int npenta=0; E_Int npyra=0; E_Int ntri=0; E_Int nquad=0;
  E_Int off=0;
  int dimids[2];

  // Create dims and vars
  for (size_t i = 0; i < connects.size(); i++)
  {
    E_Int* cv = K_PYTREE::getValueAI(connects[i], hook);
    PyObject* er = K_PYTREE::getNodeFromName1(connects[i], "ElementRange");
    E_Int* erv = K_PYTREE::getValueAI(er, hook);
    E_Int eltType = cv[0];
    E_Int size0 = erv[1]-erv[0]+1;
    
    // count elements
    if (eltType == 17) // HEXA
    {
      nhexa += size0;
    }
    else if (eltType == 10) // TETRA
    {
      ntetra += size0;
    }
    else if (eltType == 14) // PENTA
    {
      npenta += size0;
    }
    else if (eltType == 12) // PYRA
    {
      npyra += size0;
    }
    else if (eltType == 5) // TRI
    {
      ntri += size0;
    }
    else if (eltType == 7) // QUADS
    {
      nquad += size0;
    }
  }

  // create dim id and vars
  if (nhexa > 0) // HEXA
  {
    nc_def_dim(ncid, "no_of_hexaeders", nhexa, &ni_hexa);
    nc_def_dim(ncid, "points_per_hexaeder", 8, &nj_hexa);
    dimids[0] = ni_hexa; dimids[1] = nj_hexa;
    nc_def_var(ncid, "points_of_hexaeders", NC_INT, 2, dimids, &var_hexa);
  }
  if (ntetra > 0) // TETRA
  {
    nc_def_dim(ncid, "no_of_tetraeders", ntetra, &ni_tetra);
    nc_def_dim(ncid, "points_per_tetraeder", 4, &nj_tetra);
    dimids[0] = ni_tetra; dimids[1] = nj_tetra;
    nc_def_var(ncid, "points_of_tetraeders", NC_INT, 2, dimids, &var_tetra);
  }
  if (npenta > 0) // PENTA
  {
    nc_def_dim(ncid, "no_of_prisms", npenta, &ni_penta);
    nc_def_dim(ncid, "points_per_prism", 6, &nj_penta);
    dimids[0] = ni_penta; dimids[1] = nj_penta;
    nc_def_var(ncid, "points_of_prisms", NC_INT, 2, dimids, &var_penta);
  }
  if (npyra > 0) // PYRA
  {
    nc_def_dim(ncid, "no_of_pyramids", npyra, &ni_pyra);
    nc_def_dim(ncid, "points_per_pyramids", 5, &nj_pyra);
    dimids[0] = ni_pyra; dimids[1] = nj_pyra;
    nc_def_var(ncid, "points_of_pyramids", NC_INT, 2, dimids, &var_pyra);
  }
  if (ntri > 0) // TRI
  {
    nc_def_dim(ncid, "no_of_surfacetriangles", ntri, &ni_tri);
    nc_def_dim(ncid, "points_per_surfacetriangle", 3, &nj_tri);
    dimids[0] = ni_tri; dimids[1] = nj_tri;
    nc_def_var(ncid, "points_of_surfacetriangles", NC_INT, 2, dimids, &var_tri);
  }
  if (nquad > 0) // QUAD
  {      
    nc_def_dim(ncid, "no_of_surfacequadrilaterals", nquad, &ni_quad);
    nc_def_dim(ncid, "points_per_surfacequadrilateral", 4, &nj_quad);
    dimids[0] = ni_quad; dimids[1] = nj_quad;
    nc_def_var(ncid, "points_of_surfacequadrilaterals", NC_INT, 2, dimids, &var_quad);
  }

  nc_def_dim(ncid, "no_of_elements", nhexa+ntetra+npenta+npyra, &no_of_elements);
  nc_def_dim(ncid, "no_of_surfaceelements", ntri+nquad, &no_of_surfaceelements);
  
  E_Int zoneType, nvertex, ncells;
  K_PYTREE::getZoneDim(zone, zoneType, nvertex, ncells, hook);
    
  nc_def_dim(ncid, "no_of_points", nvertex, &no_of_points);
  dimids[0] = no_of_points;
  nc_def_var(ncid, "points_xc", NC_DOUBLE, 1, dimids, &var_xc);
  nc_def_var(ncid, "points_yc", NC_DOUBLE, 1, dimids, &var_yc);
  nc_def_var(ncid, "points_zc", NC_DOUBLE, 1, dimids, &var_zc);

  PyObject* zbc = K_PYTREE::getNodeFromName1(zone, "ZoneBC");
  std::vector<PyObject*> BCs;
  if (zbc != NULL) K_PYTREE::getNodesFromType1(zbc, "BC_t", BCs);
  nc_def_dim(ncid, "no_of_markers", BCs.size(), &no_of_markers);
  
  dimids[0] = no_of_surfaceelements;
  nc_def_var(ncid, "boundarymarker_of_surfaces", NC_INT, 1, dimids, &var_markers);

  nc_enddef(ncid); // end of define mode

  // Write arrays, connectivities are concatenated
  std::vector<E_Int> oldStart(connects.size()); // old start of connect
  std::vector<E_Int> oldEnd(connects.size()); // old end of connect
  std::vector<E_Int> newStart(connects.size()); // new start of connect
  int shexa=0, stetra=0, spenta=0, spyra=0, stri=0, squad=0;
  int phexa=1, ptetra=nhexa+1, ppenta=nhexa+ntetra+1, ppyra=nhexa+ntetra+npenta+1;
  int ptri=nhexa+ntetra+npenta+npyra+1, pquad=nhexa+ntetra+npenta+npyra+ntri+1;
  size_t start[2]; size_t scount[2];

  for (size_t i = 0; i < connects.size(); i++)
  {
    E_Int* cv = K_PYTREE::getValueAI(connects[i], hook);
    PyObject* ec = K_PYTREE::getNodeFromName1(connects[i], "ElementConnectivity");
    E_Int* ecv = K_PYTREE::getValueAI(ec, hook);
    PyObject* er = K_PYTREE::getNodeFromName1(connects[i], "ElementRange");
    E_Int* erv = K_PYTREE::getValueAI(er, hook);
    E_Int eltType = cv[0];
    E_Int size = erv[1]-erv[0]+1;
    E_Int size0 = size;

    if (eltType == 17) // HEXA
    { size *= 8; }
    else if (eltType == 10) // TETRA
    { size *= 4; }
    else if (eltType == 14) // PENTA
    { size *= 6; }
    else if (eltType == 12) // PYRA
    { size *= 5; }
    else if (eltType == 5) // TRI
    { size *= 3; }
    else if (eltType == 7) // QUADS
    { size *= 4; }

    int32_t* ecv2 = new int32_t [size];
    for (E_Int j = 0; j < size; j++) { ecv2[j] = ecv[j]-1; }

    if (eltType == 17) // HEXA
    {
      //nc_put_var_int(ncid, var_hexa, ecv2);
      start[0]=shexa; start[1]=0;
      scount[0]=size0; scount[1]=8;
      nc_put_vara_int(ncid, var_hexa, start, scount, ecv2);
      shexa += size0;
      oldStart[i] = erv[0]; oldEnd[i] = erv[1];
      newStart[i] = phexa;
      phexa += size0;
    }
    else if (eltType == 10) // TETRA
    {
      //nc_put_var_int(ncid, var_tetra, ecv2);
      start[0]=stetra; start[1]=0;
      scount[0]=size0; scount[1]=4;
      nc_put_vara_int(ncid, var_tetra, start, scount, ecv2);
      stetra += size0;
      oldStart[i] = erv[0]; oldEnd[i] = erv[1];
      newStart[i] = ptetra;
      ptetra += size0;
    }
    else if (eltType == 14) // PENTA
    {
      //nc_put_var_int(ncid, var_penta, ecv2);
      start[0]=spenta; start[1]=0;
      scount[0]=size0; scount[1]=6;
      nc_put_vara_int(ncid, var_penta, start, scount, ecv2);
      spenta += size0;
      oldStart[i] = erv[0]; oldEnd[i] = erv[1];
      newStart[i] = ppenta;
      ppenta += size0;
    }
    else if (eltType == 12) // PYRA
    {
      //nc_put_var_int(ncid, var_pyra, ecv2);
      start[0]=spyra; start[1]=0;
      scount[0]=size0; scount[1]=5;
      nc_put_vara_int(ncid, var_pyra, start, scount, ecv2);
      spyra += size0;
      oldStart[i] = erv[0]; oldEnd[i] = erv[1];
      newStart[i] = ppyra;
      ppyra += size0;
    }
    else if (eltType == 5) // TRI
    {
      //nc_put_var_int(ncid, var_tri, ecv2);
      start[0]=stri; start[1]=0;
      scount[0]=size0; scount[1]=3;
      nc_put_vara_int(ncid, var_tri, start, scount, ecv2);
      stri += size0;
      oldStart[i] = erv[0]; oldEnd[i] = erv[1];
      newStart[i] = ptri;
      ptri += size0;
    }
    else if (eltType == 7) // QUADS
    {
      //nc_put_var_int(ncid, var_quad, ecv2);
      start[0]=squad; start[1]=0;
      scount[0]=size0; scount[1]=4;
      nc_put_vara_int(ncid, var_quad, start, scount, ecv2);
      squad += size0;
      oldStart[i] = erv[0]; oldEnd[i] = erv[1];
      newStart[i] = pquad;
      pquad += size0;
    }
    delete [] ecv2;
  }

  E_Int possurf = nhexa+ntetra+npenta+npyra;
  //printf("possurf=%d\n", possurf);
  //for (size_t i = 0; i < oldStart.size(); i++)
  //{
  //  printf("oldStart[%d] = %d | ", i, oldStart[i]);
  //  printf("newStart[%d] = %d\n", i, newStart[i]);
  //}
  
  // Write coordinates
  PyObject* gc = K_PYTREE::getNodeFromType1(zone, "GridCoordinates_t");
    
  PyObject* xc = K_PYTREE::getNodeFromName1(gc, "CoordinateX");
  E_Float* xcv = K_PYTREE::getValueAF(xc, hook);
  nc_put_var_double(ncid, var_xc, xcv);

  PyObject* yc = K_PYTREE::getNodeFromName1(gc, "CoordinateY");
  E_Float* ycv = K_PYTREE::getValueAF(yc, hook);
  nc_put_var_double(ncid, var_yc, ycv);

  PyObject* zc = K_PYTREE::getNodeFromName1(gc, "CoordinateZ");
  E_Float* zcv = K_PYTREE::getValueAF(zc, hook);
  nc_put_var_double(ncid, var_zc, zcv);
  
  // Write boundary markers
  E_Int nmarkers = ntri+nquad;
  //E_Int nvol = ntetra+nhexa+npyra+npenta;
  int32_t* markers = new int32_t [nmarkers];
  for (E_Int i = 0; i < nmarkers; i++) markers[i] = 0;

  E_Int count = 1; // BC count (=marker count)
  for (size_t i = 0; i < BCs.size(); i++)
  {
    // PointRange
    PyObject* er = K_PYTREE::getNodeFromName1(BCs[i], "ElementRange");
    if (er == NULL) er = K_PYTREE::getNodeFromName1(BCs[i], "PointRange");
    if (er != NULL)
    {
      E_Int* erv = K_PYTREE::getValueAI(er, hook);
      //printf("BC %d: %d %d\n", count, erv[0], erv[1]);
      for (E_Int j = erv[0]; j <= erv[1]; j++)
      {
        //off = j-1-nvol;
        off = newIndex(j, oldStart, oldEnd, newStart)-possurf;
        //if (off != j-1-nvol) 
        //printf("off=%d %d -> count=%d\n", off, j-1-nvol, count);
        off = std::min(off, nmarkers-1);
        off = std::max(off, E_Int(0));
        markers[off] = count;
      }
    }

    // PointList
    PyObject* pl = K_PYTREE::getNodeFromName1(BCs[i], "PointList");
    if (pl != NULL)
    {
      E_Int ni, nj;
      E_Int* plv = K_PYTREE::getValueAI(pl, ni, nj, hook);
      for (E_Int j = 0; j < ni*nj; j++)
      {
        //off = plv[j]-1-nvol;
        off = newIndex(plv[j], oldStart, oldEnd, newStart)-possurf;
        //printf("off=%d %d\n", off, plv[j]-1-nvol);
        off = std::min(off, nmarkers-1);
        off = std::max(off, E_Int(0));
        markers[off] = count;
      }
    }
    count += 1;
  }
  nc_put_var_int(ncid, var_markers, markers);
  delete [] markers;

  // close file
  nc_close(ncid);

  // Decref
  for (size_t i = 0; i < hook.size(); i++) Py_DECREF(hook[i]);

  return 0;
}
