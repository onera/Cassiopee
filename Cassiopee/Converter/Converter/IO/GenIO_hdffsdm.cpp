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

// Binary HDF FSDM (.h5) file support
# include "GenIO.h"
# include "hdf5.h"
# include "GenIO_hdfcgns.h"

/* ------------------------------------------------------------------------- */
#if H5_VERSION_LE(1,13,9)
static herr_t feed_children_names(hid_t id, const char* name,
                                  const H5L_info_t* linfo, void* names)
#else
static herr_t feed_children_names(hid_t id, const char* name,
                                  const H5L_info2_t* linfo, void* names)
#endif
{
  /* skip names starting with a <space> */
  if (name && (name[0] == ' ')) return 0;
  int n = 0;
  char** gnames = (char**)names;
  while (gnames[n] != NULL) n++;
  gnames[n] = new char [256];
  strcpy(gnames[n], name);
  return 0;
}

/* ------------------------------------------------------------------------- */
E_Int createGridElements(hid_t id, E_Int eltType, char* name, E_Int istart, PyObject*& GE, E_Int& ncells, E_Int*& bct)
{
  // Get ncells from id
  hid_t aid = H5Aopen_by_name(id, ".", "NumberOfCells", H5P_DEFAULT, H5P_DEFAULT);
  H5Aread(aid, H5T_NATIVE_INT, &ncells); // care when long
  printf("ncells=" SF_D_ "\n", ncells);

  // GridElements
  std::vector<npy_intp> npy_dim_vals(1);
  npy_dim_vals[0] = 2;
#ifdef E_DOUBLEINT
  PyArrayObject* r1 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT64, 1);
  int64_t* pp1 = (int64_t*)PyArray_DATA(r1);
#else
  PyArrayObject* r1 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);
  int32_t* pp1 = (int32_t*)PyArray_DATA(r1);
#endif
  pp1[0] = eltType; pp1[1] = 0;
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
  PyList_Append(children1, er); Py_INCREF(er);

  // Element connectivity
  hid_t did = H5Dopen2(id, "Cell2Node", H5P_DEFAULT);
  hid_t sid = H5Dget_space(did);
  E_Int ndims = H5Sget_simple_extent_ndims(sid);
  hsize_t dims[3];
  H5Sget_simple_extent_dims(sid, dims, NULL);
  E_Int size = 1;
  for (E_Int i = 0; i < ndims; i++) size *= dims[i];
  //printf("ndims=" SF_D_ "\n", ndims);
  //printf("size=" SF_D_ "\n", size);
  npy_dim_vals[0] = size;
  hid_t tid = H5Tcopy(H5T_NATIVE_INT); H5Tset_precision(tid, 32);
  hid_t yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  hid_t mid = H5S_ALL;
  PyArrayObject* r3 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);  
  int32_t* pp3 = (int32_t*)PyArray_DATA(r3);
  H5Dread(did, yid, mid, sid, H5P_DEFAULT, pp3);
  for (E_Int i = 0; i < size; i++) pp3[i] += 1;

  PyObject* children3 = PyList_New(0);
  PyObject* ec = Py_BuildValue("[sOOs]", "ElementConnectivity", r3, children3, "DataArray_t");
  PyList_Append(children1, ec); Py_INCREF(ec);

  // Get BC tag if any
  hid_t gid = H5Lexists(id, "CellAttributes", H5P_DEFAULT);
  bct = NULL;
  if (gid > 0)
  {
    gid = H5Gopen(id, "CellAttributes", H5P_DEFAULT);
    did = H5Dopen2(gid, "CADGroupID", H5P_DEFAULT);
    sid = H5Dget_space(did);
    ndims = H5Sget_simple_extent_ndims(sid);
    H5Sget_simple_extent_dims(sid, dims, NULL);
    E_Int size = 1;
    for (E_Int i = 0; i < ndims; i++) size *= dims[i];
    //printf("ndims=" SF_D_ "\n", ndims);
    //printf("size=" SF_D_ "\n", size);
    npy_dim_vals[0] = size;
    tid = H5Tcopy(H5T_NATIVE_INT); H5Tset_precision(tid, 32);
    yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
    mid = H5S_ALL;
    int32_t* bct = new int32_t [size];
    H5Dread(did, yid, mid, sid, H5P_DEFAULT, bct);
    // count tags
    std::map<E_Int, E_Int> tagmap; // tag -> nbre elements
    for (E_Int i = 0; i < size; i++) 
    {
      if (tagmap.find(bct[i]) == tagmap.end()) tagmap[bct[i]] = 0;
      else tagmap[bct[i]] += 1;
    }
    for (const auto& pair : tagmap)
    {
      E_Int tag = pair.first;
      E_Int nfaces = pair.second; // nbre de faces pour ce tag
      printf("tag " SF_D_ " is set " SF_D_ " times.\n", tag, nfaces);
    }
  }
  return 0;
}

//=============================================================================
/* 
   hdffsdmread
*/
//=============================================================================
E_Int K_IO::GenIO::hdffsdmread(char* file, PyObject*& tree)
{
  /* Open file */
  hid_t fapl, fid;
  fapl = H5Pcreate(H5P_FILE_ACCESS);
  fid = H5Fopen(file, H5F_ACC_RDONLY, fapl);
  H5Pclose(fapl);
  if (fid < 0)
  {
    printf("Warning: hdffsdmread: can not open file %s.\n", file);
    return 1;
  }

  GenIOHdf HDF;

  // Open FS:Mesh node
  hid_t node = H5Gopen(fid, "/FS:Mesh", H5P_DEFAULT);
  if (node < 0)
  {
    printf("Warning: hdffsdmread: no FS:Mesh node.\n");
    return 1;
  }
  H5Gclose(node);

  hid_t uc = H5Gopen(fid, "/FS:Mesh/UnstructuredCells", H5P_DEFAULT);
  if (uc < 0)
  {
    printf("Warning: hdffsdmread: no FS:Mesh/UnstructuredCells node.\n");
    return 1;
  }

  // Create tree
  PyObject* children1 = PyList_New(0);
  tree = Py_BuildValue("[sOOs]", "tree", Py_None, children1, "CGNSTree");

  // Create version
  std::vector<npy_intp> npy_dim_vals(1);
  npy_dim_vals[0] = 1;
  PyObject* children2 = PyList_New(0);
  PyArrayObject* r2 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_FLOAT32, 1);
  float* pf2 = (float*)PyArray_DATA(r2);
  pf2[0] = 4.0;
  PyObject* version = Py_BuildValue("[sOOs]", "CGNSLibraryVersion", (PyObject*)r2, children2, "CGNSLibraryVersion_t");
  PyList_Append(children1, version); Py_INCREF(version);

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
  PyList_Append(children1, base); Py_INCREF(base);

  // Create Zone
  npy_dim_vals.reserve(2);
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
  PyList_Append(children3, zone); Py_INCREF(zone);

  // Create ZoneType
  PyObject* children9 = PyList_New(0);
  npy_dim_vals[0] = 12;
  PyArrayObject* r9 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_STRING, 1);
  char* pp9 = (char*)PyArray_DATA(r9);
  strcpy(pp9, "Unstructured");
  PyObject* zoneType = Py_BuildValue("[sOOs]", "ZoneType", r9, children9, "ZoneType_t");
  PyList_Append(children4, zoneType); Py_INCREF(zoneType);

  // Get nvertex
  E_Int nvertex = 0; // nbre de vertex
  node = H5Gopen(uc, "Node", H5P_DEFAULT);
  if (node < 0)
  {
    printf("Warning: hdffsdmread: no FS:Mesh/UnstructuredCells/Node node.\n");
    return 1;
  }
  hid_t aid = H5Aopen_by_name(node, ".", "NumberOfCells", H5P_DEFAULT, H5P_DEFAULT);
  hid_t tid = H5Aget_type(aid);
  H5Aread(aid, H5T_NATIVE_INT, &nvertex); // care when long
  pp4[0] = nvertex;
  printf("nvertex=" SF_D_ "\n", nvertex);

  // Create GridCoordinates
  PyObject* children5 = PyList_New(0);
  PyObject* GC = Py_BuildValue("[sOOs]", "GridCoordinates", Py_None, children5, "GridCoordinates_t");
  PyList_Append(children4, GC); Py_INCREF(GC);

  node = H5Gopen(fid, "/FS:Mesh/UnstructuredCells/Datasets", H5P_DEFAULT);
  if (node < 0)
  {
    printf("Warning: hdffsdmread: no FS:Mesh/UnstructuredCells/Datasets node.\n");
    return 1;
  }
  node = H5Gopen(fid, "/FS:Mesh/UnstructuredCells/Datasets/Coordinates", H5P_DEFAULT);
  if (node < 0)
  {
    printf("Warning: hdffsdmread: no FS:Mesh/UnstructuredCells/Datasets/Coordinates node.\n");
    return 1;
  }
  hid_t did = H5Dopen2(node, "Values", H5P_DEFAULT);
  hid_t sid = H5Dget_space(did);
  E_Int ndims = H5Sget_simple_extent_ndims(sid);
  hsize_t dims[5];
  H5Sget_simple_extent_dims(sid, dims, NULL);

  E_Int size = 1;
  for (E_Int i = 0; i < ndims; i++) size *= dims[i];
  printf("ndims=" SF_D_ "\n", ndims);
  printf("size=" SF_D_ "\n", size);
  tid = H5Tcopy(H5T_NATIVE_DOUBLE); H5Tset_precision(tid, 64);
  hid_t yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  hid_t mid = H5S_ALL;
  double* r = new double [size];
  H5Dread(did, yid, mid, sid, H5P_DEFAULT, r);

  // CoordinateX
  npy_dim_vals[0] = nvertex;
  PyArrayObject* xc = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_FLOAT64, 1);
  E_Float* pxc = (E_Float*)PyArray_DATA(xc);
  for (E_Int i = 0; i < nvertex; i++) pxc[i] = r[3*i];
  PyObject* children6 = PyList_New(0);
  PyObject* nxc = Py_BuildValue("[sOOs]", "CoordinateX", xc, children6, "DataArray_t");
  PyList_Append(children5, nxc); Py_INCREF(nxc);

  // CoordinateY
  PyArrayObject* yc = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_FLOAT64, 1);
  E_Float* pyc = (E_Float*)PyArray_DATA(yc);
  for (E_Int i = 0; i < nvertex; i++) pyc[i] = r[3*i+1];
  PyObject* children7 = PyList_New(0);
  PyObject* nyc = Py_BuildValue("[sOOs]", "CoordinateY", yc, children7, "DataArray_t");
  PyList_Append(children5, nyc); Py_INCREF(nyc);

  // CoordinateZ
  PyArrayObject* zc = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_FLOAT64, 1);
  E_Float* pzc = (E_Float*)PyArray_DATA(zc);
  for (E_Int i = 0; i < nvertex; i++) pzc[i] = r[3*i+2];
  PyObject* children8 = PyList_New(0);
  PyObject* nzc = Py_BuildValue("[sOOs]", "CoordinateZ", zc, children8, "DataArray_t");
  PyList_Append(children5, nzc); Py_INCREF(nzc);
  delete [] r;

  // Create GridElements
  E_Int nhexa=0; E_Int ntetra=0; E_Int npenta=0; E_Int npyra=0; E_Int ntri=0; E_Int nquad=0;
  hid_t* chids = HDF.getChildren(uc);
  E_Int nchildren = 0;
  while (chids[nchildren] != -1) { nchildren++; }
  char** names = new char* [nchildren+1];
  for (E_Int i = 0; i < nchildren+1; i++) names[i] = NULL;
#if H5_VERSION_LE(1,11,9)
  H5Literate(uc, H5_INDEX_NAME, H5_ITER_INC, NULL, feed_children_names, (void*)names);
#else
  H5Literate2(uc, H5_INDEX_NAME, H5_ITER_INC, NULL, feed_children_names, (void*)names);
#endif
  hid_t id = 0; E_Int c = 0; PyObject* GE;
  E_Int ncells = 0; E_Int istart=1;
  E_Int* bct; // bc tags
  std::vector<E_Int*> allbct;
  while (id != -1)
  {
    id = chids[c]; char* name = names[c]; c++;
    if (name)
    {
      if (strcmp(name, "Hexa8") == 0) 
      { 
        createGridElements(id, 17, name, istart, GE, ncells, bct);
        PyList_Append(children4, GE); Py_INCREF(GE);
        istart += ncells;
        nhexa = ncells;
        delete [] bct;
      }
      else if (strcmp(name, "Prism6") == 0) 
      { 
        createGridElements(id, 14, name, istart, GE, ncells, bct);
        PyList_Append(children4, GE); Py_INCREF(GE);
        istart += ncells;
        npenta = ncells;
        delete [] bct;
      }
      else if (strcmp(name, "Tetra4") == 0) 
      { 
        createGridElements(id, 10, name, istart, GE, ncells, bct);
        PyList_Append(children4, GE); Py_INCREF(GE);
        istart += ncells;
        ntetra = ncells;
        delete [] bct;
      }
      else if (strcmp(name, "Pyra5") == 0) 
      { 
        createGridElements(id, 12, name, istart, GE, ncells, bct);
        PyList_Append(children4, GE); Py_INCREF(GE);
        istart += ncells;
        npyra = ncells;
        delete [] bct;
      }
      else if (strcmp(name, "Tri3") == 0)
      { 
        createGridElements(id, 5, name, istart, GE, ncells, bct);
        PyList_Append(children4, GE); Py_INCREF(GE);
        istart += ncells;
        ntri = ncells;
        allbct.push_back(bct);
      }
      else if (strcmp(name, "Quad4") == 0) 
      { 
        createGridElements(id, 7, name, istart, GE, ncells, bct);
        PyList_Append(children4, GE); Py_INCREF(GE);
        istart += ncells;
        nquad = ncells;
        allbct.push_back(bct);
      }
    }
  }
  pp4[1] = istart-1;

  // Add zoneBC
  if (allbct.size() > 0)
  {
    // Create ZoneBC
    PyObject* children9 = PyList_New(0);
    PyObject* nzbc = Py_BuildValue("[sOOs]", "ZoneBC", Py_None, children9, "ZoneBC_t");
    PyList_Append(children4, nzbc); Py_INCREF(nzbc);
    
  }

  // Close file
  H5Fclose(fid);
  return 0;
}

//=============================================================================
/* 
   hdffsdmwrite
*/
//=============================================================================
E_Int K_IO::GenIO::hdffsdmwrite(char* file, PyObject* tree)
{
  printf("Error: Converter has been installed without FSDM/HDF support.\n");
  printf("Error: please install libhdf5 first for FSDM/HDF support.\n");
  return 0;
}
