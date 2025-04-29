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
// Create GridElement node for basic elements connectivity
/* ------------------------------------------------------------------------- */
E_Int createGridElements(hid_t id, E_Int eltType, char* name, E_Int istart, PyObject*& GE, E_Int& ncells, 
  int32_t*& bct, E_Int& size, std::map<E_Int, E_Int>& tagmap)
{
  // Get ncells from id
  hid_t aid = H5Aopen_by_name(id, ".", "NumberOfCells", H5P_DEFAULT, H5P_DEFAULT);
  H5Aread(aid, H5T_NATIVE_INT, &ncells); // care when long
  //printf("ncells=" SF_D_ "\n", ncells);

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
  PyList_Append(children1, er); Py_INCREF(er);

  // Element connectivity
  hid_t did = H5Dopen2(id, "Cell2Node", H5P_DEFAULT);
  hid_t sid = H5Dget_space(did);
  E_Int ndims = H5Sget_simple_extent_ndims(sid);
  hsize_t dims[3];
  H5Sget_simple_extent_dims(sid, dims, NULL);
  E_Int size2 = 1;
  for (E_Int i = 0; i < ndims; i++) size2 *= dims[i];
  //printf("ndims=" SF_D_ "\n", ndims);
  //printf("size=" SF_D_ "\n", size);
  npy_dim_vals[0] = size2;
  hid_t tid = H5Tcopy(H5T_NATIVE_INT); H5Tset_precision(tid, 32);
  hid_t yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  hid_t mid = H5S_ALL;
  PyArrayObject* r3 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);  
  int32_t* pp3 = (int32_t*)PyArray_DATA(r3);
  H5Dread(did, yid, mid, sid, H5P_DEFAULT, pp3);
  for (E_Int i = 0; i < size2; i++) pp3[i] += 1;
  PyObject* children3 = PyList_New(0);
  PyObject* ec = Py_BuildValue("[sOOs]", "ElementConnectivity", r3, children3, "DataArray_t");
  PyList_Append(children1, ec); Py_INCREF(ec);

  // Get BC tag if any
  hid_t gid = H5Lexists(id, "CellAttributes", H5P_DEFAULT);
  if (gid > 0 && (eltType == 7 || eltType == 5)) // read BC for Tri and Quad
  {
    gid = H5Gopen(id, "CellAttributes", H5P_DEFAULT);
    did = H5Dopen2(gid, "CADGroupID", H5P_DEFAULT);
    sid = H5Dget_space(did);
    ndims = H5Sget_simple_extent_ndims(sid);
    H5Sget_simple_extent_dims(sid, dims, NULL);
    E_Int size2 = 1;
    for (E_Int i = 0; i < ndims; i++) size2 *= dims[i];
    //printf("ndims=" SF_D_ "\n", ndims);
    //printf("size=" SF_D_ "\n", size2);
    npy_dim_vals[0] = size2;
    tid = H5Tcopy(H5T_NATIVE_INT); H5Tset_precision(tid, 32);
    yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
    mid = H5S_ALL;
    int32_t* bct2 = new int32_t [size2];
    H5Dread(did, yid, mid, sid, H5P_DEFAULT, bct2);
    // count tags
    for (E_Int i = 0; i < size2; i++) 
    {
      if (tagmap.find(bct2[i]) == tagmap.end()) tagmap[bct2[i]] = 1;
      else tagmap[bct2[i]] += 1;
    }
    //for (const auto& pair : tagmap) // for each tag
    //{
    //  E_Int tag = pair.first; // tag
    //  E_Int nfaces = pair.second; // nbre de faces pour ce tag
    //  printf("tag " SF_D_ " is set " SF_D_ " times.\n", tag, nfaces);
    //}
    if (bct == NULL) { bct = bct2; size = size2; }
    else 
    {
      // merge in connectivities order
      int32_t* bct3 = new int32_t [size+size2];
      for (E_Int i = 0; i < size; i++)
      {
        bct3[i] = bct[i];
      }
      for (E_Int i = 0; i < size2; i++)
      {
        bct3[i+size] = bct2[i];
      }
      delete [] bct; delete [] bct2;
      bct = bct3; size = size+size2;
    }
  }
  return 0;
}

/* ------------------------------------------------------------------------- */
// Create GridElement node for polyedral (create NGON and NFACES)
/* ------------------------------------------------------------------------- */
E_Int createGridElementsNGon(hid_t id, E_Int istart, E_Int nvertex, 
  PyObject*& NGON, PyObject*& NFACE, E_Int& nfaces, E_Int& ncells)
{
  E_Int size;
  hsize_t dims[3];
  std::vector<npy_intp> npy_dim_vals(1);

  // read Cell2NodeCounts (cell -> nodes)
  hid_t did = H5Dopen2(id, "Cell2NodeCounts", H5P_DEFAULT);
  hid_t sid = H5Dget_space(did);
  E_Int ndims = H5Sget_simple_extent_ndims(sid);
  H5Sget_simple_extent_dims(sid, dims, NULL);
  ncells = 1;
  for (E_Int i = 0; i < ndims; i++) ncells *= dims[i];
  npy_dim_vals[0] = ncells;
  hid_t tid = H5Tcopy(H5T_NATIVE_INT); H5Tset_precision(tid, 32);
  hid_t yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  hid_t mid = H5S_ALL;
  int32_t* c2nc = new int32_t [ncells];
  H5Dread(did, yid, mid, sid, H5P_DEFAULT, c2nc);
  //printf("I found %d cells.\n", ncells);
  
  // read cell2node (cell -> nodes)  
  did = H5Dopen2(id, "Cell2NodeList", H5P_DEFAULT);
  sid = H5Dget_space(did);
  ndims = H5Sget_simple_extent_ndims(sid);
  H5Sget_simple_extent_dims(sid, dims, NULL);
  size = 1;
  for (E_Int i = 0; i < ndims; i++) size *= dims[i];
  npy_dim_vals[0] = size;
  tid = H5Tcopy(H5T_NATIVE_INT); H5Tset_precision(tid, 32);
  yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  mid = H5S_ALL;
  int32_t* c2n = new int32_t [size];
  H5Dread(did, yid, mid, sid, H5P_DEFAULT, c2n);

  // read cell2facecount (nbre de face par cellule)  
  did = H5Dopen2(id, "Cell2FaceCounts", H5P_DEFAULT);
  sid = H5Dget_space(did);
  ndims = H5Sget_simple_extent_ndims(sid);
  H5Sget_simple_extent_dims(sid, dims, NULL);
  size = 1;
  for (E_Int i = 0; i < ndims; i++) size *= dims[i];
  npy_dim_vals[0] = size;
  tid = H5Tcopy(H5T_NATIVE_INT); H5Tset_precision(tid, 32);
  yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  mid = H5S_ALL;
  int32_t* c2fc = new int32_t [size];
  H5Dread(did, yid, mid, sid, H5P_DEFAULT, c2fc);

  // read face2nodecounts (face -> node cell indirect)  
  did = H5Dopen2(id, "Face2NodeCounts", H5P_DEFAULT);
  sid = H5Dget_space(did);
  ndims = H5Sget_simple_extent_ndims(sid);
  H5Sget_simple_extent_dims(sid, dims, NULL);
  nfaces = 1;
  for (E_Int i = 0; i < ndims; i++) nfaces *= dims[i];
  npy_dim_vals[0] = nfaces;
  tid = H5Tcopy(H5T_NATIVE_INT); H5Tset_precision(tid, 32);
  yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  mid = H5S_ALL;
  int32_t* f2nc = new int32_t [nfaces];
  H5Dread(did, yid, mid, sid, H5P_DEFAULT, f2nc);

  // read face2node (face -> node cell indirect)  
  did = H5Dopen2(id, "Face2NodeList", H5P_DEFAULT);
  sid = H5Dget_space(did);
  ndims = H5Sget_simple_extent_ndims(sid);
  H5Sget_simple_extent_dims(sid, dims, NULL);
  size = 1;
  for (E_Int i = 0; i < ndims; i++) size *= dims[i];
  npy_dim_vals[0] = size;
  tid = H5Tcopy(H5T_NATIVE_INT); H5Tset_precision(tid, 32);
  yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  mid = H5S_ALL;
  int32_t* f2n = new int32_t [size];
  H5Dread(did, yid, mid, sid, H5P_DEFAULT, f2n);

  // Create NGON
  npy_dim_vals[0] = 2;
#ifdef E_DOUBLEINT
  PyArrayObject* r1 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT64, 1);
  int64_t* pp1 = (int64_t*)PyArray_DATA(r1);
#else
  PyArrayObject* r1 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);
  int32_t* pp1 = (int32_t*)PyArray_DATA(r1);
#endif
  pp1[0] = 22; pp1[1] = 0;
  PyObject* children1 = PyList_New(0);
  NGON = Py_BuildValue("[sOOs]", "NGON", r1, children1, "Elements_t");

  // Element range
#ifdef E_DOUBLEINT
  PyArrayObject* r2 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT64, 1);
  int64_t* pp2 = (int64_t*)PyArray_DATA(r2);
#else
  PyArrayObject* r2 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);
  int32_t* pp2 = (int32_t*)PyArray_DATA(r2);
#endif
  pp2[0] = istart; pp2[1] = istart+nfaces-1;
  PyObject* children2 = PyList_New(0);
  PyObject* er = Py_BuildValue("[sOOs]", "ElementRange", r2, children2, "IndexRange_t");
  PyList_Append(children1, er); Py_INCREF(er);

  // Element connectivity
  npy_dim_vals[0] = size;
  PyArrayObject* r3 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);  
  int32_t* pp3 = (int32_t*)PyArray_DATA(r3);
  H5Dread(did, yid, mid, sid, H5P_DEFAULT, pp3);
  PyObject* children3 = PyList_New(0);
  PyObject* ec = Py_BuildValue("[sOOs]", "ElementConnectivity", r3, children3, "DataArray_t");
  PyList_Append(children1, ec); Py_INCREF(ec);

  E_Int cell = 0; // current cell
  E_Int lcount = 0; // local cell count
  E_Int ps = 0; // pointer faces
  E_Int pcell = 0; // pointer cells
  for (E_Int i = 0; i < nfaces; i++) // pour chaque face
  {
    for (E_Int j = 0; j < f2nc[i]; j++) // pour tous les noeuds de la face
    {
      E_Int lnode = f2n[ps]; // local node shift dans la cellule
      E_Int node = c2n[pcell+lnode]; // reel node
      pp3[ps] = node+1;
      ps++;
    }
    lcount++;
    if (lcount >= c2fc[cell]) { pcell += c2nc[cell]; cell++; lcount = 0; }
    //printf("face=%d cell=%d\n", i, cell);
  }
  //printf("=> %d %d\n", ps, size);

  // Element start offset
  npy_dim_vals[0] = nfaces+1;
  PyArrayObject* r4 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);  
  int32_t* pp4 = (int32_t*)PyArray_DATA(r4);
  pp4[0] = 0;
  for (E_Int i = 0; i < nfaces; i++) pp4[i+1] = pp4[i]+f2nc[i];
  PyObject* children4 = PyList_New(0);
  PyObject* eso = Py_BuildValue("[sOOs]", "ElementStartOffset", r4, children4, "DataArray_t");
  PyList_Append(children1, eso); Py_INCREF(eso);

  // Create NFACE
  npy_dim_vals[0] = 2;
#ifdef E_DOUBLEINT
  PyArrayObject* r5 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT64, 1);
  int64_t* pp5 = (int64_t*)PyArray_DATA(r5);
#else
  PyArrayObject* r5 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);
  int32_t* pp5 = (int32_t*)PyArray_DATA(r5);
#endif
  pp5[0] = 23; pp5[1] = 0;
  PyObject* children5 = PyList_New(0);
  NFACE = Py_BuildValue("[sOOs]", "NFACE", r5, children5, "Elements_t");
  
  // Element range
#ifdef E_DOUBLEINT
  PyArrayObject* r8 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT64, 1);
  int64_t* pp8 = (int64_t*)PyArray_DATA(r8);
#else
  PyArrayObject* r8 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);
  int32_t* pp8 = (int32_t*)PyArray_DATA(r8);
#endif
  pp8[0] = istart+nfaces; pp2[1] = istart+nfaces+ncells-1;
  PyObject* children8 = PyList_New(0);
  PyObject* er2 = Py_BuildValue("[sOOs]", "ElementRange", r8, children8, "IndexRange_t");
  PyList_Append(children5, er2); Py_INCREF(er2);

  // rebuild Element start offset
  npy_dim_vals[0] = ncells+1;
  PyArrayObject* r6 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);  
  int32_t* pp6 = (int32_t*)PyArray_DATA(r6);
  PyObject* children6 = PyList_New(0);
  PyObject* eso2 = Py_BuildValue("[sOOs]", "ElementStartOffset", r6, children6, "DataArray_t");
  PyList_Append(children5, eso2); Py_INCREF(eso2);
  pp6[0] = 0;
  for (E_Int i = 0; i < ncells; i++) pp6[i+1] = pp6[i]+c2fc[i];
  
  // rebuild Element connectivity
  size = pp6[ncells];
  npy_dim_vals[0] = size;
  PyArrayObject* r7 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);  
  int32_t* pp7 = (int32_t*)PyArray_DATA(r7);
  PyObject* children7 = PyList_New(0);
  PyObject* ec2 = Py_BuildValue("[sOOs]", "ElementConnectivity", r7, children7, "DataArray_t");
  PyList_Append(children5, ec2); Py_INCREF(ec2);
  E_Int c = 1;
  for (E_Int i = 0; i < ncells; i++)
    for (E_Int j = 0; j < c2fc[i]; j++)
      { pp7[i] = 1; }

  delete [] c2nc; delete [] c2n; delete [] c2fc; 
  delete [] f2nc; delete [] f2n; 
  return 0;

  /*
  // Invert Face2Node in Node2Face
  std::vector< std::vector<E_Int> > node2Face(nvertex);
  ps = 0;
  for (E_Int i = 0; i < nfaces; i++)
  {
    for (E_Int j = 0; j < fc[i]; i++)
    {
      //printf("face=%d, node=%d\n", i, pp3[ps+j]);
      node2Face[pp3[ps+j]-1].push_back(i);
    }
    ps += fc[i];
  }

  delete [] fc;

  // read Cell2NodeCounts
  did = H5Dopen2(id, "Face2NodeCounts", H5P_DEFAULT);
  sid = H5Dget_space(did);
  ndims = H5Sget_simple_extent_ndims(sid);
  H5Sget_simple_extent_dims(sid, dims, NULL);
  ncells = 1;
  for (E_Int i = 0; i < ndims; i++) ncells *= dims[i];
  npy_dim_vals[0] = nfaces;
  tid = H5Tcopy(H5T_NATIVE_INT); H5Tset_precision(tid, 32);
  yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  mid = H5S_ALL;
  int32_t* cc = new int32_t [ncells];
  H5Dread(did, yid, mid, sid, H5P_DEFAULT, cc);
  
  // read Cell2NodeList
  did = H5Dopen2(id, "Face2NodeList", H5P_DEFAULT);
  sid = H5Dget_space(did);
  ndims = H5Sget_simple_extent_ndims(sid);
  H5Sget_simple_extent_dims(sid, dims, NULL);
  size = 1;
  for (E_Int i = 0; i < ndims; i++) size *= dims[i];
  npy_dim_vals[0] = size;
  tid = H5Tcopy(H5T_NATIVE_INT); H5Tset_precision(tid, 32);
  yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  mid = H5S_ALL;
  int32_t* cn = new int32_t [size];
  H5Dread(did, yid, mid, sid, H5P_DEFAULT, cn);

  // get cell2Face
  ps = 0;
  std::vector< std::vector<E_Int> > cell2Face(ncells);
  for (E_Int i = 0; i < ncells; i++) // pour chaque cellule
  {
    std::map<E_Int, E_Int> facemap; // faces intervenants
    for (E_Int j = 0; j < cc[i]; j++) // pour chaque noeud
    {
      E_Int node = cn[ps+j];
      for (size_t k = 0; k < node2Face[node].size(); k++)
      {
        E_Int face = node2Face[node][k];
        if (facemap.find(face) == facemap.end()) facemap[face] = 1;
        else facemap[face] += 1;
      }
    }
    ps += cc[i];
    for (const auto& pair : facemap) // for each tag
    {
      E_Int face = pair.first;
      E_Int nn = pair.second;
      if (nn >= 3) cell2Face[i].push_back(face);
    }
  }

  delete [] cc; delete [] cn;  
  return 0;
  */
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
  std::vector<npy_intp> npy_dim_vals(2);
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
  H5Gclose(node);
  //printf("nvertex=" SF_D_ "\n", nvertex);

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
  H5Gclose(node);
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
  //printf("ndims=" SF_D_ "\n", ndims);
  //printf("size=" SF_D_ "\n", size);
  tid = H5Tcopy(H5T_NATIVE_DOUBLE); H5Tset_precision(tid, 64);
  hid_t yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  hid_t mid = H5S_ALL;
  double* r = new double [size];
  H5Dread(did, yid, mid, sid, H5P_DEFAULT, r);
  H5Dclose(did);
  H5Gclose(node);

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
  E_Int nhexa=0; E_Int ntetra=0; E_Int npenta=0; E_Int npyra=0; 
  E_Int ntri=0; E_Int nquad=0;
  E_Int npoly3d=0; E_Int npoly2d=0;
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
  E_Int ncells=0; E_Int istart=1; E_Int n3dcells=0; 
  std::map<E_Int, E_Int> tagmap; // map of bc tags
  int32_t* bct = NULL; size = 0; // bc tag merged
  while (id != -1)
  {
    id = chids[c]; char* name = names[c]; c++;
    if (name)
    {
      if (strcmp(name, "Hexa8") == 0) 
      { 
        createGridElements(id, 17, name, istart, GE, ncells, bct, size, tagmap);
        PyList_Append(children4, GE); Py_INCREF(GE);
        istart += ncells;
        nhexa = ncells;
        n3dcells += ncells;
      }
      else if (strcmp(name, "Prism6") == 0) 
      { 
        createGridElements(id, 14, name, istart, GE, ncells, bct, size, tagmap);
        PyList_Append(children4, GE); Py_INCREF(GE);
        istart += ncells;
        npenta = ncells;
        n3dcells += ncells;
      }
      else if (strcmp(name, "Tetra4") == 0) 
      { 
        createGridElements(id, 10, name, istart, GE, ncells, bct, size, tagmap);
        PyList_Append(children4, GE); Py_INCREF(GE);
        istart += ncells;
        ntetra = ncells;
        n3dcells += ncells;
      }
      else if (strcmp(name, "Pyra5") == 0) 
      { 
        createGridElements(id, 12, name, istart, GE, ncells, bct, size, tagmap);
        PyList_Append(children4, GE); Py_INCREF(GE);
        istart += ncells;
        npyra = ncells;
        n3dcells += ncells;
      }
      else if (strcmp(name, "Tri3") == 0)
      { 
        createGridElements(id, 5, name, istart, GE, ncells, bct, size, tagmap);
        PyList_Append(children4, GE); Py_INCREF(GE);
        istart += ncells;
        ntri = ncells;
      }
      else if (strcmp(name, "Quad4") == 0) 
      { 
        createGridElements(id, 7, name, istart, GE, ncells, bct, size, tagmap);
        PyList_Append(children4, GE); Py_INCREF(GE);
        istart += ncells;
        nquad = ncells;
      }
      else if (strcmp(name, "Poly3D") == 0) 
      { 
        PyObject* NGON; PyObject* NFACE; E_Int nfaces;
        createGridElementsNGon(id, istart, nvertex, NGON, NFACE, nfaces, ncells);
        PyList_Append(children4, NGON); Py_INCREF(NGON);
        //PyList_Append(children4, NFACE); Py_INCREF(NFACE);
        istart += ncells+nfaces;
        npoly3d = ncells;
        n3dcells += ncells;
      }
    }
  }
  pp4[1] = n3dcells; // pas bon

  // Add zoneBC
  if (size > 0)
  {
    // Create ZoneBC
    PyObject* children9 = PyList_New(0);
    PyObject* nzbc = Py_BuildValue("[sOOs]", "ZoneBC", Py_None, children9, "ZoneBC_t");
    PyList_Append(children4, nzbc); Py_INCREF(nzbc);
    
    // Build BCs
    char bcname[256]; char bctype[256];
    size_t offset = ntetra + nhexa + npenta + npyra;
    for (const auto& pair : tagmap) // for each tag
    {
      E_Int tag = pair.first; // tag
      E_Int nfaces = pair.second; // nbre de faces pour ce tag
      //printf("tag " SF_D_ " is set " SF_D_ " times.\n", tag, nfaces);
      if (nfaces == 0) continue;
      // Create BC_t for tag
      PyObject* children10 = PyList_New(0);
      strcpy(bctype, "FamilySpecified");
      sprintf(bcname, "BC%d", tag);
      npy_dim_vals[0] = strlen(bctype);
      PyArrayObject* r10 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_STRING, 1);
      char* pp10 = (char*)PyArray_DATA(r10);
      strcpy(pp10, bctype);
      PyObject* bc = Py_BuildValue("[sOOs]", bcname, r10, children10, "BC_t");
      PyList_Append(children9, bc); Py_INCREF(bc);
      // Create point list
      npy_dim_vals[0] = 1;
      npy_dim_vals[1] = nfaces;
      PyArrayObject* r11 = (PyArrayObject*)PyArray_EMPTY(2, &npy_dim_vals[0], NPY_INT, 1);
      int32_t* pp11 = (int32_t*)PyArray_DATA(r11);
      E_Int c = 0;
      for (E_Int i = 0; i < size; i++)
      {
        if (bct[i] == tag) { pp11[c] = i+1+offset; c++; }
      }
      PyObject* children11 = PyList_New(0);
      PyObject* pl = Py_BuildValue("[sOOs]", "PointList", r11, children11, "IndexArray_t");
      PyList_Append(children10, pl); Py_INCREF(pl);
      // Create GridLocation
      npy_dim_vals[0] = 10;
      PyArrayObject* r12 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_STRING, 1);
      char* pp12 = (char*)PyArray_DATA(r12);
      strcpy(pp12, "FaceCenter");
      PyObject* children12 = PyList_New(0);
      PyObject* gl = Py_BuildValue("[sOOs]", "GridLocation", r12, children12, "GridLocation_t");
      PyList_Append(children10, gl); Py_INCREF(gl);
      // Create BC FamilyName
      PyObject* children13 = PyList_New(0);
      sprintf(bctype, "BCType%d", tag);
      npy_dim_vals[0] = strlen(bctype);
      PyArrayObject* r13 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_STRING, 1);
      char* pp13 = (char*)PyArray_DATA(r13);
      strcpy(pp13, bctype);
      PyObject* famName = Py_BuildValue("[sOOs]", "FamilyName", r13, children13, "FamilyName_t");
      PyList_Append(children10, famName); Py_INCREF(famName);
    }
    delete [] bct;
  }

  // Close nodes
  H5Gclose(uc);

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
