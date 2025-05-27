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
// return new index when input is global orig index in connectivities
/* ------------------------------------------------------------------------- */
E_Int newIndex2(E_Int index, 
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
E_Int createGridElements(hid_t id, E_Int eltType, const char* name, E_Int istart, PyObject*& GE, E_Int& ncells, 
  int32_t*& bct, E_Int& size, std::map<E_Int, E_Int>& tagmap)
{
  // Get ncells from id
  hid_t aid = H5Aopen_by_name(id, ".", "NumberOfCells", H5P_DEFAULT, H5P_DEFAULT);
  H5Aread(aid, H5T_NATIVE_INT, &ncells); // care when long
  H5Aclose(aid);
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
  PyList_Append(children1, er); Py_DECREF(er);

  // Element connectivity
  hid_t did = H5Dopen2(id, "Cell2Node", H5P_DEFAULT);
  hid_t sid = H5Dget_space(did);
  E_Int ndims = H5Sget_simple_extent_ndims(sid);
  hsize_t dims[3];
  H5Sget_simple_extent_dims(sid, dims, NULL);
  E_Int size2 = 1;
  for (E_Int i = 0; i < ndims; i++) size2 *= dims[i];
  npy_dim_vals[0] = size2;
  hid_t mid = H5S_ALL;
#ifdef E_DOUBLEINT
  PyArrayObject* r3 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT64, 1);  
  int64_t* pp3 = (int64_t*)PyArray_DATA(r3);
  hid_t tid = H5Tcopy(H5T_NATIVE_INT64); H5Tset_precision(tid, 64); // mem data type
#else
  PyArrayObject* r3 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);  
  int32_t* pp3 = (int32_t*)PyArray_DATA(r3);
  hid_t tid = H5Tcopy(H5T_NATIVE_INT); H5Tset_precision(tid, 32);
#endif
  H5Dread(did, tid, mid, sid, H5P_DEFAULT, pp3);
  H5Dclose(did);
  for (E_Int i = 0; i < size2; i++) pp3[i] += 1;
  PyObject* children3 = PyList_New(0);
  PyObject* ec = Py_BuildValue("[sOOs]", "ElementConnectivity", r3, children3, "DataArray_t");
  PyList_Append(children1, ec); Py_DECREF(ec);

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
    
    npy_dim_vals[0] = size2;
    tid = H5Tcopy(H5T_NATIVE_INT); H5Tset_precision(tid, 32);
    mid = H5S_ALL;
    // bct always int_32
    int32_t* bct2 = new int32_t [size2];
    H5Dread(did, tid, mid, sid, H5P_DEFAULT, bct2);
    H5Dclose(did); H5Gclose(gid);
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
    if (bct == NULL) { bct = bct2; size = size2; } // first pass
    else // second pass
    {
      // merge in fixed order (TRI then QUAD)
      int32_t* bct3 = new int32_t [size+size2];
      if (eltType == 5) // if second pass is TRI
      {
        for (E_Int i = 0; i < size2; i++) bct3[i] = bct[i];
        for (E_Int i = 0; i < size; i++) bct3[i+size2] = bct2[i];
      }
      else
      {
        for (E_Int i = 0; i < size; i++) bct3[i] = bct[i];
        for (E_Int i = 0; i < size2; i++) bct3[i+size] = bct2[i];
      }
      delete [] bct; delete [] bct2;
      bct = bct3; size = size+size2;
    }
  }
  return 0;
}

/* ------------------------------------------------------------------------- */
// Create GridElement node for polyedral connectivity (create NGON and NFACES)
/* ------------------------------------------------------------------------- */
E_Int createGridElementsNGon(hid_t id, E_Int istart, E_Int nvertex, 
  PyObject*& NGON, PyObject*& NFACE, E_Int& nfaces, E_Int& ncells)
{
  E_Int size;
  hsize_t dims[3];
  std::vector<npy_intp> npy_dim_vals(1);
  hid_t tid;

  // read Cell2NodeCounts (cell -> nodes)
  hid_t did = H5Dopen2(id, "Cell2NodeCounts", H5P_DEFAULT);
  hid_t sid = H5Dget_space(did);
  E_Int ndims = H5Sget_simple_extent_ndims(sid);
  H5Sget_simple_extent_dims(sid, dims, NULL);
  ncells = 1;
  for (E_Int i = 0; i < ndims; i++) ncells *= dims[i];
  npy_dim_vals[0] = ncells;
#ifdef E_DOUBLEINT
  tid = H5Tcopy(H5T_NATIVE_INT64); H5Tset_precision(tid, 64);
#else
  tid = H5Tcopy(H5T_NATIVE_INT); H5Tset_precision(tid, 32);
#endif
  hid_t mid = H5S_ALL;
  E_Int* c2nc = new E_Int [ncells];
  H5Dread(did, tid, mid, sid, H5P_DEFAULT, c2nc);
  H5Dclose(did);

  // read cell2node (cell -> nodes)  
  did = H5Dopen2(id, "Cell2NodeList", H5P_DEFAULT);
  sid = H5Dget_space(did);
  ndims = H5Sget_simple_extent_ndims(sid);
  H5Sget_simple_extent_dims(sid, dims, NULL);
  size = 1;
  for (E_Int i = 0; i < ndims; i++) size *= dims[i];
  npy_dim_vals[0] = size;
#ifdef E_DOUBLEINT  
  tid = H5Tcopy(H5T_NATIVE_INT64); H5Tset_precision(tid, 64);
#else
  tid = H5Tcopy(H5T_NATIVE_INT); H5Tset_precision(tid, 32);
#endif
  mid = H5S_ALL;
  E_Int* c2n = new E_Int [size];
  H5Dread(did, tid, mid, sid, H5P_DEFAULT, c2n);
  H5Dclose(did);

  // read cell2facecount (nbre de faces par cellule)  
  did = H5Dopen2(id, "Cell2FaceCounts", H5P_DEFAULT);
  sid = H5Dget_space(did);
  ndims = H5Sget_simple_extent_ndims(sid);
  H5Sget_simple_extent_dims(sid, dims, NULL);
  size = 1;
  for (E_Int i = 0; i < ndims; i++) size *= dims[i];
  npy_dim_vals[0] = size;
#ifdef E_DOUBLEINT  
  tid = H5Tcopy(H5T_NATIVE_INT64); H5Tset_precision(tid, 64);
#else
  tid = H5Tcopy(H5T_NATIVE_INT); H5Tset_precision(tid, 32);
#endif
  mid = H5S_ALL;
  E_Int* c2fc = new E_Int [size];
  H5Dread(did, tid, mid, sid, H5P_DEFAULT, c2fc);
  H5Dclose(did);

  // Check total number of faces and number of faces per element * number of elts
  E_Int size1 = 0;
  for (E_Int i = 0; i < size; i++) size1 += c2fc[i];
  //printf("total number of faces by elements: " SF_D_ "\n", size1);

  // read face2nodecounts (face -> node cell indirect)  
  did = H5Dopen2(id, "Face2NodeCounts", H5P_DEFAULT);
  sid = H5Dget_space(did);
  ndims = H5Sget_simple_extent_ndims(sid);
  H5Sget_simple_extent_dims(sid, dims, NULL);
  nfaces = 1;
  for (E_Int i = 0; i < ndims; i++) nfaces *= dims[i];
  npy_dim_vals[0] = nfaces;
#ifdef E_DOUBLEINT  
  tid = H5Tcopy(H5T_NATIVE_INT64); H5Tset_precision(tid, 64);
#else
  tid = H5Tcopy(H5T_NATIVE_INT); H5Tset_precision(tid, 32);
#endif
  mid = H5S_ALL;
  E_Int* f2nc = new E_Int [nfaces];
  H5Dread(did, tid, mid, sid, H5P_DEFAULT, f2nc);
  H5Dclose(did);
  //printf("number of faces (duplicated): " SF_D_ "\n", nfaces);

  // read face2node (face -> node cell indirect from elt)
  did = H5Dopen2(id, "Face2NodeList", H5P_DEFAULT);
  sid = H5Dget_space(did);
  ndims = H5Sget_simple_extent_ndims(sid);
  H5Sget_simple_extent_dims(sid, dims, NULL);
  size = 1;
  for (E_Int i = 0; i < ndims; i++) size *= dims[i];
  npy_dim_vals[0] = size;
#ifdef E_DOUBLEINT  
  tid = H5Tcopy(H5T_NATIVE_INT64); H5Tset_precision(tid, 64);
#else
  tid = H5Tcopy(H5T_NATIVE_INT); H5Tset_precision(tid, 32);
#endif
  mid = H5S_ALL;
  E_Int* f2n = new E_Int [size];
  H5Dread(did, tid, mid, sid, H5P_DEFAULT, f2n);
  H5Dclose(did);

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

  // Element range du NGON
  npy_dim_vals[0] = 2;
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
  PyList_Append(children1, er); Py_DECREF(er);

  // Element connectivity du NGON
  npy_dim_vals[0] = size;
#ifdef E_DOUBLEINT
  PyArrayObject* r3 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT64, 1);  
  int64_t* pp3 = (int64_t*)PyArray_DATA(r3);
#else
  PyArrayObject* r3 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);  
  int32_t* pp3 = (int32_t*)PyArray_DATA(r3);
#endif
  //H5Dread(did, tid, mid, sid, H5P_DEFAULT, pp3);
  for (E_Int i = 0; i < size; i++) pp3[i] = f2n[i];
  PyObject* children3 = PyList_New(0);
  PyObject* ec = Py_BuildValue("[sOOs]", "ElementConnectivity", r3, children3, "DataArray_t");
  PyList_Append(children1, ec); Py_DECREF(ec);

  E_Int cell = 0; // current cell
  E_Int lcount = 0; // local cell count
  E_Int ps = 0; // pointer faces pp3
  E_Int pcell = 0; // pointer cells
  for (E_Int i = 0; i < nfaces; i++) // pour chaque face dup
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
  }

  // Element start offset du NGON
  npy_dim_vals[0] = nfaces+1;
#ifdef E_DOUBLEINT
  PyArrayObject* r4 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT64, 1);  
  int64_t* pp4 = (int64_t*)PyArray_DATA(r4);
#else
  PyArrayObject* r4 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);  
  int32_t* pp4 = (int32_t*)PyArray_DATA(r4);
#endif
  pp4[0] = 0;
  for (E_Int i = 0; i < nfaces; i++) pp4[i+1] = pp4[i]+f2nc[i];
  PyObject* children4 = PyList_New(0);
  PyObject* eso = Py_BuildValue("[sOOs]", "ElementStartOffset", r4, children4, "DataArray_t");
  PyList_Append(children1, eso); Py_DECREF(eso);

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
  
  // Element range du NFACE
  npy_dim_vals[0] = 2;
#ifdef E_DOUBLEINT
  PyArrayObject* r8 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT64, 1);
  int64_t* pp8 = (int64_t*)PyArray_DATA(r8);
#else
  PyArrayObject* r8 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);
  int32_t* pp8 = (int32_t*)PyArray_DATA(r8);
#endif
  pp8[0] = istart+nfaces; pp8[1] = istart+nfaces+ncells-1;
  PyObject* children8 = PyList_New(0);
  PyObject* er2 = Py_BuildValue("[sOOs]", "ElementRange", r8, children8, "IndexRange_t");
  PyList_Append(children5, er2); Py_DECREF(er2);

  // rebuild Element start offset du NFACE
  npy_dim_vals[0] = ncells+1;
#ifdef E_DOUBLEINT
  PyArrayObject* r6 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT64, 1);  
  int64_t* pp6 = (int64_t*)PyArray_DATA(r6);
#else
  PyArrayObject* r6 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);  
  int32_t* pp6 = (int32_t*)PyArray_DATA(r6);
#endif
  PyObject* children6 = PyList_New(0);
  PyObject* eso2 = Py_BuildValue("[sOOs]", "ElementStartOffset", r6, children6, "DataArray_t");
  PyList_Append(children5, eso2); Py_DECREF(eso2);
  pp6[0] = 0;
  for (E_Int i = 0; i < ncells; i++) pp6[i+1] = pp6[i]+c2fc[i];
  
  // rebuild Element connectivity du NFACE
  size = pp6[ncells];
  npy_dim_vals[0] = size;
#ifdef E_DOUBLEINT
  PyArrayObject* r7 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT64, 1);
  int64_t* pp7 = (int64_t*)PyArray_DATA(r7);
#else
  PyArrayObject* r7 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_INT32, 1);  
  int32_t* pp7 = (int32_t*)PyArray_DATA(r7);
#endif
  PyObject* children7 = PyList_New(0);
  PyObject* ec2 = Py_BuildValue("[sOOs]", "ElementConnectivity", r7, children7, "DataArray_t");
  PyList_Append(children5, ec2); Py_DECREF(ec2);
  E_Int c = 0;
  for (E_Int i = 0; i < ncells; i++)
    for (E_Int j = 0; j < c2fc[i]; j++)
      { pp7[c] = c+1; c++; }

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

  //GenIOHdf HDF;

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
  H5Aclose(aid); H5Gclose(node);
  //printf("nvertex=" SF_D_ "\n", nvertex);

  // Create GridCoordinates
  PyObject* children5 = PyList_New(0);
  PyObject* GC = Py_BuildValue("[sOOs]", "GridCoordinates", Py_None, children5, "GridCoordinates_t");
  PyList_Append(children4, GC); Py_DECREF(GC);

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
  tid = H5Tcopy(H5T_NATIVE_DOUBLE); H5Tset_precision(tid, 64);
  hid_t mid = H5S_ALL;
  double* r = new double [size];
  H5Dread(did, tid, mid, sid, H5P_DEFAULT, r);
  H5Dclose(did); H5Gclose(node);

  // CoordinateX
  npy_dim_vals[0] = nvertex;
  PyArrayObject* xc = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_FLOAT64, 1);
  E_Float* pxc = (E_Float*)PyArray_DATA(xc);
  for (E_Int i = 0; i < nvertex; i++) pxc[i] = r[3*i];
  PyObject* children6 = PyList_New(0);
  PyObject* nxc = Py_BuildValue("[sOOs]", "CoordinateX", xc, children6, "DataArray_t");
  PyList_Append(children5, nxc); Py_DECREF(nxc);

  // CoordinateY
  PyArrayObject* yc = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_FLOAT64, 1);
  E_Float* pyc = (E_Float*)PyArray_DATA(yc);
  for (E_Int i = 0; i < nvertex; i++) pyc[i] = r[3*i+1];
  PyObject* children7 = PyList_New(0);
  PyObject* nyc = Py_BuildValue("[sOOs]", "CoordinateY", yc, children7, "DataArray_t");
  PyList_Append(children5, nyc); Py_DECREF(nyc);

  // CoordinateZ
  PyArrayObject* zc = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_FLOAT64, 1);
  E_Float* pzc = (E_Float*)PyArray_DATA(zc);
  for (E_Int i = 0; i < nvertex; i++) pzc[i] = r[3*i+2];
  PyObject* children8 = PyList_New(0);
  PyObject* nzc = Py_BuildValue("[sOOs]", "CoordinateZ", zc, children8, "DataArray_t");
  PyList_Append(children5, nzc); Py_DECREF(nzc);
  delete [] r;

  // Create FlowSolution
  if (H5Lexists(fid, "/FS:Mesh/UnstructuredCells/Datasets/AugState", H5P_DEFAULT) > 0)
  {
    node = H5Gopen(fid, "/FS:Mesh/UnstructuredCells/Datasets/AugState", H5P_DEFAULT);
    hid_t aid = H5Aopen_by_name(node, ".", "NumberOfVariables", H5P_DEFAULT, H5P_DEFAULT);
    E_Int nvars = 0;
    H5Aread(aid, H5T_NATIVE_INT, &nvars);
    H5Aclose(aid);
    aid = H5Aopen_by_name(node, ".", "NumberOfCells", H5P_DEFAULT, H5P_DEFAULT);
    E_Int ncells = 0;
    H5Aread(aid, H5T_NATIVE_INT, &ncells);
    H5Aclose(aid);
    //printf("nvars = %d, ncells=%d\n", nvars, ncells);

    dims[0] = ncells; dims[1] = nvars;
    did = H5Screate_simple(2, dims, NULL);
    tid = H5Tcopy(H5T_NATIVE_DOUBLE); H5Tset_precision(tid, 64);
    hid_t vid = H5Dopen2(node, "Values", H5P_DEFAULT);
    
    PyObject* childrenFS = PyList_New(0);
    PyObject* FS = Py_BuildValue("[sOOs]", "FlowSolution#Centers", Py_None, childrenFS, "FlowSolution_t");
    PyList_Append(children4, FS); Py_DECREF(FS);

    npy_dim_vals[0] = 10;
    PyObject* children14 = PyList_New(0);
    PyArrayObject* r14 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_STRING, 1);
    char* pp14 = (char*)PyArray_DATA(r14);
    K_STRING::cpy(pp14, "CellCenter", npy_dim_vals[0], false);
    PyObject* gl = Py_BuildValue("[sOOs]", "GridLocation", r14, children14, "GridLocation_t");
    PyList_Append(childrenFS, gl); Py_DECREF(gl);

    hsize_t start[2]; hsize_t scount[2];
    char name[35]; char name2[76];
    for (E_Int n = 0; n < nvars; n++)
    {
      // Read name
      sprintf(name, "Variable" SF_D_, n);
      hid_t gid2 = H5Gopen(node, name, H5P_DEFAULT);
      aid = H5Aopen_by_name(gid2, ".", "Name", H5P_DEFAULT, H5P_DEFAULT);
      hid_t atype = H5Aget_type(aid);
      size = H5Tget_size(atype);
      H5Aread(aid, atype, name2);
      name2[size] = '\0';
      //H5Aread(aid, H5T_C_S1, name2);
      H5Aclose(aid); H5Gclose(gid2);

      // Read values (partial)
      npy_dim_vals[0] = ncells;
      hid_t sid = H5Scopy(did);
      start[0] = 0; start[1] = n;
      scount[0] = ncells; scount[1] = 1;
      H5Sselect_hyperslab(sid, H5S_SELECT_SET, start, NULL, scount, NULL);
      hid_t memspace = H5Screate_simple(2, scount, NULL);
      PyArrayObject* f = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_FLOAT64, 1);
      E_Float* fp = (E_Float*)PyArray_DATA(f);
      H5Dread(vid, tid, memspace, sid, H5P_DEFAULT, fp);
      PyObject* fc = PyList_New(0);
      PyObject* fl = Py_BuildValue("[sOOs]", name2, f, fc, "DataArray_t");
      PyList_Append(childrenFS, fl); Py_DECREF(fl);
    }
    H5Dclose(vid); H5Sclose(did);
    H5Gclose(node);
  }

  // Read MappingCellType2Index
  if (H5Lexists(fid, "/FS:Mesh/UnstructuredCells/MappingCellType2Index", H5P_DEFAULT) > 0)
  {
    char** names = new char* [20+1];
    node = H5Gopen(fid, "/FS:Mesh/UnstructuredCells/MappingCellType2Index", H5P_DEFAULT);
    // Hard coded order: Tetra, Pyra, Penta, Hexa, Tri, Quad.
    H5Gclose(node);
    delete [] names;
  }

  // Create GridElements
  E_Int nhexa=0; E_Int ntetra=0; E_Int npenta=0; E_Int npyra=0; 
  E_Int ntri=0; E_Int nquad=0;
  E_Int npoly3d=0; E_Int npoly2d=0;
  /*
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
  */
  hid_t id = 0; PyObject* GE;
  E_Int ncells=0; E_Int istart=1; E_Int n3dcells=0; 
  std::map<E_Int, E_Int> tagmap; // map of bc tags
  int32_t* bct = NULL; size = 0; // bc tag merged

  if (H5Lexists(fid, "/FS:Mesh/UnstructuredCells/Tetra4", H5P_DEFAULT) > 0)
  {
    id = H5Gopen(fid, "/FS:Mesh/UnstructuredCells/Tetra4", H5P_DEFAULT);
    createGridElements(id, 10, "Tetra4", istart, GE, ncells, bct, size, tagmap);
    PyList_Append(children4, GE); Py_DECREF(GE);
    istart += ncells;
    ntetra = ncells;
    n3dcells += ncells;
    H5Gclose(id);
  }

  if (H5Lexists(fid, "/FS:Mesh/UnstructuredCells/Pyra5", H5P_DEFAULT) > 0)
  {
    id = H5Gopen(fid, "/FS:Mesh/UnstructuredCells/Pyra5", H5P_DEFAULT);
    createGridElements(id, 12, "Pyra5", istart, GE, ncells, bct, size, tagmap);
    PyList_Append(children4, GE); Py_DECREF(GE);
    istart += ncells;
    npyra = ncells;
    n3dcells += ncells;
    H5Gclose(id);
  }

  if (H5Lexists(fid, "/FS:Mesh/UnstructuredCells/Prism6", H5P_DEFAULT) > 0)
  {
    id = H5Gopen(fid, "/FS:Mesh/UnstructuredCells/Prism6", H5P_DEFAULT);
    createGridElements(id, 14, "Prism6", istart, GE, ncells, bct, size, tagmap);
    PyList_Append(children4, GE); Py_DECREF(GE);
    istart += ncells;
    npenta = ncells;
    n3dcells += ncells;
    H5Gclose(id);
  }

  if (H5Lexists(fid, "/FS:Mesh/UnstructuredCells/Hexa8", H5P_DEFAULT) > 0)
  {
    id = H5Gopen(fid, "/FS:Mesh/UnstructuredCells/Hexa8", H5P_DEFAULT);
    createGridElements(id, 17, "Hexa8", istart, GE, ncells, bct, size, tagmap);
    PyList_Append(children4, GE); Py_DECREF(GE);
    istart += ncells;
    nhexa = ncells;
    n3dcells += ncells;
    H5Gclose(id);
  }

  if (H5Lexists(fid, "/FS:Mesh/UnstructuredCells/Tri3", H5P_DEFAULT) > 0)
  {
    id = H5Gopen(fid, "/FS:Mesh/UnstructuredCells/Tri3", H5P_DEFAULT);
    createGridElements(id, 5, "Tri3", istart, GE, ncells, bct, size, tagmap);
    PyList_Append(children4, GE); Py_DECREF(GE);
    istart += ncells;
    ntri = ncells;
    H5Gclose(id);
  }

  if (H5Lexists(fid, "/FS:Mesh/UnstructuredCells/Quad4", H5P_DEFAULT) > 0)
  {
    id = H5Gopen(fid, "/FS:Mesh/UnstructuredCells/Quad4", H5P_DEFAULT);
    createGridElements(id, 7, "Quad4", istart, GE, ncells, bct, size, tagmap);
    PyList_Append(children4, GE); Py_DECREF(GE);
    istart += ncells;
    nquad = ncells;
    H5Gclose(id);
  }

  if (H5Lexists(fid, "/FS:Mesh/UnstructuredCells/Poly3D", H5P_DEFAULT) > 0)
  {
    id = H5Gopen(fid, "/FS:Mesh/UnstructuredCells/Poly3D", H5P_DEFAULT);
    PyObject* NGON; PyObject* NFACE; E_Int nfaces;
    createGridElementsNGon(id, istart, nvertex, NGON, NFACE, nfaces, ncells);
    PyList_Append(children4, NGON); Py_DECREF(NGON);
    PyList_Append(children4, NFACE); Py_DECREF(NFACE);
    istart += ncells+nfaces;
    npoly3d = ncells;
    n3dcells += ncells;
    H5Gclose(id);
  }

  /*
  while (id != -1)
  {
    id = chids[c]; char* name = names[c]; c++;
    if (name)
    {
      if (strcmp(name, "Hexa8") == 0) 
      { 
        createGridElements(id, 17, name, istart, GE, ncells, bct, size, tagmap);
        PyList_Append(children4, GE); Py_DECREF(GE);
        istart += ncells;
        nhexa = ncells;
        n3dcells += ncells;
      }
      else if (strcmp(name, "Prism6") == 0) 
      { 
        createGridElements(id, 14, name, istart, GE, ncells, bct, size, tagmap);
        PyList_Append(children4, GE); Py_DECREF(GE);
        istart += ncells;
        npenta = ncells;
        n3dcells += ncells;
      }
      else if (strcmp(name, "Tetra4") == 0) 
      { 
        createGridElements(id, 10, name, istart, GE, ncells, bct, size, tagmap);
        PyList_Append(children4, GE); Py_DECREF(GE);
        istart += ncells;
        ntetra = ncells;
        n3dcells += ncells;
      }
      else if (strcmp(name, "Pyra5") == 0) 
      { 
        createGridElements(id, 12, name, istart, GE, ncells, bct, size, tagmap);
        PyList_Append(children4, GE); Py_DECREF(GE);
        istart += ncells;
        npyra = ncells;
        n3dcells += ncells;
      }
      else if (strcmp(name, "Tri3") == 0)
      { 
        createGridElements(id, 5, name, istart, GE, ncells, bct, size, tagmap);
        PyList_Append(children4, GE); Py_DECREF(GE);
        istart += ncells;
        ntri = ncells;
      }
      else if (strcmp(name, "Quad4") == 0) 
      { 
        createGridElements(id, 7, name, istart, GE, ncells, bct, size, tagmap);
        PyList_Append(children4, GE); Py_DECREF(GE);
        istart += ncells;
        nquad = ncells;
      }
      else if (strcmp(name, "Poly3D") == 0) 
      { 
        PyObject* NGON; PyObject* NFACE; E_Int nfaces;
        createGridElementsNGon(id, istart, nvertex, NGON, NFACE, nfaces, ncells);
        PyList_Append(children4, NGON); Py_DECREF(NGON);
        PyList_Append(children4, NFACE); Py_DECREF(NFACE);
        istart += ncells+nfaces;
        npoly3d = ncells;
        n3dcells += ncells;
      }
    }
  }
  */

  pp4[1] = n3dcells;

  // Add zoneBC
  if (size > 0)
  {
    // Create ZoneBC
    PyObject* children9 = PyList_New(0);
    PyObject* nzbc = Py_BuildValue("[sOOs]", "ZoneBC", Py_None, children9, "ZoneBC_t");
    PyList_Append(children4, nzbc); Py_DECREF(nzbc);
    
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
      sprintf(bcname, "BC" SF_D_, tag);
      npy_dim_vals[0] = strlen(bctype);
      PyArrayObject* r10 = (PyArrayObject*)PyArray_EMPTY(1, &npy_dim_vals[0], NPY_STRING, 1);
      char* pp10 = (char*)PyArray_DATA(r10);
      K_STRING::cpy(pp10, bctype, npy_dim_vals[0], false);
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
      for (E_Int i = 0; i < size; i++)
      {
        if (bct[i] == tag) { pp11[c] = i+1+offset; c++; }
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
      K_STRING::cpy(pp13, bctype, npy_dim_vals[0], false);
      PyObject* famName = Py_BuildValue("[sOOs]", "FamilyName", r13, children13, "FamilyName_t");
      PyList_Append(children10, famName); Py_DECREF(famName);
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
  if (tree == Py_None)
  {
    // nothing to write
    return 1;
  }

  /* Ouverture du fichier pour l'ecriture */
  hid_t fapl, fid, capl;
  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);  
  H5Pset_libver_bounds(fapl, KHDFVERSION, KHDFVERSION);
  capl = H5Pcreate(H5P_FILE_CREATE);
  H5Pset_link_creation_order(capl, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);
  fid = H5Fcreate(file, H5F_ACC_TRUNC, capl, fapl);
  H5Pclose(fapl); H5Pclose(capl);

  if (fid < 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "hdffsdmwrite: can not open file.");
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
 
  // Get ncells, nvertex
  std::vector<PyArrayObject*> hook;
  E_Int nvertex = K_PYTREE::getNumberOfPointsOfZone(zone, hook);

  // Write
  hid_t gid, gid2, aid, did, tid, vid;
  hsize_t dims[2];
  char* name = new char [35]; 
  
  // write FS:Mesh group and version
  gid = H5Gcreate(fid, "FS:Mesh", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dims[0] = 1;
  did = H5Screate_simple(1, dims, NULL);
  aid = H5Acreate(gid, "Version", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT);
  int32_t version = 2;
  H5Awrite(aid, H5T_NATIVE_INT, &version);
  H5Aclose(aid); H5Sclose(did);

  // Write UnstructuredCells
  hid_t uc = H5Gcreate(gid, "UnstructuredCells", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(gid);

  // Write DataSets
  hid_t ds = H5Gcreate(uc, "Datasets", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Write Coordinates
  hid_t coord = H5Gcreate(ds, "Coordinates", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dims[0] = 1;
  did = H5Screate_simple(1, dims, NULL);
#ifdef E_DOUBLEINT
  aid = H5Acreate(coord, "NumberOfCellTypes", H5T_NATIVE_INT64, did, H5P_DEFAULT, H5P_DEFAULT);
  int64_t numberOfCellTypes = 1;
  H5Awrite(aid, H5T_NATIVE_INT64, &numberOfCellTypes);
#else
  aid = H5Acreate(coord, "NumberOfCellTypes", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT);
  int32_t numberOfCellTypes = 1;
  H5Awrite(aid, H5T_NATIVE_INT, &numberOfCellTypes);
#endif
  H5Aclose(aid); H5Sclose(did);
  dims[0] = 1;
  did = H5Screate_simple(1, dims, NULL);
#ifdef E_DOUBLEINT
  aid = H5Acreate(coord, "NumberOfCells", H5T_NATIVE_INT64, did, H5P_DEFAULT, H5P_DEFAULT);
  int64_t numberOfCells = nvertex;
  H5Awrite(aid, H5T_NATIVE_INT64, &numberOfCells);
#else
  aid = H5Acreate(coord, "NumberOfCells", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT);
  int32_t numberOfCells = nvertex;
  H5Awrite(aid, H5T_NATIVE_INT, &numberOfCells);
#endif
  H5Aclose(aid); H5Sclose(did);
  dims[0] = 1;
  did = H5Screate_simple(1, dims, NULL);
#ifdef E_DOUBLEINT
  aid = H5Acreate(coord, "NumberOfVariables", H5T_NATIVE_INT64, did, H5P_DEFAULT, H5P_DEFAULT);
  int64_t numberOfVariables = 3;
  H5Awrite(aid, H5T_NATIVE_INT64, &numberOfVariables);
#else
  aid = H5Acreate(coord, "NumberOfVariables", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT);
  int32_t numberOfVariables = 3;
  H5Awrite(aid, H5T_NATIVE_INT, &numberOfVariables);
#endif
  H5Aclose(aid); H5Sclose(did);
  dims[0] = 1;
  did = H5Screate_simple(1, dims, NULL);
#ifdef E_DOUBLEINT
  aid = H5Acreate(coord, "SpansAllCells", H5T_NATIVE_INT64, did, H5P_DEFAULT, H5P_DEFAULT);
  int64_t span = 1;
  H5Awrite(aid, H5T_NATIVE_INT64, &span);
#else
  aid = H5Acreate(coord, "SpansAllCells", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT);
  int32_t span = 1;
  H5Awrite(aid, H5T_NATIVE_INT, &span);
#endif
  H5Aclose(aid); H5Sclose(did);

  // Create CellType0
  gid = H5Gcreate(coord, "CellType0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  did = H5Screate(H5S_SCALAR);
  tid = H5Tcopy(H5T_C_S1); H5Tset_size(tid, 5);
  aid = H5Acreate(gid, "Name", tid, did, H5P_DEFAULT, H5P_DEFAULT);
  strcpy(name, "Node"); H5Awrite(aid, tid, name);
  H5Aclose(aid); H5Sclose(did); H5Tclose(tid);
  H5Gclose(gid);

  // Create Variables
  gid = H5Gcreate(coord, "Variable0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  did = H5Screate(H5S_SCALAR);
  tid = H5Tcopy(H5T_C_S1); H5Tset_size(tid, 12);
  //H5Tset_size(tid, H5T_VARIABLE);
  aid = H5Acreate(gid, "Name", tid, did, H5P_DEFAULT, H5P_DEFAULT);
  strcpy(name, "CoordinateX"); H5Awrite(aid, tid, name);
  H5Aclose(aid); H5Sclose(did); H5Tclose(tid);
  gid2 = H5Gcreate(gid, "DataSpecification", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  did = H5Screate(H5S_SCALAR); 
  tid = H5Tcopy(H5T_C_S1); H5Tset_size(tid, 11);
  aid = H5Acreate(gid2, "InfoString", tid, did, H5P_DEFAULT, H5P_DEFAULT);
  strcpy(name, "[m]"); H5Awrite(aid, tid, name);
  H5Aclose(aid); H5Sclose(did); H5Tclose(tid);
  did = H5Screate(H5S_SCALAR); 
  tid = H5Tcopy(H5T_C_S1); H5Tset_size(tid, 12);
  aid = H5Acreate(gid2, "Type", tid, did, H5P_DEFAULT, H5P_DEFAULT);
  strcpy(name, "Dimensional"); H5Awrite(aid, tid, name);
  H5Aclose(aid); H5Sclose(did); H5Tclose(tid); 
  H5Gclose(gid2); H5Gclose(gid);

  gid = H5Gcreate(coord, "Variable1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  did = H5Screate(H5S_SCALAR);
  tid = H5Tcopy(H5T_C_S1); H5Tset_size(tid, 12);
  aid = H5Acreate(gid, "Name", tid, did, H5P_DEFAULT, H5P_DEFAULT);
  strcpy(name, "CoordinateY"); H5Awrite(aid, tid, name);
  H5Aclose(aid); H5Sclose(did); H5Tclose(tid);
  gid2 = H5Gcreate(gid, "DataSpecification", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  did = H5Screate(H5S_SCALAR); 
  tid = H5Tcopy(H5T_C_S1); H5Tset_size(tid, 4);
  aid = H5Acreate(gid2, "InfoString", tid, did, H5P_DEFAULT, H5P_DEFAULT);
  strcpy(name, "[m]"); H5Awrite(aid, tid, name);
  H5Aclose(aid); H5Sclose(did); H5Tclose(tid);
  did = H5Screate(H5S_SCALAR); 
  tid = H5Tcopy(H5T_C_S1); H5Tset_size(tid, 12);
  aid = H5Acreate(gid2, "Type", tid, did, H5P_DEFAULT, H5P_DEFAULT);
  strcpy(name, "Dimensional"); H5Awrite(aid, tid, name);
  H5Aclose(aid); H5Sclose(did); H5Tclose(tid);
  H5Gclose(gid2); H5Gclose(gid);
  
  gid = H5Gcreate(coord, "Variable2", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  did = H5Screate(H5S_SCALAR);
  tid = H5Tcopy(H5T_C_S1); H5Tset_size(tid, 12);
  aid = H5Acreate(gid, "Name", tid, did, H5P_DEFAULT, H5P_DEFAULT);
  strcpy(name, "CoordinateZ"); H5Awrite(aid, tid, name);
  H5Aclose(aid); H5Sclose(did); H5Tclose(tid);
  gid2 = H5Gcreate(gid, "DataSpecification", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  did = H5Screate(H5S_SCALAR); 
  tid = H5Tcopy(H5T_C_S1); H5Tset_size(tid, 11);
  aid = H5Acreate(gid2, "InfoString", tid, did, H5P_DEFAULT, H5P_DEFAULT);
  strcpy(name, "[m]"); H5Awrite(aid, tid, name);
  H5Aclose(aid); H5Sclose(did); H5Tclose(tid);
  did = H5Screate(H5S_SCALAR); 
  tid = H5Tcopy(H5T_C_S1); H5Tset_size(tid, 12);
  aid = H5Acreate(gid2, "Type", tid, did, H5P_DEFAULT, H5P_DEFAULT);
  strcpy(name, "Dimensional"); H5Awrite(aid, tid, name);
  H5Aclose(aid); H5Sclose(did); H5Tclose(tid);
  H5Gclose(gid2); H5Gclose(gid);

  // Write Coordinates/Values (compact)
  PyObject* gc = K_PYTREE::getNodeFromName1(zone, "GridCoordinates");
    
  PyObject* xc = K_PYTREE::getNodeFromName1(gc, "CoordinateX");
  E_Float* xcv = K_PYTREE::getValueAF(xc, hook);
  PyObject* yc = K_PYTREE::getNodeFromName1(gc, "CoordinateY");
  E_Float* ycv = K_PYTREE::getValueAF(yc, hook);
  PyObject* zc = K_PYTREE::getNodeFromName1(gc, "CoordinateZ");
  E_Float* zcv = K_PYTREE::getValueAF(zc, hook);
  double* data = new double [nvertex*3];
  for (E_Int i = 0; i < nvertex; i++)
  {
    data[3*i] = xcv[i];
    data[3*i+1] = ycv[i];
    data[3*i+2] = zcv[i];
  }
  dims[0] = nvertex; dims[1] = 3;
  did = H5Screate_simple(2, dims, NULL);
  vid = H5Dcreate(coord, "Values", H5T_NATIVE_DOUBLE, did, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(vid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Dclose(vid); H5Sclose(did);
  delete [] data;
  H5Gclose(coord); H5Gclose(ds);

  // Write AugmentedState
  PyObject* FS = K_PYTREE::getNodeFromName1(zone, "FlowSolution#Centers");
  if (FS != NULL)
  { /* a finir */

  }

  // Write Node
  gid = H5Gcreate(uc, "Node", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dims[0] = 1;
  did = H5Screate_simple(1, dims, NULL);
  aid = H5Acreate(gid, "NumberOfCells", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT);
  numberOfCells = nvertex;
  H5Awrite(aid, H5T_NATIVE_INT, &numberOfCells);
  H5Aclose(aid); H5Sclose(did);
  //gid2 = H5Gcreate(gid, "CellAttributes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  //dims[0] = nvertex;
  //did = H5Screate_simple(1, dims, NULL);
  //vid = H5Dcreate(gid2, "GlobalNumber", H5T_NATIVE_DOUBLE, did, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  //int32_t* numbers = new int32_t [nvertex];
  //for (E_Int i = 0; i < nvertex; i++) numbers[i] = i;
  //H5Dwrite(vid, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, numbers);
  //delete [] numbers;
  //H5Dclose(vid); H5Sclose(did); H5Gclose(gid2); 
  H5Gclose(gid);

  // Write and merge connectivities
  std::vector<PyObject*> connects;
  K_PYTREE::getNodesFromType1(zone, "Elements_t", connects);

  E_Int nhexa=0; E_Int ntetra=0; E_Int npenta=0; E_Int npyra=0; E_Int ntri=0; E_Int nquad=0;

  // count elements
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

  // Create connectivities
  hid_t vhexa=0, vtetra=0, vpenta=0, vpyra=0, vtri=0, vquad=0; // connect dataset
  hid_t dhexa=0, dtetra=0, dpenta=0, dpyra=0, dtri=0, dquad=0; // connect dims
  hid_t mtri=0, mquad=0; // markers dataset

  if (ntetra > 0) // TETRA
  {
    gid = H5Gcreate(uc, "Tetra4", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dims[0] = 1;
    did = H5Screate_simple(1, dims, NULL);
#ifdef E_DOUBLEINT
    aid = H5Acreate(gid, "NumberOfCells", H5T_NATIVE_INT64, did, H5P_DEFAULT, H5P_DEFAULT);
    int64_t numberOfCells = ntetra;
    H5Awrite(aid, H5T_NATIVE_INT64, &numberOfCells);
#else
    aid = H5Acreate(gid, "NumberOfCells", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT);
    int32_t numberOfCells = ntetra;
    H5Awrite(aid, H5T_NATIVE_INT, &numberOfCells);
#endif
    H5Aclose(aid); H5Sclose(did);
    dims[0] = ntetra; dims[1] = 4;
    dtetra = H5Screate_simple(2, dims, NULL);
#ifdef E_DOUBLEINT
    vtetra = H5Dcreate(gid, "Cell2Node", H5T_NATIVE_INT64, dtetra, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
    vtetra = H5Dcreate(gid, "Cell2Node", H5T_NATIVE_INT, dtetra, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
    H5Gclose(gid);
  }
  if (npyra > 0) // PYRA
  {
    gid = H5Gcreate(uc, "Pyra5", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dims[0] = 1;
    did = H5Screate_simple(1, dims, NULL);
#ifdef E_DOUBLEINT
    aid = H5Acreate(gid, "NumberOfCells", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT);
    int64_t numberOfCells = npyra;
    H5Awrite(aid, H5T_NATIVE_INT64, &numberOfCells);
#else
    aid = H5Acreate(gid, "NumberOfCells", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT);
    int32_t numberOfCells = npyra;
    H5Awrite(aid, H5T_NATIVE_INT, &numberOfCells);
#endif
    H5Aclose(aid); H5Sclose(did);
    dims[0] = npyra; dims[1] = 5;
    dpyra = H5Screate_simple(2, dims, NULL);
#ifdef E_DOUBLEINT
    vpyra = H5Dcreate(gid, "Cell2Node", H5T_NATIVE_INT64, dpyra, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
    vpyra = H5Dcreate(gid, "Cell2Node", H5T_NATIVE_INT, dpyra, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
    H5Gclose(gid);
  }
  if (npenta > 0) // PENTA
  {
    gid = H5Gcreate(uc, "Prism6", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dims[0] = 1;
    did = H5Screate_simple(1, dims, NULL);
#ifdef E_DOUBLEINT
    aid = H5Acreate(gid, "NumberOfCells", H5T_NATIVE_INT64, did, H5P_DEFAULT, H5P_DEFAULT);
    int64_t numberOfCells = npenta;
    H5Awrite(aid, H5T_NATIVE_INT64, &numberOfCells);
#else
    aid = H5Acreate(gid, "NumberOfCells", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT);
    int32_t numberOfCells = npenta;
    H5Awrite(aid, H5T_NATIVE_INT, &numberOfCells);
#endif
    H5Aclose(aid); H5Sclose(did);
    dims[0] = npenta; dims[1] = 6;
    dpenta = H5Screate_simple(2, dims, NULL);
#ifdef E_DOUBLEINT
    vpenta = H5Dcreate(gid, "Cell2Node", H5T_NATIVE_INT64, dpenta, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
    vpenta = H5Dcreate(gid, "Cell2Node", H5T_NATIVE_INT, dpenta, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
    H5Gclose(gid);
  }
  if (nhexa > 0) // HEXA
  {
    gid = H5Gcreate(uc, "Hexa8", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dims[0] = 1;
    did = H5Screate_simple(1, dims, NULL);
#ifdef E_DOUBLEINT
    aid = H5Acreate(gid, "NumberOfCells", H5T_NATIVE_INT64, did, H5P_DEFAULT, H5P_DEFAULT);
    int64_t numberOfCells = nhexa;
    H5Awrite(aid, H5T_NATIVE_INT64, &numberOfCells);
#else
    aid = H5Acreate(gid, "NumberOfCells", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT);
    int32_t numberOfCells = nhexa;
    H5Awrite(aid, H5T_NATIVE_INT, &numberOfCells);
#endif
    H5Aclose(aid); H5Sclose(did);
    dims[0] = nhexa; dims[1] = 8;
    dhexa = H5Screate_simple(2, dims, NULL);
#ifdef E_DOUBLEINT
    vhexa = H5Dcreate(gid, "Cell2Node", H5T_NATIVE_INT64, dhexa, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
    vhexa = H5Dcreate(gid, "Cell2Node", H5T_NATIVE_INT, dhexa, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
    H5Gclose(gid);
  }
  if (ntri > 9) // TRI
  {
    gid = H5Gcreate(uc, "Tri3", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dims[0] = 1;
    did = H5Screate_simple(1, dims, NULL);
#ifdef E_DOUBLEINT
    aid = H5Acreate(gid, "NumberOfCells", H5T_NATIVE_INT64, did, H5P_DEFAULT, H5P_DEFAULT);
    int64_t numberOfCells = ntri;
    H5Awrite(aid, H5T_NATIVE_INT64, &numberOfCells);
#else
    aid = H5Acreate(gid, "NumberOfCells", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT);
    int32_t numberOfCells = ntri;
    H5Awrite(aid, H5T_NATIVE_INT, &numberOfCells);
#endif
    H5Aclose(aid); H5Sclose(did);
    dims[0] = ntri; dims[1] = 3;
    dtri = H5Screate_simple(2, dims, NULL);
#ifdef E_DOUBLEINT
    vtri = H5Dcreate(gid, "Cell2Node", H5T_NATIVE_INT64, dtri, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
    vtri = H5Dcreate(gid, "Cell2Node", H5T_NATIVE_INT, dtri, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
    gid2 = H5Gcreate(gid, "CellAttributes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dims[0] = ntri;
    did = H5Screate_simple(1, dims, NULL);
    mtri = H5Dcreate(gid2, "CADGroupID", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(gid); H5Gclose(gid2); H5Sclose(did);
  } 
  if (nquad > 0) // QUAD
  {
    gid = H5Gcreate(uc, "Quad4", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dims[0] = 1;
    did = H5Screate_simple(1, dims, NULL);
#ifdef E_DOUBLEINT
    aid = H5Acreate(gid, "NumberOfCells", H5T_NATIVE_INT64, did, H5P_DEFAULT, H5P_DEFAULT);
    int64_t numberOfCells = nquad;
    H5Awrite(aid, H5T_NATIVE_INT64, &numberOfCells);
#else
    aid = H5Acreate(gid, "NumberOfCells", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT);
    int32_t numberOfCells = nquad;
    H5Awrite(aid, H5T_NATIVE_INT, &numberOfCells);
#endif
    H5Aclose(aid); H5Sclose(did);
    dims[0] = nquad; dims[1] = 4;
    dquad = H5Screate_simple(2, dims, NULL);
#ifdef E_DOUBLEINT
    vquad = H5Dcreate(gid, "Cell2Node", H5T_NATIVE_INT64, dquad, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
    vquad = H5Dcreate(gid, "Cell2Node", H5T_NATIVE_INT, dquad, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
    hid_t gid2 = H5Gcreate(gid, "CellAttributes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dims[0] = nquad;
    did = H5Screate_simple(1, dims, NULL);
    mquad = H5Dcreate(gid2, "CADGroupID", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(gid); H5Gclose(gid2); H5Sclose(did);
  }

  std::vector<E_Int> oldStart(connects.size()); // old start of connect
  std::vector<E_Int> oldEnd(connects.size()); // old end of connect
  std::vector<E_Int> newStart(connects.size()); // new start of connect
  E_Int shexa=0, stetra=0, spenta=0, spyra=0, stri=0, squad=0;
  E_Int phexa=1, ptetra=nhexa+1, ppenta=nhexa+ntetra+1, ppyra=nhexa+ntetra+npenta+1;
  E_Int ptri=nhexa+ntetra+npenta+npyra+1, pquad=nhexa+ntetra+npenta+npyra+ntri+1;
  hsize_t start[2]; hsize_t scount[2];
  E_Int* ngon = NULL; E_Int* ngonOffset = NULL; E_Int npolyFaces = -1;
  E_Int* nface = NULL; E_Int* nfaceOffset = NULL; E_Int npolyCells = -1;

  for (size_t i = 0; i < connects.size(); i++)
  {
    E_Int* cv = K_PYTREE::getValueAI(connects[i], hook);
    PyObject* er = K_PYTREE::getNodeFromName1(connects[i], "ElementRange");
    E_Int* erv = K_PYTREE::getValueAI(er, hook);
    PyObject* ec = K_PYTREE::getNodeFromName1(connects[i], "ElementConnectivity");
    E_Int* ecv = K_PYTREE::getValueAI(ec, hook);
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
    else if (eltType == 22) // NGON
    { size = 0; }
    else if (eltType == 23) // NFACE
    { size = 0; }
    
    E_Int* ecv2 = new E_Int [size];
    for (E_Int j = 0; j < size; j++) ecv2[j] = ecv[j]-1;

    if (eltType == 17) // HEXA
    {
      //H5Dwrite(vhexa, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ecv2);
      start[0]=shexa; start[1]=0;
      scount[0]=size0; scount[1]=8; 
      hid_t sid = H5Scopy(dhexa);     
      H5Sselect_hyperslab(sid, H5S_SELECT_SET, start, NULL, scount, NULL);
      hid_t memspace = H5Screate_simple(2, scount, NULL);
#ifdef E_DOUBLEINT
      H5Dwrite(vhexa, H5T_NATIVE_INT64, memspace, sid, H5P_DEFAULT, ecv2);
#else
      H5Dwrite(vhexa, H5T_NATIVE_INT, memspace, sid, H5P_DEFAULT, ecv2);
#endif
      shexa += size0;
      oldStart[i] = erv[0]; oldEnd[i] = erv[1];
      newStart[i] = phexa;
      phexa += size0;
    }
    else if (eltType == 10) // TETRA
    {
      //H5Dwrite(vtetra, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ecv2);
      start[0]=stetra; start[1]=0;
      scount[0]=size0; scount[1]=4;
      hid_t sid = H5Scopy(dtetra);     
      H5Sselect_hyperslab(sid, H5S_SELECT_SET, start, NULL, scount, NULL);
      hid_t memspace = H5Screate_simple(2, scount, NULL);
#ifdef E_DOUBLEINT
      H5Dwrite(vtetra, H5T_NATIVE_INT64, memspace, sid, H5P_DEFAULT, ecv2);
#else
      H5Dwrite(vtetra, H5T_NATIVE_INT, memspace, sid, H5P_DEFAULT, ecv2);
#endif
      stetra += size0;
      oldStart[i] = erv[0]; oldEnd[i] = erv[1];
      newStart[i] = ptetra;
      ptetra += size0;
    }
    else if (eltType == 14) // PENTA
    {
      //H5Dwrite(vpenta, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ecv2);
      start[0]=spenta; start[1]=0;
      scount[0]=size0; scount[1]=6;
      hid_t sid = H5Scopy(dpenta);     
      H5Sselect_hyperslab(sid, H5S_SELECT_SET, start, NULL, scount, NULL);
      hid_t memspace = H5Screate_simple(2, scount, NULL);
#ifdef E_DOUBLEINT
      H5Dwrite(vpenta, H5T_NATIVE_INT64, memspace, sid, H5P_DEFAULT, ecv2);
#else
      H5Dwrite(vpenta, H5T_NATIVE_INT, memspace, sid, H5P_DEFAULT, ecv2);
#endif
      spenta += size0;
      oldStart[i] = erv[0]; oldEnd[i] = erv[1];
      newStart[i] = ppenta;
      ppenta += size0;
    }
    else if (eltType == 12) // PYRA
    {
      //H5Dwrite(vpyra, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ecv2);
      start[0]=spyra; start[1]=0;
      scount[0]=size0; scount[1]=5;
      hid_t sid = H5Scopy(dpyra);     
      H5Sselect_hyperslab(sid, H5S_SELECT_SET, start, NULL, scount, NULL);
      hid_t memspace = H5Screate_simple(2, scount, NULL);
#ifdef E_DOUBLEINT
      H5Dwrite(vpyra, H5T_NATIVE_INT64, memspace, sid, H5P_DEFAULT, ecv2);
#else
      H5Dwrite(vpyra, H5T_NATIVE_INT, memspace, sid, H5P_DEFAULT, ecv2);
#endif
      spyra += size0;
      oldStart[i] = erv[0]; oldEnd[i] = erv[1];
      newStart[i] = ppyra;
      ppyra += size0;
    }
    else if (eltType == 5) // TRI
    {
      //H5Dwrite(vtri, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ecv2);
      start[0]=stri; start[1]=0;
      scount[0]=size0; scount[1]=3;
      hid_t sid = H5Scopy(dtri);     
      H5Sselect_hyperslab(sid, H5S_SELECT_SET, start, NULL, scount, NULL);
      hid_t memspace = H5Screate_simple(2, scount, NULL);
#ifdef E_DOUBLEINT
      H5Dwrite(vtri, H5T_NATIVE_INT64, memspace, sid, H5P_DEFAULT, ecv2);
#else
      H5Dwrite(vtri, H5T_NATIVE_INT, memspace, sid, H5P_DEFAULT, ecv2);
#endif
      stri += size0;
      oldStart[i] = erv[0]; oldEnd[i] = erv[1];
      newStart[i] = ptri;
      ptri += size0;
    }
    else if (eltType == 7) // QUAD
    {
      //H5Dwrite(vquad, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ecv2);
      start[0]=squad; start[1]=0;
      scount[0]=size0; scount[1]=4;
      hid_t sid = H5Scopy(dquad);
      H5Sselect_none(sid);
      H5Sselect_hyperslab(sid, H5S_SELECT_SET, start, NULL, scount, NULL);
      hid_t memspace = H5Screate_simple(2, scount, NULL);
#ifdef E_DOUBLEINT
      H5Dwrite(vquad, H5T_NATIVE_INT64, memspace, sid, H5P_DEFAULT, ecv2);
#else
      H5Dwrite(vquad, H5T_NATIVE_INT, memspace, sid, H5P_DEFAULT, ecv2);
#endif
      squad += size0;
      oldStart[i] = erv[0]; oldEnd[i] = erv[1];
      newStart[i] = pquad;
      pquad += size0;
    }
    else if (eltType == 22) // NGON
    {
      ngon = ecv;
      npolyFaces = size0;
      PyObject* off = K_PYTREE::getNodeFromName1(connects[i], "ElementStartOffset");
      if (off != NULL)
      {
        E_Int* offv = K_PYTREE::getValueAI(off, hook);
        ngonOffset = offv;
      }
    }
    else if (eltType == 23)
    {
      npolyCells = size0;
      nface = ecv;
      PyObject* off = K_PYTREE::getNodeFromName1(connects[i], "ElementStartOffset"); // must be v4
      if (off != NULL)
      {
        E_Int* offv = K_PYTREE::getValueAI(off, hook);
        nfaceOffset = offv;
      }
    }
    delete [] ecv2;
  }

  E_Int possurf = nhexa+ntetra+npenta+npyra;
  if (nhexa > 0) { H5Dclose(vhexa); H5Sclose(dhexa); }
  if (ntetra > 0) { H5Dclose(vtetra); H5Sclose(dtetra); }
  if (npenta > 0) { H5Dclose(vpenta); H5Sclose(dpenta); }
  if (npyra > 0) { H5Dclose(vpyra); H5Sclose(dpyra); }
  if (ntri > 0) { H5Dclose(vtri); H5Sclose(dtri); }
  if (nquad > 0) { H5Dclose(vquad); H5Sclose(dquad); }

  // Boundary conditions
  PyObject* zbc = K_PYTREE::getNodeFromName1(zone, "ZoneBC");
  std::vector<PyObject*> BCs;
  if (zbc != NULL) K_PYTREE::getNodesFromType1(zbc, "BC_t", BCs);

  // Boundary conditions for ME
  if (ntri > 0 || nquad > 0)
  {
    int32_t* markersTri = new int32_t [ntri];
    int32_t* markersQuad = new int32_t [nquad];
    E_Int off = 0;
    for (E_Int i = 0; i < ntri; i++) markersTri[i] = 0;
    for (E_Int i = 0; i < nquad; i++) markersQuad[i] = 0;

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
          off = newIndex2(j, oldStart, oldEnd, newStart)-possurf;
          if (off < ntri) 
          {
            off = std::min(off, ntri-1);
            off = std::max(off, E_Int(0));
            markersTri[off] = count;
          }
          else 
          {
            off = off-ntri;
            off = std::min(off, nquad-1);
            off = std::max(off, E_Int(0));
            markersQuad[off] = count;
          }
          //if (off != j-1-nvol) 
          //printf("off=%d %d -> count=%d\n", off, j-1-nvol, count);
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
          off = newIndex2(plv[j], oldStart, oldEnd, newStart)-possurf;
          if (off < ntri) 
          { 
            off = std::min(off, ntri-1);
            off = std::max(off, E_Int(0));
            markersTri[off] = count;
          }
          else 
          {
            off = off-ntri;
            off = std::min(off, nquad-1);
            off = std::max(off, E_Int(0));
            markersQuad[off] = count;
          }
          //printf("off=%d %d\n", off, plv[j]-1-nvol);
        }
      }
      count += 1;
    }

    if (ntri > 0)
    {
      H5Dwrite(mtri, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, markersTri);
      H5Dclose(mtri);
      delete [] markersTri;
    }
    if (nquad) 
    {
      H5Dwrite(mquad, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, markersQuad);
      H5Dclose(mquad);
      delete [] markersQuad;
    }
  }

  // Create specific connectivities for NGON
  if (ngon != NULL && ngonOffset != NULL && nface != NULL && nfaceOffset != NULL)
  {
    // Build cell2Node connectivity
    E_Int* cell2NodeCount = new E_Int [npolyCells];
    E_Int* cell2FaceCount = new E_Int [npolyCells];
    E_Int* cell2NodeOffset = new E_Int [npolyCells+1]; // perso
    cell2NodeOffset[0] = 0;
    std::set<E_Int> nodes;

    E_Int size1 = 0; E_Int size0 = 0;
    E_Int c = 0;
    for (E_Int i = 0; i < npolyCells; i++)
    {
      nodes.clear();
      for (E_Int j = nfaceOffset[i]; j < nfaceOffset[i+1]; j++)
      {
        E_Int f = std::abs(nface[j])-1;
        for (E_Int n = ngonOffset[f]; n < ngonOffset[f+1]; n++)
        {
          E_Int node = ngon[n]-1;
          nodes.insert(node);
        }
      }
      size0 = nodes.size();
      cell2NodeCount[i] = size0;
      size1 += size0;
      cell2FaceCount[i] = nfaceOffset[i+1]-nfaceOffset[i];
      cell2NodeOffset[i+1] = cell2NodeOffset[i]+size0;
    }

    E_Int* cell2Node = new E_Int [size1];
    c = 0;
    for (E_Int i = 0; i < npolyCells; i++)
    {
      nodes.clear();
      for (E_Int j = nfaceOffset[i]; j < nfaceOffset[i+1]; j++)
      {
        E_Int f = std::abs(nface[j])-1;
        for (E_Int n = ngonOffset[f]; n < ngonOffset[f+1]; n++)
        {
          E_Int node = ngon[n]-1;
          nodes.insert(node);
        }
      }
      for (auto n : nodes)
      { 
        cell2Node[c] = n; c++;
      }
    }
    
    // Build face2node connectivity. No duplicated faces.
    /*
    int32_t* face2NodeCount = new int32_t [npolyFaces];
    E_Int size2 = ngonOffset[npolyFaces];
    int32_t* face2Node = new int32_t [size2];
    c = 0;
    E_Int loc = 0; E_Int elt = 0; // element containing face
    for (E_Int i = 0; i < npolyFaces; i++)
    {
      for (E_Int j = ngonOffset[i]; j < ngonOffset[i+1]; j++)
      {
        E_Int node = ngon[j]-1; // node index
        E_Int off = 0; // shift of node in cell
        // Cherche le noeud dans elt
        for (E_Int k = cell2NodeOffset[elt]; k < cell2NodeOffset[elt+1]; k++)
        {
          E_Int n = cell2Node[k];
          if (n == node) break;
          off++;
        }
        face2Node[c] = off; c++; // must be an offset on elt2Node
      }
      face2NodeCount[i] = ngonOffset[i+1] - ngonOffset[i];
      loc++;
      if (loc > cell2FaceCount[elt]) { elt++; loc = 0; }
    }
    */

    // Build face2node connectivity. Faces are duplicated by elements. They must be ordered.
    E_Int npolyFacesD = 0;
    for (E_Int i = 0; i < npolyCells; i++) npolyFacesD += cell2FaceCount[i];
    E_Int* face2NodeCount = new E_Int [npolyFacesD];

    E_Int size2 = 0; // size of face2Node
    for (E_Int i = 0; i < npolyCells; i++) // pour chaque cellule
    {
      for (E_Int f = nfaceOffset[i]; f < nfaceOffset[i+1]; f++) // pour toutes les faces
      {
        // face orig dans le pytree
        E_Int forig = std::abs(nface[f])-1;
        size2 += ngonOffset[forig+1] - ngonOffset[forig]; 
      }
    }
    E_Int* face2Node = new E_Int [size2];
    
    E_Int floc = 0; // numerotation par element des cellules fsdm
    c = 0;
    for (E_Int i = 0; i < npolyCells; i++) // pour chaque cellule
    {
      //E_Int nf = cell2FaceCount[i]; // nbre de faces de la cellule
      for (E_Int f = nfaceOffset[i]; f < nfaceOffset[i+1]; f++) // pour toutes les faces
      {
        // face orig dans le pytree
        E_Int forig = std::abs(nface[f])-1;

        if (nface[f] > 0)
        {
          for (E_Int j = ngonOffset[forig]; j < ngonOffset[forig+1]; j++) // pour tous les noeuds
          {
            E_Int node = ngon[j]-1; // node index
            E_Int off = 0; // shift of node in cell
            // Cherche le noeud dans elt
            for (E_Int k = cell2NodeOffset[i]; k < cell2NodeOffset[i+1]; k++)
            {
              E_Int n = cell2Node[k];
              if (n == node) break;
              off++;
            }
            face2Node[c] = off; c++; // must be an ordered offset on elt2Node 
          }
        }
        else
        {
          for (E_Int j = ngonOffset[forig+1]-1; j >= ngonOffset[forig]; j--) // pour tous les noeuds
          {
            E_Int node = ngon[j]-1; // node index
            E_Int off = 0; // shift of node in cell
            // Cherche le noeud dans elt
            for (E_Int k = cell2NodeOffset[i]; k < cell2NodeOffset[i+1]; k++)
            {
              E_Int n = cell2Node[k];
              if (n == node) break;
              off++;
            }
            face2Node[c] = off; c++; // must be an ordered offset on elt2Node 
          }
        }
        face2NodeCount[floc] = ngonOffset[forig+1] - ngonOffset[forig]; floc++;
      }
    }

    gid = H5Gcreate(uc, "Poly3D", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dims[0] = 1;
    did = H5Screate_simple(1, dims, NULL);
#ifdef E_DOUBLEINT
    aid = H5Acreate(gid, "NumberOfCells", H5T_NATIVE_INT64, did, H5P_DEFAULT, H5P_DEFAULT);
    int64_t numberOfCells = npolyCells;
    H5Awrite(aid, H5T_NATIVE_INT64, &numberOfCells);
#else
    aid = H5Acreate(gid, "NumberOfCells", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT);
    int32_t numberOfCells = npolyCells;
    H5Awrite(aid, H5T_NATIVE_INT, &numberOfCells);
#endif
    H5Aclose(aid); H5Sclose(did);
    dims[0] = npolyCells;
    did = H5Screate_simple(1, dims, NULL);
#ifdef E_DOUBLEINT
    vid = H5Dcreate(gid, "Cell2NodeCounts", H5T_NATIVE_INT64, did, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(vid, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell2NodeCount);
#else
    vid = H5Dcreate(gid, "Cell2NodeCounts", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(vid, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell2NodeCount);
#endif
    H5Dclose(vid); H5Sclose(did);
    dims[0] = size1;
    did = H5Screate_simple(1, dims, NULL);
#ifdef E_DOUBLEINT
    vid = H5Dcreate(gid, "Cell2NodeList", H5T_NATIVE_INT64, did, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(vid, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell2Node);
#else
    vid = H5Dcreate(gid, "Cell2NodeList", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(vid, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell2Node);
#endif
    H5Dclose(vid); H5Sclose(did);
    dims[0] = npolyCells;
    did = H5Screate_simple(1, dims, NULL);
#ifdef E_DOUBLEINT
    vid = H5Dcreate(gid, "Cell2FaceCounts", H5T_NATIVE_INT64, did, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(vid, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell2FaceCount);
#else
    vid = H5Dcreate(gid, "Cell2FaceCounts", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(vid, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell2FaceCount);
#endif
    H5Dclose(vid); H5Sclose(did);
    dims[0] = npolyFacesD;
    did = H5Screate_simple(1, dims, NULL);
#ifdef E_DOUBLEINT
    vid = H5Dcreate(gid, "Face2NodeCounts", H5T_NATIVE_INT64, did, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(vid, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, face2NodeCount);
#else
    vid = H5Dcreate(gid, "Face2NodeCounts", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(vid, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, face2NodeCount);
#endif
    H5Dclose(vid); H5Sclose(did);
    dims[0] = size2;
    did = H5Screate_simple(1, dims, NULL);
#ifdef E_DOUBLEINT
    vid = H5Dcreate(gid, "Face2NodeList", H5T_NATIVE_INT64, did, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(vid, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, face2Node);
#else
    vid = H5Dcreate(gid, "Face2NodeList", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(vid, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, face2Node);
#endif
    H5Dclose(vid); H5Sclose(did);
    H5Gclose(gid);

    delete [] cell2NodeCount;
    delete [] cell2FaceCount;
    delete [] cell2NodeOffset;
    delete [] cell2Node;
    delete [] face2NodeCount;
    delete [] face2Node;
  }

  // Boundary conditions for NGON
  if (ngon != NULL && ngonOffset != NULL && nface != NULL && nfaceOffset != NULL)
  {
    // dimensionnement Poly2d
    E_Int size = 0; E_Int nc = 0;
    for (size_t i = 0; i < BCs.size(); i++)
    {
      PyObject* pl = K_PYTREE::getNodeFromName1(BCs[i], "PointList"); // face centers
      if (pl != NULL)
      {
        E_Int s0, s1;
        E_Int* p = K_PYTREE::getValueAI(pl, s0, s1, hook);
        for (E_Int f = 0; f < s0*s1; f++)
        {
          size += ngonOffset[p[f]]-ngonOffset[p[f]-1];
        }
        nc += s0*s1;
      }
      PyObject* pr = K_PYTREE::getNodeFromName1(BCs[i], "PointRange"); // face range
      if (pr != NULL)
      {
        E_Int* p = K_PYTREE::getValueAI(pr, hook);
        for (E_Int f = p[0]; f < p[1]; f++) size += ngonOffset[f]-ngonOffset[f-1];
        nc += p[1]-p[0]+1;
      }
    }
    
    int32_t* poly2dCount = new int32_t [nc];
    int32_t* poly2d = new int32_t [size];
    int32_t* marker = new int32_t [nc];

    // Remplissage poly2d
    E_Int bcCount = 1; E_Int faceCount = 0; E_Int ptr = 0;
    for (size_t i = 0; i < BCs.size(); i++)
    {
      PyObject* pl = K_PYTREE::getNodeFromName1(BCs[i], "PointList"); // face centers
      if (pl != NULL)
      {
        E_Int s0, s1;
        E_Int* p = K_PYTREE::getValueAI(pl, s0, s1, hook);
        for (E_Int f = 0; f < s0*s1; f++) // pour chaque face
        {
          E_Int size0 = ngonOffset[p[f]]-ngonOffset[p[f]-1];
          poly2dCount[faceCount] = size0;
          marker[faceCount] = bcCount;
          for (E_Int j = ngonOffset[p[f]-1]; j < ngonOffset[p[f]]; j++) { poly2d[ptr] = ngon[j]-1; ptr++; }
          faceCount++;
        }
      }
      PyObject* pr = K_PYTREE::getNodeFromName1(BCs[i], "PointRange"); // face range
      if (pr != NULL)
      {
        E_Int* p = K_PYTREE::getValueAI(pr, hook);
        for (E_Int f = p[0]; f < p[1]; f++) // pour chaque face 
        {
          E_Int size0 = ngonOffset[f+1]-ngonOffset[f];
          poly2dCount[faceCount] = size0;
          marker[faceCount] = bcCount;
          for (E_Int j = ngonOffset[f-1]; j < ngonOffset[f]; j++) { poly2d[ptr] = ngon[j]-1; ptr++; }
          faceCount++;
        }
      }
      bcCount += 1;
    }

    if (nc > 0)
    {
      gid = H5Gcreate(uc, "Poly2D", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      dims[0] = 1;
      did = H5Screate_simple(1, dims, NULL);
#ifdef E_DOUBLEINT
      aid = H5Acreate(gid, "NumberOfCells", H5T_NATIVE_INT64, did, H5P_DEFAULT, H5P_DEFAULT);
      int64_t numberOfCells = nc;
      H5Awrite(aid, H5T_NATIVE_INT64, &numberOfCells);
#else
      aid = H5Acreate(gid, "NumberOfCells", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT);
      int32_t numberOfCells = nc;
      H5Awrite(aid, H5T_NATIVE_INT, &numberOfCells);
#endif
      H5Aclose(aid); H5Sclose(did);
      dims[0] = nc;
      did = H5Screate_simple(1, dims, NULL);
#ifdef E_DOUBLEINT
      vid = H5Dcreate(gid, "Cell2NodeCounts", H5T_NATIVE_INT64, did, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(vid, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, poly2dCount);
#else
      vid = H5Dcreate(gid, "Cell2NodeCounts", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(vid, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, poly2dCount);
#endif
      H5Dclose(vid); H5Sclose(did);
      dims[0] = size;
      did = H5Screate_simple(1, dims, NULL);
#ifdef E_DOUBLEINT
      vid = H5Dcreate(gid, "Cell2NodeList", H5T_NATIVE_INT64, did, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(vid, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, poly2d);
#else
      vid = H5Dcreate(gid, "Cell2NodeList", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(vid, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, poly2d);
#endif
      H5Dclose(vid); H5Sclose(did);
      gid2 = H5Gcreate(gid, "CellAttributes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      dims[0] = nc;
      did = H5Screate_simple(1, dims, NULL);
      vid = H5Dcreate(gid2, "CADGroupID", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(vid, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, marker);
      H5Dclose(vid); H5Sclose(did);
      H5Gclose(gid2);
      H5Gclose(gid);
    }
    delete [] poly2dCount;
    delete [] poly2d;
    delete [] marker;
  }

  // Mapping cellType to index
  /*
  gid = H5Gcreate(uc, "MappingCellType2Index", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dims[0] = 1;
  did = H5Screate_simple(1, dims, NULL);
  aid = H5Acreate(gid, "Node", H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT);
  int32_t value = 0;
  H5Awrite(aid, H5T_NATIVE_INT, &value);
  H5Aclose(aid); H5Sclose(did);

  for (size_t i = 0; i < connects.size(); i++)
  {
    E_Int* cv = K_PYTREE::getValueAI(connects[i], hook);
    E_Int eltType = cv[0];
    if (eltType == 17) // HEXA
    { strcpy(name, "Hexa8"); }
    else if (eltType == 10) // TETRA
    { strcpy(name, "Tetra4"); }
    else if (eltType == 14) // PENTA
    { strcpy(name, "Prism6"); }
    else if (eltType == 12) // PYRA
    { strcpy(name, "Pyra5"); }
    else if (eltType == 5) // TRI
    { strcpy(name, "Tri3"); }
    else if (eltType == 7) // QUADS
    { strcpy(name, "Quad7"); }
    else if (eltType == 22) // NGON
    { strcpy(name, "Poly3D"); }
    else if (eltType == 23) // NFACE
    { strcpy(name, "Poly2D"); }
    else strcpy(name, "Unknown");

    dims[0] = 1;
    did = H5Screate_simple(1, dims, NULL);
    aid = H5Acreate(gid, name, H5T_NATIVE_INT, did, H5P_DEFAULT, H5P_DEFAULT);
    int32_t value = i;
    H5Awrite(aid, H5T_NATIVE_INT, &value);
    H5Aclose(aid); H5Sclose(did);
  }
  H5Gclose(gid);
  */
 
  // Name of cell attribute values
  gid = H5Gcreate(uc, "NamesOfCellAttributeValues", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  gid2 = H5Gcreate(gid, "CADGroupID", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  for (size_t i = 0; i < BCs.size(); i++)
  {
    char* bcname = K_PYTREE::getNodeName(BCs[i], hook);
    did = H5Screate(H5S_SCALAR); 
    tid = H5Tcopy(H5T_C_S1); H5Tset_size(tid, strlen(bcname)+1);
    sprintf(name, "%zu", i+1);
    aid = H5Acreate(gid2, name, tid, did, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(aid, tid, bcname);
    H5Aclose(aid); H5Sclose(did); H5Tclose(tid);
  }
  H5Gclose(gid2); H5Gclose(gid);
  delete [] name;

  // close
  H5Gclose(uc);
  H5Fclose(fid);

  // Decref
  for (size_t i = 0; i < hook.size(); i++) Py_DECREF(hook[i]);

  return 0;
}
