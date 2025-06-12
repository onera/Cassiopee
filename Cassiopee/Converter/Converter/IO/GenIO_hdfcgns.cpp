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

// Binary HDF CGNS file support
#ifdef _MPI
#if defined(_WIN64)
# define __int64 long long
#endif
# include "mpi.h"
# include "mpi4py/mpi4py.h"
#else
#define MPI_Comm void
#endif

// to be suppressed in next release
#define FORCEPERIODICR4 0

# include "GenIO.h"
# include "kcore.h"
# include "hdf5.h"
# include "GenIO_hdfcgns.h"
using namespace std;
using namespace K_IO;

/* ------------------------------------------------------------------------- */
E_Int checkCompressionFilters()
{
  unsigned int filter_info;
  E_Int avail = H5Zfilter_avail(H5Z_FILTER_SZIP);
  if (!avail) 
  {
    printf ("Warning: hdf: szip filter not available.\n");
    return 1;
  }
  H5Zget_filter_info (H5Z_FILTER_SZIP, &filter_info);
  if ( !(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED) ||
      !(filter_info & H5Z_FILTER_CONFIG_DECODE_ENABLED) ) 
  {
    printf ("Warning: hdf: szip filter not available for encoding and decoding.\n");
    return 1;
  }

  // Example of creating compressed dataSet
  /*
  dcpl = H5Pcreate(H5P_DATASET_CREATE);
  status = H5Pset_szip (dcpl, H5_SZIP_NN_OPTION_MASK, 8);
  status = H5Pset_chunk(dcpl, 2, chunk);

  dset = H5Dcreate (file, DATASET, H5T_STD_I32LE, space, H5P_DEFAULT, dcpl,
                    H5P_DEFAULT);
  status = H5Dwrite (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                wdata[0]);
  */
  return 1;
}

/* ------------------------------------------------------------------------- */
hid_t K_IO::GenIOHdf::ADF_to_HDF_datatype(const char *tp)
{
  if (!strcmp(tp, L3T_B1)) return H5T_NATIVE_UCHAR;
  else if (!strcmp(tp, L3T_C1)) return H5T_NATIVE_CHAR;
  else if (!strcmp(tp, L3T_I4)) return _NATIVE_INT;
  else if (!strcmp(tp, L3T_R8)) return _NATIVE_DOUBLE;
  else if (!strcmp(tp, L3T_I1)) return H5T_NATIVE_INT8;
  else if (!strcmp(tp, L3T_I8)) return _NATIVE_LONG;
  else if (!strcmp(tp, L3T_U4)) return H5T_NATIVE_UINT32;
  else if (!strcmp(tp, L3T_U8)) return H5T_NATIVE_UINT64;
  else if (!strcmp(tp, L3T_R4)) return _NATIVE_FLOAT;
  else return 0;
}
/* ------------------------------------------------------------------------- */
int HDF_Get_Attribute_As_Integer(hid_t nodeid, const char *name, int *value)
{
  hid_t aid; herr_t status;
  aid = H5Aopen_by_name(nodeid, ".", name, H5P_DEFAULT, H5P_DEFAULT);

  if (aid < 0) return 0;
  status = H5Aread(aid, H5T_NATIVE_INT, value);
  H5Aclose(aid);
  if (status < 0) return 0;
  return *value;
}

/* ------------------------------------------------------------------------- */
char* HDF_Get_Attribute_As_String(hid_t nodeid, const char* name, char* value)
{
  hid_t aid, tid;
  value[0] = '\0';
  aid = H5Aopen_by_name(nodeid, ".", name, H5P_DEFAULT, H5P_DEFAULT);

  if (aid > 0)
  {
    tid = H5Aget_type(aid);
    if (tid < 0)
    {
      H5Aclose(aid);
      return 0;
    }
    H5Aread(aid, tid, value);
    H5Tclose(tid); H5Aclose(aid);
  }
  return value;
}
/* ------------------------------------------------------------------------- */
int HDF_Set_Attribute_As_String(hid_t nodeid, const char *name, char *value)
{
  hid_t aid, tid;
  aid = H5Aopen_by_name(nodeid, ".", name, H5P_DEFAULT, H5P_DEFAULT);

  if (aid > 0)
  {
    tid = H5Aget_type(aid);
    if (tid < 0)
    {
      H5Aclose(aid);
      return 0;
    }
    H5Awrite(aid, tid, value);
    H5Tclose(tid); H5Aclose(aid);
  }
  return 1; // OK
}

/* ------------------------------------------------------------------------- */
char *HDF_Get_Attribute_As_Data(hid_t nodeid, const char *name, char *value)
{
  hid_t did;
  value[0] = '\0';
  did = H5Dopen2(nodeid, name, H5P_DEFAULT);
  if (did > 0)
  {
    
#if defined(_MPI) && defined(H5_HAVE_PARALLEL)
    // hid_t xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    // hid_t ret        = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_INDEPENDENT);
    // H5Dread(did, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, xfer_plist, value);
    H5Dread(did, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, value);
#else
    H5Dread(did, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, value);
#endif
    H5Dclose(did);
  }
  return value;
}

/* ------------------------------------------------------------------------- */
#if H5_VERSION_LE(1,11,9)
static herr_t count_children(hid_t id, const char* name, void* count)
#elif H5_VERSION_LE(1,13,9)
static herr_t count_children(hid_t id, const char* name, const H5L_info_t* linfo, void* count)
#else
static herr_t count_children(hid_t id, const char* name, const H5L_info2_t* linfo, void* count)
#endif
{
  if (name && (name[0] != ' '))
  {
    (*((int*)count))++;
  }
  return 0;
}
/* ------------------------------------------------------------------------- */
#if H5_VERSION_LE(1,11,9)
static herr_t gfind_name(hid_t id, const char *nm, void *snm)
#elif H5_VERSION_LE(1,13,9)
static herr_t gfind_name(hid_t id, const char *nm, const H5L_info_t* linfo, void *snm)
#else
static herr_t gfind_name(hid_t id, const char *nm, const H5L_info2_t* linfo, void *snm)
#endif
{
  /*  printf("GFIND [%s][%s]\n",nm,snm);fflush(stdout); */
  if (!strcmp(nm, (char*)snm)) return 1;
  return 0;
}
#if H5_VERSION_LE(1,11,9)
#define has_child(ID,NAME) H5Giterate(ID, ".", NULL, gfind_name, (void*)NAME)
#else
#define has_child(ID,NAME) H5Lexists(ID, NAME, H5P_DEFAULT)
//#define has_child(ID,NAME) H5Literate2(ID, H5_INDEX_CRT_ORDER, H5_ITER_NATIVE, NULL, gfind_name, (void*)NAME)
#endif
/* ------------------------------------------------------------------------- */
#if H5_VERSION_LE(1,11,9)
static herr_t delete_children(hid_t id, const char* name, void* data)
#elif H5_VERSION_LE(1,13,9)
static herr_t delete_children(hid_t id, const char* name, const H5L_info_t* linfo, void* data)
#else
static herr_t delete_children(hid_t id, const char* name, const H5L_info2_t* linfo, void* data)
#endif
{
  /* do not change link id with actual here, stop deletion at link node */
  if (name && (name[0] == ' ')) /* leaf node */
  {
    H5Ldelete(id, name, H5P_DEFAULT);
  }
  else 
  {
    /* should use H5Literate */
#if H5_VERSION_LE(1,11,9)
    H5Giterate(id, name, NULL, delete_children, data);
#else
    H5Literate_by_name2(id, name, H5_INDEX_CRT_ORDER, H5_ITER_INC, NULL, delete_children, data, H5P_DEFAULT);
#endif
    H5Ldelete(id, name, H5P_DEFAULT);
  }
  return 0;
}

/* ------------------------------------------------------------------------- */
#if H5_VERSION_LE(1,13,9)
static herr_t feed_children_ids(hid_t id, const char* name,
                                const H5L_info_t* linfo, void* idlist)
#else
static herr_t feed_children_ids(hid_t id, const char* name,
                                const H5L_info2_t* linfo, void* idlist)
#endif
{
  hid_t cid; int n;

  /* skip names starting with a <space> */
  if (name && (name[0] == ' ')) return 0;
  cid = H5Gopen(id, name, H5P_DEFAULT); // group is closed in loadOne
  /* set id  */
  n = 0;
  while (((hid_t*)idlist)[n] != -1) n++;
  ((hid_t*)idlist)[n] = cid;
  return 0;
}

/* ------------------------------------------------------------------------- */
// Retourne 1 si path est dans links
E_Int checkPathInLinks(const char* path, PyObject* links)
{
  E_Int li = PyList_Size(links);
  for (E_Int i = 0; i < li; i++)
  {
    PyObject* tp = PyList_GetItem(links, i);
    PyObject* d = PyList_GetItem(tp, 3);
    char* s;
    if (PyString_Check(d)) s = PyString_AsString(d);
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(d)) s = (char*)PyUnicode_AsUTF8(d);
#endif
    else s = NULL;
    //printf("path=%s links=%s\n", path, s);
    if (strcmp(s, path) == 0) return 1;
  }
  return 0;
}

/* ------------------------------------------------------------------------- */
hid_t* K_IO::GenIOHdf::getChildren(hid_t nodeid)
{
  hid_t* idlist;
  hid_t gpl;
  int nchildren, n;
  unsigned order=0;

  nchildren = 0;
#if H5_VERSION_LE(1,11,9)
  H5Giterate(nodeid, ".", NULL, count_children, (void*)&nchildren);
#else
  //H5Literate2(nodeid, H5_INDEX_CRT_ORDER, H5_ITER_INC, NULL, count_children, (void*)&nchildren);
  H5Literate2(nodeid, H5_INDEX_NAME, H5_ITER_INC, NULL, count_children, (void*)&nchildren);
#endif
  if (!nchildren) return NULL;

  idlist = (hid_t*)malloc(sizeof(hid_t)*(nchildren+1));
  /* use last -1 as sentinel */
  for (n = 0; n <= nchildren; n++) {idlist[n] = (hid_t)-1;}

  // le fichier a-t'il ete enregistre avec order?
  if (H5Iget_type(nodeid) == H5I_GROUP)
  {
    gpl = H5Gget_create_plist(nodeid);
    if (H5Iis_valid(gpl)) H5Pget_link_creation_order(gpl, &order);
    H5Pclose(gpl);
  }
  else order = 0;

  if ((order & H5_INDEX_CRT_ORDER) == H5_INDEX_CRT_ORDER)
  {
    // avec order
#if H5_VERSION_LE(1,11,9)
    H5Literate(nodeid, H5_INDEX_CRT_ORDER, H5_ITER_INC,
               NULL, feed_children_ids, (void*)idlist);
#else
    H5Literate2(nodeid, H5_INDEX_CRT_ORDER, H5_ITER_INC,
               NULL, feed_children_ids, (void*)idlist);
#endif
  }
  else
  {
    // sans order
#if H5_VERSION_LE(1,11,9)
    H5Literate(nodeid, H5_INDEX_NAME, H5_ITER_INC,
               NULL, feed_children_ids, (void*)idlist);
#else
    H5Literate2(nodeid, H5_INDEX_NAME, H5_ITER_INC,
                NULL, feed_children_ids, (void*)idlist);
#endif
  }

  // Cette fonction ne marche que si les noeuds ont ete cree avec
  // la property creation_order
  //int status;
  //status = H5Literate(nodeid, H5_INDEX_CRT_ORDER, H5_ITER_INC,
  //                    NULL, feed_children_ids, (void *)idlist);
  // Celle-ci utilise l'ordre natif
  //status = H5Literate(nodeid, H5_INDEX_NAME, H5_ITER_NATIVE,
  //                    NULL, feed_children_ids, (void *)idlist);
  // Celle-ci lit les enfants dans l'ordre alphabetique (comme cgnsview)
  //H5Giterate(nodeid, ".",
  //           NULL, feed_children_ids_list, (void *)idlist);

  // Remplace les enfants dont les noeuds sont des liens
  /*
  H5G_stat_t sb;
  char _dtype[CGNSMAXLABEL+1];
  for (n = 0; n < nchildren; n++)
  {
    hid_t id = idlist[n];
    HDF_Get_Attribute_As_String(id, L3S_DTYPE, dtype);
    if (strcmp(dtype, "LK") == 0)
    {
      herr_t herr = H5Gget_objinfo(id, L3S_LINK, (hbool_t)0, &sb);
      if (herr < 0)
      {
        printf("Error: hdfcgnsread: error opening link file.\n");
        // supprimer le noeud?
      }
      else
      {
        hid_t lid = H5Gopen2(id, L3S_LINK, H5P_DEFAULT);
        idlist[n] = lid; // du coup le noeud de lien apparait avec le nouveau nom du noeud
        H5Gclose(id); //H5Gclose(herr);
      }
    }
  }
  */

  return idlist;
}

//=============================================================================
int HDF_Get_DataDimensions(hid_t nid, hsize_t *dims)
{
  int n, ndims;
  hsize_t int_dim_vals[CGNSMAXDIM];
  hid_t did, sid;

  L3M_CLEARDIMS(dims);
  did = H5Dopen2(nid, L3S_DATA, H5P_DEFAULT);
  sid = H5Dget_space(did);
  ndims = H5Sget_simple_extent_ndims(sid);
  H5Sget_simple_extent_dims(sid, int_dim_vals, NULL);

  for (n = 0; n < ndims; n++){ dims[n] = int_dim_vals[n]; }

  H5Dclose(did); H5Sclose(sid);

  return 1;
}

/* ------------------------------------------------------------------------- */
int HDF_Add_Attribute_As_String(hid_t nodeid, const char *name,
                                const char *value)
{
  hid_t sid, tid, aid;
  herr_t status;
  hsize_t dim;
  char buff[L3C_MAX_ATTRIB_SIZE+1];

  if (!strcmp(name, L3S_DTYPE)) dim = (hsize_t)(L3C_MAX_DTYPE+1);
  else dim = (hsize_t)(L3C_MAX_ATTRIB_SIZE+1);

  sid = H5Screate(H5S_SCALAR);
  if (sid < 0) return 0;

  tid = H5Tcopy(H5T_C_S1);
  if (tid < 0)
  {
    H5Sclose(sid); return 0;
  }
  if (H5Tset_size(tid, dim) < 0)
  {
    H5Tclose(tid); H5Sclose(sid);
    return 0;
  }
  aid = H5Acreate(nodeid, name, tid, sid, H5P_DEFAULT, H5P_DEFAULT);
  if (aid < 0)
  {
    H5Tclose(tid); H5Sclose(sid); return 0;
  }
  if (strlen(value) > L3C_MAX_ATTRIB_SIZE)
    printf("Warning: hdfcgnswrite: %s node name has been truncated.\n", value);
  memset(buff, 0, dim); strcpy(buff, value);
  status = H5Awrite(aid, tid, buff);

  H5Aclose(aid); H5Tclose(tid); H5Sclose(sid);

  /* autre facon -
  sid = H5Screate_simple(1, &dim, NULL);
  aid = H5Acreate(nodeid, name, H5T_NATIVE_CHAR, sid, H5P_DEFAULT, H5P_DEFAULT);
  memset(buff, 0, dim); strcpy(buff, value);
  status = H5Awrite(aid, H5T_NATIVE_CHAR, buff);
  */

  if (status < 0) return 0;
  return 1;
}
/* ------------------------------------------------------------------------- */
int HDF_Add_Attribute_As_Integer(hid_t nodeid, const char *name, int value)
{
  hid_t sid, aid;
  herr_t status;
  hsize_t dim;

  dim = 1;
  sid = H5Screate_simple(1, &dim, NULL);
  if (sid < 0) return 0;

  aid = H5Acreate(nodeid, name, H5T_NATIVE_INT, sid, H5P_DEFAULT, H5P_DEFAULT);
  if (aid < 0)
  {
    H5Sclose(sid); return 0;
  }

  status = H5Awrite(aid, H5T_NATIVE_INT, &value);
  H5Aclose(aid); H5Sclose(sid);

  if (status < 0) return 0;
  return 1;
}

/* ------------------------------------------------------------------------- */
int HDF_Add_Attribute_As_Data(hid_t id, const char *name,
                              const char *value,int size)
{
  hid_t   sid,did;
  hsize_t dim;
  herr_t  status;

  dim = (hsize_t)(size+1);
  sid = H5Screate_simple(1,&dim,NULL);
  if (sid < 0)
  {
    return 0;
  }
  did = H5Dcreate2(id, name, H5T_NATIVE_CHAR, sid,
                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (did < 0)
  {
    H5Sclose(sid);
    return 0;
  }
  status = H5Dwrite(did,H5T_NATIVE_CHAR,H5S_ALL,H5S_ALL,H5P_DEFAULT,value);
  H5Dclose(did);
  H5Sclose(sid);
  
  // char buff[L3C_MAX_ATTRIB_SIZE+1];
  // memset(buff, 0, dim); strcpy(buff, value);
  // status = H5Awrite(aid, tid, buff);

  if (status < 0) return 0;
  return 1;
}
//=============================================================================
int K_IO::GenIOHdf::getSingleI4(hid_t node, hid_t tid)
{
  int r;
  hid_t did, yid;
  did = H5Dopen2(node, L3S_DATA, H5P_DEFAULT);
  yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  H5Dread(did, yid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &r);
  H5Dclose(did); H5Tclose(yid);
  return r;
}

//=============================================================================
E_LONG K_IO::GenIOHdf::getSingleI8(hid_t node, hid_t tid)
{
  E_LONG r;
  hid_t did, yid;
  did = H5Dopen2(node, L3S_DATA, H5P_DEFAULT);
  yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  H5Dread(did, yid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &r);
  H5Dclose(did); H5Tclose(yid);
  return r;
}
//=============================================================================
float K_IO::GenIOHdf::getSingleR4(hid_t node, hid_t tid)
{
  float r;
  hid_t did,yid;
  did = H5Dopen2(node, L3S_DATA, H5P_DEFAULT);
  yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  H5Dread(did, yid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &r);
  H5Dclose(did); H5Tclose(yid);
  return r;
}

//=============================================================================
double K_IO::GenIOHdf::getSingleR8(hid_t node, hid_t tid)
{
  double r;
  hid_t did, yid;
  did = H5Dopen2(node, L3S_DATA, H5P_DEFAULT);
  yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  H5Dread(did, yid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &r);
  H5Dclose(did); H5Tclose(yid);
  return r;
}

//=============================================================================
PyObject* K_IO::GenIOHdf::getArrayI1(hid_t node, hid_t tid,
                                     int dim, hsize_t* dims,
                                     hid_t mid,            /* mem_space_id */
                                     hid_t sid)
{
  IMPORTNUMPY;
  PyArrayObject* r = NULL;

  // Create numpy: en INT8
  vector<npy_intp> npy_dim_vals(dim);
  for (E_Int nn = 0; nn < dim; nn++) npy_dim_vals[nn] = dims[nn];
  r = (PyArrayObject*)PyArray_EMPTY(dim, &npy_dim_vals[0], NPY_BYTE, 1);

  hid_t did, yid;
  did = H5Dopen2(node, L3S_DATA, H5P_DEFAULT);
  yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  
#if defined(_MPI) && defined(H5_HAVE_PARALLEL)
  if (_ismpi == 1)    /** HDF is executed in parallel context and compiled in MPI **/
  {
    hid_t xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    hid_t ret        = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    H5Dread(did, yid, mid, sid, xfer_plist, PyArray_DATA(r));
  }
  else /** HDF is executed in sequential context and compiled in MPI **/
  {
    H5Dread(did, yid, mid, sid, H5P_DEFAULT, PyArray_DATA(r));
  }
#else  /** HDF is executed in sequential context and compiled in sequential **/
  H5Dread(did, yid, mid, sid, H5P_DEFAULT, PyArray_DATA(r));
#endif
  H5Tclose(yid); H5Dclose(did);

  return (PyObject*)r;
}

//=============================================================================
PyObject* K_IO::GenIOHdf::getArrayI4Raw(hid_t node, hid_t tid,
                                        int dim, hsize_t* dims,
                                        hid_t mid,
                                        hid_t sid)
{
  IMPORTNUMPY;
  PyArrayObject* r = NULL;
  // Create numpy en NPY_INT
  vector<npy_intp> npy_dim_vals(dim);
  for (E_Int nn = 0; nn < dim; nn++) npy_dim_vals[nn] = dims[nn];
  r = (PyArrayObject*)PyArray_EMPTY(dim, &npy_dim_vals[0], NPY_INT, 1);

  hid_t did, yid;
  did = H5Dopen2(node, L3S_DATA, H5P_DEFAULT);
  yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  
#if defined(_MPI) && defined(H5_HAVE_PARALLEL)
  if (_ismpi == 1)    /** HDF is executed in parallel context and compiled in MPI **/
  {
    hid_t xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    hid_t ret        = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    H5Dread(did, yid, mid, sid, xfer_plist, PyArray_DATA(r));
  }
  else /** HDF is executed in sequential context and compiled in MPI **/
  {
    H5Dread(did, yid, mid, sid, H5P_DEFAULT, PyArray_DATA(r));
  }
#else  /** HDF is executed in sequential context and compiled in sequential **/
  H5Dread(did, yid, mid, sid, H5P_DEFAULT, PyArray_DATA(r));
#endif
  H5Tclose(yid); H5Dclose(did);

  return (PyObject*)r;
}

//=============================================================================
// Get arrayi4, le convertit en i8
PyObject* K_IO::GenIOHdf::getArrayI42I8(hid_t node, hid_t tid,
                                        int dim, hsize_t* dims,
                                        hid_t mid,
                                        hid_t sid)
{
  IMPORTNUMPY;
  int s, sizem;
  PyArrayObject* r = NULL;
  sizem = 1;
  for (s = 0; s < dim; s++) sizem = sizem*dims[s];
  int* ptr = (int*)::malloc(sizem*sizeof(int));
  if (!ptr) { Py_INCREF(Py_None); return Py_None; }

  hid_t did,yid;
  did = H5Dopen2(node, L3S_DATA, H5P_DEFAULT);
  yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  
#if defined(_MPI) && defined(H5_HAVE_PARALLEL)
  if (_ismpi == 1)    // HDF is executed in parallel context and compiled in MPI
  {
    hid_t xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    hid_t ret        = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    H5Dread(did, yid, mid, sid, xfer_plist, ptr);
  }
  else // HDF is executed in sequential context and compiled in MPI
  {
    H5Dread(did, yid, mid, sid, H5P_DEFAULT, ptr);
  }
#else // HDF is executed in sequential context and compiled in sequential
  H5Dread(did, yid, mid, sid, H5P_DEFAULT, ptr);
#endif
  
  H5Tclose(yid); H5Dclose(did);

  int64_t* ptr2 = (int64_t*)::malloc(sizem*sizeof(int64_t));
  for (int n = 0; n < sizem; n++)
  {
    ptr2[n] = (int64_t)(ptr[n]);
  }
  free(ptr);

  // Create numpy
  vector<npy_intp> npy_dim_vals(dim);
  for (E_Int nn = 0; nn < dim; nn++)
  {
    npy_dim_vals[nn] = dims[nn];
  }
  r = (PyArrayObject*)PyArray_EMPTY(dim, &npy_dim_vals[0], NPY_INT64, 1);
  memcpy(PyArray_DATA(r), ptr2, sizem*sizeof(int64_t));
  free(ptr2);

  return (PyObject*)r;
}

//=============================================================================
PyObject* K_IO::GenIOHdf::getArrayI8Raw(hid_t node, hid_t tid,
                                       int dim, hsize_t* dims,
                                       hid_t mid,  /* mem_space_id */
                                       hid_t sid)
{
  IMPORTNUMPY;
  PyArrayObject* r = NULL;
  vector<npy_intp> npy_dim_vals(dim);
  for (E_Int nn = 0; nn < dim; nn++) npy_dim_vals[nn] = dims[nn];
  r = (PyArrayObject*)PyArray_EMPTY(dim, &npy_dim_vals[0], NPY_INT64, 1);

  hid_t did, yid;
  did = H5Dopen2(node, L3S_DATA, H5P_DEFAULT);
  yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  
#if defined(_MPI) && defined(H5_HAVE_PARALLEL)
  if (_ismpi == 1)    /** HDF is executed in parallel context and compiled in MPI **/
  {
    hid_t xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    hid_t ret        = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    H5Dread(did, yid, mid, sid, xfer_plist, PyArray_DATA(r));
  }
  else /** HDF is executed in sequential context and compiled in MPI **/
  {
    H5Dread(did, yid, mid, sid, H5P_DEFAULT, PyArray_DATA(r));
  }
#else  /** HDF is executed in sequential context and compiled in sequential **/
  H5Dread(did, yid, mid, sid, H5P_DEFAULT, PyArray_DATA(r));
#endif
  H5Tclose(yid); H5Dclose(did);

  return (PyObject*)r;
}

//=============================================================================
// Get arrayi8, le convertit en i4
PyObject* K_IO::GenIOHdf::getArrayI82I4(hid_t node, hid_t tid,
                                        int dim, hsize_t* dims,
                                        hid_t mid,
                                        hid_t sid)
{
  IMPORTNUMPY;
  int s, sizem;
  PyArrayObject* r = NULL;
  sizem = 1;
  for (s = 0; s < dim; s++) sizem = sizem*dims[s];
  E_LONG* ptr = (E_LONG*)::malloc(sizem*sizeof(E_LONG));
  if (!ptr) { Py_INCREF(Py_None); return Py_None; }

  hid_t did,yid;
  did = H5Dopen2(node, L3S_DATA, H5P_DEFAULT);
  yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  
#if defined(_MPI) && defined(H5_HAVE_PARALLEL)  
  if (_ismpi == 1)    // HDF is executed in parallel context and compiled in MPI
  {
    hid_t xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    hid_t ret        = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    H5Dread(did, yid, mid, sid, xfer_plist, ptr);
  }
  else // HDF is executed in sequential context and compiled in MPI
  {
    H5Dread(did, yid, mid, sid, H5P_DEFAULT, ptr);
  }
#else // HDF is executed in sequential context and compiled in sequential
  H5Dread(did, yid, mid, sid, H5P_DEFAULT, ptr);
#endif
  
  H5Tclose(yid); H5Dclose(did);

  int* ptr2 = (int*)::malloc(sizem*sizeof(int));
  for (int n = 0; n < sizem; n++)
  {
    ptr2[n] = (int)(ptr[n]);
  }
  free(ptr);

  // Create numpy : toujours en INT
  vector<npy_intp> npy_dim_vals(dim);
  for (E_Int nn = 0; nn < dim; nn++)
  {
    npy_dim_vals[nn] = dims[nn];
  }
  r = (PyArrayObject*)PyArray_EMPTY(dim, &npy_dim_vals[0], NPY_INT, 1);
  memcpy(PyArray_DATA(r), ptr2, sizem*sizeof(int));
  free(ptr2);

  return (PyObject*)r;
}

//=============================================================================
// Get arrayi8, le convertit en i4 avec warning check
PyObject* K_IO::GenIOHdf::getArrayI82I4C(hid_t node, hid_t tid,
                                         int dim, hsize_t* dims,
                                         hid_t mid,
                                         hid_t sid)
{
  IMPORTNUMPY;
  int s, sizem;
  PyArrayObject* r = NULL;
  sizem = 1;
  for (s = 0; s < dim; s++) sizem = sizem*dims[s];
  E_LONG* ptr = (E_LONG*)::malloc(sizem*sizeof(E_LONG));
  if (!ptr) { Py_INCREF(Py_None); return Py_None; }

  hid_t did,yid;
  did = H5Dopen2(node, L3S_DATA, H5P_DEFAULT);
  yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  
#if defined(_MPI) && defined(H5_HAVE_PARALLEL)  
  if (_ismpi == 1)    // HDF is executed in parallel context and compiled in MPI
  {
    hid_t xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    hid_t ret        = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    H5Dread(did, yid, mid, sid, xfer_plist, ptr);
  }
  else // HDF is executed in sequential context and compiled in MPI
  {
    H5Dread(did, yid, mid, sid, H5P_DEFAULT, ptr);
  }
#else // HDF is executed in sequential context and compiled in sequential
  H5Dread(did, yid, mid, sid, H5P_DEFAULT, ptr);
#endif
  
  H5Tclose(yid); H5Dclose(did);

  int* ptr2 = (int*)::malloc(sizem*sizeof(int));
  E_Boolean exceed = false;
  int64_t val;
  int64_t maxInt = 2;
  maxInt = (maxInt<<31)-2;
  for (int64_t n = 0; n < sizem; n++)
  {
    val = ptr[n];
    if (val > maxInt) exceed = true;
    ptr2[n] = (int)(val);
  }
  free(ptr);

  if (exceed == true) printf("Warning: int64 value badly converted to int32.\n");

  // Create numpy : toujours en INT
  vector<npy_intp> npy_dim_vals(dim);
  for (E_Int nn = 0; nn < dim; nn++)
  {
    npy_dim_vals[nn] = dims[nn];
  }
  r = (PyArrayObject*)PyArray_EMPTY(dim, &npy_dim_vals[0], NPY_INT, 1);
  memcpy(PyArray_DATA(r), ptr2, sizem*sizeof(int));
  free(ptr2);

  return (PyObject*)r;
}

//=============================================================================
// drivers
PyObject* K_IO::GenIOHdf::getArrayI8(hid_t node, hid_t tid,
                                     int dim, hsize_t* dims,
                                     hid_t mid,
                                     hid_t sid)
{
#ifdef E_DOUBLEINT
  if (_readIntMode == 0) return getArrayI8Raw(node, tid, dim, dims, mid, sid);
  else return getArrayI8Raw(node, tid, dim, dims, mid, sid);
#else  
  if (_readIntMode == 0) return getArrayI82I4C(node, tid, dim, dims, mid, sid);
  else return getArrayI8Raw(node, tid, dim, dims, mid, sid);
#endif
}

PyObject* K_IO::GenIOHdf::getArrayI4(hid_t node, hid_t tid,
                                     int dim, hsize_t* dims,
                                     hid_t mid,
                                     hid_t sid)
{
#ifdef E_DOUBLEINT
  if (_readIntMode == 0) return getArrayI42I8(node, tid, dim, dims, mid, sid);
  else return getArrayI4Raw(node, tid, dim, dims, mid, sid);
#else
  if (_readIntMode == 0) return getArrayI4Raw(node, tid, dim, dims, mid, sid);
  else return getArrayI4Raw(node, tid, dim, dims, mid, sid);
#endif
}

//=============================================================================
PyObject* K_IO::GenIOHdf::getArrayR4Skel(hid_t node, hid_t tid,
                                         int dim, hsize_t* dims)
{
  if (_maxFloatSize == 0) { Py_INCREF(Py_None); return Py_None; }
  int  s, sizem;
  sizem = 1;
  for (s = 0; s < dim; s++) sizem = sizem*dims[s];
  if (sizem < _maxFloatSize) return getArrayR4(node, tid, dim, dims);
  else { Py_INCREF(Py_None); return Py_None; }
}

//=============================================================================
// avec conversion des R4 en R8
PyObject* K_IO::GenIOHdf::getArrayR42R8(hid_t node, hid_t tid,
                                        int dim, hsize_t* dims,
                                        hid_t mid,
                                        hid_t sid)
{
  IMPORTNUMPY;
  int s, sizem;
  PyArrayObject* r = NULL;

  sizem = 1;
  for (s = 0; s < dim; s++) sizem = sizem*dims[s];
  float* ptr = (float*)::malloc(sizem*sizeof(float));
  if (!ptr) { Py_INCREF(Py_None); return Py_None; }

  hid_t did, yid;
  did = H5Dopen2(node, L3S_DATA, H5P_DEFAULT);
  yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  
#if defined(_MPI) && defined(H5_HAVE_PARALLEL)  
  if (_ismpi == 1)    // HDF is executed in parallel context and compiled in MPI
  {
    hid_t xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    hid_t ret        = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    H5Dread(did, yid, mid, sid, xfer_plist, ptr);
  }
  else // HDF is executed in sequential context and compiled in MPI
  {
    H5Dread(did, yid, mid, sid, H5P_DEFAULT, ptr);
  }
#else  // HDF is executed in sequential context and compiled in sequential
  H5Dread(did, yid, mid, sid, H5P_DEFAULT, ptr);
#endif
  H5Tclose(yid); H5Dclose(did);

  double* ptr2 = (double*)::malloc(sizem*sizeof(double));
  for (int n = 0; n < sizem; n++)
  {
    ptr2[n] = (double)(ptr[n]);
  }
  free(ptr);

  // Create numpy: toujours en DOUBLE
  vector<npy_intp> npy_dim_vals(dim);
  for (E_Int nn = 0; nn < dim; nn++)
  {
    npy_dim_vals[nn] = dims[nn];
  }
  r = (PyArrayObject*)PyArray_EMPTY(dim, &npy_dim_vals[0], NPY_DOUBLE, 1);
  memcpy(PyArray_DATA(r), ptr2, sizem*sizeof(double));
  free(ptr2);

  return (PyObject*)r;
}
//=============================================================================
PyObject* K_IO::GenIOHdf::getArrayR4Raw(hid_t node, hid_t tid, int dim, hsize_t* dims,
                                        hid_t mid,  /* mem_space_id */
                                        hid_t sid)
{
  IMPORTNUMPY;
  PyArrayObject* r = NULL;
  vector<npy_intp> npy_dim_vals(dim);
  hid_t did, yid;
  for (E_Int nn = 0; nn < dim; nn++) npy_dim_vals[nn] = dims[nn];
  r = (PyArrayObject*)PyArray_EMPTY(dim, &npy_dim_vals[0], NPY_FLOAT, 1);

  did = H5Dopen2(node, L3S_DATA, H5P_DEFAULT);
  yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);

#if defined(_MPI) && defined(H5_HAVE_PARALLEL)
  if (_ismpi == 1)    /** HDF is executed in parallel context and compiled in MPI **/
  {
    hid_t xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    hid_t ret        = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    //hid_t ret        = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_INDEPENDENT);
    H5Dread(did, yid, mid, sid, xfer_plist, PyArray_DATA(r));
  }     
  else  /** HDF is executed in sequential context and compiled in MPI **/
  {
    H5Dread(did, yid, mid, sid, H5P_DEFAULT, PyArray_DATA(r));
  }
#else /** HDF is executed in sequential context and compiled in sequential **/
  H5Dread(did, yid, mid, sid, H5P_DEFAULT, PyArray_DATA(r));
#endif
  H5Tclose(yid); H5Dclose(did);
  
  return (PyObject*)r;
}

//=============================================================================
PyObject* K_IO::GenIOHdf::getArrayR8Skel(hid_t node, hid_t tid,
                                         int dim, hsize_t* dims)
{
  if (_maxFloatSize == 0) { Py_INCREF(Py_None); return Py_None; }
  int  s;
  hsize_t sizem;
  sizem = 1;
  for (s = 0; s < dim; s++) sizem = sizem*dims[s];
  if (sizem < (hsize_t)_maxFloatSize) return getArrayR8(node, tid, dim, dims);
  else { Py_INCREF(Py_None); return Py_None; }
}

//=============================================================================
PyObject* K_IO::GenIOHdf::getArrayI1Skel(hid_t node, hid_t tid,
                                         int dim, hsize_t* dims)
{
  if (_maxFloatSize == 0) { Py_INCREF(Py_None); return Py_None; }
  int  s;
  hsize_t sizem;
  sizem = 1;
  for (s = 0; s < dim; s++) sizem = sizem*dims[s];
  if (sizem < (hsize_t)_maxFloatSize) return getArrayI1(node, tid, dim, dims);
  else { Py_INCREF(Py_None); return Py_None; }
}

//=============================================================================
PyObject* K_IO::GenIOHdf::getArrayI8Skel(hid_t node, hid_t tid,
                                         int dim, hsize_t* dims)
{
  if (_maxFloatSize == 0) { Py_INCREF(Py_None); return Py_None; }
  int  s;
  hsize_t sizem;
  sizem = 1;
  for (s = 0; s < dim; s++) sizem = sizem*dims[s];
  if (sizem < (hsize_t)_maxFloatSize) return getArrayI8(node, tid, dim, dims);
  else { Py_INCREF(Py_None); return Py_None; }
}

//=============================================================================
PyObject* K_IO::GenIOHdf::getArrayI4Skel(hid_t node, hid_t tid,
                                         int dim, hsize_t* dims)
{
  if (_maxFloatSize == 0) { Py_INCREF(Py_None); return Py_None; }
  int  s;
  hsize_t sizem;
  sizem = 1;
  for (s = 0; s < dim; s++) sizem = sizem*dims[s];
  if (sizem < (hsize_t)_maxFloatSize) return getArrayI4(node, tid, dim, dims);
  else { Py_INCREF(Py_None); return Py_None; }
}

//=============================================================================
PyObject* K_IO::GenIOHdf::getArrayR8(hid_t node, hid_t tid, int dim, hsize_t* dims,
                                     hid_t mid,  /* mem_space_id */
                                     hid_t sid)
{
  IMPORTNUMPY;
  PyArrayObject* r = NULL;
  vector<npy_intp> npy_dim_vals(dim);
  hid_t did, yid;
  // for (E_Int nn = 0; nn < dim; nn++) printf("getArrayR8 of dims :: %d \n", dims[nn]);
  // printf("K_IO::GenIOHdf::getArrayR8 \n ");
  // Create numpy: toujours en DOUBLE
  for (E_Int nn = 0; nn < dim; nn++) npy_dim_vals[nn] = dims[nn];
  r = (PyArrayObject*)PyArray_EMPTY(dim, &npy_dim_vals[0], NPY_DOUBLE, 1);

  did = H5Dopen2(node, L3S_DATA, H5P_DEFAULT);
  yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);

#if defined(_MPI) && defined(H5_HAVE_PARALLEL)
  if (_ismpi == 1)    /** HDF is executed in parallel context and compiled in MPI **/
  {
    hid_t xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    hid_t ret        = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    //hid_t ret        = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_INDEPENDENT);
    H5Dread(did, yid, mid, sid, xfer_plist, PyArray_DATA(r));
  }     
  else  /** HDF is executed in sequential context and compiled in MPI **/
  {
    H5Dread(did, yid, mid, sid, H5P_DEFAULT, PyArray_DATA(r));
  }
#else /** HDF is executed in sequential context and compiled in sequential **/
  H5Dread(did, yid, mid, sid, H5P_DEFAULT, PyArray_DATA(r));
#endif
  H5Tclose(yid); H5Dclose(did);
  
  return (PyObject*)r;
}

//=============================================================================
char* K_IO::GenIOHdf::getArrayC1(hid_t node, hid_t tid, int dim, hsize_t* dims)
{
  IMPORTNUMPY;
  int s, sizem;

  sizem = 1;
  for (s = 0; s < dim; s++) sizem = sizem*dims[s];
  char* ptr = (char*)::malloc((sizem+1)*sizeof(char));

  hid_t did, yid;
  did = H5Dopen2(node, L3S_DATA, H5P_DEFAULT);
  yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);

#if defined(_MPI) && defined(H5_HAVE_PARALLEL)
  if (_ismpi == 1)    /** HDF is executed in parallel context and compiled in MPI **/
  {
    hid_t xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    hid_t ret        = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    H5Dread(did, yid, H5S_ALL, H5S_ALL, xfer_plist, ptr);
  }     /** HDF is executed in sequential context and compiled in MPI **/
  else
  {
    H5Dread(did, yid, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptr);
  }
#else /** HDF is executed in sequential context and compiled in sequential **/
   H5Dread(did, yid, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptr);
#endif
  H5Tclose(yid); H5Dclose(did);
  ptr[sizem] = '\0';
  return ptr;
}

// ===============================================================================
// Generic function to get data in HDF - NO cast - Allows contiguous or Interlaced 
// ===============================================================================
PyObject* K_IO::GenIOHdf::getArrayContigous(hid_t     node, 
                                            hid_t     tid, 
                                            int       dim,
                                            hsize_t*  dims,
                                            int       NPYtype,
                                            hid_t     mid,  /* mem_space_id */
                                            hid_t     sid, 
                                            PyObject* data)
{
  IMPORTNUMPY;
  vector<npy_intp> npy_dim_vals(dim);
  hid_t            did, yid;

  // Create numpy: No cast 
  if (data == NULL)
  {
    // printf("Allocate data  ... \n");
    // for (E_Int nn = 0; nn < dim; nn++) printf("getArrayR8 of dims :: %d \n", dims[nn]);
    for (E_Int nn = 0; nn < dim; nn++) npy_dim_vals[nn] = dims[nn];
    data = PyArray_EMPTY(dim, &npy_dim_vals[0], NPYtype, 1);
  }

  did = H5Dopen2(node, L3S_DATA, H5P_DEFAULT);
  yid = H5Tget_native_type(tid, H5T_DIR_ASCEND);

#if defined(_MPI) && defined(H5_HAVE_PARALLEL)
  if (_ismpi == 1) /** HDF is executed in parallel context and compiled in MPI **/
  {
   // printf("getArrayContigous H5_HAVE_PARALLEL / _ismpi ON \n ");
   hid_t xfer_plist = H5Pcreate(H5P_DATASET_XFER);
   hid_t ret        = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
   H5Dread(did, yid, mid, sid, xfer_plist, PyArray_DATA( (PyArrayObject*) data));
 }
 else
 {
   H5Dread(did, yid, mid, sid, H5P_DEFAULT, PyArray_DATA( (PyArrayObject*) data));
 }
#else
  H5Dread(did, yid, mid, sid, H5P_DEFAULT, PyArray_DATA((PyArrayObject*) data));
#endif
  H5Tclose(yid); H5Dclose(did);

  return data;
}

//=============================================================================
/*
   hdfcgnsread
   IN: file: file name
   OUT: tree, arbre CGNS/python charge
   OUT: dataShape: si dataShape != NULL, retourne les chemins de noeuds + shape
   OUT: links: si links != NULL, retourne les chemins des noeuds de link
   IN: skeleton: 0 (full), 1 (only skeleton loaded)
   IN: maxFloatSize: si skeleton=1, load si shape < maxFloatSize
   IN: maxDepth: profondeur max de load
   IN: readIntMode: 0: convert int to Cassiopee compile type, 1: read as in file
   IN: skipTypes: types to skip
*/
//=============================================================================
E_Int K_IO::GenIO::hdfcgnsread(char* file, PyObject*& tree, PyObject* dataShape, PyObject* links, 
                               int skeleton, int maxFloatSize, int maxDepth, int readIntMode,
                               PyObject* skipTypes)
{
  tree = PyList_New(4);

  /* Open file */
  hid_t fapl, fid;
  fapl = H5Pcreate(H5P_FILE_ACCESS);
  fid = H5Fopen(file, H5F_ACC_RDONLY, fapl);
  //printf("open Avant =%d\n",H5Fget_obj_count(fid, H5F_OBJ_ALL));
  // printf("open=%d\n",H5Fget_obj_count(fid, H5F_OBJ_ALL));
  H5Pclose(fapl);
    
  if (fid < 0)
  {
    printf("Warning: hdfcgnsread: can not open file %s.\n", file);
    return 1;
  }
  hid_t gid = H5Gopen(fid, "/", H5P_DEFAULT);
  GenIOHdf HDF;
  HDF._readIntMode = readIntMode;

  /* Prepare skip types */
  if (skipTypes != NULL)
  {
    E_Int skipSize = PyList_Size(skipTypes);
    for (E_Int i = 0; i < skipSize; i++)
    {
      char* typeToSkip = NULL;
      char* nameToSkip = NULL;
      PyObject* l = PyList_GetItem(skipTypes, i);
      if (PyString_Check(l)) typeToSkip = PyString_AsString(l);
#if PY_VERSION_HEX >= 0x03000000
      else if (PyUnicode_Check(l)) typeToSkip = (char*)PyUnicode_AsUTF8(l); 
#endif
      else if (PyTuple_Check(l)) 
      {
        /* Get name */
        PyObject* m1 = PyTuple_GetItem(l, 0);
        if (PyString_Check(m1)) nameToSkip = PyString_AsString(m1);
#if PY_VERSION_HEX >= 0x03000000
        else if (PyUnicode_Check(m1)) nameToSkip = (char*)PyUnicode_AsUTF8(m1);
#endif
        /* Get type */
        PyObject* m2 = PyTuple_GetItem(l, 1);
        if (PyString_Check(m2)) typeToSkip = PyString_AsString(m2);
#if PY_VERSION_HEX >= 0x03000000
        else if (PyUnicode_Check(m2)) typeToSkip = (char*)PyUnicode_AsUTF8(m2);
#endif
      }
      /* Create map to skip */
      if (nameToSkip != NULL) 
        HDF._skipNameAndTypes[make_pair(string(nameToSkip), string(typeToSkip))] = true;
      else HDF._skipTypes[string(typeToSkip)] = true;
    }
  }

  /* Recursive load */
  HDF._skeleton = skeleton;
  HDF._ismpi = 0;
  HDF._maxFloatSize = maxFloatSize;
  if (maxDepth == -1) HDF._maxDepth = 1e6;
  else HDF._maxDepth = maxDepth;
  HDF._fatherStack.push_front(gid);
  HDF._stringStack.push_front("");
  HDF.loadOne(tree, 0, dataShape, links);
  HDF._fatherStack.pop_front();
  HDF._stringStack.pop_front();

  /* Close */
  HDF._fatherStack.clear();

  // DBX
  /*
  E_Int opened = H5Fget_obj_count(fid, H5F_OBJ_ALL);
  printf("opened=%d\n", opened);
  if (opened > 0)
  {
    hid_t* objIds = new hid_t [opened];
    E_Int howmany = H5Fget_obj_ids(fid, H5F_OBJ_ALL, opened, objIds);
    char name[1024];
    for (E_Int i = 0; i < howmany; i++) 
    {
      hid_t anobj = *objIds++;
      H5I_type_t ot = H5Iget_type(anobj);
      herr_t status = H5Iget_name(anobj, name, 1024);
      printf(" %d: type %d, name %s\n",i,ot,name);
    }  
    delete [] objIds;
  }
  */
  // END DBX

  H5Fclose(fid);  
  return 0;
}

//==========================================================================
// Retourne le groupe qui correspond a path a partir d'un groupe start
// path: /Base/cart12
//==========================================================================
hid_t K_IO::GenIOHdf::openGroupWithLinks(hid_t start, char* path)
{
  //printf("%s\n", path); fflush(stdout);
  char _dtype[CGNSMAXLABEL+1]; char _name[CGNSMAXLABEL+1];

  hid_t gid = H5Gopen(start, "/", H5P_DEFAULT);
  
  char local[L3C_MAX_ATTRIB_SIZE+1]; E_Int c = 0;
  char* p = path;
  while (*p != '\0')
  {
    if (*p == '/')
    {
      local[c] = '\0';
      if (strlen(local) > 0)
      {
        start = gid;
        HDF_Get_Attribute_As_String(start, L3S_NAME, _name);
        HDF_Get_Attribute_As_String(start, L3S_DTYPE, _dtype);
        //printf("current: %s %s\n", _name, _dtype);
        //printf("trying to open : %s\n", local);

        if (strcmp(_dtype, "LK") != 0)
        {
          gid = H5Gopen(start, local, H5P_DEFAULT);
          H5Gclose(start);
        }
        else
        {
          H5L_info_t sb;
          herr_t herr = H5Lget_info(start, L3S_LINK, &sb, H5P_DEFAULT);

          if (herr < 0)
          {
            printf("Error: hdfcgnsread: error opening link file.\n");
          }
          hid_t lid = H5Gopen2(start, L3S_LINK, H5P_DEFAULT);
          H5Gclose(start); //H5Gclose(herr);
          gid = H5Gopen(lid, local, H5P_DEFAULT);
        }
      }
      c = 0;
    }
    else
    {
      local[c] = *p; c++;
    }
    p++;
  }
  // last
  local[c] = '\0';
  if (strlen(local) > 0)
  {
    start = gid;
    HDF_Get_Attribute_As_String(start, L3S_NAME, _name);
    HDF_Get_Attribute_As_String(start, L3S_DTYPE, _dtype);
    //printf("current: %s %s\n", _name, _dtype);
    //printf("trying to open : %s\n", local);

    if (strcmp(_dtype, "LK") != 0)
    {
      gid = H5Gopen(start, local, H5P_DEFAULT);
      H5Gclose(start);
    }
    else
    {
      H5L_info_t sb; 
      herr_t herr = H5Lget_info(start, L3S_LINK, &sb, H5P_DEFAULT);

      if (herr < 0)
      {
        printf("Error: hdfcgnsread: error opening link file.\n");
      }
      hid_t lid = H5Gopen2(start, L3S_LINK, H5P_DEFAULT);
      H5Gclose(start); //H5Gclose(herr);
      gid = H5Gopen(lid, local, H5P_DEFAULT);
    }    
  }
  return gid;
}

//=============================================================================
// Lit les paths specifies dans le fichier file.
// Retourne une liste d'objets pythons contenant les noeuds pointes par les
// chemins
//=============================================================================
PyObject* K_IO::GenIO::hdfcgnsReadFromPaths(char* file, PyObject* paths,
                                            E_Int maxFloatSize, E_Int maxDepth,
                                            E_Int readIntMode, 
                                            PyObject* dataShape,
                                            PyObject* skipTypes,
                                            PyObject* mpi4pyCom)
{
  if (PyList_Check(paths) == false)
  {
    PyErr_SetString(PyExc_TypeError,
                    "hdfread: paths must be a list of strings.");
    return NULL;
  }
  E_Int size = PyList_Size(paths);

  /* Open file */
  hid_t fapl, fid;
  fapl = H5Pcreate(H5P_FILE_ACCESS);
  
  PyObject* ret = PyList_New(0);
  GenIOHdf HDF;
  HDF._readIntMode = readIntMode;
  HDF._ismpi = 0;
  HDF._skeleton = 0;
  PyObject* node;

  /* Prepare skip types */
  if (skipTypes != NULL)
  {
    E_Int skipSize = PyList_Size(skipTypes);
    for (E_Int i = 0; i < skipSize; i++)
    {
      char* typeToSkip = NULL;
      char* nameToSkip = NULL;
      PyObject* l = PyList_GetItem(skipTypes, i);
      if (PyString_Check(l)) typeToSkip = PyString_AsString(l);
#if PY_VERSION_HEX >= 0x03000000
      else if (PyUnicode_Check(l)) typeToSkip = (char*)PyUnicode_AsUTF8(l); 
#endif
      else if (PyTuple_Check(l))
      {
        /* Get name */
        PyObject* m1 = PyTuple_GetItem(l, 0);
        if (PyString_Check(m1)) nameToSkip = PyString_AsString(m1);
#if PY_VERSION_HEX >= 0x03000000
        else if (PyUnicode_Check(m1)) nameToSkip = (char*)PyUnicode_AsUTF8(m1);
#endif
        /* Get type */
        PyObject* m2 = PyTuple_GetItem(l, 1);
        if (PyString_Check(m2)) typeToSkip = PyString_AsString(m2);
#if PY_VERSION_HEX >= 0x03000000
        else if (PyUnicode_Check(m2)) typeToSkip = (char*)PyUnicode_AsUTF8(m2);
#endif
      }
      /* Create map to skip */
      if (nameToSkip != NULL)
        HDF._skipNameAndTypes[make_pair(string(nameToSkip), string(typeToSkip))] = true;
      else HDF._skipTypes[string(typeToSkip)] = true;
    }
  }

#if defined(_MPI) && defined(H5_HAVE_PARALLEL)
  if (mpi4pyCom != Py_None && mpi4pyCom != NULL) HDF._ismpi = 1;
  else HDF._ismpi = 0;
#else
  HDF._ismpi = 0;
#endif

#if defined(_MPI) && defined(H5_HAVE_PARALLEL)
  if (HDF._ismpi == 1)
  {
    void* pt_comm = (void*)&(((PyMPICommObject*)mpi4pyCom)->ob_mpi);
    MPI_Comm comm = *((MPI_Comm*) pt_comm);
    MPI_Info info   = MPI_INFO_NULL;
    H5Pset_fapl_mpio(fapl, comm, info);
   }
#endif

  fid = H5Fopen(file, H5F_ACC_RDONLY, fapl);
  if (fid < 0)
  {
    PyErr_SetString(PyExc_TypeError, "hdfread: cannot open file.");
    return NULL;
  }
  H5Pclose(fapl);

  for (E_Int i = 0; i < size; i++)
  {
    char* path; PyObject* l;
    l = PyList_GetItem(paths, i);
    if (PyString_Check(l)) path = PyString_AsString(l);
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(l)) path = (char*)PyUnicode_AsUTF8(l);
#endif
    else
    {
      PyErr_SetString(PyExc_TypeError, "hdfcgnsread: paths must be strings.");
      return NULL;
    }
    /* Open group in HDF corresponding to path */
    hid_t gid = HDF.openGroupWithLinks(fid, path);
    if (gid < 0)
    { 
      printf("Warning: hdfcgnsread: cannot find this path.\n");
    }
    else if (maxDepth == 0)
    {
      std::string shortPath = path;
      size_t pos = shortPath.find_last_of("/");
      if (pos == std::string::npos) pos = 0; 
      shortPath = shortPath.erase(pos);
      HDF._stringStack.push_front(shortPath); // short path      
      node = HDF.createNode(gid, dataShape, NULL);
      HDF._stringStack.pop_front();

      PyObject* children = PyList_New(0);
      PyList_SetItem(node, 2, children);
      PyList_Append(ret, node); Py_DECREF(node);
    }
    else
    {
      if (maxFloatSize >= 0)
      { HDF._maxFloatSize = maxFloatSize; HDF._skeleton = 1; }
      if (maxDepth >= 0) HDF._maxDepth = maxDepth;
      else HDF._maxDepth = 1e6;

      std::string shortPath = path;
      size_t pos = shortPath.find_last_of("/");
      if (pos == std::string::npos) pos = 0;
      shortPath = shortPath.erase(pos);
      
      HDF._stringStack.push_front(shortPath); // short path            
      node = HDF.createNode(gid, dataShape, NULL);
      HDF._stringStack.pop_front();

      HDF._stringStack.push_front(path);
      HDF._fatherStack.push_front(gid);
      HDF.loadOne(node, 0, dataShape);
      HDF._fatherStack.pop_front();
      HDF._stringStack.pop_front();
      HDF._fatherStack.clear();
      HDF._stringStack.clear();
      PyList_Append(ret, node); Py_DECREF(node);
    }
  }
  H5Fclose(fid);
  return ret;
}

//=============================================================================
PyObject* K_IO::GenIOHdf::loadOne(PyObject* tree, int depth, 
                                  PyObject* dataShape, PyObject* links)
{
  hid_t* sonList;
  PyObject* node;
  PyObject* l = PyList_New(0);
  PyList_SetItem(tree, 2, l);
  hid_t father = _fatherStack.front();
  if (depth >= _maxDepth) {H5Gclose(father); return tree;}
  sonList = getChildren(father);
  int c = 0;
  if (_skipTypes.size() == 0  && _skipNameAndTypes.size() == 0)
  {
    while (sonList != NULL && sonList[c] != (hid_t)-1)
    {
      node = createNode(sonList[c], dataShape, links);
      PyList_Append(l, node); Py_DECREF(node);
      _stringStack.push_front(_currentPath);
      _fatherStack.push_front(sonList[c]);
      loadOne(node, depth+1, dataShape, links);
      _fatherStack.pop_front();
      _stringStack.pop_front();
      c++;
    }
  }
  else
  {
    /* with Skip */
    while (sonList != NULL && sonList[c] != (hid_t)-1)
    {
      node = createNode(sonList[c], dataShape, links);
      if (isAnodeToSkip())
      {
        Py_DECREF(node);
      }
      else
      {
        PyList_Append(l, node); Py_DECREF(node);
        _stringStack.push_front(_currentPath);
        _fatherStack.push_front(sonList[c]);
        loadOne(node, depth+1, dataShape, links);
        _fatherStack.pop_front();
        _stringStack.pop_front();
      }
      c++;
    }
  }

  free(sonList);
  H5Gclose(father);
  return tree;
}

//=============================================================================
bool K_IO::GenIOHdf::isAnodeToSkip()
{
  //return (_skipTypes.find(_type) != _skipTypes.end());
  bool result = false;
  if (_skipTypes.size() != 0) result = (_skipTypes.find(_type) != _skipTypes.end());
  if (_skipNameAndTypes.size() != 0)
    result = (_skipNameAndTypes.find(make_pair(_name, _type)) != _skipNameAndTypes.end());
  return result;
}

//=============================================================================
// Cree un noeud du pyTree
// node peut etre modifie dans le cas d'un lien
// Remplit dataShape si non NULL
// Remplit links si non null
//=============================================================================
PyObject* K_IO::GenIOHdf::createNode(hid_t& node, PyObject* dataShape, PyObject* links)
{
  // Nom du noeud de l'arbre
  HDF_Get_Attribute_As_String(node, L3S_NAME, _name);
  if (_stringStack.size() > 0)
  {
    std::string path = _stringStack.front();
    _currentPath = path+"/"+string(_name);
  }
  
  // Remplace les liens eventuellement (change node par le node cible)
  HDF_Get_Attribute_As_String(node, L3S_DTYPE, _dtype);
  //printf("create node %s %s\n", _name, _dtype); fflush(stdout);
  
  if (strcmp(_dtype, "LK") == 0)
  {
    // store link data
    if (links != NULL)
    {
      H5L_info_t lk;
      char querybuff[512];
      const char *file; const char *path;
      char rfile[512]; char rpath[512];
      H5Lget_info(node, L3S_LINK, &lk, H5P_DEFAULT);
      H5Lget_val(node, L3S_LINK, querybuff, sizeof(querybuff), H5P_DEFAULT);
      if (lk.type == H5L_TYPE_EXTERNAL)
      {
        H5Lunpack_elink_val(querybuff, lk.u.val_size, NULL, &file, &path);
        strcpy(rfile, file); strcpy(rpath, path);
      }
      else
      {
        strcpy(rpath, querybuff);
        rfile[0] = '\0';
      }

      //printf("link %s\n", _name);
      //printf("link file: %s\n", rfile);
      //printf("link current path: %s\n", _currentPath.c_str());
      //printf("link target path: %s\n", rpath);
      
      char lksearch[128]; strcpy(lksearch, ".");
      // list of ['targetdirectory', 'targetfilename', 'targetpath', 'currentpath']
      PyObject* v = Py_BuildValue("[s,s,s,s]", lksearch, rfile, rpath, _currentPath.c_str());
      PyList_Append(links, v); Py_DECREF(v);
    }

    // follow link
    H5L_info_t sb;
    herr_t herr = H5Lget_info(node, L3S_LINK, &sb, H5P_DEFAULT);
    
    if (herr < 0)
    {
      printf("Error: hdfcgnsread: error opening file referenced in links.\n");
      PyObject* s = Py_BuildValue("[sOOs]", _name, Py_None, Py_None, _type);
      return s;
    }
    else
    {
      hid_t lid = H5Gopen2(node, L3S_LINK, H5P_DEFAULT);
      if (lid < 0)
      {
        printf("Error: node %s referenced by link was not found.\n", _name);
        PyObject* s = Py_BuildValue("[sOOs]", _name, Py_None, Py_None, _type);
        return s;
      }
      H5Gclose(node); //H5Gclose(herr);
      node = lid;
    }
  }
  /* fin des liens */

  // Type du noeud
  HDF_Get_Attribute_As_String(node, L3S_LABEL, _type);

  // Valeur du noeud
  int dim = 0;
  //npy_intp npy_dim_vals[1]; npy_dim_vals[0] = 1;
  npy_intp npy_dim_vals2[1];
  int d;
  HDF_Get_Attribute_As_String(node, L3S_DTYPE, _dtype);

  hid_t tid;
  tid = ADF_to_HDF_datatype(_dtype);

  if (strcmp(_dtype, L3T_MT) != 0)
  {
    HDF_Get_DataDimensions(node, _dims2);
    for (d = 0; d < CGNSMAXDIM; d++)
    { if (_dims2[d] == -1) break; }
    dim = d;

    // inverse les dimensions (indexation fortran)
    for (d = 0; d < dim; d++) _dims[d] = _dims2[dim-d-1];
    //for (d = 0; d < dim; d++)  printf("%d \n", dims[d]);
    //printf("\n");
  }

  PyObject* v = NULL;
  if (strcmp(_dtype, L3T_I4) == 0)
  {
    if (_skeleton == 1 && (strcmp(_type, "IndexArray_t") == 0 || strcmp(_type, "DataArray_t") == 0)) v = getArrayI4Skel(node, tid, dim, _dims);
    else v = getArrayI4(node, tid, dim, _dims);
  }
  else if (strcmp(_dtype, L3T_R8) == 0)
  {
    if (_skeleton == 1 && strcmp(_type, "DataArray_t") == 0) v = getArrayR8Skel(node, tid, dim, _dims);
    else v = getArrayR8(node, tid, dim, _dims);
  }
  else if (strcmp(_dtype, L3T_I8) == 0)
  {
    if (_skeleton == 1 && (strcmp(_type, "IndexArray_t") == 0 || strcmp(_type, "DataArray_t") == 0)) v = getArrayI8Skel(node, tid, dim, _dims);
    else v = getArrayI8(node, tid, dim, _dims);
  }
  else if (strcmp(_dtype, L3T_R4) == 0)
  {
    if (_skeleton == 1 && strcmp(_type, "DataArray_t") == 0) v = getArrayR4Skel(node, tid, dim, _dims);
    else v = getArrayR4(node, tid, dim, _dims);
  }
  else if (strcmp(_dtype, L3T_I1) == 0)
  {
    if (_skeleton == 1 && (strcmp(_type, "IndexArray_t") == 0 || strcmp(_type, "DataArray_t") == 0)) v = getArrayI1Skel(node, tid, dim, _dims);
    else v = getArrayI1(node, tid, dim, _dims);
  }
  else if (strcmp(_dtype, L3T_C1) == 0)
  {  
    IMPORTNUMPY;
    char* s = getArrayC1(node, tid, dim, _dims);
    if (dim == 1)
    {
      E_Int l = strlen(s); npy_dim_vals2[0] = l;
      v = PyArray_EMPTY(1, npy_dim_vals2, NPY_STRING, 1);
      memcpy(PyArray_DATA((PyArrayObject*)v), s, l*sizeof(char));
      free(s);
    }
    else 
    {
      npy_intp* npy_dim_vals = new npy_intp[dim];
      E_Int l = 1;
      for (E_Int i = 0; i < dim; i++) 
      { npy_dim_vals[i] = _dims[i]; l = l*_dims[i]; } 
      v = PyArray_EMPTY(dim, npy_dim_vals, NPY_STRING, 1);
      memcpy(PyArray_DATA((PyArrayObject*)v), s, l*sizeof(char));
      free(s); delete [] npy_dim_vals;
    }
  }
  else if (strcmp(_dtype, L3T_MT) == 0)
  {
    v = Py_None; Py_INCREF(Py_None);
  }
  else
  {
    printf("Warning: hdfcgnsread: unknown type of node: %s.\n", _dtype);
  }

  if (dataShape != NULL)
  {
    // printf("Ajoute la shape au chemin \n");
    PyObject* result = PyTuple_New(3);
    PyObject* stid   = Py_BuildValue("h", 1);
    PyObject* type   = Py_BuildValue("s", _dtype);
        
    PyTuple_SetItem(result, 0, stid);
    PyTuple_SetItem(result, 1, type); 

    PyObject* shape = PyTuple_New(dim);
    for (int i = 0; i < dim; i++) 
    {
      // PyObject* temp = Py_BuildValue("%d", _dims[i]);
      PyObject* temp = Py_BuildValue("h", _dims[i]);
      PyTuple_SetItem(shape, i, temp);
    }
    PyTuple_SetItem(result, 2, shape); 
    PyDict_SetItemString(dataShape, _currentPath.c_str(), result); 
  }

  // if (tid != 0) H5Tclose(tid);
  PyObject* s = Py_BuildValue("[sOOs]", _name, v, Py_None, _type);
  Py_DECREF(v);
  return s;
}

//=============================================================================
/*
   hdfcgnswrite
*/
//=============================================================================
E_Int K_IO::GenIO::hdfcgnswrite(char* file, PyObject* tree, PyObject* links,
                                int writeIntMode, int writeRealMode)
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
  //H5Pset_libver_bounds(fapl, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
  H5Pset_libver_bounds(fapl, KHDFVERSION, KHDFVERSION);

  capl = H5Pcreate(H5P_FILE_CREATE);
  H5Pset_link_creation_order(capl, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);

  //lapl = H5Pcreate(H5P_LINK_ACCESS);
  //H5Pset_nlinks(lapl, ADF_MAXIMUM_LINK_DEPTH);

  fid = H5Fcreate(file, H5F_ACC_TRUNC, capl, fapl);
  H5Pclose(fapl); H5Pclose(capl);

  if (fid < 0)
  {
    printf("Warning: hdfcgnswrite: can not open file %s.\n", file);
    return 1;
  }

  GenIOHdf HDF;
  HDF._ismpi = 0;
  HDF._maxDepth = 1e6;
  HDF._writeIntMode = writeIntMode;
  HDF._writeRealMode = writeRealMode;
  //printf("writeIntMode=%d\n", writeIntMode);
  //printf("writeRealMode=%d\n", writeRealMode);
  
  // Ajout version... au root node
  hid_t gid = H5Gopen2(fid, "/", H5P_DEFAULT);
  HDF_Add_Attribute_As_String(gid, L3S_NAME, L3S_ROOTNODENAME);
  HDF_Add_Attribute_As_String(gid, L3S_LABEL, L3S_ROOTNODETYPE);
  HDF_Add_Attribute_As_String(gid, L3S_DTYPE, L3T_MT);
  char format[CGNSMAXLABEL];
  hid_t type = H5Tcopy(H5T_NATIVE_FLOAT);
  if (H5Tequal(type, H5T_IEEE_F32BE))
    strcpy(format, "IEEE_BIG_32");
  else if (H5Tequal(type, H5T_IEEE_F32LE))
    strcpy(format, "IEEE_LITTLE_32");
  else if (H5Tequal(type, H5T_IEEE_F64BE))
    strcpy(format, "IEEE_BIG_64");
  else if (H5Tequal(type, H5T_IEEE_F64LE))
    strcpy(format, "IEEE_LITTLE_64");
  else
    sprintf(format, "NATIVE_%d", (int)H5Tget_precision(type));
  H5Tclose(type);
  HDF.setArrayC1(gid, format, (char*)L3S_FORMAT);
  unsigned int maj, min, rel;
  H5get_libversion(&maj, &min, &rel);
  char version[CGNSMAXLABEL+1];
  memset(version, 0, CGNSMAXLABEL+1);
  sprintf(version, "HDF5 Version %u.%u.%u", maj, min, rel);
  HDF.setArrayC1(gid, version, (char*)L3S_VERSION);

  PyObject* o;
  int listsize = PyList_Size(tree);
  for (int n = 0; n < listsize; n++) // pour chaque Base
  {
    o = PyList_GetItem(tree, n);
    HDF._stringStack.push_front("");
    HDF._fatherStack.push_front(gid);
    HDF.dumpOne(o, 0, links);
    HDF._fatherStack.pop_front();
    HDF._stringStack.pop_front();
  }
  
  //H5Gclose(gid);

  /** Manage links (a l'ecriture) */
  /* List of ['targetdirectory', 'targetfilename', 'targetpath', 'currentpath',0] */
  E_Int size;
  if (links == NULL) size = 0;
  else size = PyList_Size(links);
  
  for (E_Int i = 0; i < size; i++)
  {
    PyObject* llink  = PyList_GetItem(links, i);
    char* dir_file; char* tgt_file; char* tgt_path; char* cur_path;
    PyObject* l = PyList_GetItem(llink, 0);
    if (PyString_Check(l)) dir_file = PyString_AsString(l);
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(l)) dir_file = (char*)PyUnicode_AsUTF8(l);
#endif
    else dir_file = NULL;
    l = PyList_GetItem(llink, 1);
    if (PyString_Check(l)) tgt_file = PyString_AsString(l);
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(l)) tgt_file = (char*)PyUnicode_AsUTF8(l);
#endif
    else tgt_file = NULL;
    char* tgt_file_all = NULL;
    if (tgt_file != NULL && dir_file != NULL)
    {
      E_Int size = strlen(tgt_file)+strlen(dir_file)+2;
      tgt_file_all = new char [size];
      strcpy(tgt_file_all, dir_file);
      strcat(tgt_file_all, "/");
      strcat(tgt_file_all, tgt_file);
    }
    else if (tgt_file != NULL)
    {
      E_Int size = strlen(tgt_file)+1;
      tgt_file_all = new char [size];
      strcpy(tgt_file_all, tgt_file);
    }
    //printf("all: %s\n", tgt_file_all);

    l = PyList_GetItem(llink, 2);
    if (PyString_Check(l)) tgt_path = PyString_AsString(l);
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(l)) tgt_path = (char*)PyUnicode_AsUTF8(l);
#endif
    else tgt_path = NULL;
    l = PyList_GetItem(llink, 3);
    if (PyString_Check(l)) cur_path = PyString_AsString(l);
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(l)) cur_path = (char*)PyUnicode_AsUTF8(l);
#endif
    else cur_path = NULL;
    
    // > Dramatic verbose
    //printf(" tgt_file_all : %s \n", tgt_file_all);
    //printf(" tgt_path : %s \n", tgt_path);
    //printf(" cur_path : %s \n", cur_path);
    
    /* Rip end of path to get parent */
    char* startPath; char* name;
    ripEndOfPath(cur_path, startPath);
    getEndOfPath(cur_path, name);
    hid_t gidp = H5Gopen(fid, startPath, H5P_DEFAULT);
    delete [] startPath;
    
    /* node existe deja */
    if (has_child(gidp, name))
    {
#if H5_VERSION_LE(1,11,9)
      H5Giterate(gidp, name, NULL, delete_children, NULL);
      H5Gunlink(gidp, name);
#else
      H5Literate_by_name2(gidp, name, H5_INDEX_CRT_ORDER, H5_ITER_INC, NULL, delete_children, NULL, H5P_DEFAULT);
      H5Ldelete(gidp, name, H5P_DEFAULT);
#endif
    }

    /* Create link node */
    hid_t nid = H5Gcreate2(gidp, name, H5P_DEFAULT, HDF._group, H5P_DEFAULT);
    if (nid < 0) {printf("Error: nid is invalid.\n");}
    
    HDF_Add_Attribute_As_String(nid, L3S_NAME, name);
    HDF_Add_Attribute_As_String(nid, L3S_DTYPE, L3T_LK);
    HDF_Add_Attribute_As_String(nid, L3S_LABEL, "");
  
    delete [] name;
  
    /** Make the link effective **/
    //HDF_Add_Attribute_As_Data(nid, L3S_PATH, cur_path, strlen(cur_path));
    HDF_Add_Attribute_As_Data(nid, L3S_PATH, tgt_path, strlen(tgt_path));
    
    if (strcmp(tgt_file_all, "./") == 0 || strcmp(tgt_file_all, "/") == 0)
    {
      // soft
      H5Lcreate_soft(tgt_path, nid, L3S_LINK, H5P_DEFAULT, H5P_DEFAULT);
    }
    else
    {
      // external
      H5Lcreate_external(tgt_file_all, tgt_path, nid, L3S_LINK, H5P_DEFAULT, H5P_DEFAULT);
      HDF_Add_Attribute_As_Data(nid, L3S_FILE, tgt_file, strlen(tgt_file));
    }

    H5Gclose(gidp); H5Gclose(nid);
    delete [] tgt_file_all;
      
  } /* end link */

  // DBX
  /*
  E_Int opened = H5Fget_obj_count(fid, H5F_OBJ_ALL);
  printf("opened=" SF_D_ "\n", opened);
  if (opened > 0)
  {
    hid_t* objIds = new hid_t [opened];
    E_Int howmany = H5Fget_obj_ids(fid, H5F_OBJ_ALL, opened, objIds);
    char name[1024];
    for (E_Int i = 0; i < howmany; i++ ) 
    {
      hid_t anobj = *objIds++;
      H5I_type_t ot = H5Iget_type(anobj);
      herr_t status = H5Iget_name(anobj, name, 1024);
      printf(" %d: type %d, name %s\n", i, ot, name);
    }  
    delete [] objIds;
  }
  */
  // END DBX

  H5Fclose(fid);
  
  return 0;
}
//=============================================================================
// soit un chemin /A/B/C retourne /A/B
//=============================================================================
void K_IO::GenIO::ripEndOfPath(char* path, char*& startPath)
{
  E_Int i, j;
  E_Int l = strlen(path);
  startPath = new char [l+1];
  for (i = l-1; i >= 0; i--)
  { if (path[i] == '/') break; }
  for (j = 0; j < i; j++) startPath[j] = path[j];
  startPath[i] = '\0';
}

//=============================================================================
// soit un chemin /A/B/C retourne C
//=============================================================================
void K_IO::GenIO::getEndOfPath(char* path, char*& EndPath)
{
  E_Int i, j;
  E_Int l = strlen(path);
  EndPath = new char [l+1];
  for (i = l-1; i >= 0; i--)
  { if (path[i] == '/') break; }
  for (j = 0; j < l-i; j++) EndPath[j] = path[j+i+1];
  EndPath[l-i] = '\0';
}

//=============================================================================
// Ecrit seulement les chemins specifies de l'arbre
//=============================================================================
E_Int K_IO::GenIO::hdfcgnsWritePaths(char* file, PyObject* treeList,
                                     PyObject* paths, PyObject* links, 
                                     E_Int maxDepth, E_Int mode)
{
  if (PyList_Check(paths) == false)
  {
    PyErr_SetString(PyExc_TypeError,
                    "hdfcgnswrite: paths must be a list of strings.");
    return 1;
  }
  E_Int size = PyList_Size(paths);

  /* Ouverture du fichier pour l'ecriture */
  hid_t fapl, fid;
  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);
  fid = H5Fopen(file, H5F_ACC_RDWR, fapl);
  H5Pclose(fapl);

  if (fid < 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "hdfcgnswrite: can not open file.");
    return 1;
  }

  GenIOHdf HDF;
  if (maxDepth >= 0) HDF._maxDepth = maxDepth;
  else HDF._maxDepth = 1e6;

  for (E_Int i = 0; i < size; i++)
  {
    char* path; PyObject* l;
    l = PyList_GetItem(paths, i);
    if (PyString_Check(l)) path = PyString_AsString(l);
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(l)) path = (char*)PyUnicode_AsUTF8(l);
#endif
    else
    {
      PyErr_SetString(PyExc_TypeError, "hdfwrite: paths must be strings.");
      return 1;
    }
    
    PyObject* node = PyList_GetItem(treeList, i);
    hid_t gidp = H5Gopen(fid, path, H5P_DEFAULT);
    
    if (gidp < 0)
      printf("Warning: hdfcgnswrite: cannot write this path %s.\n", path);
    else
    {
      if (mode == 0) // append mixed
      {
        // Find node[0]
        char* nodeName; PyObject* l;
        l = PyList_GetItem(node,0);
        if (PyString_Check(l)) nodeName = PyString_AsString(l);
#if PY_VERSION_HEX >= 0x03000000
        else if (PyUnicode_Check(l)) nodeName = (char*)PyUnicode_AsUTF8(l);
#endif
        else nodeName = NULL;
        char spath[512];
        strcpy(spath, path);
        strcat(spath, "/");
        strcat(spath, nodeName);

        herr_t status = H5Lexists(gidp, nodeName, H5P_DEFAULT);
        if (status > 0)
        {
          H5Ldelete(fid, spath, H5P_DEFAULT);
        }
        
        HDF._fatherStack.push_front(gidp);
        HDF.dumpOne(node, 0, links);
        HDF._fatherStack.pop_front();
        H5Gclose(gidp);
      }
      else if (mode == 2) // append pur
      {
        HDF._fatherStack.push_front(gidp);
        HDF.dumpOne(node, 0, links);
        HDF._fatherStack.pop_front();
        H5Gclose(gidp);
      }
      else if (mode == 1) // replace
      {
        H5Gclose(gidp);
        vector<char*> pelts;
        getABFromPath(path, pelts);
        char* pp = new char [strlen(path)];
        strcpy(pp, "/");
        for (size_t i = 0; i < pelts.size()-1; i++) { strcat(pp, pelts[i]); strcat(pp, "/"); }
        for (size_t i = 0; i < pelts.size(); i++) delete [] pelts[i];
        gidp = H5Gopen(fid, pp, H5P_DEFAULT); // parent node
        delete [] pp;
        
        H5Ldelete(fid, path, H5P_DEFAULT); // gain space for hdf > 1.10
        HDF._fatherStack.push_front(gidp);
        HDF.dumpOne(node, 0, links);
        HDF._fatherStack.pop_front();
      }
    }
  }

  //printf("open=%d\n", H5Fget_obj_count(fid, H5F_OBJ_ALL));
  H5Fclose(fid);
  
  return 0;
}

//=============================================================================
PyObject* K_IO::GenIOHdf::dumpOne(PyObject* tree, int depth, PyObject* links)
{
  // ecrit le noeud courant
  hid_t node = _fatherStack.front();

  if (depth > _maxDepth) return tree;

  // Check if node is a link, set _currentPath
  if (_stringStack.size() > 0)
  {
    PyObject* pname = PyList_GetItem(tree, 0);
    char* name;
    if (PyString_Check(pname)) name = PyString_AsString(pname);
    #if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(pname)) name = (char*)PyUnicode_AsUTF8(pname);
#endif
    else name = NULL;
    std::string path = _stringStack.front();
    _currentPath = path+"/"+string(name);
    //printf("write node: writing %s\n", _currentPath.c_str());
    if (links != NULL && checkPathInLinks(_currentPath.c_str(), links) == 1) return tree;
  }
  
  // Write node
  hid_t child = writeNode(node, tree);
  
  std::string parentPath = _currentPath;

  // Dump les enfants
  if (PyList_Check(tree) == true && PyList_Size(tree) > 3)
  {
    PyObject* children = PyList_GetItem(tree, 2);
    if (PyList_Check(children) == true)
    {
      int nChildren = PyList_Size(children);
      for (E_Int i = 0; i < nChildren; i++)
      {
        _stringStack.push_front(parentPath);
        _fatherStack.push_front(child);
        dumpOne(PyList_GetItem(children, i), depth+1, links);
        _fatherStack.pop_front();
        _stringStack.pop_front();
      }
    }
  }
  H5Gclose(child);
  return tree;
}
//=============================================================================
hid_t K_IO::GenIOHdf::writeNode(hid_t node, PyObject* tree)
{
  IMPORTNUMPY;
  char s1[CGNSMAXLABEL+1];
  char s2[CGNSMAXLABEL+1];
  PyObject* pname = PyList_GetItem(tree, 0);
  PyObject* plabel = PyList_GetItem(tree, 3);
  char* name; char* label;
  if (PyString_Check(pname)) name = PyString_AsString(pname);
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(pname)) name = (char*)PyUnicode_AsUTF8(pname);
#endif
  else name = NULL;
  if (PyString_Check(plabel)) label = PyString_AsString(plabel);
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(plabel)) label = (char*)PyUnicode_AsUTF8(plabel);
#endif
  else label = NULL;
  strcpy(s1, name); strcpy(s2, label);

  // Creation du noeud, prop. est requis par cgnslib
  hid_t child;
  herr_t status = H5Lexists(node, s1, H5P_DEFAULT);
  if (status == true) // group already exist, strange
  {
    //printf("duplicated node\n");
    child = H5Gopen2(node, s1, H5P_DEFAULT);
    return child;
  }
  else child = H5Gcreate2(node, s1, H5P_DEFAULT, _group, H5P_DEFAULT);

  HDF_Add_Attribute_As_String(child, L3S_NAME, s1);
  HDF_Add_Attribute_As_String(child, L3S_LABEL, s2);
  HDF_Add_Attribute_As_Integer(child, L3S_FLAGS, 1);

  // Ecriture de la valeur
  PyObject* v = PyList_GetItem(tree, 1);

  if (v == Py_None)
  {
    // direct dans l'attribut
    HDF_Add_Attribute_As_String(child, L3S_DTYPE, L3T_MT);
  }
  else if (PyString_Check(v))
  {
    setArrayC1(child, PyString_AsString(v));
    HDF_Add_Attribute_As_String(child, L3S_DTYPE, L3T_C1);
  }
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(v))
  {
    setArrayC1(child, (char*)PyUnicode_AsUTF8(v));
    HDF_Add_Attribute_As_String(child, L3S_DTYPE, L3T_C1);
  }
#endif
  else if (PyInt_Check(v))
  {
    setSingleI4(child, PyInt_AsLong(v));
  }
  else if (PyFloat_Check(v))
  {
    if (strcmp(name, "CGNSLibraryVersion") == 0)
      setSingleR4(child, PyFloat_AsDouble(v));
    else
      setSingleR8(child, PyFloat_AsDouble(v));
  }
  else if (PyArray_Check(v))
  {
    PyArrayObject* ar = (PyArrayObject*)v;
    int dim = PyArray_NDIM(ar);
    hsize_t* dims = new hsize_t [dim];
    E_Int typeNum = PyArray_TYPE(ar);
    E_Int elSize = PyArray_ITEMSIZE(ar);
    for (E_Int n = 0; n < dim; n++)
    {
      // Fortran indexing
      dims[n] = PyArray_DIMS(ar)[dim-n-1];
    }

    if (dim == 1 && dims[0] == 1) // valeur simple
    {
      if (typeNum == NPY_DOUBLE)
      {
        if (strcmp(name, "CGNSLibraryVersion") == 0)
        {
          double* ptr = (double*)PyArray_DATA(ar);
          setSingleR4(child, (float)ptr[0]);
        }
        else
        {
          double* ptr = (double*)PyArray_DATA(ar);
          setSingleR8(child, ptr[0]);
        }
      }
      else if (typeNum == NPY_FLOAT)
      {
        float* ptr = (float*)PyArray_DATA(ar);
        setSingleR4(child, ptr[0]);
      }
      else if (typeNum == NPY_INT || typeNum == NPY_LONG || typeNum == NPY_INT64)
      {
        if (elSize == 4)
        {
          int* ptr = (int*)PyArray_DATA(ar);
          setSingleI4(child, ptr[0]);
        }
        else
        {
          E_LONG* ptr = (E_LONG*)PyArray_DATA(ar);
          setSingleI8(child, ptr[0]);
        }
      }
      else if (typeNum == NPY_STRING ||
               typeNum == NPY_BYTE ||
               //typeNum == NPY_SBYTE ||
               typeNum == NPY_UBYTE)
      {
        E_Int diml = PyArray_DIMS(ar)[0];
        char* buf = new char [diml+1];
        char* pt = (char*)PyArray_DATA(ar);
        for (E_Int i = 0; i < diml; i++) buf[i] = pt[i]; 
        //strncpy(buf, (char*)PyArray_DATA(ar), diml); // pb align
        buf[diml] = '\0';
        setArrayC1(child, buf);
        HDF_Add_Attribute_As_String(child, L3S_DTYPE, L3T_C1);
        delete [] buf;
      }
    }
    else // tableau
    {
      if (typeNum == NPY_DOUBLE)
      {
#if (FORCEPERIODICR4 == 1)
        // patch pour l'ancienne norme CGNS
        if (strcmp(name, "RotationCenter") == 0 ||
            strcmp(name, "RotationAngle") == 0 ||
            strcmp(name, "RotationRateVector") == 0 ||
            strcmp(name, "Translation") == 0)
        {
          E_Int s = PyArray_SIZE(ar);
          float* buf = new float [s];
          double* ptr = (double*)PyArray_DATA(ar);
          for (E_Int i = 0; i < s; i++) buf[i] = ptr[i];
          setArrayR4(child, buf, dim, dims);
          delete [] buf;
        }
        else
        {
          setArrayR8(child, (double*)PyArray_DATA(ar), dim, dims);
        }
#else
        setArrayR8(child, (double*)PyArray_DATA(ar), dim, dims);
#endif
      }
      else if (typeNum == NPY_INT || typeNum == NPY_INT64 || typeNum == NPY_LONG)
      {
       if (elSize == 4)
       {
         setArrayI4(child, (int*)PyArray_DATA(ar), dim, dims);
       }
       else if (elSize == 1)
       {
         setArrayI1(child, (char*)PyArray_DATA(ar), dim, dims);
       }
       else if (elSize == 8)
       {
         if (strcmp(label, "CGNSBase_t") == 0 ||
             strcmp(label, "Elements_t") == 0) // to comply with cgns 4 norm
         {
          // convert to i4
          E_Int s = PyArray_SIZE(ar);
          int* buf = new int [s];
          E_Int* ptr = (E_Int*)PyArray_DATA(ar);
          for (E_Int i = 0; i < s; i++) buf[i] = (int)ptr[i];
          setArrayI4(child, buf, dim, dims); // raw
          delete [] buf;
         }
         else
         {
          setArrayI8(child, (E_LONG*)PyArray_DATA(ar), dim, dims); // driver best/raw
         }
       }
      }
      else if (typeNum == NPY_BYTE)
      {
         setArrayI1(child, (char*)PyArray_DATA(ar), dim, dims);
      }
      else if (typeNum == NPY_STRING ||
               //typeNum == NPY_BYTE ||
               //typeNum == NPY_SBYTE ||
               typeNum == NPY_UBYTE)
      {
        if (dim == 1)
        {
          E_Int diml = PyArray_DIMS(ar)[0];
          char* buf = new char [diml+1];
          //strncpy(buf, (char*)PyArray_DATA(ar), diml); // pb align
          char* pt = (char*)PyArray_DATA(ar);
          for (E_Int i = 0; i < diml; i++) buf[i] = pt[i];
          buf[diml] = '\0';
          setArrayC1(child, buf);
          HDF_Add_Attribute_As_String(child, L3S_DTYPE, L3T_C1);
          delete [] buf;
        }
        else
        { 
          setArrayC1(child, (char*)PyArray_DATA(ar), dim, dims);
        }
      }
      else if (typeNum == NPY_FLOAT)
      {
        //E_Int s = PyArray_Size(v);
        //double* buf = new double [s];
        //for (int i = 0; i < s; i++) buf[i] = (double)PyArray_DATA(ar)[i];
        //setArrayR8(child, buf, dim, dims);
        //delete [] buf;
        setArrayR4(child, (float*)PyArray_DATA(ar), dim, dims);
      }
    }
    delete [] dims;
  }
  return child;
}

//=============================================================================
// modify name, type and value from pytree node
hid_t K_IO::GenIOHdf::modifyNode(hid_t node, PyObject* tree)
{
  IMPORTNUMPY;
  char s1[CGNSMAXLABEL+1];
  char s2[CGNSMAXLABEL+1];
  PyObject* pname = PyList_GetItem(tree, 0);
  
  char* name = NULL;
  if (PyString_Check(pname)) name = PyString_AsString(pname);
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(pname)) name = (char*)PyUnicode_AsUTF8(pname);
#endif
  PyObject* plabel = PyList_GetItem(tree, 3);
  char* label = NULL;
  if (PyString_Check(plabel)) label = PyString_AsString(plabel);
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(plabel)) label = (char*)PyUnicode_AsUTF8(plabel);
#endif
  strcpy(s1, name); strcpy(s2, label);

  HDF_Set_Attribute_As_String(node, L3S_NAME, s1);
  HDF_Set_Attribute_As_String(node, L3S_LABEL, s2);
  //HDF_Set_Attribute_As_Integer(node, L3S_FLAGS, 1);
  
  // Missing replace value!!
  hid_t child = H5Gopen2(node, s1, H5P_DEFAULT);

  return child;
}
//=============================================================================
hid_t K_IO::GenIOHdf::setSingleR4(hid_t node, float data)
{
  hsize_t dim; hsize_t dims[1];
  dim = 1; dims[0] = 1;

  // data type
  hid_t tid = H5Tcopy(H5T_NATIVE_FLOAT); H5Tset_precision(tid, 32);
  // Create dataspace
  hid_t sid = H5Screate_simple(dim, dims, NULL);
  // Create dataset
  hid_t did = H5Dcreate(node, L3S_DATA, tid, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hid_t mid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  H5Dwrite(did, mid, H5S_ALL, sid, H5P_DEFAULT, &data);
  H5Tclose(tid); H5Dclose(did); H5Sclose(sid); H5Tclose(mid);
  HDF_Add_Attribute_As_String(node, L3S_DTYPE, L3T_R4);
  return node;
}

//=============================================================================
hid_t K_IO::GenIOHdf::setSingleR8(hid_t node, double data)
{
  hsize_t dim; hsize_t dims[1];
  dim = 1; dims[0] = 1;

  // data type
  hid_t tid = H5Tcopy(H5T_NATIVE_DOUBLE); H5Tset_precision(tid, 64);
  // Create dataspace
  hid_t sid = H5Screate_simple(dim, dims, NULL);
  // Create dataset
  hid_t did = H5Dcreate(node, L3S_DATA, tid, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hid_t mid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  H5Dwrite(did, mid, H5S_ALL, sid, H5P_DEFAULT, &data);
  H5Tclose(tid); H5Dclose(did); H5Sclose(sid); H5Tclose(mid);
  HDF_Add_Attribute_As_String(node, L3S_DTYPE, L3T_R8);
  return node;
}

//=============================================================================
hid_t K_IO::GenIOHdf::setSingleI4(hid_t node, int data)
{
  hsize_t dim; hsize_t dims[1];
  dim = 1; dims[0] = 1;

  // data type
  hid_t tid = H5Tcopy(H5T_NATIVE_INT); H5Tset_precision(tid, 32);
  // Create dataspace
  hid_t sid = H5Screate_simple(dim, dims, NULL);
  // Create dataset
  hid_t did = H5Dcreate(node, L3S_DATA, tid, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hid_t mid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  H5Dwrite(did, mid, H5S_ALL, sid, H5P_DEFAULT, &data);
  H5Tclose(tid); H5Dclose(did); H5Sclose(sid); H5Tclose(mid);
  HDF_Add_Attribute_As_String(node, L3S_DTYPE, L3T_I4);
  return node;
}

//=============================================================================
hid_t K_IO::GenIOHdf::setSingleI8(hid_t node, E_LONG data)
{
  hsize_t dim; hsize_t dims[1];
  dim = 1; dims[0] = 1;

  // data type
  hid_t tid = H5Tcopy(H5T_NATIVE_LONG); H5Tset_precision(tid, 64);
  // Create dataspace
  hid_t sid = H5Screate_simple(dim, dims, NULL);
  // Create dataset
  hid_t did = H5Dcreate(node, L3S_DATA, tid, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  
  hid_t mid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  H5Dwrite(did, mid, H5S_ALL, sid, H5P_DEFAULT, &data);
  H5Tclose(tid); H5Dclose(did); H5Sclose(sid); H5Tclose(mid);
  HDF_Add_Attribute_As_String(node, L3S_DTYPE, L3T_I8);
  return node;
}

//=============================================================================
hid_t K_IO::GenIOHdf::setArrayC1(hid_t node, char* data, char* label)
{
  hsize_t dim; hsize_t dims[1];

  dim = 1;
  if (strcmp(label, L3S_VERSION) == 0) dims[0] = 33;
  else dims[0] = strlen(data);

  // data type
  hid_t tid = H5Tcopy(H5T_NATIVE_CHAR);
  // Create dataspace
  hid_t sid = H5Screate_simple(dim, dims, NULL);
  // Create dataset
  hid_t did = H5Dcreate(node, label, tid, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hid_t mid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  H5Dwrite(did, mid, H5S_ALL, sid, H5P_DEFAULT, data);
  H5Tclose(tid); H5Dclose(did); H5Sclose(sid); H5Tclose(mid);
  return node;
}

//=============================================================================
hid_t K_IO::GenIOHdf::setArrayC1(hid_t node, char* data, int idim, hsize_t* idims)
{
  hsize_t dim; hsize_t* dims;
  dim = idim; dims = (hsize_t*)malloc(sizeof(hsize_t)*dim);
  for (E_Int i = 0; i < idim; i++) dims[i] = idims[i];

  // data type
  hid_t tid = H5Tcopy(H5T_NATIVE_CHAR);
  // Create dataspace
  hid_t sid = H5Screate_simple(dim, dims, NULL);
  // Create dataset
  hid_t did = H5Dcreate(node, L3S_DATA, tid, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hid_t mid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  H5Dwrite(did, mid, H5S_ALL, sid, H5P_DEFAULT, data);
  H5Tclose(tid); H5Dclose(did); H5Sclose(sid); H5Tclose(mid);
  HDF_Add_Attribute_As_String(node, L3S_DTYPE, L3T_C1);
  free(dims);
  return node;
}

//=============================================================================
hid_t K_IO::GenIOHdf::setArrayI1(hid_t node, char* data, int idim, hsize_t* idims)
{
  hsize_t dim; hsize_t* dims;
  dim = idim; dims = (hsize_t*)malloc(sizeof(hsize_t)*dim);
  for (E_Int i = 0; i < idim; i++) dims[i] = idims[i];
  // data type
  hid_t tid = H5Tcopy(H5T_NATIVE_INT8); H5Tset_precision(tid, 8);
  // Create dataspace
  hid_t sid = H5Screate_simple(dim, dims, NULL);
  // Create dataset
  hid_t did = H5Dcreate(node, L3S_DATA, tid, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hid_t mid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  H5Dwrite(did, mid, H5S_ALL, sid, H5P_DEFAULT, data);
  H5Tclose(tid); H5Dclose(did); H5Sclose(sid); H5Tclose(mid);
  HDF_Add_Attribute_As_String(node, L3S_DTYPE, L3T_I1);
  free(dims);
  return node;
}

//=============================================================================
hid_t K_IO::GenIOHdf::setArrayI4(hid_t node, int* data, int idim, hsize_t* idims)
{
  hsize_t dim; hsize_t* dims;
  dim = idim; dims = (hsize_t*)malloc(sizeof(hsize_t)*dim);
  for (E_Int i = 0; i < idim; i++) dims[i] = idims[i];

  // data type
  hid_t tid = H5Tcopy(H5T_NATIVE_INT); H5Tset_precision(tid, 32);
  // Create dataspace
  hid_t sid = H5Screate_simple(dim, dims, NULL);
  // Create dataset
  hid_t did = H5Dcreate(node, L3S_DATA, tid, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hid_t mid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  H5Dwrite(did, mid, H5S_ALL, sid, H5P_DEFAULT, data);
  H5Tclose(tid); H5Dclose(did); H5Sclose(sid); H5Tclose(mid);
  HDF_Add_Attribute_As_String(node, L3S_DTYPE, L3T_I4);
  free(dims);
  return node;
}

//=============================================================================
// write a i8 array raw or best.
hid_t K_IO::GenIOHdf::setArrayI8(hid_t node, E_LONG* data, int idim, hsize_t* idims)
{
  if (_writeIntMode == 0) return setArrayI8Raw(node, data, idim, idims);
  else return setArrayI8B(node, data, idim, idims);
}

//=============================================================================
// export as i8 with i8 in memory
//=============================================================================
hid_t K_IO::GenIOHdf::setArrayI8Raw(hid_t node, E_LONG* data, int idim, hsize_t* idims)
{
  hsize_t dim; hsize_t* dims;
  dim = idim; dims = (hsize_t*)malloc(sizeof(hsize_t)*dim);
  for (E_Int i = 0; i < idim; i++) dims[i] = idims[i];

  // data type
  hid_t tid = H5Tcopy(H5T_NATIVE_LONG); H5Tset_precision(tid, 64);
  // Create dataspace
  hid_t sid = H5Screate_simple(dim, dims, NULL);
  // Create dataset
  hid_t did = H5Dcreate(node, L3S_DATA, tid, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hid_t mid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  H5Dwrite(did, mid, H5S_ALL, sid, H5P_DEFAULT, data);
  H5Tclose(tid); H5Dclose(did); H5Sclose(sid); H5Tclose(mid);
  HDF_Add_Attribute_As_String(node, L3S_DTYPE, L3T_I8);
  free(dims);
  return node;
}

//=============================================================================
// Best : if possible exports without loss to i4 else export to i8
//=============================================================================
hid_t K_IO::GenIOHdf::setArrayI8B(hid_t node, E_LONG* data, int idim, hsize_t* idims)
{
  hsize_t dim; hsize_t* dims;
  dim = idim; dims = (hsize_t*)malloc(sizeof(hsize_t)*dim);
  int64_t size = 1;
  for (E_Int i = 0; i < idim; i++) { dims[i] = idims[i]; size = size* dims[i]; }

  // check int size in array
  int64_t maxInt = 2;
  maxInt = (maxInt<<31)-2;
  E_Boolean exceed = false;
  for (int64_t i = 0; i < size; i++)
  {
    if (data[i] > maxInt) { exceed = true; break; }
  }

  if (exceed)
  { // write in i8
    hid_t tid = H5Tcopy(H5T_NATIVE_LONG); H5Tset_precision(tid, 64);
    hid_t sid = H5Screate_simple(dim, dims, NULL);
    hid_t did = H5Dcreate(node, L3S_DATA, tid, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t mid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
    H5Dwrite(did, mid, H5S_ALL, sid, H5P_DEFAULT, data);
    H5Tclose(tid); H5Dclose(did); H5Sclose(sid); H5Tclose(mid);
    HDF_Add_Attribute_As_String(node, L3S_DTYPE, L3T_I8);
  }
  else
  { // write in i4
    hid_t tid = H5Tcopy(H5T_NATIVE_INT); H5Tset_precision(tid, 32);
    hid_t sid = H5Screate_simple(dim, dims, NULL);
    hid_t did = H5Dcreate(node, L3S_DATA, tid, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t mid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
    H5Dwrite(did, mid, H5S_ALL, sid, H5P_DEFAULT, data);
    H5Tclose(tid); H5Dclose(did); H5Sclose(sid); H5Tclose(mid);
    HDF_Add_Attribute_As_String(node, L3S_DTYPE, L3T_I4);
  }
  free(dims);
  return node;
}

//=============================================================================
// export as R4 with R4 in memory
//=============================================================================
hid_t K_IO::GenIOHdf::setArrayR4(hid_t node, float* data,
                                 int idim, hsize_t* idims)
{
  hsize_t dim; hsize_t* dims;
  dim = idim; dims = (hsize_t*)malloc(sizeof(hsize_t)*dim);
  for (E_Int i = 0; i < idim; i++) dims[i] = idims[i];

  // data type
  hid_t tid = H5Tcopy(H5T_NATIVE_FLOAT); H5Tset_precision(tid, 32);
  // Create dataspace
  hid_t sid = H5Screate_simple(dim, dims, NULL);
  // Create dataset
  hid_t did = H5Dcreate(node, L3S_DATA, tid, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hid_t mid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  H5Dwrite(did, mid, H5S_ALL, sid, H5P_DEFAULT, data);
  H5Tclose(tid); H5Dclose(did); H5Sclose(sid); H5Tclose(mid);
  HDF_Add_Attribute_As_String(node, L3S_DTYPE, L3T_R4);
  free(dims);
  return node;
}

//=============================================================================
hid_t K_IO::GenIOHdf::setArrayR8(hid_t node, double* data,
                                 int idim, hsize_t* idims)
{
  if (_writeRealMode == 0) return setArrayR8Raw(node, data, idim, idims);
  else return setArrayR82R4(node, data, idim, idims);
}

//=============================================================================
// export as R8 with R8 in memory
//=============================================================================
hid_t K_IO::GenIOHdf::setArrayR8Raw(hid_t node, double* data,
                                    int idim, hsize_t* idims)
{
  hsize_t dim; hsize_t* dims;
  dim = idim; dims = (hsize_t*)malloc(sizeof(hsize_t)*dim);
  for (E_Int i = 0; i < idim; i++) dims[i] = idims[i];

  // data type
  hid_t tid = H5Tcopy(H5T_NATIVE_DOUBLE); H5Tset_precision(tid, 64);
  // Create dataspace
  hid_t sid = H5Screate_simple(dim, dims, NULL);
  // Create dataset
  hid_t did = H5Dcreate(node, L3S_DATA, tid, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hid_t mid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  H5Dwrite(did, mid, H5S_ALL, sid, H5P_DEFAULT, data);
  H5Tclose(tid); H5Dclose(did); H5Sclose(sid); H5Tclose(mid);
  HDF_Add_Attribute_As_String(node, L3S_DTYPE, L3T_R8);
  free(dims);
  return node;
}

//=============================================================================
// export as R4 with R8 in memory
//=============================================================================
hid_t K_IO::GenIOHdf::setArrayR82R4(hid_t node, double* data,
                                    int idim, hsize_t* idims)
{
  hsize_t dim; hsize_t* dims;
  dim = idim; dims = (hsize_t*)malloc(sizeof(hsize_t)*dim);
  for (E_Int i = 0; i < idim; i++) dims[i] = idims[i];

  // data type
  hid_t tid = H5Tcopy(H5T_NATIVE_FLOAT); H5Tset_precision(tid, 32);
  // Create dataspace
  hid_t sid = H5Screate_simple(dim, dims, NULL);
  // Create dataset
  hid_t did = H5Dcreate(node, L3S_DATA, tid, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hid_t mid = H5Tcopy(H5T_NATIVE_DOUBLE); H5Tset_precision(tid, 64);
  H5Dwrite(did, mid, H5S_ALL, sid, H5P_DEFAULT, data);
  H5Tclose(tid); H5Dclose(did); H5Sclose(sid); H5Tclose(mid);
  HDF_Add_Attribute_As_String(node, L3S_DTYPE, L3T_R4);
  free(dims);
  return node;
}

//=============================================================================
// Delete les chemins specifies du fichier (et tous les sous-noeuds attaches)
// Attention ne libere pas de place dans le fichier
//=============================================================================
E_Int K_IO::GenIO::hdfcgnsDeletePaths(char* file,
                                      PyObject* paths)
{
  if (PyList_Check(paths) == false)
  {
    PyErr_SetString(PyExc_TypeError,
                    "hdfdelete: paths must be a list of strings.");
    return 1;
  }
  E_Int size = PyList_Size(paths);
  for (E_Int i = 0; i < size; i++)
  {
    PyObject* o = PyList_GetItem(paths, i);
    E_Int isString = 0;
    if (PyString_Check(o)) isString = 1;
#if PY_VERSION_HEX >= 0x03000000
    if (PyUnicode_Check(o)) isString = 1;
#endif
    if (isString == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "hdfdelete: paths must be a list of strings.");
      return 1;
    }
  }

  /* Ouverture du fichier pour l'ecriture */
  hid_t fapl, fid;
  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);
  fid = H5Fopen(file, H5F_ACC_RDWR, fapl);
  H5Pclose(fapl);

  if (fid < 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "hdfdelete: can not open file.");
    return 1;
  }

  GenIOHdf HDF;
  for (E_Int i = 0; i < size; i++)
  {
    PyObject* o = PyList_GetItem(paths, i);
    char* path = NULL;
    if (PyString_Check(o)) path = PyString_AsString(o);
#if PY_VERSION_HEX >= 0x03000000
      else if (PyUnicode_Check(o)) path = (char*)PyUnicode_AsUTF8(o);
#endif
    hid_t gidp = H5Gopen(fid, path, H5P_DEFAULT);
    
    if (gidp < 0)
      printf("Warning: hdfcgnsdelete: cannot write this path %s.\n", path);
    else
    {
      H5Gclose(gidp);
      H5Ldelete(fid, path, H5P_DEFAULT); // gain space for hdf > 1.10
    }
  }

  H5Fclose(fid);
  return 0;
}

#include "GenIO_hdfcgns_partialMPI.cpp"
