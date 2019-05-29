/*
    Copyright 2013-2019 Onera.

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

// Lecture partielle des noeuds decrits dans Filter (chemins et filtre)
PyObject* K_IO::GenIO::hdfcgnsReadFromPathsPartial(char* file,
                                                   PyObject* Filter,
                                                   void* comm)
{
  hid_t fapl, fid, ret;

  PyObject *key, *DataSpaceDIM;
  Py_ssize_t pos = 0;

  // > Out
  PyObject* NodesList;
  
  /* Check */
  if (PyDict_Check(Filter) == false)
  {
    PyErr_SetString(PyExc_TypeError, "hdfread: Filter must be a dict of paths.");
    return NULL;
  }

  /* Begin */
  // int OutputType = 0;  // List
  int OutputType = 1;  // Dict
  if (OutputType == 0){NodesList = PyList_New(0);}
  else                {NodesList = PyDict_New( );}

  /* Open file */
  fapl = H5Pcreate(H5P_FILE_ACCESS);

  // H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);
  /* Group Level */
  GenIOHdf HDF;
  HDF._skeleton = 0;

#if defined(_MPI) && defined(H5_HAVE_PARALLEL)
  HDF._ismpi = 1;
#else
  HDF._ismpi = 0;
#endif

#if defined(_MPI) && defined(H5_HAVE_PARALLEL)
  if (HDF._ismpi == 1){
     /* Mpi context */
     MPI_Comm* comm2 = (MPI_Comm*)comm;
     MPI_Info info   = MPI_INFO_NULL;
     ret             = H5Pset_fapl_mpio(fapl, *comm2, info);
   }
#endif  

  /* File Level */
  fid = H5Fopen(file, H5F_ACC_RDONLY, fapl);

  if (fid < 0)
  {
    printf("Warning: hdfcgnsReadFromPathsPartial: can not open file %s.\n", file);
    return Py_None;
  }
  
  H5Pclose(fapl);
    
#if defined(_MPI) && defined(H5_HAVE_PARALLEL)
  int nRank, myRank;
  MPI_Comm* comm2 = (MPI_Comm*)comm;
  MPI_Comm_size(*comm2, &nRank);
  MPI_Comm_rank(*comm2, &myRank);
  // printf("[%d] - open Avant =%d\n",myRank, H5Fget_obj_count(fid, H5F_OBJ_ALL));
#endif
  
  PyObject* node;

  while (PyDict_Next(Filter, &pos, &key, &DataSpaceDIM))
  {
    // Multiple path or Not ?
    E_Boolean isKeyString = false;
    if (PyString_Check(key)) isKeyString = true;
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(key)) isKeyString = true; 
#endif
    if (isKeyString)
    {
      E_Int FilterSize = PyList_Size(DataSpaceDIM);
      /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */
      /* Check if key is String values */
      
      /* Check if key is String values */
      if (PyList_Check(DataSpaceDIM) == false)
      {
        PyErr_SetString(PyExc_TypeError, "hdfread: DataSpaceDIM must be a list of numbers.");
        return NULL;
      }
      if (FilterSize < 9)  /** Dans le cas particulier Contigous with only one path **/
      {
        printf("FilterSize: %d \n", FilterSize);
        PyErr_SetString(PyExc_TypeError, "hdfread: FilterSize must be a list of 9 numbers.");
        return NULL;
      }
      /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */

      /* Get path */
      char* path = NULL;
      if (PyString_Check(key)) path = PyString_AsString(key);
#if PY_VERSION_HEX >= 0x03000000
      else if (PyUnicode_Check(key)) path = PyBytes_AsString(PyUnicode_AsUTF8String(key)); 
#endif
      // printf("path 1 ...  %s\n", path);

      HDF._path = path;

      /* Open group in HDF corresponding to path */
      hid_t gid = HDF.openGroupWithLinks(fid, path);  

      /* Fill data space */
      HDF.fillDataSpaceWithFilter(DataSpaceDIM);
      node = HDF.createNodePartial(gid);

      /* Close */
      H5Gclose(gid);

      if (OutputType == 0) {PyList_Append(NodesList, node); Py_DECREF(node);}
      else
      {
        PyDict_SetItemString(NodesList, HDF._path, PyList_GetItem(node,1));
        Py_DECREF(node);
      }
    }
    else if (PyTuple_Check(key) == true && PyList_Check(DataSpaceDIM) == true)  /** Contigous or Interlaced of field of same size **/
    {
      E_Int FilterSize = PyList_Size(DataSpaceDIM);
      /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */
      // if (PyList_Check(DataSpaceDIM) == false)
      // {
      //   PyErr_SetString(PyExc_TypeError, "hdfread: DataSpaceDIM must be a list of numbers.");
      //   return NULL;
      // }
      if (FilterSize != 10)
      {
        PyErr_SetString(PyExc_TypeError, "hdfread: FilterSize must be a list of 9 numbers.");
        return NULL;
      }
      /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */

      /** Number of Field to load **/
      int nField = PyTuple_Size(key);

      PyObject* lpath;
      PyObject* data = NULL;

      /** Loop over List of Path to load **/
      for (int iField = 0; iField < nField; iField++)
      {
        lpath      = PyTuple_GetItem(key, iField);
        char* path = NULL;
        if (PyString_Check(lpath)) path = PyString_AsString(lpath);
#if PY_VERSION_HEX >= 0x03000000
        else if (PyUnicode_Check(lpath)) path = PyBytes_AsString(PyUnicode_AsUTF8String(lpath));
#endif

        /** Verbose **/
        // printf("path contigous ...  %s\n", path);

        /** Store the first path **/
        if (iField == 0) {HDF._path = path;}

        /* Open group in HDF corresponding to path */              
        hid_t gid = HDF.openGroupWithLinks(fid, path);  
        
        /* Fill data space */
        HDF.fillDataSpaceWithFilter(DataSpaceDIM);
        data = HDF.createNodePartialContigous(gid, iField, nField, data);

        /** Verbose **/
        // double* ptr = (double *) PyArray_DATA((PyArrayObject * ) data);
        // for(int n=0; n<24; n++)
        // {
        //   printf("data[%d] : %f \n", n , ptr[n]);
        // }

        /* Close */
        H5Gclose(gid);

      }

      /** Build PyTree Node **/
      if (OutputType == 0)
      {
        PyObject* node = Py_BuildValue("[sOOs]", HDF._path, data, Py_None, HDF._type);
        PyList_Append(NodesList, node);
        Py_DECREF(node);
      }
      else
      {
        PyDict_SetItemString(NodesList, HDF._path, data);
      }

      /** Free ref **/
      Py_DECREF(data);
      // Py_DECREF(lpath); core ...
    }
    /* ----------------------------------------------------------------------------- */
    else if (PyTuple_Check(key) == true && PyTuple_Check(DataSpaceDIM) == true) /** Contigous or Interlaced of field of different size **/
    {
      /** DataSpaceGlob contains the global size of the desired array **/
      /** In this case DataSpaceDim is a list of DataSpace ...        **/
      /** Each DataSpace is assign to a current path                  **/

      PyObject* lpath;
      PyObject* LocalDataSetDim=NULL;
      PyObject* data = NULL;

      int nField = PyTuple_Size(key);
      int oField = 0;

      /* On parse tous les chemins du tuple de chemins */
      assert(nField > 0 && "Number of field must be greater than zero");
      for (int iField = 0; iField < nField; iField++)
      {
        /** Get current path **/
        lpath      = PyTuple_GetItem(key, iField);
        char* path = NULL;
        if (PyString_Check(lpath)) path = PyString_AsString(lpath);
#if PY_VERSION_HEX >= 0x03000000
        else if (PyUnicode_Check(lpath)) path = PyBytes_AsString(PyUnicode_AsUTF8String(lpath)); 
#endif
        /** Store the first path **/
        if (iField == 0){HDF._path = path;}

        /* Open group in HDF corresponding to path */
        hid_t gid = HDF.openGroupWithLinks(fid, path);  

        /** Get the current dataSpace to Load ( in file sens ) **/
        LocalDataSetDim = PyTuple_GetItem(DataSpaceDIM, iField);

        /** LocalDataSetDim is a list normaly **/
        assert(PyList_Size(LocalDataSetDim) == 10 && "Wrong size for LocalDataSetDim");
        //assert(lsize == 10);

        /** Verbose **/
        // printf("path ...  %s -> %d / %d \n", path, lsize, oField);

        /** Fill data space **/
        HDF.fillDataSpaceWithFilter(LocalDataSetDim);

        /** Read the data **/
        data = HDF.createNodePartialContigous(gid, iField, oField, data);        
        
        /** Verbose **/
        // int* ptr = (int *) PyArray_DATA((PyArrayObject * ) data);
        // for(int n=0; n<24; n++)
        // {
        //   printf("data[%d] : %d \n", n , ptr[n]);
        // }

        /* Close */
        H5Gclose(gid);
      }

      /** Build PyTree Node **/
      if (OutputType == 0)
      {
        PyObject* node = Py_BuildValue("[sOOs]", HDF._path, data, Py_None, HDF._type);
        PyList_Append(NodesList, node);
        Py_DECREF(node);
      }
      else
      {
        PyDict_SetItemString(NodesList, HDF._path, data);
      }

      /** Free ref **/
      Py_DECREF(data);
      if (LocalDataSetDim != NULL) Py_DECREF(LocalDataSetDim);
    }
    /* ----------------------------------------------------------------------------- */
    else
    {
      PyErr_SetString(PyExc_TypeError, "hdfread: paths must be a list of strings.");
      return NULL;
    }
  }

  /* Close */
  // printf("open=%d\n",H5Fget_obj_count(fid, H5F_OBJ_ALL));
  /*
#if defined(_MPI) && defined(H5_HAVE_PARALLEL)
  if (HDF._ismpi == 1)
  {
    // printf("[%d] - open=%d\n",myRank, H5Fget_obj_count(fid, H5F_OBJ_ALL));
  }
#else
  printf("open=%ld\n",H5Fget_obj_count(fid, H5F_OBJ_ALL));
#endif
  */
  H5Fclose(fid);

  return NodesList;
}


//=============================================================================
// Cree un noeud du pyTree
// node peut etre modifie dans le cas d'un lien
//=============================================================================
PyObject* K_IO::GenIOHdf::createNodePartial(hid_t& node)
{
  /* ***************************************************** */
  /* Declaration */
  int d;
  int dim = 0;
  /* ***************************************************** */

  /* Nom du noeud de l'arbre */
  HDF_Get_Attribute_As_String(node, L3S_NAME, _name);

  /* Remplace les liens eventuellement (change node par le node cible) */
  HDF_Get_Attribute_As_String(node, L3S_DTYPE, _dtype);
  if (strcmp(_dtype, "LK") == 0)
  {
    H5G_stat_t sb; /* Object information */
    herr_t herr = H5Gget_objinfo(node, L3S_LINK, (hbool_t)0, &sb);
    if (herr < 0)
    {
      printf("Error: hdfcgnsread: error opening link file.\n");
      // supprimer le noeud?
    }
    else
    {
      hid_t lid = H5Gopen2(node, L3S_LINK, H5P_DEFAULT);
      H5Gclose(node); //H5Gclose(herr);
      node = lid;
    }
  }
  /* fin des liens */

  /* Type du noeud */
  HDF_Get_Attribute_As_String(node, L3S_LABEL, _type);

  // Valeur du noeud
  dim = 0;
  //npy_intp npy_dim_vals[1]; npy_dim_vals[0] = 1;
  npy_intp npy_dim_vals2[1];
  HDF_Get_Attribute_As_String(node, L3S_DTYPE, _dtype);

  hid_t tid;
  tid = ADF_to_HDF_datatype(_dtype);

  if (strcmp(_dtype, L3T_MT) != 0)
  {
    HDF_Get_DataDimensionsPartial(node, _dims2, DataSpace.Src_Offset,
                                                DataSpace.Src_Stride,
                                                DataSpace.Src_Count ,
                                                DataSpace.Src_Block);
    for (d = 0; d < CGNSMAXDIM; d++)
    { if (_dims2[d] == -1) break; }
    dim = d;

    // inverse les dimensions (indexation fortran)
    L3M_CLEARDIMS(_dims);
    // To doux : Flags for array order C or Fortran ?
    // C'est Faux -> A voir avec Christophe car dans ce cas on recupere la dim de l'utilisateur pas via le HDF...
    for (d = 0; d < dim; d++) _dims[d] = _dims2[dim-d-1];
    // for (d = 0; d < dim; d++) _dims[d] = _dims2[d];
    // for (d = 0; d < dim; d++)  printf("%d \n", _dims[d]);
    // printf("\n");
  }

  /* Init capsule */
  PyObject* v = NULL;

  /* Prepare entry DataSpace */
  hid_t sid = createDataSpaceEntry(node, DataSpace.Dst_Offset,
                                         DataSpace.Dst_Stride,
                                         DataSpace.Dst_Count ,
                                         DataSpace.Dst_Block);

  /* Prepare output DataSpace */
  hid_t mid = createDataSpaceOutput(node, _dims2, DataSpace.Src_Offset,
                                                  DataSpace.Src_Stride,
                                                  DataSpace.Src_Count ,
                                                  DataSpace.Src_Block);

  /* Read the data */
  if (strcmp(_dtype, L3T_I4) == 0)
  {
    v = getArrayI4(node, tid, dim, _dims, mid, sid);
  }
  else if (strcmp(_dtype, L3T_I8) == 0)
  {
    v = getArrayI8(node, tid, dim, _dims, mid, sid);
  }
  else if (strcmp(_dtype, L3T_R4) == 0)
  {
    if (_skeleton == 1) v = getArrayR4Skel(node, tid, dim, _dims);
    else v = getArrayR4(node, tid, dim, _dims, mid, sid);
  }
  else if (strcmp(_dtype, L3T_R8) == 0)
  {
    if (_skeleton == 1) v = getArrayR8Skel(node, tid, dim, _dims);
    else v = getArrayR8(node, tid, dim, _dims, mid, sid);
  }
  else if (strcmp(_dtype, L3T_C1) == 0)
  {
    IMPORTNUMPY;
    char* s = getArrayC1(node, tid, dim, _dims);
    if (dim == 1)
    {
      E_Int l = strlen(s); npy_dim_vals2[0] = l;
      v = PyArray_EMPTY(1, npy_dim_vals2, NPY_CHAR, 1);
      memcpy(PyArray_DATA((PyArrayObject*)v), s, l*sizeof(char));
      free(s);
    }
    else
    {
      npy_intp* npy_dim_vals = new npy_intp[dim];
      E_Int l = 1;
      for (E_Int i = 0; i < dim; i++) 
      { npy_dim_vals[i] = _dims[i]; l = l*_dims[i]; } 
      v = PyArray_EMPTY(dim, npy_dim_vals, NPY_CHAR, 1);
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

  /* Close */
  H5Sclose(mid);
  H5Sclose(sid);

  /** Verbose **/
  // int* ptr = (int*) PyArray_DATA((PyArrayObject * ) v);
  // for(int n=0; n<24; n++)
  // {
  //   printf("data[%d] : %d \n", n , ptr[n]);
  // }

  //if (tid != 0) H5Tclose(tid);
  // PyObject* s = Py_BuildValue("[sOOs]", _name, v, Py_None, _type);
  PyObject* s = Py_BuildValue("[sOOs]", _path, v, Py_None, _type);
  Py_DECREF(v);
  return s;
}

//=============================================================================
// Cree un noeud du pyTree
// node peut etre modifie dans le cas d'un lien
//=============================================================================
PyObject* K_IO::GenIOHdf::createNodePartialContigous(hid_t&    node,
                                                     int       iField,  /* Current field to read */
                                                     int&      nField,  /* Total number of field to read */
                                                     PyObject* data)
{
  /* ***************************************************** */
  /* Declaration */
  int d;
  int dim = 0;

  /* Hdf id */
  hid_t      tid;
  H5G_stat_t sb; /* Object information */
  herr_t     herr;
  /* ***************************************************** */

  /* Nom du noeud de l'arbre */
  HDF_Get_Attribute_As_String(node, L3S_NAME, _name);

  // printf("K_IO::GenIOHdf::createNodePartialContigous \n");
  /* Remplace les liens eventuellement (change node par le node cible) */
  HDF_Get_Attribute_As_String(node, L3S_DTYPE, _dtype);
  // printf("K_IO::GenIOHdf::createNodePartialContigous lk=[%s] \n", _dtype);
  if (strcmp(_dtype, "LK") == 0)
  {
    herr = H5Gget_objinfo(node, L3S_LINK, (hbool_t)0, &sb);
    if (herr < 0)
    {
      printf("Error: hdfcgnsread: error opening link file.\n");
      // supprimer le noeud?
    }
    else
    {
      hid_t lid = H5Gopen2(node, L3S_LINK, H5P_DEFAULT);
      H5Gclose(node); //H5Gclose(herr);
      node = lid;
    }
  }
  /* Fin des liens */

  /* Type du noeud */
  HDF_Get_Attribute_As_String(node, L3S_LABEL, _type);

  // Valeur du noeud
  dim = 0;
  //npy_intp npy_dim_vals[1]; npy_dim_vals[0] = 1;
  HDF_Get_Attribute_As_String(node, L3S_DTYPE, _dtype);

  tid = ADF_to_HDF_datatype(_dtype);

  if (strcmp(_dtype, L3T_MT) != 0)
  {
    HDF_Get_DataDimensionsPartial(node, _dims2, DataSpace.Src_Offset,
                                                DataSpace.Src_Stride,
                                                DataSpace.Src_Count ,
                                                DataSpace.Src_Block);
    for (d = 0; d < CGNSMAXDIM; d++)
    { if (_dims2[d] == -1) break; }
    dim = d;

    // inverse les dimensions (indexation fortran)
    L3M_CLEARDIMS(_dims);
    // To doux : Flags for array order C or Fortran ?
    // for (d = 0; d < dim; d++) _dims[d] = _dims2[dim-d-1]*nField;

    if(DataSpace.Flags[0] == 0)  /** Contigous **/
    {
      /** Contigous Fortran == Interlaced C
       *   So we put at the end the nField to made contigous (in fortran sens )
       *   Finalement le tableau numpy sera (i,j,k,nField)
       */
      dim += 1;
      _dims2[dim-1] = nField;
      for (d = 0; d < dim; d++) _dims[d] = _dims2[dim-d-1];

      /** ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: **/
      /** On rajoute une dimension au dataSpace à la faim (Car HDF est C order -_- )**/
      DataSpace_t DataSpaceSave = DataSpace;

      DataSpace.Src_Offset[0] = (hsize_t) iField;
      DataSpace.Src_Count[0]  = 1;
      DataSpace.Src_Stride[0] = (hsize_t) nField;
      DataSpace.Src_Block[0]  = 1;

      for(int n=1;n<L3C_MAX_DIMS-1;n++)
      {
        DataSpace.Src_Offset[n] = DataSpaceSave.Src_Offset[n-1];
        DataSpace.Src_Count[n]  = DataSpaceSave.Src_Count[n-1] ;
        DataSpace.Src_Stride[n] = DataSpaceSave.Src_Stride[n-1];
        DataSpace.Src_Block[n]  = DataSpaceSave.Src_Block[n-1] ;
      }
      /** ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: **/

    }
    else if(DataSpace.Flags[0] == 1) /** Interlaced  **/
    {

      dim += 1;
      _dims2[dim-1] = nField;
      for (d = 0; d < dim; d++) _dims[d] = _dims2[dim-d-1];

      /** You know the hat trick ? **/
      for (d = 0; d < CGNSMAXDIM; d++) _tmp[d]   = _dims2[d];
      for (d = 0; d < CGNSMAXDIM; d++) _dims2[d] = _dims[d];
      for (d = 0; d < CGNSMAXDIM; d++) _dims[d]  = _tmp[d];

      /** ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: **/
      /** On rajoute une dimension au dataSpace au début (Car HDF est C order -_- ) **/
      DataSpace.Src_Offset[dim-1] = (hsize_t) iField;
      DataSpace.Src_Count[dim-1]  = 1;
      DataSpace.Src_Stride[dim-1] = (hsize_t) nField;
      DataSpace.Src_Block[dim-1]  = 1;
      /** ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: **/
    }
    else if(DataSpace.Flags[0] == 2)  /** Contigous Multiple path **/
    {
      /** 
       *  Just modify the dataSpace ( Memory ) in order to append 
       *  to the current data -> For now only 1D array 
       */
      assert(dim==1);

      /* Verbose */
      // printf("dim      : %d \n", dim);
      // printf("dim Glob : %d \n", DataSpace.GlobDataSetDim[0]);

      _dims2[0] = DataSpace.GlobDataSetDim[0];               /* Assume 1D */
      for (d = 0; d < dim; d++) _dims[d] = _dims2[dim-d-1];
        
      /** ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: **/
      DataSpace.Src_Offset[0] = (hsize_t) nField;
      DataSpace.Src_Count[0]  = (hsize_t) DataSpace.Dst_Count[0];
      DataSpace.Src_Stride[0] = 1;
      DataSpace.Src_Block[0]  = 1;
      
      /** ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: **/
      /* Go to next offset **/
      nField += (int) DataSpace.Dst_Count[0];

    }
    else
    {
      printf("Flags not recognized : 0 = Contigous // 1 = Interlaced \n");
    }

  }

  /* Prepare entry DataSpace */
  hid_t sid = createDataSpaceEntry(node, DataSpace.Dst_Offset,
                                         DataSpace.Dst_Stride,
                                         DataSpace.Dst_Count ,
                                         DataSpace.Dst_Block);

  /* Il faut creer un dim out je pense ? -> Watch out _dims is overwrite here ... */
  hid_t mid = createDataSpaceOutput(node, _dims, DataSpace.Src_Offset,
                                                 DataSpace.Src_Stride,
                                                 DataSpace.Src_Count ,
                                                 DataSpace.Src_Block);

  /* Read the data */
  if (strcmp(_dtype, L3T_I4) == 0)
  {
    data = getArrayContigous(node, tid, dim, _dims2, NPY_INT, mid, sid, data);
  }
  else if (strcmp(_dtype, L3T_I8) == 0)
  {
    data = getArrayContigous(node, tid, dim, _dims2, NPY_LONG, mid, sid, data);
  }
  else if (strcmp(_dtype, L3T_R4) == 0)
  {
    data = getArrayContigous(node, tid, dim, _dims2, NPY_FLOAT, mid, sid, data);
  }
  else if (strcmp(_dtype, L3T_R8) == 0)
  {
    data = getArrayContigous(node, tid, dim, _dims2, NPY_DOUBLE, mid, sid, data);
  }
  else if (strcmp(_dtype, L3T_C1) == 0)
  {
    printf("Contigous load af char * not allowed ...\n");
    exit(1);
  }
  else if (strcmp(_dtype, L3T_MT) == 0)
  {
    data = Py_None; Py_INCREF(Py_None);
  }
  else
  {
    printf("Warning: hdfcgnsread: unknown type of node: %s.\n", _dtype);
  }

  /* Close */
  H5Sclose(mid);
  H5Sclose(sid);

  /** Return new ptr on data **/
  return data;
}

//=============================================================================
/*
   hdfcgnswrite
*/
//=============================================================================
E_Int K_IO::GenIO::hdfcgnsWritePathsPartial(char* file, PyObject* tree,
                                                        PyObject* Filter,
                                                        int skeleton,
                                                        void* comm)
{
  /* ***************************************************** */
  /* Declaration */
  hid_t fapl, fid, ret/*, capl*/;

  PyObject   *key, *DataSpaceDIM;
  Py_ssize_t  pos = 0;
  /* ***************************************************** */

  /* Check */
  if (PyDict_Check(Filter) == false)
  {
    PyErr_SetString(PyExc_TypeError, "hdfread: Filter must be a dict of paths.");
    return 0;
  }

  /* Ouverture du fichier pour l'ecriture */
  fapl = H5Pcreate(H5P_FILE_ACCESS);
  // H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);

  // capl = H5Pcreate(H5P_FILE_CREATE);
  // H5Pset_link_creation_order(capl, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);

  // /* Create the file collectively */
  // fid = H5Fcreate(file, H5F_ACC_TRUNC, capl, fapl);
  // H5Pclose(fapl); H5Pclose(capl);

  GenIOHdf HDF;
  HDF._skeleton = skeleton;

#if defined(_MPI) && defined(H5_HAVE_PARALLEL)
  HDF._ismpi = 1;
#else
  HDF._ismpi = 0;
#endif

#if defined(_MPI) && defined(H5_HAVE_PARALLEL)
  if (HDF._ismpi == 1)
  {
     /* Mpi context */
     MPI_Comm* comm2 = (MPI_Comm*)comm;
     MPI_Info info   = MPI_INFO_NULL;
     ret             = H5Pset_fapl_mpio(fapl, *comm2, info);
  }
#endif  
   
  /* Access to the file collectively */
  fid = H5Fopen(file, H5F_ACC_RDWR, fapl); 
  if (fid < 0)
  {
    PyErr_SetString(PyExc_TypeError, "hdfwritepartial: can not open file.");
    return 1;
  }
  H5Pclose(fapl);

#if defined(_MPI) && defined(H5_HAVE_PARALLEL)
  int nRank, myRank;
  MPI_Comm* comm2 = (MPI_Comm*)comm;
  MPI_Comm_size(*comm2, &nRank);
  MPI_Comm_rank(*comm2, &myRank);
  // printf("[%d] - open Avant =%d\n",myRank, H5Fget_obj_count(fid, H5F_OBJ_ALL));
#endif
  

  while(PyDict_Next(Filter, &pos, &key, &DataSpaceDIM))
  {
    E_Int FilterSize = PyList_Size(DataSpaceDIM);
    /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */
    /* Check if key is String values */
    E_Boolean isKeyString = false;
    if (PyString_Check(key)) isKeyString = true;
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(key)) isKeyString = true; 
#endif
    if (isKeyString == false)
    {
      PyErr_SetString(PyExc_TypeError, "hdfread: paths must be a list of strings.");
      return 0;
    }
    /* Check if key is String values */
    if (PyList_Check(DataSpaceDIM) == false)
    {
      PyErr_SetString(PyExc_TypeError, "hdfread: DataSpaceDIM must be a list of numbers.");
      return 0;
    }
    if(FilterSize != 9)
    {
      PyErr_SetString(PyExc_TypeError, "hdfread: FilterSize must be a list of 9 numbers.");
      return 0;
    }
    /* ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo */

    /* Get path */
    char* path = NULL;
    if (PyString_Check(key)) path = PyString_AsString(key);
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(key)) path = PyBytes_AsString(PyUnicode_AsUTF8String(key));
#endif
    // printf("path to write ...  %s\n", path);

    /* Open group in HDF corresponding to path */
    hid_t gid = HDF.openGroupWithLinks(fid, path);  

    /* GetNodeByPath */
    PyObject* node = K_PYTREE::getNodeFromPath(tree, path);

    /* Fill data space */
    HDF.fillDataSpaceWithFilter(DataSpaceDIM);

    hid_t did = HDF.writeNodePartial(gid, node);

    H5Gclose(gid); // did is gidp
  }

  // printf("open=%d\n", H5Fget_obj_count(fid, H5F_OBJ_ALL));
  H5Fclose(fid);

  // printf("hdfcgnswriteFromPathPartial End \n");
  return 1;
}

//=============================================================================
hid_t K_IO::GenIOHdf::writeNodePartial(hid_t     node,
                                       PyObject *tree)
{
  IMPORTNUMPY;
  /* ***************************************************** */
  /* Declaration */
  char  s1[CGNSMAXLABEL+1];
  char  s2[CGNSMAXLABEL+1];
  // hid_t child;
  /* ***************************************************** */
  PyObject* pname  = PyList_GetItem(tree, 0);
  char* name = NULL;
  if (PyString_Check(pname)) name = PyString_AsString(pname);
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(pname)) name = PyBytes_AsString(PyUnicode_AsUTF8String(pname));
#endif
  PyObject* plabel = PyList_GetItem(tree, 3);
  char* label = NULL;
  if (PyString_Check(plabel)) label = PyString_AsString(plabel);
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(plabel)) label = PyBytes_AsString(PyUnicode_AsUTF8String(plabel));
#endif
  strcpy(s1, name); strcpy(s2, label);

  // printf("writeNodePartial -> tre::name  %s \n", name);
  // printf("writeNodePartial -> tre::label %s \n", label);

  // Creation du noeud, prop. est requis par cgnslib
  // child = H5Gcreate2(node, s1, H5P_DEFAULT, _group, H5P_DEFAULT);
  // HDF_Add_Attribute_As_String(child, L3S_NAME, s1);
  // HDF_Add_Attribute_As_String(child, L3S_LABEL, s2);
  // HDF_Add_Attribute_As_Integer(child, L3S_FLAGS, 1);

  // Ecriture de la valeur
  PyObject* v = PyList_GetItem(tree, 1);

  if (v == Py_None)
  {
    // direct dans l'attribut
    HDF_Add_Attribute_As_String(node, L3S_DTYPE, L3T_MT);
  }
  else if (PyString_Check(v))
  {
    setArrayC1(node, PyString_AsString(v));
    HDF_Add_Attribute_As_String(node, L3S_DTYPE, L3T_C1);
  }
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(v))
  {
    setArrayC1(node, PyBytes_AsString(PyUnicode_AsUTF8String(v)));
    HDF_Add_Attribute_As_String(node, L3S_DTYPE, L3T_C1);
  }
#endif
  else if (PyInt_Check(v))
  {
    setSingleI4(node, PyInt_AsLong(v));
  }
  else if (PyFloat_Check(v))
  {
    if (strcmp(name, "CGNSLibraryVersion") == 0)
      setSingleR4(node, PyFloat_AsDouble(v));
    else
      setSingleR8(node, PyFloat_AsDouble(v));
  }
  else if (PyArray_Check(v))
  {
    PyArrayObject* ar = (PyArrayObject*)v;
    int dim = PyArray_NDIM(ar);
    int* dims = new int [dim];
    //int typeNum = ar->descr->type_num;
    //int elSize = ar->descr->elsize;
    int typeNum = PyArray_TYPE(ar);
    int elSize = PyArray_ITEMSIZE(ar);
    for (int n = 0; n < dim; n++)
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
          setSingleR4(node, (float)ptr[0]);
        }
        else
        {
          double* ptr = (double*)PyArray_DATA(ar);
          setSingleR8(node, ptr[0]);
        }
      }
      else if (typeNum == NPY_INT)
      {
        if (elSize == 4)
        {
          int* ptr = (int*)PyArray_DATA(ar);
          setSingleI4(node, ptr[0]);
        }
        else
        {
          E_LONG* ptr = (E_LONG*)PyArray_DATA(ar);
          setSingleI8(node, ptr[0]);
        }
      }
      else if (typeNum == NPY_CHAR ||
               typeNum == NPY_STRING ||
               typeNum == NPY_BYTE ||
               //typeNum == NPY_SBYTE ||
               typeNum == NPY_UBYTE )
      {
        E_Int diml = PyArray_DIMS(ar)[0];
        char* buf = new char [diml+1];
        strncpy(buf, (char*)PyArray_DATA(ar), diml);
        buf[diml] = '\0';
        setArrayC1(node, buf);
        HDF_Add_Attribute_As_String(node, L3S_DTYPE, L3T_C1);
        delete [] buf;
      }
      else if (typeNum == NPY_LONG)
      {
        if (elSize == 4)
        {
          int* ptr = (int*)PyArray_DATA(ar);
          setSingleI4(node, ptr[0]);
        }
        else
        {
          E_LONG* ptr = (E_LONG*)PyArray_DATA(ar);
          setSingleI8(node, ptr[0]);
        }
      }
      else if (typeNum == NPY_FLOAT)
      {
        float* ptr = (float*)PyArray_DATA(ar);
        setSingleR4(node, ptr[0]);
      }
    }
    else // tableau
    {
      if (typeNum == NPY_DOUBLE)
      {
        // patch pour la norme ADF
        if (strcmp(name, "RotationCenter") == 0 ||
            strcmp(name, "RotationAngle") == 0 ||
            strcmp(name, "RotationRateVector") == 0 ||
            strcmp(name, "Translation") == 0)
        {
          E_Int s = PyArray_Size(v);
          float* buf = new float [s];
          double* ptr = (double*)PyArray_DATA(ar);
          for (int i = 0; i < s; i++) buf[i] = ptr[i];
          setArrayR4(node, buf, dim, dims);
          delete [] buf;
        }
        else
        {
          setArrayPartial(node, (void*)PyArray_DATA(ar), dim, dims,
                          _NATIVE_DOUBLE, (char*)L3T_R8);
        }
      }
      else if (typeNum == NPY_INT)
      {
        if (elSize == 4)
        {
          setArrayPartial(node, (void*)PyArray_DATA(ar), dim, dims,
                                 _NATIVE_INT, (char*)L3T_I4);
        }
        else
        {
          setArrayPartial(node, (void*)PyArray_DATA(ar), dim, dims,
                                 _NATIVE_LONG, (char*)L3T_I8);
        }
      }
      else if (typeNum == NPY_CHAR ||
               typeNum == NPY_STRING ||
               typeNum == NPY_BYTE ||
               //typeNum == NPY_SBYTE ||
               typeNum == NPY_UBYTE )
      {
        E_Int diml = PyArray_DIMS(ar)[0];
        char* buf = new char [diml+1];
        strncpy(buf, (char*)PyArray_DATA(ar), diml);
        buf[diml] = '\0';
        setArrayC1(node, buf);
        HDF_Add_Attribute_As_String(node, L3S_DTYPE, L3T_C1);
        delete [] buf;
      }
      else if (typeNum == NPY_LONG)
      {
        if (elSize == 4)
        {
          setArrayPartial(node, (void*)PyArray_DATA(ar), dim, dims,
                                 _NATIVE_INT, L3T_I4);
        }
        else
        {
          setArrayPartial(node, (void*)PyArray_DATA(ar), dim, dims,
                                 _NATIVE_LONG, L3T_I8);
        }
        //E_Int s = PyArray_Size(v);
        //int* buf = new int [s];
        //for (int i = 0; i < s; i++) buf[i] = (int)PyArray_DATA(ar)[i];
        //setArrayI4(node, buf, dim, dims);
        //delete [] buf;
      }
      else if (typeNum == NPY_FLOAT)
      {
        //E_Int s = PyArray_Size(v);
        //double* buf = new double [s];
        //for (int i = 0; i < s; i++) buf[i] = (double)PyArray_DATA(ar)[i];
        //setArrayR8(node, buf, dim, dims);
        //delete [] buf;
        setArrayPartial(node, (void*)PyArray_DATA(ar), dim, dims,
                               _NATIVE_FLOAT, L3T_R4);
      }
    }
    delete [] dims;
  }
  //H5Pclose(gapl);
  // Py_DECREF(v);
  return node;
}


//=============================================================================
hid_t K_IO::GenIOHdf::setArrayPartial(hid_t node, void* data, int idim, int* idims,
                                      hid_t DataType, char *CGNSType)
{
  //hid_t    acc_tpl1;              /* File access templates */
  hid_t    xfer_plist;            /* Dataset transfer properties list */
  hid_t    sid;                   /* Dataspace ID */
  hid_t    file_dataspace;        /* File dataspace ID */
  hid_t    mem_dataspace;         /* memory dataspace ID */
  hid_t    dataset;               /* Dataset ID */
  hsize_t  dim;
  hsize_t* dims;
  herr_t   ret;                   /* Generic return value */
  /* ***************************************************** */
  /* > Begin */
  // printf("setArrayPartial\n");
  /* TODO : */
  /* A nettoyer idims et LocalDataSetDim c'est pareil */

  dim = idim; dims = (hsize_t*)malloc(sizeof(hsize_t)*dim);
  for (E_Int i = 0; i < idim; i++) {dims[i] = idims[i];}

  /* Manage data type */
  // hid_t tid = H5Tcopy(H5T_NATIVE_DOUBLE); H5Tset_precision(tid, 64);
  hid_t tid = DataType; //H5Tcopy(H5T_NATIVE_DOUBLE); H5Tset_precision(tid, 64);

  /* -------------------------------------------------------------- */
  /* Creation du dataset FICHIER collectif */
  sid = H5Screate_simple(dim, DataSpace.GlobDataSetDim, NULL);
  if (sid < 0) {printf("Fail in setArrayPartial::H5Screate_simple\n");}

#if defined(_MPI) && defined(H5_HAVE_PARALLEL)
  /* Create a dataset collectively */
  dataset   = H5Dcreate2(node, L3S_DATA, tid, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset < 0) {printf("Fail in setArrayPartial::H5Dcreate2\n");}
#else
  // printf("setArrayPartial Sequential \n");
  /* Create a dataset at skeleton write */
  if (_skeleton == 1)
  {
    // printf("setArrayPartial Sequential Skeleton \n");
    dataset   = H5Dcreate2(node, L3S_DATA, tid, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset < 0) {printf("Fail in setArrayPartial::H5Dcreate2\n");}
    ret = H5Dclose(dataset);
    H5Sclose(sid);
    return node;
  }
  else   /* No Skeleton */
  {
    /* Open DataSet previously create with Skeleton */
    // printf("setArrayPartial Sequential Open \n");
    dataset   = H5Dopen2(node, L3S_DATA, H5P_DEFAULT);
  }

#endif
  hid_t mid = H5Tget_native_type(tid, H5T_DIR_ASCEND);
  /* -------------------------------------------------------------- */

  /* -------------------------------------------------------------- */
  /* Create a file dataspace independently */
  file_dataspace = H5Dget_space(dataset);
  if (file_dataspace < 0) {printf("Fail in setArrayPartial::H5Dget_space\n");}

  ret = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET,
                            DataSpace.Dst_Offset, DataSpace.Dst_Stride,
                            DataSpace.Dst_Count, DataSpace.Dst_Block);
  /* -------------------------------------------------------------- */

  /* -------------------------------------------------------------- */
  /* create a memory dataspace independently */
  mem_dataspace = H5Screate_simple(dim, dims, NULL);
  if (mem_dataspace < 0) {printf("Fail in setArrayPartial::H5Screate_simple\n");}

  ret = H5Sselect_hyperslab(mem_dataspace, H5S_SELECT_SET,
                            DataSpace.Src_Offset, DataSpace.Src_Stride,
                            DataSpace.Src_Count, DataSpace.Src_Block);
  /* -------------------------------------------------------------- */

  /* -------------------------------------------------------------- */
  /* Set up the collective transfer properties list */
  xfer_plist = H5Pcreate (H5P_DATASET_XFER);
  if (xfer_plist < 0) {printf("Fail in setArrayPartial::H5Pcreate\n");}

#if defined(_MPI) && defined(H5_HAVE_PARALLEL)
  ret = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
  // ret = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_INDEPENDENT);
  if (ret < 0) {printf("Fail in setArrayPartial::H5Pset_dxpl_mpio\n");}
#endif
  /* -------------------------------------------------------------- */

  /* -------------------------------------------------------------- */
  /* Write data collectively */
  ret = H5Dwrite(dataset, mid, mem_dataspace, file_dataspace, xfer_plist, data);
  if (ret < 0) {printf("Fail in setArrayPartial::H5Dwrite\n");}
  /* -------------------------------------------------------------- */

  /* -------------------------------------------------------------- */
  /* MISC */
  HDF_Set_Attribute_As_String(node, L3S_DTYPE, CGNSType);
  /* -------------------------------------------------------------- */

  /* -------------------------------------------------------------- */
  /* Close */
  H5Sclose(file_dataspace);
  H5Sclose(mem_dataspace);
  H5Pclose(xfer_plist);
  ret = H5Dclose(dataset);
  H5Sclose(sid);
  // H5Tclose(tid);
  H5Tclose(mid);
  /* -------------------------------------------------------------- */

  free(dims);
  return node;
}

/* ------------------------------------------------------------------------- */
// Recopie le filtre (python) en HDF (DataSpace)
void K_IO::GenIOHdf::fillDataSpaceWithFilter(PyObject* Filter)
{
  /* ***************************************************** */
  // C'est une copie ici ou pas ?
  /* ***************************************************** */
  /* Fill Source DataSpace */
  fillArrayLongWithList(Filter, 0, DataSpace.Src_Offset);
  fillArrayLongWithList(Filter, 1, DataSpace.Src_Stride);
  fillArrayLongWithList(Filter, 2, DataSpace.Src_Count );
  fillArrayLongWithList(Filter, 3, DataSpace.Src_Block );

  /* Fill Destination DataSpace */
  fillArrayLongWithList(Filter, 4, DataSpace.Dst_Offset);
  fillArrayLongWithList(Filter, 5, DataSpace.Dst_Stride);
  fillArrayLongWithList(Filter, 6, DataSpace.Dst_Count );
  fillArrayLongWithList(Filter, 7, DataSpace.Dst_Block );

  /* Fill global array */
  fillArrayLongWithList(Filter, 8, DataSpace.GlobDataSetDim);

  /* Fill flags */
  if(PyList_Size(Filter) > 9)
  {
    // printf("Fill additional ... \n");
    fillArrayLongWithList(Filter, 9, DataSpace.Flags);
  }

  // printf(" DataSpace->Dst_Block  %d\n",  (int)DataSpace.GlobDataSetDim[0]);
  // printf(" DataSpace->Dst_Block  %d\n",  (int)DataSpace.Src_Offset[3]);

  /** Permute DataSpace **/
  /** THIS STEP IS MANDATORY BECAUSE DATA SPACE /HYPERSLAB IN HDF IS C ORDER AND CGNS IS FORTRAN ORDER **/
  /** Better make reverse that in interface **/
  DataSpace_t DataSpaceSave = DataSpace;

  int i=0;
  for (int n=0; n<L3C_MAX_DIMS; n++)
  {
    DataSpace.Src_Offset[n] = -1;
    DataSpace.Src_Count[n]  = -1;
    DataSpace.Src_Stride[n] = -1;
    DataSpace.Src_Block[n]  = -1;
    DataSpace.GlobDataSetDim[n] = -1;
    if (DataSpaceSave.Src_Offset[L3C_MAX_DIMS-n-1] != -1)
    {
      DataSpace.Src_Offset[i] = DataSpaceSave.Src_Offset[L3C_MAX_DIMS-n-1];
      DataSpace.Src_Count[i]  = DataSpaceSave.Src_Count[L3C_MAX_DIMS-n-1];
      DataSpace.Src_Stride[i] = DataSpaceSave.Src_Stride[L3C_MAX_DIMS-n-1];
      DataSpace.Src_Block[i]  = DataSpaceSave.Src_Block[L3C_MAX_DIMS-n-1];
      DataSpace.GlobDataSetDim[i]  = DataSpaceSave.GlobDataSetDim[L3C_MAX_DIMS-n-1];
      i++;
    }
  }

  i=0;
  for (int n=0; n<L3C_MAX_DIMS; n++)
  {
    DataSpace.Dst_Offset[n] = -1;
    DataSpace.Dst_Count[n]  = -1;
    DataSpace.Dst_Stride[n] = -1;
    DataSpace.Dst_Block[n]  = -1;
    if (DataSpaceSave.Dst_Offset[L3C_MAX_DIMS-n-1] != -1)
    {
      DataSpace.Dst_Offset[i] = DataSpaceSave.Dst_Offset[L3C_MAX_DIMS-n-1];
      DataSpace.Dst_Count[i]  = DataSpaceSave.Dst_Count[L3C_MAX_DIMS-n-1];
      DataSpace.Dst_Stride[i] = DataSpaceSave.Dst_Stride[L3C_MAX_DIMS-n-1];
      DataSpace.Dst_Block[i]  = DataSpaceSave.Dst_Block[L3C_MAX_DIMS-n-1];
      i++;
    }
  }

  // for(int n=0; n<3; n++)
  // {
  //   printf("----------------------- \n");
  //   printf("DataSpace.Dst_Offset     [%d] : %d \n", n, DataSpace.Dst_Offset[n]     );
  //   printf("DataSpace.Dst_Count      [%d] : %d \n", n, DataSpace.Dst_Count[n]      );
  //   printf("DataSpace.Dst_Stride     [%d] : %d \n", n, DataSpace.Dst_Stride[n]     );
  //   printf("DataSpace.Dst_Block      [%d] : %d \n", n, DataSpace.Dst_Block[n]      );
  //   printf("DataSpace.Src_Offset     [%d] : %d \n", n, DataSpace.Src_Offset[n]     );
  //   printf("DataSpace.Src_Count      [%d] : %d \n", n, DataSpace.Src_Count[n]      );
  //   printf("DataSpace.Src_Stride     [%d] : %d \n", n, DataSpace.Src_Stride[n]     );
  //   printf("DataSpace.Src_Block      [%d] : %d \n", n, DataSpace.Src_Block[n]      );
  //   printf("DataSpace.GlobDataSetDim [%d] : %d \n", n, DataSpace.GlobDataSetDim[n] );
  // }
}

// IN: obj: liste de listes python
// OUT: recopie le no item dans val
void fillArrayLongWithList(PyObject* obj, int item, hsize_t *val)
{
  PyObject* SubList;
  int       n;
  int       s;

  SubList = PyList_GetItem(obj, item);
  for (n = 0; n < PyList_Size(SubList); n++)
  {
    val[n] = PyInt_AsLong(PyList_GetItem(SubList,n));
  }
  s = n;
  for (n = s; n < L3C_MAX_DIMS; n++)
  {
    val[n] = -1;
  }
  // Py_DECREF(SubList); // Variable local donc pas besoin ? -> Surtout pas je Plante ...
}

/* ------------------------------------------------------------------------- */
// Retourne la dimension du tableau partiel dans dims
int HDF_Get_DataDimensionsPartial(hid_t nid, int *dims,
                                  hsize_t *dst_offset,
                                  hsize_t *dst_stride,
                                  hsize_t *dst_count,
                                  hsize_t *dst_block)
{
  int   n;
  int   ndims;
  L3M_CLEARDIMS(dims);
  ndims = 0;
  for (n = 0; dst_count[n] != -1; n++)
  {
    ndims += 1;
    dims[n] = dst_count[n]*dst_block[n];
  }
  dims[ndims] = -1; // sentinelle
  return 1;
}

//=============================================================================
// Cree un dataspace HDF - Entry
//=============================================================================
hid_t createDataSpaceEntry(hid_t nid, hsize_t *src_offset,
                                      hsize_t *src_stride,
                                      hsize_t *src_count,
                                      hsize_t *src_block)
{
  int       /*n, dst_ndims,*/src_ndims;
  hsize_t   src_dim_vals[L3C_MAX_DIMS];
  hid_t     did,sid;
  herr_t    stat;
  
  did = H5Dopen2(nid,L3S_DATA,H5P_DEFAULT);
  sid = H5Dget_space(did);

  src_ndims = H5Sget_simple_extent_ndims(sid);
  H5Sget_simple_extent_dims(sid, src_dim_vals, NULL);

  src_dim_vals[src_ndims]=-1;

  stat = H5Sselect_hyperslab(sid, H5S_SELECT_SET,
                             src_offset, src_stride, src_count, src_block);

  /* Fermeture */
  H5Dclose(did);

  /* CAUTION - sid need to be close after */
  return sid;
}

//=============================================================================
// Cree un dataspace HDF - Output
//=============================================================================
hid_t createDataSpaceOutput(hid_t nid, int     *dst_dims,
                                       hsize_t *dst_offset,
                                       hsize_t *dst_stride,
                                       hsize_t *dst_count,
                                       hsize_t *dst_block)
{
  /* ***************************************************** */
  /* Declaration */
  int       n, dst_ndims;
  hsize_t   dst_dim_vals[L3C_MAX_DIMS];
  hid_t     /*yid,tid,*/mid;
  hssize_t  dst_size;
  herr_t    stat;
  
  dst_size  = 1;
  dst_ndims = 0;

  for (n = 0; dst_dims[n] != -1; n++)
  {
    dst_size *= dst_count[n];
    dst_size *= dst_block[n];
    dst_ndims += 1;
    dst_dim_vals[n] = dst_dims[n];
  }
  dst_dim_vals[dst_ndims]=-1;

  /* Create DataSpace */
  mid  = H5Screate_simple(dst_ndims, dst_dim_vals, NULL);
  stat = H5Sselect_hyperslab(mid, H5S_SELECT_SET,
                             dst_offset, dst_stride, dst_count, dst_block);

  /* CAUTION - mid need to be closed after */
  return mid;
}
