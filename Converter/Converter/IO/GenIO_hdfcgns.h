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

#include "GenIO.h"
#include "kcore.h"
#include "hdf5.h"
#include <map>

// For now, force output of v1.8 of HDF
#if H5_VERSION_GE(1,10,2)
#define KHDFVERSION H5F_LIBVER_V18
#else
#define KHDFVERSION H5F_LIBVER_LATEST
#endif

//# define CGNSMAXLABEL 33
# define CGNSMAXLABEL 128
# define CGNSMAXDIM 20

#define L3S_NAME           "name"
#define L3S_LABEL          "label"
#define L3S_DTYPE          "type"
#define L3S_FLAGS          "flags"
#define L3S_VERSION        " hdf5version"
#define L3S_FORMAT         " format"
#define L3S_DATA           " data"
#define L3S_FILE           " file"
#define L3S_PATH           " path"
#define L3S_LINK           " link"

#define L3S_ROOTNODENAME  "HDF5 MotherNode"
#define L3S_ROOTNODETYPE  "Root Node of HDF5 File"

#define L3C_MAX_DIMS        12
#define L3C_MAX_LINK_DEPTH  100
#define L3C_MAX_SIZE_BUFFER 4096
// Mettre a 128 cette valeur, des que cgnsview peut prendre en compte les labels > 32
#define L3C_MAX_ATTRIB_SIZE 32
#define L3C_MAX_DTYPE       2

#define L3T_MT "MT"
#define L3T_LK "LK"
#define L3T_B1 "B1"
#define L3T_C1 "C1"
#define L3T_I4 "I4"
#define L3T_I8 "I8"
#define L3T_U4 "U4"
#define L3T_U8 "U8"
#define L3T_R4 "R4"
#define L3T_R8 "R8"

// si 0 check la propriete order pour iterate, sinon prend l'ordre natif
#define ITERATE 0

// Usefull inline function
#define L3M_CLEARDIMS(dims) \
  {int __nn; for (__nn=0; __nn < L3C_MAX_DIMS; __nn++){dims[__nn]=-1;};}


/* ------------------------------------------------------------------------- */
typedef struct DataSpace_t
{
  hsize_t Src_Offset[L3C_MAX_DIMS];  /* Begin                 */
  hsize_t Src_Count[L3C_MAX_DIMS];   /* Number of Entry       */
  hsize_t Src_Stride[L3C_MAX_DIMS];  /* Stride Block to Block */
  hsize_t Src_Block[L3C_MAX_DIMS];   /* Number of Block Entry */

  hsize_t Dst_Offset[L3C_MAX_DIMS];  /* Begin                 */
  hsize_t Dst_Count[L3C_MAX_DIMS];   /* Number of Entry       */
  hsize_t Dst_Stride[L3C_MAX_DIMS];  /* Stride Block to Block */
  hsize_t Dst_Block[L3C_MAX_DIMS];   /* Number of Block Entry */

  hsize_t GlobDataSetDim[L3C_MAX_DIMS];   /* Number of Block Entry */

  hsize_t Flags[L3C_MAX_DIMS];   /* Number of Block Entry */

} DataSpace_t;

// --------------------------------------------------------------------------------
namespace K_IO
{
class GenIOHdf
{
  public:
    /* Load one with depth deep */
    PyObject* loadOne(PyObject* tree, int depth, 
                      PyObject* dataShape=NULL, PyObject* links=NULL);

    /* get* of a current HDFObject */
    hid_t* getChildren(hid_t);
    char*  getName(double node);
    char*  getLabel(double node);
    void   getType(double node, char* type, int dim, int* dims);
    bool   isAnodeToSkip();

    /* Create HDF Method */
    PyObject* createNode(hid_t& node, PyObject* dataShape=NULL, PyObject* links=NULL);
    PyObject* createNodePartial(hid_t& node);
    PyObject* createNodePartialContigous(hid_t& node, int iField, int &nField, PyObject* data);

    /* Write HDF Method */
    hid_t writeNode(hid_t node, PyObject* tree);
    hid_t writeNodePartial(hid_t node, PyObject* tree);
    hid_t modifyNode(hid_t node, PyObject* tree);
    hid_t openGroupWithLinks(hid_t start, char* path);

    /* Method to getSingle in HDF */
    int    getSingleI4(hid_t node, hid_t tid);
    E_LONG getSingleI8(hid_t node, hid_t tid);
    float  getSingleR4(hid_t node, hid_t tid);
    double getSingleR8(hid_t node, hid_t tid);

    /* Method to getArray in HDF */
    PyObject* getArrayI8Skel(hid_t node, hid_t tid, int dim, int* dims);
    PyObject* getArrayI4Skel(hid_t node, hid_t tid, int dim, int* dims);
    PyObject* getArrayR8Skel(hid_t node, hid_t tid, int dim, int* dims);
    PyObject* getArrayR4Skel(hid_t node, hid_t tid, int dim, int* dims);

    PyObject* getArrayR8(hid_t node, hid_t tid, int dim, int* dims, hid_t mid=H5S_ALL, hid_t sid=H5S_ALL);
    PyObject* getArrayR4(hid_t node, hid_t tid, int dim, int* dims, hid_t mid=H5S_ALL, hid_t sid=H5S_ALL);
    PyObject* getArrayI4(hid_t node, hid_t tid, int dim, int* dims, hid_t mid=H5S_ALL, hid_t sid=H5S_ALL);
    PyObject* getArrayI8(hid_t node, hid_t tid, int dim, int* dims, hid_t mid=H5S_ALL, hid_t sid=H5S_ALL);
    char* getArrayC1(hid_t node, hid_t tid, int dim, int* dims);

    /* Method for contiguous array **/
    PyObject* getArrayContigous(hid_t node, hid_t tid, int dim, int* dims, int NPYtype, hid_t mid=H5S_ALL, hid_t sid=H5S_ALL, PyObject* r = NULL);

    /* Method to setSingle in HDF */
    hid_t setSingleR4(hid_t node, float  data);
    hid_t setSingleR8(hid_t node, double data);
    hid_t setSingleI4(hid_t node, int    data);
    hid_t setSingleI8(hid_t node, E_LONG data);

    /* Method to setArray in HDF */
    hid_t setArrayR4(hid_t node, float*  data, int dim, int *dims);
    hid_t setArrayR8(hid_t node, double* data, int dim, int *dims);
    hid_t setArrayI4(hid_t node, int*    data, int dim, int *dims);
    hid_t setArrayI8(hid_t node, E_LONG* data, int dim, int *dims);
    hid_t setArrayC1(hid_t node, char*   data, char* label=(char*)L3S_DATA);
    hid_t setArrayC1(hid_t node, char*   data, int dim, int *dims);

    /* Method to setPartialArray in HDF */
    hid_t setArrayPartial(hid_t node, void* data, int idim, int* idims,
                          hid_t DataType, char *CGNSType);

    /* DataSpace Fill */
    void fillDataSpaceWithFilter(PyObject* Filter);

    /* Full dump of a tree */
    PyObject* dumpOne(PyObject* tree, int depth, PyObject* links=NULL);
    hid_t ADF_to_HDF_datatype(const char *tp);

  /* Constructor */
  GenIOHdf()
  {
    /* Fill some attributes */
    _skeleton=0; _maxFloatSize=1e6; _maxDepth=1e6;

    /* Create basic data types used everywhere */
    _NATIVE_FLOAT  = H5Tcopy(H5T_NATIVE_FLOAT ); H5Tset_precision(_NATIVE_FLOAT , 32);
    _NATIVE_DOUBLE = H5Tcopy(H5T_NATIVE_DOUBLE); H5Tset_precision(_NATIVE_DOUBLE, 64);
    _NATIVE_INT    = H5Tcopy(H5T_NATIVE_INT   ); H5Tset_precision(_NATIVE_INT   , 32);
    _NATIVE_LONG   = H5Tcopy(H5T_NATIVE_LONG  ); H5Tset_precision(_NATIVE_LONG  , 64);

    /* Group creation */
    _group = H5Pcreate(H5P_GROUP_CREATE);
    H5Pset_link_creation_order(_group, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);

    /* Prepare DataSpace for partial load  */

  }

  /* Destructor  */
  ~GenIOHdf()
  {
    /* Close all object */
    H5Tclose(_NATIVE_FLOAT);
    H5Tclose(_NATIVE_DOUBLE);
    H5Tclose(_NATIVE_INT);
    H5Tclose(_NATIVE_LONG);
    H5Pclose(_group);
  }

  /* Public attributes   */
  public:
    std::list<hid_t> _fatherStack;
    std::list<std::string> _stringStack;
    std::map<std::string, bool> _skipTypes;
    hid_t _NATIVE_FLOAT;
    hid_t _NATIVE_DOUBLE;
    hid_t _NATIVE_INT;
    hid_t _NATIVE_LONG;
    hid_t _group;
    int   _skeleton;
    int   _maxFloatSize;
    int   _maxDepth;
    char  _type[CGNSMAXLABEL+1];
    char  _name[CGNSMAXLABEL+1];
    char  _dtype[CGNSMAXLABEL+1];
    int   _dims[CGNSMAXDIM];
    int   _dims2[CGNSMAXDIM];
    int   _tmp[CGNSMAXDIM];  /* Just use to swap dims and dims2 */
    /* Store CGNS Path */
    char  *_path;
    std::string _currentPath;

    /* Store the communicator if MPI */
    int  _ismpi;

    /* DataSpace for partial load */
    DataSpace_t DataSpace;
};
}

// --------------------------------------------------------------------------------


void fillArrayLongWithList( PyObject* obj, int item, hsize_t *val);
// void fillDataSpaceWithList( PyObject* obj, DataSpace_t *DataSpace);

int HDF_Get_DataDimensionsPartial(hid_t nid, int     *dims,
                                             hsize_t *dst_offset,
                                             hsize_t *dst_stride,
                                             hsize_t *dst_count,
                                             hsize_t *dst_block);

hid_t createDataSpaceEntry(hid_t nid, hsize_t *src_offset,
                                      hsize_t *src_stride,
                                      hsize_t *src_count,
                                      hsize_t *src_block);

hid_t createDataSpaceOutput(hid_t nid, int     *dst_dims,
                                       hsize_t *dst_offset,
                                       hsize_t *dst_stride,
                                       hsize_t *dst_count,
                                       hsize_t *dst_block);
